#!/usr/bin/env python3
"""draws.py — census/draws.csv: every RNG draw site in model/ + frozen-contract
ordering sites (plan §3 a–g), the machine map the `// parity:` annotations and
phase-2..9 drift checks run against.

Kinds:
  rng-def    distribution primitives in util/random.h (define values-consumed)
  rng-call   call of one of those methods (receiver.classified stream)
  master     master_RNG seed/checkpoint touches + MasterRng(0,0) fixed streams
  gsl-direct direct gsl_cdf_*/gsl_ran_* use outside random.h
  order      non-draw frozen-ordering contract sites (stable_sorts, FP product
             order, validation-before-seed, measure-id order)

contract column = plan §3 letter when the site belongs to a–g, else ''.
Usage: python3 util/census/draws.py   (from repo root; repo-wide greps)
"""
import os, re, csv, collections

METHODS = ('uniform_01|gauss|gamma|log_normal|max_multi_log_normal|beta|'
           'betaWithMean|poisson|bernoulli|uniform|weibull')
CALL = re.compile(r'\b(\w+)(?:\.|->)(?:rng\.|rng->)?(' + METHODS + r')\s*\(')
DEF = re.compile(r'^\s*(?:inline\s+)?(?:double|bool|int)\s+(' + METHODS + r')\s*\(')
GSL = re.compile(r'\b(gsl_(?:ran|cdf)_\w+)\s*\(')

# plan §3 frozen contracts, by file + line-window (curated from §3 + code read)
def contract_for(f, ln, call, kind):
    if f == 'model/util/Deployments.h' or f.endswith('/Deployments.h'):
        return 'a' if kind == 'rng-call' or call.startswith('gsl_cdf') else ''
    if f.endswith('CMDecisionTree.cpp') or f.endswith('DecisionTree5Day.cpp'):
        return 'd'
    if f.endswith('ITN.cpp') and call in ('poisson', 'bernoulli', 'uniform_01'):
        return 'e'
    if f.endswith('ImportedInfections.h'):
        return 'e'
    if f.endswith('main.cpp') and kind == 'master':
        return 'f'
    return ''

NOTE = {  # random.h primitives: what the definition freezes
    'uniform_01': 'one gen_double per call = one stream position; basis of every contract',
    'bernoulli': 'consumes exactly ONE uniform_01 (prob compare) — basis of contracts a/d/e',
    'uniform': 'consumes exactly ONE uniform_01 (scaled int) — do not reorder vs bernoulli',
    'max_multi_log_normal': 'exactly ONE uniform_01 + one ugaussian_Pinv; CDF-inversion identity frozen',
    'gauss': 'GSL gsl_ran_gaussian on wrapped stream — value sequence depends on GSL version (F2)',
    'gamma': 'GSL gsl_ran_gamma on wrapped stream (F2)',
    'log_normal': 'GSL gsl_ran_lognormal on wrapped stream (F2)',
    'beta': 'GSL gsl_ran_beta on wrapped stream (F2)',
    'betaWithMean': 'wraps beta() — one beta draw; a=m*b/(1-m) algebra frozen',
    'poisson': 'GSL gsl_ran_poisson on wrapped stream (F2)',
    'weibull': 'GSL gsl_ran_weibull on wrapped stream (F2)',
}

def consumer(f, ln, fcache):
    if fcache.get('f') != f:
        fcache.clear(); fcache['f'] = f
        fcache['ctx'] = st.scopes(f)
    ctx = fcache['ctx'].get(ln, ())
    owner = next((n for k, n in reversed(ctx) if k in ('class', 'struct')), '')
    fn = next((n for k, n in reversed(ctx) if k == 'function'), '')
    return (owner + '::' if owner else '') + (fn or '(file-scope)')

import importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("st", HERE + '/statics.py')
st = importlib.util.module_from_spec(spec); spec.loader.exec_module(st)

rows, fcache = [], {}
src = {}
for root, _, names in os.walk('model'):
    for n in names:
        if n.endswith(('.cpp', '.h')):
            p = os.path.join(root, n)
            src[p] = open(p, encoding='utf-8', errors='replace').read().splitlines()

for p, lines in sorted(src.items()):
    for i, txt in enumerate(lines, 1):
        c = txt.find('//')
        base = txt if c < 0 else txt[:c]
        if base.strip().startswith(('*', '///')):
            continue
        kind = call = receiver = ''
        if p == 'model/util/random.h':
            m = DEF.match(base)
            if m:
                kind, call = 'rng-def', m.group(1)
                receiver = 'self'
        if not kind and 'master_RNG' in base:
            if re.search(r'\bmaster_RNG\s*\.\s*(seed|checkpoint|gen_seed)', base):
                kind, call = 'master', re.search(r'master_RNG\s*\.\s*(\w+)', base).group(1)
                receiver = 'master_RNG'
            elif 'master_RNG' in base and re.search(r'\bmaster_RNG\b', base):
                kind, call = 'master', 'master_RNG'
                receiver = 'master_RNG'
        if not kind and re.search(r'\bMasterRng\s+\w+\s*\(', base):
            kind, call = 'master', 'MasterRng-fixed-seed'
            receiver = re.search(r'MasterRng\s+(\w+)', base).group(1)
        if not kind:
            m = CALL.search(base)
            if m and 'scnXml' not in base:
                kind, receiver, call = 'rng-call', m.group(1), m.group(2)
        if not kind:
            m = GSL.search(base)
            if m and p != 'model/util/random.h':
                kind, call = 'gsl-direct', m.group(1)
                receiver = 'gsl'
        if not kind:
            continue
        stream = ('master' if kind == 'master' and receiver == 'master_RNG' else
                  'fixed-seed' if kind == 'master' else
                  'host-local' if receiver in ('human', 'hostData') or '.rng' in base else
                  'local' if kind in ('rng-call', 'rng-def') else
                  'gsl-shared' if kind == 'gsl-direct' else 'mixed')
        rows.append([p, i, kind, call, receiver, stream,
                     consumer(p, i, fcache) if kind != 'rng-def' else 'util::RNG',
                     contract_for(p, i, call, kind),
                     NOTE.get(call, '') if kind in ('rng-def',) else ''])

# ---- non-draw ordering contract sites (curated; plan §3 b, c, f, g) ----
ORDER = [  # (file, line-or-None+marker, letter, note, marker)
    ('model/interventions/InterventionManager.cpp', None, 'b',
     'stable_sort order = reproducibility contract; comparator ties resolve by insertion order',
     r'stable_sort'),
    ('model/Transmission/PerHost.cpp', None, 'c',
     'history-ordered activeComponents FP product — multiply in vector order, '
     'not reassociated', r'^double PerHost::entoAvailabilityHetVecItv'),
    ('model/main.cpp', None, 'f',
     'ALL validation (XMLChecker + option rules) completes BEFORE master_RNG.seed() '
     '(load->check->seed); changing validation timing changes observable error order',
     r'master_RNG\s*\.\s*seed'),
    ('model/mon/OutMeasures.h', None, 'g',
     'defs[] index == internal Measure id; outId order == output column order — '
     'string->enum migration must preserve', r'defs'),
    ('model/Transmission/PerHost.h', None, '',
     'MasterRng(0,0) availability-percentile seed is load-bearing for goldens — do not "clean up"',
     r'MasterRng\s+seed\s*\(\s*0\s*,\s*0'),
    ('model/Clinical/CMDecisionTree.cpp', None, 'd',
     'cumulative keys summed in XSD outcome order (map<double> numeric key order)',
     r'cum_p \+= outcome'),
    ('model/Clinical/CMDecisionTree.cpp', None, 'd',
     'child chain order == XSD sequence of decision-tree element choices',
     r'self->children\.push_back'),
]
for f, ln, letter, note, marker in ORDER:
    if ln is None:  # locate by marker scan (comment lines excluded: markers
        for i, txt in enumerate(src[f], 1):  # must hit code, not own insertions)
            if txt.lstrip().startswith('//'):
                continue
            if re.search(marker, txt):
                ln = i; break
    rows.append([f, ln, 'order', '(ordering)', '', 'n/a',
                 consumer(f, ln, fcache), letter, note])

rows.sort(key=lambda r: (r[0], int(r[1]) if r[1] else 0))
w = csv.writer(open('census/draws.csv', 'w', encoding='utf-8'))
w.writerow(['file', 'line', 'kind', 'call', 'receiver', 'stream', 'consumer',
            'contract', 'note'])
for r in rows:
    w.writerow(r)
print(collections.Counter(r[2] for r in rows), 'total', len(rows),
      '-> census/draws.csv')
