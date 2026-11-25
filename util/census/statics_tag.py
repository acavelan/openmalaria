#!/usr/bin/env python3
"""statics_tag.py — tag census/statics.raw.csv rows {init|runtime|checkpointed}
(+ a test-reset column) + run_context + one-line rationale -> census/statics.csv.

Machine evidence per variable: every write site (assignment, compound-assign,
push_back/reset/resize/clear/insert/erase, ++/--, `& stream` serialization) of
the variable name, attributed to its enclosing function via the scope walk
from statics.py:
  * unittest/ file (only for variables a test can even reach)   -> test-reset
  * function matches /checkpoint|restore/ or `& stream`         -> checkpointed
  * function matches /init|load|setup|read|parse|set|main/      -> init
  * anything else                                               -> runtime
Cross-file unqualified mentions are ignored (only qualified Class::name or the
same directory counts) so anonymous-namespace name collisions can't lie.
Tag precedence: checkpointed > runtime > init; const tables with no writes are
'init' (static-initialised). Declarator lines (= init value at the definition
line) are not mutations. run_context=yes is the frozen plan-§5 list
(RUN_CONTEXT below). Ambiguous rows land at confidence=needs-review.

Run: util/census/statics.sh build-p1 && python3 util/census/statics_tag.py
"""
import csv, re, subprocess, collections, importlib.util, os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("st", HERE + '/statics.py')
st = importlib.util.module_from_spec(spec); spec.loader.exec_module(st)

INIT_FN = re.compile(r'\b(init\w*|load\w*|setup|read\w*|parse\w*|set\w*|main)\b')
CKPT_FN = re.compile(r'\b(\w*[Cc]heckpoint\w*|restore\w*)\b')
EXTERN = {'master_RNG', 'runtime'}  # the two header-extern globals

def WRITE(name):
    # name, optionally followed by [subscript]/.member chains, then a mutation:
    # assignment/compound-assign/++/-- or a mutating method call, or `& stream`.
    head = (r'\b' + re.escape(name) + r'\b(?:\[[^\]]*\]|\.\w+)*\s*(?:'
            r'=(?!=)|[-+*/&^%|]=|\+=|-=|\*=|/=|%=|\+\+(?!\+)|--(?!-)|&\s*stream|'
            r'\.(?:push_back|pop_back|reset|resize|clear|emplace(?:_back)?|assign|'
            r'insert|erase|set\w*|loadFromXML|makeConstant|open|close|write|seekp|'
            r'seekg|flush|swap|add\w*|register\w*|save\w*)\s*\()')
    # stream-extraction into the variable: `in >> x[i]`
    ext = r'>>\s*(?:[\w:]+\.)*' + re.escape(name) + r'\b'
    return re.compile(head + '|' + ext)

# run_context=yes — plan §5: the ~55 that become RunContext members
RUN_CONTEXT = [
    (r'sim\.h', r'interval|steps_per_year|years_per_step', 'SimClock/SimData'),
    (r'sim\.h', r's_start|s_end|s_max_human_age|s_t0|s_t1|s_interv', 'SimClock/SimData'),
    (r'sim\.h', r'in_update', 'SimClock in_update assert'),
    (r'random\.(h|cpp)', r'master_RNG', 'master RNG'),
    (r'ModelOptions\.h', r'options', 'ModelOptions u64/PAIRS'),
    (r'mon/Monitoring\.(h|cpp)', r'runtime', 'monitoring sink (MeasureStore)'),
    (r'mon/Continuous\.cpp', r'ctsOStream|streamOffset|registered|toReport|ctsPeriod|duringInit', 'ctsout sink'),
    (r'Genotypes\.(h|cpp)', r'cum_initial_freqs|alleleCodes|nextAlleleCode|genotypes|current_mode|interv_mode|N_genotypes', 'Genotypes registry'),
    (r'Diagnostic\.(h|cpp)', r'diagnostic_set|monitoring_diagnostic', 'Diagnostic registry'),
    (r'Treatments\.h', r'treatments', 'Treatments registry'),
    (r'CMDecisionTree\.cpp', r'decision_library', 'decision-tree registry'),
    (r'LSTMDrugType\.cpp', r'drugTypes|drugTypeNames|drugsInUse', 'drug registry'),
    (r'LSTMTreatments\.cpp', r'schedules|scheduleNames|dosages|dosagesNames', 'drug schedule registry'),
    (r'PopulationAgeStructure\.h', r'ageGroupBounds|ageGroupPercent|cumAgeProp', 'age-structure registry'),
    (r'InterventionManager\.(h|cpp)', r'identifierMap|humanComponents|continuous|timed|nextTimed|importedInfections', 'intervention registries'),
    (r'interventions/(GVI|IRS|ITN)\.h', r'componentsByIndex', 'component type registries'),
    (r'Vaccine\.h', r'params|reportComponent', 'vaccine component registry'),
    (r'CommonWithinHost\.h', r'createInfection|checkpointedInfection', 'within-host factory registry'),
    (r'PerHost\.h', r'params|entoAvailabilityPercentiles', 'transmission scenario params (PerHost statics)'),
    (r'AgeGroupInterpolation\.cpp', r'singleton', 'AgeGroupInterpolation default'),
]

def run_ctx(file, name):
    if name == 'MAX_INFECTIONS':
        return '', 'constexpr clamp -> SlotAllocator in phase 6, not RunContext'
    for rx, names, role in RUN_CONTEXT:
        if re.search(rx, file) and re.fullmatch(names, name):
            s = SURVIVOR.get((file, name))
            return 'yes', 'RunContext: ' + role + (' — ' + s if s else '')
    return 'no', ''

# phase-3 survivors: planned for RunContext, stayed static, legality
# verified (TwoRunsSuite green + T1 bit-exact with them shared):
SURVIVOR = {
    ('model/util/AgeGroupInterpolation.cpp', 'singleton'):
        'stateless dummy (every method throws), never written — no per-run state',
    ('model/util/ModelOptions.h', 'options'):
        'fully reassigned per run by ModelOptions::init (xml/util.cpp:771); argv bits are process-constant',
}

# hand-verified fixes the write-regex cannot see (GSL mutates the workspace
# through the pointer; comparators/dummy are stateless => keep init, just high)
MANUAL = {
    ('model/PkPd/Drug/LSTMDrugConversion.cpp', 134, 'gsl_intgr_conv_wksp'):
        ('runtime', 'high', 'GSL integration workspace allocated at static-init, mutated by gsl_integration_qag on every call'),
    ('model/PkPd/Drug/LSTMDrugThreeComp.cpp', 128, 'gsl_intgr_wksp'):
        ('runtime', 'high', 'GSL integration workspace allocated at static-init, mutated by gsl_integration_qag on every call'),
    ('model/Host/NeonatalMortality.cpp', 56, 'ageLb'):
        ('init', 'high', 'static-init constant (fromYearsI uses constexpr DAYS_IN_YEAR)'),
    ('model/Host/NeonatalMortality.cpp', 56, 'ageUb'):
        ('init', 'high', 'static-init constant (fromYearsI uses constexpr DAYS_IN_YEAR)'),
    ('model/Host/WithinHost/CommonWithinHost.cpp', 291, 'infGenotypeSorter'):
        ('init', 'high', 'stateless comparator functor'),
    ('model/interventions/Deployments.h', 702, 'byDeployTime'):
        ('init', 'high', 'stateless comparator functor'),
    ('model/util/AgeGroupInterpolation.cpp', 56, 'singleton'):
        ('init', 'high', 'never-written dummy interpolator default'),
}

def fn_of(fp, ln, cache):
    if fp not in cache: cache[fp] = st.scopes(fp)
    for k, n in reversed(cache[fp].get(ln, ())):
        if k == 'function': return n
    return '(file-scope)'

def scan_writes(name, f0, ln0, scope, linkage, src, fcache):
    is_extern = linkage == 'external' and name in EXTERN
    testable = scope in ('namespace-global', 'static-class-member')
    w = {'init': [], 'runtime': [], 'checkpointed': [], 'test-reset': []}
    for sf, lines in src.items():
        same_dir = sf.startswith(f0.rsplit('/', 1)[0] + '/')
        qual = (not is_extern) and sf != f0 and not same_dir and not sf.startswith('unittest/')
        wr = WRITE(name)
        qrx = re.compile(r'\b\w+::%s\b' % re.escape(name)) if qual else None
        for i, txt in enumerate(lines, 1):
            if sf == f0 and i == ln0: continue  # declarator/initializer, not mutation
            if sf.startswith('unittest') and not testable: continue
            if name not in txt: continue
            c = txt.find('//')
            if 0 <= c < txt.find(name): continue
            if not wr.search(txt): continue
            if qual and not qrx.search(txt): continue
            if sf.startswith('unittest'):
                w['test-reset'].append((sf, i)); continue
            fn = fn_of(sf, i, fcache)
            if CKPT_FN.search(fn) or re.search(r'\b%s\b\s*&\s*stream' % re.escape(name), txt):
                w['checkpointed'].append((sf, i))
            elif INIT_FN.search(fn): w['init'].append((sf, i))
            else: w['runtime'].append((sf, i))
    return w

def main():
    rows = list(csv.DictReader(open('census/statics.raw.csv', encoding='utf-8')))
    src = {}
    for root in ('model', 'unittest'):
        for dirpath, _, names in os.walk(root):
            for n in names:
                if n.endswith(('.cpp', '.h', '.hpp')):
                    p = os.path.join(dirpath, n)
                    src[p] = open(p, encoding='utf-8', errors='replace').read().splitlines()
    fcache = {}
    out = open('census/statics.csv', 'w', encoding='utf-8')
    out.write('file,line,name,type,scope,owner_class,namespace,linkage,definition,init,'
              'other_sites,tag,test_reset,run_context,rationale,confidence\n')
    stats = collections.Counter()
    needs_review = []
    for r in rows:
        f, ln, name = r['file'], int(r['line']), r['name']
        const = r['type'].split() and r['type'].split()[0] in ('const', 'constexpr', 'inline') \
            or 'const' in r['type'].split()[:2]
        if r['scope'] == 'not-in-build':
            tag, conf = 'n/a', 'needs-review'
            rationale = ('file not compiled (dead source) — see census/latent_bugs.md'
                         if 'DescriptiveIPT' in f else 'debug-build-only (#ifndef NDEBUG)')
            w = {k: [] for k in ('init', 'runtime', 'checkpointed', 'test-reset')}
        elif name in EXTERN and f.split('/')[-1].startswith('random'):
            tag, conf = 'checkpointed', 'high'
            rationale = 'seeded main.cpp:158; drawn everywhere; state via model/checkpoint.cpp:84,128'
            w = {k: [] for k in ('init', 'runtime', 'checkpointed', 'test-reset')}
            w['test-reset'] = [('unittest', 0)]
        else:
            w = scan_writes(name, f, ln, r['scope'], r['linkage'], src, fcache)
            kinds = {k for k, v in w.items() if v and k != 'test-reset'}
            if 'checkpointed' in kinds: tag = 'checkpointed'
            elif 'runtime' in kinds: tag = 'runtime'
            elif kinds == {'init'}: tag = 'init'
            elif not kinds: tag = 'init'
            else: tag = 'runtime'
            conf = 'high' if (kinds or const or w['test-reset']) else 'needs-review'
            if len(kinds) > 1 and 'init' in kinds and 'runtime' in kinds:
                # set-at-init AND mutated: fine (init-time params mutated later
                # are exactly runtime state), still high if both classified
                conf = 'high'
            bits = []
            for k in ('init', 'runtime', 'checkpointed', 'test-reset'):
                if w[k]:
                    loc = ';'.join('%s:%d' % t for t in w[k][:2])
                    bits.append('%s-write@%s' % (k, loc))
            if not bits:
                bits.append('const/static-init table (no writes)' if const
                            else 'no writes found in model/ (read-only after construction or via alias)')
            rationale = ' '.join(bits)
        manual = MANUAL.get((f, ln, name))
        if manual: tag, conf, rationale = manual
        rc, rcnote = run_ctx(f, name)
        if rcnote: rationale = (rcnote + ' | ' + rationale)[:200]
        test_reset = 'yes' if w['test-reset'] else 'no'
        out.write(','.join('"%s"' % str(c).replace('"', "'") for c in
            [f, ln, name, r['type'], r['scope'], r['owner_class'], r['namespace'],
             r['linkage'], r['definition'], r['init'], r['other_sites'], tag,
             test_reset, rc, rationale[:200], conf]) + '\n')
        stats[tag] += 1
        stats['rc:' + (rc or 'n/a')] += 1
        if conf == 'needs-review': needs_review.append('%s:%s %s [%s] %s' % (f, ln, name, tag, rationale[:60]))
    out.close()
    print(dict(stats))
    print('needs-review:', len(needs_review))
    with open('/tmp/opencode/p1/needs_review.txt', 'w') as fh:
        fh.write('\n'.join(needs_review) + '\n')

if __name__ == '__main__':
    main()
