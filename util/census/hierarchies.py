#!/usr/bin/env python3
"""hierarchies.sh support: machine-extract polymorphic class hierarchies in model/.

Method: comment-stripped regex parse of class/struct declarators (name + base
list + body brace-matched) over model/**/*.h|cpp; a class is polymorphic if its
body (or an inherited model base's body) declares `virtual`. Roots = polymorphic
model classes whose bases are not model-polymorphic. LOC per hierarchy = sum of
cloc code lines over the files defining its classes.

Usage: python3 hierarchies.py            (run from repo root; writes census/hierarchies.csv)
"""
import os, re, json, subprocess, collections, csv

SRC = []
for root, dirs, names in os.walk('model'):
    for n in names:
        if n.endswith(('.h', '.cpp')):
            SRC.append(os.path.join(root, n))

def strip_comments(text):
    out, i, n, res = [], 0, len(text), []
    in_c = in_s = in_cmt = in_line = False
    while i < n:
        c = text[i]
        if in_cmt:
            if text.startswith('*/', i): in_cmt = False; i += 2; continue
            i += 1; continue
        if in_line:
            if c == '\n': in_line = False; res.append(c)
            i += 1; continue
        if in_s:
            res.append(c)
            if c == '\\': res.append(text[i+1]); i += 2; continue
            if c == '"': in_s = False
            i += 1; continue
        if c == '/' and text.startswith('/*', i): in_cmt = True; i += 2; continue
        if c == '/' and text.startswith('//', i): in_line = True; i += 2; continue
        if c == '"': in_s = True
        res.append(c); i += 1
    return ''.join(res)

CLS = re.compile(r'\bclass\s+([A-Za-z_]\w*)(?:\s+(?:final|OM_\w+))?\s*(?::([^;{]*))?\s*\{', re.S)
STR = re.compile(r'\bstruct\s+([A-Za-z_]\w*)\s*(?::([^;{]*))?\s*\{', re.S)

classes = {}   # name -> dict
def body_at(text, brace_idx):
    depth = 0; i = brace_idx
    while i < len(text):
        if text[i] == '{': depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0: return text[brace_idx:i+1], text[:i].count('\n') + 1
        i += 1
    return text[brace_idx:], text[:i].count('\n') + 1

for f in SRC:
    text = strip_comments(open(f, encoding='utf-8', errors='replace').read())
    for rx in (CLS, STR):
        for m in rx.finditer(text):
            name, bases_txt = m.group(1), (m.group(2) or '')
            body, _ = body_at(text, text.index('{', m.end() - 1))
            bases = []
            for b in bases_txt.split(','):
                b = re.sub(r'\b(public|private|protected|virtual)\b', '', b).strip()
                b = re.sub(r'<.*', '', b).split('::')[-1].strip()
                if b: bases.append(b)
            if name in classes and len(body) < len(classes[name]['body']):
                continue
            classes[name] = dict(name=name, file=f, bases=bases, body=body,
                                 line=text[:m.start()].count('\n') + 1)

def is_poly(c):
    return re.search(r'\b(virtual\b|override\b)', c['body']) is not None

poly = {n for n, c in classes.items() if is_poly(c)}
# namespace by path
def ns_of(c):
    m = re.search(r'namespace\s+(\w+)', c['body'][:0] or '')  # unused
    p = c['file']
    for ns in ('interventions', 'Clinical', 'Transmission', 'Host', 'mon', 'PkPd', 'util'):
        if '/%s/' % ns in p or p.endswith('/%s.h' % ns) or p.endswith('/mon/%s' % ns): pass
    return p

def ancestors(n, seen=None):
    seen = seen or set()
    out = set()
    for b in classes.get(n, {}).get('bases', []):
        if b in classes and b in poly and b not in seen:
            seen.add(b); out.add(b); out |= ancestors(b, seen)
    return out

roots = {n for n in poly if not any(b in poly and b in classes for b in classes[n]['bases'])}
# merge roots that are in the same file family? no — keep AST truth, curate below

children = collections.defaultdict(list)
for n in poly:
    for b in classes[n]['bases']:
        if b in poly: children[b].append(n)

def tree(n):
    out = [n]
    for c in children.get(n, []):
        out += tree(c)
    return out

# cloc per file
cloc = json.loads(subprocess.run(['cloc', '--json', '--quiet', '--by-file'] + SRC,
                   capture_output=True, text=True).stdout)
def fileloc(f):
    e = cloc.get(f) or cloc.get('SUM')
    return e['code'] if e and isinstance(e, dict) and 'code' in e else 0

# phase/annotation curation — plan §7 rows 4a/4b/5/7:
# clinical decision tree = 4a, interventions = 7, rest = 4b; phase-5 = cycle
# break touching interventions bases in PerHost.h; util/vtable notes are judgment.
CURATE = {
    'AgeGroupInterpolation': ('4b', 'yes', 'linear/none only (cspline dead, F6) -> enum + free fn; virtual eval is a per-human-day cost (phase 6/8 scratch cols)'),
    'AnophelesModel': ('4b', 'yes', 'SimpleMPD -> optional member + 4 guarded hooks (§5); solver flip at phase 8'),
    'CMDecisionTree': ('4a', 'yes', 'compiled decision tree + flat node-array interpreter (§5); child order = XSD sequence (contract d)'),
    'ClinicalModel': ('4a/4b', 'yes', '3-config -> 2-state collapse (§5); DT interpreter owned by 4a'),
    'DecayFunction': ('4b', 'no', 'exponential/extrap piecewise -> enum + free fn'),
    'EmergenceRateEstimator': ('4b', 'no', 'solver vs fitter -> tag+switch; ATSB secant path phase 8/9'),
    'HumanInterventionComponent': ('7', 'yes', 'dense components by type-rank + date-cursor (§5); bases cycle out of PerHost.h at phase 5'),
    'Infection': ('4b', 'yes', 'slots <=21/host behind ONE SlotAllocator; variant tag columns (§5)'),
    'InfectionIncidenceModel': ('4b', 'no', 'merged incidence 4-config -> 1 (§5)'),
    'LSTMDrug': ('4b', 'no', 'PK variants as enum + free fn formulas (PkPd kept, rewritten flat §2)'),
    'PathogenesisModel': ('4b', 'no', 'variant switch; submodels flattened'),
    'PerHostInterventionData': ('5->7', 'yes', 'interventions bases out of PerHost.h (5); Copy effect structs on sidecar (7)'),
    'Sampler': ('4b', 'no', 'gamma/lognormal -> enum dispatch; inline GSL ports phase 9 (F2 gates)'),
    'TimedDeployment': ('7', 'yes', 'flat TimedDeploymentKind enum + payload struct (§9)'),
    'TransmissionModel': ('4b', 'no', 'PerHost class deleted (5); dense per-host arrays; one vector LocalRNG stream (§9)'),
    'WHInterface': ('4b', 'yes', '69 virtuals = largest hierarchy; variants as equations over shared machinery (§5)'),
    'base_exception': ('none', 'no', 'util error types; not model state — no rewrite'),
}

rows = []
for r in sorted(roots):
    members = tree(r)
    if len(members) == 1 and r in ('base_exception',):  # exceptions: keep (tiny)
        pass
    files = sorted({classes[m]['file'] for m in members})
    loc = sum(fileloc(f) for f in files)
    virts = sum(len(re.findall(r'\bvirtual\b', classes[m]['body'])) for m in members)
    rows.append(dict(root=r, derived=len(members) - 1, members=len(members),
                     virts=virts, loc=loc,
                     files=';'.join(files),
                     depth=max(len(ancestors(m)) + 1 for m in members)))
w = csv.DictWriter(open('census/hierarchies.csv', 'w'),
                   fieldnames=['root', 'derived_count', 'members', 'virtuals',
                               'vtable_heavy', 'total_code_LOC', 'depth',
                               'planned_phase', 'config_state_note', 'files'])
w.writeheader()
for r in rows:
    phase, heavy, note = CURATE.get(r['root'], ('4b', 'no', ''))
    w.writerow(dict(root=r['root'], derived_count=r['derived'], members=r['members'],
                    virtuals=r['virts'], vtable_heavy=heavy,
                    total_code_LOC=r['loc'], depth=r['depth'],
                    planned_phase=phase, config_state_note=note, files=r['files']))
print(len(rows), 'polymorphic roots ->', 'census/hierarchies.csv')

