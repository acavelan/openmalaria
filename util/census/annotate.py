#!/usr/bin/env python3
"""annotate.py — insert `// parity:` comment lines at every draw/contract site
from census/draws.csv (plan §7 row 1; user ruling 2026-11-19: comments are
LOC-exempt). COMMENT-ONLY: inserts a whole line above the site; never touches
existing lines. Idempotent (skips if identical comment within 3 lines above).

Grouping: consecutive same-file, same-consumer, same-contract sites within
GROUP_GAP lines get ONE comment above the first (tight-loop rule).

Usage: python3 util/census/annotate.py [--dry]
"""
import csv, sys, re

GROUP_GAP = 20
DRY = '--dry' in sys.argv

# one-line annotation text per (kind/stream) and per contract letter
CONTRACT_TEXT = {
    'a': 'a — deploy short-circuit: RNG draw is the LAST test of the deploy '
         'condition; filter/age/cohort checks must precede it (plan §3a)',
    'b': 'b — stable_sort = reproducibility contract; must stay stable (plan §3b)',
    'c': 'c — history-ordered activeComponents FP product; keep vector order, '
         'no reassociation (plan §3c)',
    'd': 'd — decision-tree contract: child order = XSD sequence, cumulative '
         'keys, one uniform_01 per Random node (plan §3d)',
    'e': 'e — per-day Poisson/bernoulli draw counts frozen (plan §3e)',
    'f': 'f — validation completes BEFORE master_RNG.seed(): load->check->seed '
         '(plan §3f); never seed before checks',
    'g': 'g — defs[] index = internal Measure id; outId order = output column '
         'order (plan §3g)',
}
# site-specific overrides (substring of file, substring of consumer/call)
OVERRIDE = [
    ('Deployments.h', 'TimedHumanDeployment', 'a — per-human heterogeneity '
     'lognormal/beta draws (gsl_cdf pure quantiles, no stream cost) then '
     'bernoulli-last; per-candidate draw order/count frozen (plan §3a)'),
    ('CMDecisionTree.cpp', 'cumulative', 'd — cumulative keys summed in XSD '
     'outcome order; map<double> numeric key order; last key must be 1 (plan §3d)'),
    ('CMDecisionTree.cpp', 'child chain', 'd — child chain built in XSD '
     'element-sequence order (plan §3d)'),
    ('ITN.cpp', '', 'e — per-update hole/rip Poisson draw counts frozen '
     '(≤1+1 draws per protected human-day; skip when rates 0) (plan §3e)'),
    ('ImportedInfections.h', '', 'e — bernoulli import draw per human-day; '
     'rate->draw->importInfection(order) frozen (plan §3e)'),
    ('DecisionTree5Day.cpp', '', 'd — clinical decision fraction draw order '
     '(plan §3d)'),
    ('checkpoint.cpp', '', 'order — checkpoint serialises master RNG state; '
     'v2 stream format frozen (do not reorder with other state)'),
    ('random.cpp', 'master_RNG', 'order — the master stream: static-init here, '
     'seeded in main after validation (f); checkpointed in model/checkpoint.cpp'),
    ('random.h', 'master_RNG', 'order — single extern master stream declaration'),
    ('Human.cpp', 'master_RNG', 'order — each Human copies+seeds its LocalRng '
     'from the master at construction: human-creation order == master '
     'consumption order; never reorder Human ctor'),
    ('VectorModel.cpp', 'master_RNG', 'order — vector-model LocalRng copy-seeded '
     'from master at construction; construction order frozen'),
]
STREAM_TEXT = {
    'host-local': 'order — per-host Xoshiro stream: this call’s position in '
                  'the draw sequence fixes the host RNG state',
    'local': 'order — local RNG stream: draw call order fixes stream state',
    'gsl-shared': 'order — pure GSL CDF/quantile call (no stream consumption); '
                  'keep the identical GSL function (F2)',
    'fixed-seed': 'order — MasterRng(0,0) availability-percentile seed is '
                  'load-bearing for goldens — do not “clean up” (PerHost.h:86-89)',
    'master': 'order — master stream: seed/checkpoint only, never direct draws',
}
DEF_LETTER = {'bernoulli': 'a/d/e', 'uniform': 'a/d/e'}

def text_for(r):
    if r['call'] == 'bernoulli' and r['contract'] == 'a':
        return CONTRACT_TEXT['a']  # the last-test itself needs its own wording
    for fsig, csig, txt in OVERRIDE:
        if fsig in r['file'] and (csig == '' or csig in r['consumer'] or
                                  csig in r['note'] or csig in r['call']):
            return txt
    if r['contract'] and r['kind'] != 'order':
        return CONTRACT_TEXT[r['contract']]
    if r['kind'] == 'order':
        return (r['contract'] + ' — ' if r['contract'] else 'order — ') + r['note']
    if r['kind'] == 'rng-def':
        letter = DEF_LETTER.get(r['call'], 'F2' if 'GSL' in r['note'] else 'order')
        return letter + ' — ' + r['note']
    return STREAM_TEXT[r['stream']]

rows = list(csv.DictReader(open('census/draws.csv', encoding='utf-8')))
rows.sort(key=lambda r: (r['file'], int(r['line'])))

# group: consecutive sites same file+consumer-ish+text-window -> one insert
targets = []  # (file, line, text, consumer) merged
last = None   # (file, last-absorbed line, letter-prefix, consumer) sliding window
for r in rows:
    txt = text_for(r)
    if last and last[0] == r['file'] and int(r['line']) - last[1] <= GROUP_GAP and \
       last[2] == txt.split(' — ')[0] and last[3] == r['consumer'] and \
       r['kind'] in ('rng-call', 'gsl-direct') and r['call'] != 'bernoulli':
        last = (r['file'], int(r['line']), last[2], r['consumer'])
        continue  # covered by the comment above the group head
    if targets and targets[-1][0] == r['file'] and targets[-1][1] == int(r['line']):
        if targets[-1][2].split(' — ')[0] != txt.split(' — ')[0]:
            targets[-1][2] += ' | ' + txt  # same line, different contract letter
        continue  # same letter = same statement -> one comment suffices
    targets.append([r['file'], int(r['line']), txt, r['consumer']])
    last = (r['file'], int(r['line']), txt.split(' — ')[0], r['consumer'])

inserted = skipped = 0
byfile = {}
for f, ln, txt, _c in targets:
    byfile.setdefault(f, []).append((ln, txt))
for f, ts in byfile.items():
    lines = open(f, encoding='utf-8').read().splitlines(keepends=True)
    for ln, txt in sorted(ts, reverse=True):  # bottom-up: line numbers stay valid
        at = ln - 1
        prev = lines[max(0, at - 3):at]
        if any(l.lstrip().startswith('// parity:') and txt[:40] in l for l in prev):
            skipped += 1
            continue
        if at > 0 and lines[at - 1].rstrip().endswith('\\'):
            print('SKIP macro continuation', f, ln); skipped += 1; continue
        indent = re.match(r'\s*', lines[at]).group(0)
        eol = '\r\n' if lines[at].endswith('\r\n') else '\n'
        lines.insert(at, indent + '// parity: ' + txt + eol)
        inserted += 1
    if not DRY:
        open(f, 'w', encoding='utf-8', newline='').write(''.join(lines))
print(('DRY ' if DRY else '') + 'targets:', len(targets), 'inserted:', inserted,
      'skipped(existing/skip):', skipped)
