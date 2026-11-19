#!/usr/bin/env python3
"""xmlglue.py — census/xml_glue.csv: every generated scnXml type referenced from
model/, its read sites, the model consumers, and the planned model/xml/ target
file (plan §5: conversion glue per XSD module, XSD codegen KEPT — §8 F7).

Rows:
  kind=type    one per referenced scnXml type: module, type, target_xml_file,
               ref_count, sites (file:line[:8]), consumers (enclosing model
               classes/functions), review
  kind=include one per `#include "schema/*.h"` in model/ (type=header, sites=file:line)

Consumer attribution reuses the brace-depth scope walk from statics.py.
Usage: python3 util/census/xmlglue.py   (from repo root)
"""
import os, re, csv, glob, importlib.util, collections

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("st", HERE + '/statics.py')
st = importlib.util.module_from_spec(spec); spec.loader.exec_module(st)

# ---- generated schema: type -> module (schema header without .h) ----
decl = re.compile(r'^\s*class\s+(\w+)\s*[;:]')
type_mod = {}
for h in sorted(glob.glob('build-p1/schema/*.h')):
    mod = os.path.basename(h)[:-2]
    for line in open(h, encoding='utf-8', errors='replace'):
        m = decl.match(line)
        if m:
            type_mod.setdefault(m.group(1), mod)

# ---- model-side references ----
refs = collections.defaultdict(list)   # type -> [(file, line)]
includes = []                          # (file, line, header)
src = {}
for root, _, names in os.walk('model'):
    for n in names:
        if not n.endswith(('.cpp', '.h')): continue
        p = os.path.join(root, n)
        lines = open(p, encoding='utf-8', errors='replace').read().splitlines()
        src[p] = lines
        for i, txt in enumerate(lines, 1):
            m = re.search(r'#include\s+"schema/(\w+)\.h"', txt)
            if m:
                includes.append((p, i, m.group(1)))
                continue
            c = txt.find('//')
            for t in re.findall(r'\bscnXml::(\w+)', txt):
                if 0 <= c < txt.find('scnXml::'): break
                refs[t].append((p, i))

fcache = {}
def consumer(p, ln):
    if p not in fcache: fcache[p] = st.scopes(p)
    ctx = fcache[p].get(ln, ())
    owner = next((n for k, n in reversed(ctx) if k in ('class', 'struct')), '')
    fn = next((n for k, n in reversed(ctx) if k == 'function'), '')
    return owner or fn or '(file-scope)'

TARGET = lambda mod: 'model/xml/%s.cpp' % mod
w = csv.writer(open('census/xml_glue.csv', 'w', encoding='utf-8'))
w.writerow(['kind', 'module', 'type', 'target_xml_file', 'ref_count', 'sites',
            'consumers', 'review'])
unknown = []
DEAD_NOTE = {  # hand-curated overrides for refs the machine cannot classify
    'IPTDescription': ('?', 'DEAD: type removed from current XSD (last in scenario_30); '
                       'only referenced by dead file DescriptiveIPTInfection.*'),
}
for t in sorted(refs):
    mod = type_mod.get(t)
    kind = 'type'
    if mod is None and t.startswith('parse'):
        kind, mod = 'fn', 'scenario'  # codegen entry point, not a type
    note = ''
    if t in DEAD_NOTE:
        mod, note = DEAD_NOTE[t]
    if mod is None:
        unknown.append(t); mod = '?'
    sites = refs[t]
    cons = sorted({consumer(f, l) for f, l in sites})
    # review: glue already inside an XML-consuming member (setFromXML/loadFromXML/
    # ctor reading the subtree) vs scattered inline reads phase 2 must hoist
    inline = sum(1 for f, l in sites
                 if consumer(f, l) not in ('(file-scope)',) and not re.search(
                     r'FromXML|Xml|load|init', consumer(f, l), re.I))
    review = note or (('hoist: %d/%d refs outside *FromXML/init readers' % (inline, len(sites))) if inline else 'ok')
    w.writerow([kind, mod, t, TARGET(mod) if mod != '?' else '?', len(sites),
                ';'.join('%s:%d' % s for s in sites[:8]) + (';…' if len(sites) > 8 else ''),
                ';'.join(cons[:6]) + (';…' if len(cons) > 6 else ''), review])
for (f, l, hdr) in sorted(includes):
    w.writerow(['include', hdr, hdr + '.h', TARGET(hdr), 1, '%s:%d' % (f, l), '', ''])
n_t = sum(1 for t in refs)
print('types:', n_t, 'includes:', len(includes), 'unknown-type refs:', unknown,
      '-> census/xml_glue.csv')
