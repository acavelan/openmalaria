#!/usr/bin/env python3
"""statics.py — turn clang-query detailed-ast output into census/statics.raw.csv.

Usage: statics.py <clang-query-output> <out.csv>

Parses each `Binding for "v":` VarDecl dump line (absolute spelling location
only), dedups redecls by (file, line, name), prefers the 'definition' record,
then classifies scope (file-static / namespace-global / static-class-member /
function-static) and owner class/namespace via a brace-depth walk of the
spelling file.
"""
import re, sys, collections

VAR = re.compile(r"^VarDecl 0x\w+ (?:prev 0x\w+ )?<(/[^>:]+):(\d+):\d+,[^>]*> (?P<mid>.*)$")
TYPE = re.compile(r"(?P<types>'[^']*'(?:\s*:'[^']*')? )(?P<rest>.*)$")

def parse(path):
    seen = {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        prev = False
        for line in fh:
            if prev:
                m = VAR.match(line)
                if m:
                    t = TYPE.search(m.group('mid'))
                    if t:
                        f, ln = m.group(1), int(m.group(2))
                        name = m.group('mid')[:t.start()].split()[-1]
                        typ = t.group('types').strip()
                        if "':'" in typ:  # clang sugar form 'Alias':'Real' -> keep real
                            typ = typ.split("':'", 1)[1]
                        typ = typ.strip("'")
                        flags = t.group('rest').strip()
                        key = (f, ln, name)
                        if key not in seen or ('definition' in flags and 'definition' not in seen[key][2]):
                            seen[key] = (typ, name, flags)
            prev = line.startswith('Binding for "v":')
    return seen

DECL = re.compile(r'\b(namespace|class|struct|union)\s+(?:OM_\w+\s+)?([A-Za-z_]\w*)?\b')
FN = re.compile(r'([A-Za-z_~]\w*)\s*(?:<[^<>]*>)?\s*\([^;{}()]*\)\s*(?:const\s*)?(?:noexcept\s*)?(?:override\s*)?(?:final\s*)?$')
NOSCOPE = {'if', 'for', 'while', 'switch', 'catch', 'return', 'sizeof', 'throw', 'else', 'do', 'new', 'and', 'or', 'not', 'case'}

def scopes(path):
    """line -> tuple(scope kind, name) stack, brace-depth walk of a source file."""
    try:
        lines = open(path, encoding='utf-8', errors='replace').read().splitlines()
    except OSError:
        return {}
    out, stack, pending = {}, [], None
    in_comment = False
    for i, raw in enumerate(lines, 1):
        buf, j = [], 0
        while j < len(raw):
            if in_comment:
                k = raw.find('*/', j)
                if k < 0: j = len(raw)
                else: in_comment = False; j = k + 2
            else:
                k1, k2 = raw.find('//', j), raw.find('/*', j)
                if k1 >= 0 and (k2 < 0 or k1 < k2): buf.append(raw[j:k1]); j = len(raw)
                elif k2 >= 0: buf.append(raw[j:k2]); in_comment = True; j = k2 + 2
                else: buf.append(raw[j:]); j = len(raw)
        code = re.sub(r'"(?:\\.|[^"\\])*"', '""', ''.join(buf))
        code = re.sub(r"'(?:\\.|[^'\\])*'", "''", code)
        d = DECL.search(code)
        if d and d.group(2) is not None and not re.search(r'\b(class|struct)\b[^;{]*\b\w+\s*\(', code):
            pending = (d.group(1), d.group(2))
        elif d and d.group(1) == 'namespace':
            pending = ('namespace', d.group(2) or '(anonymous)')
        fm = FN.search(code.split('{')[0].strip())
        if fm and fm.group(1) not in NOSCOPE and not re.match(r'(extern|typedef|using|template|\breturn\b|\bif\b|\bfor\b|\bwhile\b|\bswitch\b)', code.strip()):
            # function definition opens at its '{' on this or a later line; only
            # trust it when the signature line doesn't end with ';' (declaration)
            if not code.strip().endswith(';'):
                pending = ('function', fm.group(1))
        for ch in code:
            if ch == '{':
                if pending is None:
                    # multiline function signature: look back for name(...) closing
                    # just before this brace (single-path regex, no nested star)
                    back = '\n'.join(lines[max(0, i - 12):i - 1])
                    bm = re.search(r'([A-Za-z_~]\w*)\s*(?:<[^<>]*>)?\([^;{()]*\)\s*(?:(?:const|noexcept|override|final)\s*)*$', back, re.S)
                    if bm and bm.group(1) not in NOSCOPE:
                        pending = ('function', bm.group(1))
                stack.append(pending or ('block', ''))
                pending = None
            elif ch == '}':
                if stack: stack.pop()
                pending = None
        if ';' in code: pending = None
        out[i] = tuple(stack)
    return out

def classify(f, ln, name, typ, flags, caches):
    rel = f.split('/fork/')[-1]
    if f not in caches: caches[f] = scopes(f)
    ctx = caches[f].get(ln, ())
    kinds = {k for k, _ in ctx}
    owner = next((n for k, n in reversed(ctx) if k in ('class', 'struct', 'union')), '')
    ns = next((n for k, n in reversed(ctx) if k == 'namespace'), '')
    if 'function' in kinds:
        scope = 'function-static'
    elif owner:
        scope = 'static-class-member'
    elif 'internal-linkage' in flags:
        scope = 'file-static'
    else:
        scope = 'namespace-global'
    definition = 'yes' if (not re.search(r'\bextern\b', flags) and
                           not (scope == 'static-class-member' and re.search(r'\bstatic\b', flags))) else 'no'
    return [rel, ln, name, typ, scope, owner, ns,
            'internal' if 'internal-linkage' in flags else
            'external' if 'external-linkage' in flags else 'other',
            definition,
            'yes' if any(w in flags for w in ('cinit', 'callinit', 'constinit', 'initlist')) else 'no',
            flags, '']

def main():
    seen = parse(sys.argv[1])
    caches = {}
    records = []  # (name, type, row, is_def)
    for (f, ln, name), (typ, _, flags) in sorted(seen.items()):
        if not f.split('/fork/')[-1].startswith('model/'):
            continue
        is_def = ('definition' in flags or
                  (not re.search(r'\bextern\b', flags) and not re.search(r'\bstatic\b', flags)))
        records.append((name, typ, classify(f, ln, name, typ, flags, caches), is_def))
    # merge redecls (header decl + cpp definition) sharing (name, type):
    groups = collections.defaultdict(list)
    for rec in records:
        groups[(rec[0], rec[1])].append(rec)
    rows = []
    for (name, typ), members in sorted(groups.items()):
        defs = [m for m in members if m[3]]
        if len(members) > 1 and len(defs) == 1:
            others = [m for m in members if not m[3]]
            cls = next((m for m in others if m[2][4] == 'static-class-member'), None)
            kept, rest = (cls, [m for m in members if m is not cls]) if cls else (defs[0], others)
            kept = list(kept[2])
            kept[8] = 'yes'
            kept[9] = 'yes' if any(m[2][9] == 'yes' for m in members) else 'no'
            kept[11] = ';'.join('{}:{}'.format(m[2][0], m[2][1]) for m in rest)
            rows.append(kept)
        else:
            rows.extend(m[2] for m in members)
    rows.sort()
    out = open(sys.argv[2], 'w', encoding='utf-8')
    out.write('file,line,name,type,scope,owner_class,namespace,linkage,definition,init,other_sites,clang_flags\n')
    for r in rows:
        out.write(','.join('"{}"'.format(str(c).replace('"', "'")) for c in r) + '\n')
    out.close()
    print('statics:', len(rows), '->', sys.argv[2])


if __name__ == '__main__':
    main()
