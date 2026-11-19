#!/usr/bin/env bash
# Machine-extract every static-storage-duration variable in model/ (file-static,
# namespace-scope global, static class member, function-local static).
#
# Method: clang-query AST match over every model/*.cpp TU in the build's
# compile_commands.json:
#     match varDecl(hasGlobalStorage(), unless(isExpansionInSystemHeader()))
# Scope kind + owner class/namespace are then resolved statically from the
# spelling file with a brace-depth scope walk (statics.py).
#
# Prereq (gate build + compile DB):
#   cmake -B build-p1 -DCMAKE_BUILD_TYPE=Release \
#     -DCMAKE_CXX_FLAGS="-ffp-contract=off" -DOM_BOXTEST_BINARY=ON \
#     -DCMAKE_EXPORT_COMPILE_COMMANDS=ON . && cmake --build build-p1 -j
#
# Usage:  util/census/statics.sh [build-dir]      # default build-p1
# Writes census/statics.raw.csv (machine stage). census/statics.csv is the
# reviewed/tagged census; diff the two when re-running to spot new/removed statics.
set -euo pipefail
BUILD="${1:-build-p1}"
cd "$(dirname "$0")/../.."
[ -f "$BUILD/compile_commands.json" ] || { echo "run cmake with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON on $BUILD"; exit 1; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

cat > "$TMP/statics.que" <<'EOF'
set bind-root false
set output diag
enable output detailed-ast
match varDecl(hasGlobalStorage(), unless(isExpansionInSystemHeader())).bind("v")
EOF

python3 - "$BUILD" > "$TMP/tus" <<'PY'
import json, sys
db = json.load(open(sys.argv[1] + '/compile_commands.json'))
for e in db:
    if '/model/' in e['file'] and e['file'].endswith('.cpp'):
        print(e['file'])
PY

xargs -a "$TMP/tus" -P "$(nproc)" -I{} \
    env B="$BUILD" Q="$TMP/statics.que" T="$TMP" sh -c \
    'clang-query -p "$B" -f "$Q" "$1" > "$T/$(echo "$1" | tr / _).ast" 2>/dev/null || echo "WARN $1" >&2' _ {}

cat "$TMP"/*.ast > "$TMP/all.ast"
python3 "$(dirname "$0")/statics.py" "$TMP/all.ast" census/statics.raw.csv

# Fallback regex stage: static data declarations in model/ sources that the AST
# pass cannot see (dead/uncompiled files, #ifndef NDEBUG-only members under a
# Release compile). Appended with scope/note=not-in-build so the census stays
# complete and re-runnable.
python3 - census/statics.raw.csv <<'PY'
import re, csv, sys
have = set()
rows = list(csv.reader(open(sys.argv[1], encoding='utf-8')))
for r in rows[1:]:
    have.add((r[0], int(r[1])))
out = open(sys.argv[1], 'a', encoding='utf-8')
added = 0
import subprocess
scan = subprocess.run(['grep', '-rnE',
    r'^[[:space:]]+(static|static constexpr|constexpr static)[[:space:]]+[^;(]*[[:space:]*][A-Za-z_][A-Za-z0-9_]*[[:space:]]*(=[^;]*)?;',
    '--include=*.h', '--include=*.cpp', 'model'], capture_output=True, text=True).stdout
for line in scan.splitlines():
    f, ln, txt = line.split(':', 2)
    ln = int(ln)
    if (f, ln) in have:
        continue
    txt = txt.strip()
    if txt.startswith(('//', '/*', '*')):
        continue
    after = txt.split(None, 1)[1]
    if '(' in after.split('=')[0].split(';')[0] or after.startswith('assert'):
        continue  # static member function
    name = re.findall(r'([A-Za-z_]\w*)\s*(=|;)', txt)[0][0]
    toks = after.split('=')[0].strip('; ').split()
    typ = ' '.join(toks[:-1])
    row = [f, ln, name, typ, 'not-in-build', '', '', 'external', 'no',
           'yes' if '=' in txt else 'no', '',
           'regex fallback: not visible in gate build (uncompiled or NDEBUG/ifdef)']
    out.write(','.join('"{}"'.format(c) for c in row) + '\n')
    added += 1
out.close()
print('fallback rows added:', added)
PY

