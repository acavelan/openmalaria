#!/bin/sh
# monids.sh — measure-id map + drift gate (plan §3g).
#
#   ./util/census/monids.sh                    regenerate census/mon_measure_ids.csv
#   ./util/census/monids.sh compare OLD NEW    FAIL if the two binaries' survey
#                                              measure-id streams (set AND
#                                              first-appearance order) differ
#
# Static truth: OutMeasures.h defs[] index = internal Measure id; outId == the
# 'measure' token in writeText (Monitoring.cpp) and in OMOUTBIN v2 measure
# blocks; nNMFTreatments outId=-1 is never output; isObsoleteMeasure() names
# keep their defs slot (ids never reused). ctsout order = map<string,…>
# alphabetical by registered output name (mon/Continuous.cpp).
#
# Default compare scenario: test/scenario10.xml (fast; 3 distinct ids).
# Full sensitivity: SCENARIO=util/example/example_scenario.xml (18 ids, ~50 s).
set -e
cd "$(dirname "$0")/../.."

if [ "$1" = compare ]; then
    OLD=$(readlink -f "$2"); NEW=$(readlink -f "$3")
    [ -x "$OLD" ] && [ -x "$NEW" ] || { echo "usage: $0 compare OLD NEW"; exit 2; }
    TMP=$(mktemp -d)
    # self-contained tmpdir wiring (mirrors util/parity.py): scenario + its XSD
    # alongside, resources via --resource-path (run.py does the same for ctest)
    SC=${SCENARIO:-test/scenario10.xml}
    cp "$SC" "$TMP/scenario.xml"
    XSD=$(sed -n 's/.*schemaLocation="[^"]* \([^" ]*\.xsd\)".*/\1/p' "$SC" | head -1)
    for cand in schema/$XSD test/$XSD build-p1/schema/$XSD; do
        [ -f "$cand" ] && cp "$cand" "$TMP/" && break
    done
    mkdir -p "$TMP/o" "$TMP/n"
    for x in "$TMP"/*.xsd; do cp "$x" "$TMP/o/"; cp "$x" "$TMP/n/"; done 2>/dev/null || true
    ( cd "$TMP/o" && "$OLD" --resource-path "$OLDPWD/test" --scenario "$TMP/scenario.xml" ) >/dev/null 2>/dev/null
    ( cd "$TMP/n" && "$NEW" --resource-path "$OLDPWD/test" --scenario "$TMP/scenario.xml" ) >/dev/null 2>/dev/null
    [ -s "$TMP/o/output.txt" ] && [ -s "$TMP/n/output.txt" ] || { echo "FAIL: run produced no output"; exit 3; }
    python3 - "$TMP/o/output.txt" "$TMP/n/output.txt" <<'PY'
import sys
def ids(p):
    seen, allids = [], set()
    for line in open(p):
        f = line.split('\t')
        if len(f) >= 4 and f[0].lstrip('-').isdigit() and f[2].lstrip('-').isdigit():
            i = f[2]
            allids.add(i)
            if i not in seen: seen.append(i)
    return allids, seen
ao, so = ids(sys.argv[1])
an, sn = ids(sys.argv[2])
bad = 0
if ao != an:
    print('FAIL measure-id set drift:', sorted(ao ^ an, key=int)); bad = 1
elif so != sn:
    print('FAIL measure-id first-appearance order drift at index',
          next(i for i, (x, y) in enumerate(zip(so, sn)) if x != y)); bad = 1
print('measure ids: %d distinct, %d first-seen order — %s' %
      (len(an), len(sn), 'PASS' if not bad else 'FAIL'))
sys.exit(bad)
PY
    rm -rf "$TMP"
    exit $?
fi

python3 - <<'PY'
import re, csv
src = open('model/mon/OutMeasures.h', encoding='utf-8').read()
body = src[src.index('defs[] = {'):src.index('inline bool isObsoleteMeasure')]
obsolete = set(re.findall(r'"(\w+)"', src[src.index('isObsoleteMeasure'):src.index('MeasureCount =')]))
dimnames = {'Age': 'Age', 'Cohort': 'Cohort', 'Species': 'Species',
            'Genotype': 'Genotype', 'Drug': 'Drug', 'None': 'None'}
rows = []
for i, m in enumerate(re.finditer(
        r'\{\s*"(\w+)"\s*,\s*(-?\d+)\s*,\s*(false|true)\s*,\s*([^}]+?)\s*\}', body)):
    name, out, dbl, dims = m.groups()
    dims = '+'.join(dimnames.get(d.strip(), d.strip())
                    for d in dims.replace('Dim::', '').split('|'))
    rows.append([i, out, name, 'double' if dbl == 'true' else 'int', dims,
                 'no', 'not output (outId=-1)' if out == '-1' else ''])
for name in sorted(obsolete):  # rejected by isObsoleteMeasure() at <enable>-parse
    rows.append(['none', '-1', name, '?', '?', 'yes',
                 'not in defs; isObsoleteMeasure rejects it in scenario <enable>'])
w = csv.writer(open('census/mon_measure_ids.csv', 'w', encoding='utf-8'))
w.writerow(['defs_index', 'outId', 'name', 'kind', 'dims', 'obsolete', 'note'])
for r in rows:
    w.writerow(r)
ndefs = sum(1 for r in rows if r[0] != 'none')
src_mc = re.search(r'MeasureCount = sizeof\(defs\)/sizeof\(defs\[0\]\)', src.replace(' ', ''))
print('defs rows:', ndefs, '(== MeasureCount via sizeof; verify by compile)',
      'obsolete rows:', len(rows) - ndefs,
      'outId<0:', sum(int(r[1]) < 0 for r in rows), '-> census/mon_measure_ids.csv')

# self-protection: the ordering contracts this map encodes must stay in the tree
import sys
mon = open('model/mon/init.cpp').read() + open('model/mon/Monitoring.h').read() \
    + open('model/mon/Monitoring.cpp').read()
must = [
    ('reportedMeasures sorted by outId (store/output order)',
     r'std::sort\(reportedMeasures.*lhs\.outId < rhs\.outId'),
    ('allCauseIMR special-cased out of reportedMeasures', r'reportIMR = om\.outId'),
    ('MeasureStore flat index nests nDrugs factor', r'nDrugs \*'),
]
for what, rx in must:
    ok = re.search(rx, mon, re.S) is not None
    print(('OK   ' if ok else 'FAIL ') + what)
    if not ok:
        sys.exit(1)
PY
