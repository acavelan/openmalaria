#!/bin/bash
# Benchmark harness (plan v3.3): pinned core 5, serialized one binary at a
# time, interleaved rounds, 2 warmup rounds then 5 measured, /usr/bin/time %U,
# median + MAD reported. All paths repo-relative; override the scenario with
# SCEN=/path/to.xml (KPI scenario = util/example/example_scenario.xml, pop 2000;
# scaling copies: util/bench-scen-10000.xml, or sed popSize for others).
#
# usage: util/bench-om.sh LOGFILE BIN1 [BIN2 ...]
set -u
LOG=$1; shift
BINS=("$@")
REPO="$(cd "$(dirname "$0")/.." && pwd)"
SCEN=${SCEN:-$REPO/util/example/example_scenario.xml}
RES=$REPO/test
XSD=$REPO/schema/scenario_50.xsd
CORE=${CORE:-5}
RUNS=${TMPDIR:-/tmp}/om-bench
ROUNDS=7   # first 2 = warmup, last 5 measured
for r in $(seq 1 $ROUNDS); do
  for b in "${BINS[@]}"; do
    d=$RUNS/$(basename "$(dirname "$b")")
    mkdir -p "$d"; ln -sf "$XSD" "$d/scenario_50.xsd"
    if ! ( cd "$d" && /usr/bin/time -f %U -o "$d/.time" taskset -c $CORE "$b" \
         --resource-path "$RES" --scenario "$SCEN" --output-format bin --output output.bin >/dev/null 2>"$d/.err" ); then
      echo "ABORT: $(basename "$b") failed: $(head -2 "$d/.err")"; exit 1
    fi
    echo "round=$r bin=$(basename "$(dirname "$b")") user_s=$(cat "$d/.time")" | tee -a "$LOG"
  done
done
python3 - "$LOG" "${BINS[@]}" <<'EOF'
import sys, statistics
log = sys.argv[1]; bins = sys.argv[2:]
runs = {}
for line in open(log):
    p = dict(kv.split("=", 1) for kv in line.split())
    if int(p["round"]) > 2:
        runs.setdefault(p["bin"], []).append(float(p["user_s"]))
for b in bins:
    name = b.split("/")[-2]
    v = sorted(runs.get(name, []))
    if len(v) < 5: print("%-10s INCOMPLETE %s" % (name, v)); continue
    med = statistics.median(v)
    mad = statistics.median([abs(x - med) for x in v])
    print("%-10s median=%8.2fs mad=%.3f runs=%s" % (name, med, mad, v))
EOF
