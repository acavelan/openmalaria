#!/bin/bash
# Benchmark harness (plan v3.3): pinned core 5, serialized one binary at a
# time, interleaved rounds with alternating binary order, 2 warmup rounds then
# 5 measured. User CPU is the
# primary metric (median + MAD); system/elapsed/RSS/page-fault medians are
# also reported. All paths repo-relative; override the scenario with
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
{
  printf 'scenario=%s\ncore=%s\nglibc_tunables=%s\n' "$SCEN" "$CORE" "${GLIBC_TUNABLES-}"
  for b in "${BINS[@]}"; do
    printf 'bin=%s path=%s sha256=%s\n' \
      "$(basename "$(dirname "$b")")" "$b" "$(sha256sum "$b" | cut -d' ' -f1)"
  done
} > "$LOG"
for r in $(seq 1 $ROUNDS); do
  for ((i = 0; i < ${#BINS[@]}; ++i)); do
    j=$i
    ((r % 2 == 0)) && j=$((${#BINS[@]} - 1 - i))
    b=${BINS[$j]}
    d=$RUNS/$(basename "$(dirname "$b")")
    mkdir -p "$d"; ln -sf "$XSD" "$d/scenario_50.xsd"
    if ! ( cd "$d" && /usr/bin/time -f '%U %S %e %M %F %R' -o "$d/.time" taskset -c "$CORE" "$b" \
         --resource-path "$RES" --scenario "$SCEN" --output-format bin --output output.bin >/dev/null 2>"$d/.err" ); then
      echo "ABORT: $(basename "$b") failed: $(head -2 "$d/.err")"; exit 1
    fi
    read -r user system elapsed rss major minor < "$d/.time"
    echo "round=$r bin=$(basename "$(dirname "$b")") user_s=$user system_s=$system elapsed_s=$elapsed max_rss_kb=$rss major_faults=$major minor_faults=$minor" | tee -a "$LOG"
  done
done
python3 - "$LOG" "${BINS[@]}" <<'EOF'
import sys, statistics
log = sys.argv[1]; bins = sys.argv[2:]
runs = {}
for line in open(log):
    if not line.startswith("round="):
        continue
    p = dict(kv.split("=", 1) for kv in line.split())
    if int(p["round"]) > 2:
        runs.setdefault(p["bin"], []).append({k: float(p[k]) for k in
            ("user_s", "system_s", "elapsed_s", "max_rss_kb", "major_faults", "minor_faults")})
for b in bins:
    name = b.split("/")[-2]
    rows = runs.get(name, [])
    v = sorted(row["user_s"] for row in rows)
    if len(v) < 5: print("%-10s INCOMPLETE %s" % (name, v)); continue
    med = statistics.median(v)
    mad = statistics.median([abs(x - med) for x in v])
    print("%-10s user_median=%8.2fs user_mad=%.3fs system_median=%.2fs elapsed_median=%.2fs max_rss_median=%.0fKB major_faults_median=%.1f minor_faults_median=%.1f user_runs=%s" %
          (name, med, mad,
           *(statistics.median(row[k] for row in rows) for k in
             ("system_s", "elapsed_s", "max_rss_kb", "major_faults", "minor_faults")), v))
EOF
