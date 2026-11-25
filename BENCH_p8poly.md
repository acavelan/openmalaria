# Bench: p8-poly vs p7g baseline (2026-11-20, final)

Harness: taskset -c 5, `/usr/bin/time -f %U` user-CPU, scenario
/tmp/opencode/om-bench/example_scenario.xml (popSize 2000), resource-path
/home/razlock/git/fork/test, 1 discarded warmup per arm, **10 interleaved
reps** (paired old/new per rep). Flags identical (CMakeCache:
`-O3 -DNDEBUG -ffp-contract=off`). Binaries: p7g `7120add4…`,
p8-poly@`c9732f64` `df1d5b66…`.

Window: sibling agents running other gates (load ~7 on a 24-core box; no HT
sibling on c5 — `siblings: 24`). That raises BOTH arms ~1.5% vs the historical
quiet-window baseline (p7g here 12.21 vs session baseline 12.15); the
same-window paired comparison is the valid number.

| arm | raw reps (user s) | median | MAD |
|---|---|---|---|
| p7g    | 12.51 12.59 12.22 12.19 12.14 12.17 12.16 12.22 12.23 12.19 | **12.21** | 0.03 |
| p8-poly | 12.66 12.50 12.35 12.36 12.36 12.34 12.33 12.40 12.38 12.34 | **12.36** | 0.02 |

**Δ = +0.15 s = +1.3%** (9 of 10 pairs positive, Δ range −0.09…+0.22) — inside
the ≤ +1.5% expectation for virtual dispatch over final-call dispatch.

## Retraction
An earlier draft of this file claimed +27% (15.44 vs 12.15). That was a
measurement artifact: the arms ran in *different* windows while parallel
parity gates (-j12) hammered all cores; pinned user-CPU absorbs the memory
bandwidth pressure. The "static kind dispatch = 15.23 s" attribution
experiment was measured inside the same loaded window and is likewise void —
no LTO-inlining cliff exists. p8-poly costs ~1.3%.
