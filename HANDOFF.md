# Phase 7c HANDOFF — comparison memo, ISA pricing, SoA-sweep vectorization, LOC sweep (branch `p7c-vectorize` off tag `p7b-green` = 0c895755)

Crash-recovery file. Update EVERY green commit. Work order: 7c-A memo (done),
7c-B ISA experiment (done), 7c-C vectorize the SoA sweeps (the user's explicit ask),
7c-D LOC sweep (dead code/wrappers/duplication/xml/ compaction), 7c-E final gate +
report. Full p7b handoff: `git show p7b-green:HANDOFF.md`.

## Environment
- Gate build: `build-p7c/` = Release `-DCMAKE_CXX_FLAGS=-ffp-contract=off`
  `-DOM_BOXTEST_BINARY=ON`. `cmake --build build-p7c -j16`. unittest from
  `build-p7c/unittest/` → `./unittest` (".OK!", 94 cases). ctest: `cd build-p7c && ctest -j16` (52).
- Perf build (lever baseline): `build-p7c-a2/` = exact p7b-a2 recipe (Release
  `-ffp-contract=off -flto=auto -fprofile-use -fprofile-correction -Wno-missing-profile`,
  tests OFF). Verified this session: bench 11.98 (raw 11.96–12.00) — yardstick reproduces.
- `build-p7c-a2n/` = a2 + `-march=native -mtune=native` — 7c-B record, KEEP, do not adopt.
- T1 reference: `/tmp/opencode/om-ref/A0/openMalaria` + recipe.json (hash-checked by runner).
  `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new <bin>
  --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a -j16` (53, ~17 s);
  `--tiers b --seeds 1-10 -j16` (10, ~35 s). Never run benches concurrently with parity/ctest.
- BENCH: **MUST `cd /tmp/opencode/om-bench` first** (XSD symlink lives there; bench.sh does NOT
  cd — running it elsewhere = instant validation error = "0.00" lines, trap hit once this
  session). `./bench.sh <bin> /home/razlock/git/fork/util/example/example_scenario.xml` =
  1 discard + 5× `/usr/bin/time -f %U`, taskset -c 5. Noise floor max(1.0%, 2×MAD).
- perf events work at paranoid=2 via `cpu_core/…​/u` pmu prefix (`perf stat -e
  instructions,cycles,branches,branch-misses` — hybrid atom events show "<not counted>",
  harmless). `perf record -F 999` for software sampling.

## Current yardsticks (this branch, HEAD green at p7b-green = 0c895755)
| binary | CPU-s | vs frozen baseline-A2 11.34 | note |
|---|---|---|---|
| gate A0-class (build-p7c) | 12.47 | — | p7b stamp |
| **build-p7c-a2 (lever baseline)** | **11.98 MAD .010** | **0.947** | re-measured 7c-B |
| build-p7c-a2n (native) | 12.12 MAD .010 | 0.936 | +1.17 % slower, T1 63/63 clean — NOT adopted |

## Status
- [x] branch `p7c-vectorize` off `p7b-green`; builds build-p7c / -a2 / -a2n configured+green
- [x] 7c-A: `census/antigravity_comparison.md` — staircase table (bars NOT comparable:
      wall-mix+startup vs CPU-s, ISA asymmetry, 1e-6 text vs bit-exact parity), per-lever
      adoption table, the four "only works at 1e-6" flags (hand-LU solver, tbv multiply
      reorder, hand-rolled CDFs, value-changing bug fixes), ISA pricing (from 7c-B).
- [x] 7c-B ISA: a2n bit-clean (63/63, 0 FMA emissions) but **−1.17 % SLOWER** (IPC 3.53→3.45,
      −1.05 % instr / +1.28 % cycles) ⇒ NOT proposed; plan's -march ban costs nothing here.
      PERF row 7c-B. This prices the biggest bundled ingredient of their 1.20→1.46 step.
- [x] **AUDIT (7c-B, source-verified): the p7b "A2/P-go" yardstick binary never consumed
      profile data.** `-fprofile-use` resolves `.gcda` beside each object (probe:
      `-Wmissing-profile` names it); none existed in build-p7b-a2; om-pgo holds only
      `#…​build-p7b-pgogen#…​`-encoded names; `-Wno-missing-profile` silenced everything.
      ⇒ "PGO collapsed 1.087→1.003" (row 7b-A) was never actually tested. Fix mechanism
      probe-verified: rename profile files to the new build-dir encoding + `-fprofile-dir=DIR`.
      PGO on the current tree gets re-priced properly in 7c-E (plan: build-p7c-pgogen corpus
      on HEAD → rename → build-p7c-a2p with -fprofile-dir → T1 tier-b + bench).
- [x] 7c-C: SoA sweep vectorization — **C3 KEPT** (immunitySurvivalFactor host-tick hoist,
      −1.04 % = 11.85, commit 5125841b). C1 pyrogen-hoist, C2 age-sweep hoist, C6 vectorUpdate
      two-pass staging all measured sub-floor ⇒ reverted+logged (PERF 7c-C rows). Vec-info
      evidence: hot TUs already run 126 auto-vectorized loops, default cost model reports zero
      missed; rest blocked by control flow / clobbering libm+recordStat = phase-8/9 territory.
      Pre-phase-8 bit-exact vectorization is exhausted.
- [ ] 7c-D: LOC sweep (≥-50 LOC deletions = own commit, gate + tier-a each)
- [ ] 7c-E: final gate (T1 a+b 63/63, ctest, unittest, cloc table, benches incl. a2p,
      perf top-5, PERF 7c-F scoreboard, report, HANDOFF close)

## 7c-C work notes (vectorize the sweeps — the user's explicit ask)
- Pattern target: pass-1 stage elementwise quantities into flat run-owned scratch
  (HostScratch [host] and [host×genotype] f64 slabs, capacity reused never shrunk),
  pass-2 SCALAR ordered reduce (no reassociation — phase 8 owns reductions going SIMD).
- Elementwise-only loops: contiguous scratch, __restrict/local buffers, no early return
  in the map lane, no / or libm in the hot lane. `-fopt-info-vec-optimized` must show
  the intended loops; per-loop `// vec-blocked: needs phase-8 break` where not fixable.
- Hot chain to grep in the vec-missed report: Host::update stages,
  calculateEIR_vector/vectorUpdate, probTransmissionToMosquito_*, updateKappa,
  immunitySurvivalFactor/updateImmuneStatus, Descriptive updateDensity, iv_update.
- Per-structure commits, measure on build-p7c-a2 rebuild (same recipe = 11.98 baseline),
  keep only >floor wins, revert+log the rest.

## Gates log
- (setup) 3 builds configured; gate build: ctest 52/52 ✓, unittest .OK! ✓; a2 bench 11.98
  (raw 11.96–12.00) = yardstick reproduced.
- (7c-B) a2n: objdump 0 FMA; parity tier-a 53/53; tier-b s1–10 10/10; bench 12.12 MAD .010
  (raw 12.10–12.13); perf stat −1.05 % instr, +1.28 % cycles. Verdict: not adopted, logged.
- (7c-C3) commit 5125841b: gate ctest 52/52, unittest OK, tier-a 53/53; a2 11.85 (n=20).
- (7c-C1/C2/C6) each measured sub-floor on a2 n=20 ⇒ reverted, logged in PERF 7c-C. Final
  tree = HEAD 5125841b; a2 same-session re-measure 11.835 (n=20, raw 11.81–11.88). Vec-info
  artifacts: /tmp/opencode/vec_all.txt.

## Traps (inherited, all live, + new)
- NEW: benches only from `/tmp/opencode/om-bench` cwd; never overlap bench with any other load.
  Scenario = absolute `/home/razlock/git/fork/util/example/example_scenario.xml` (bare/wrong
  name ⇒ all-0.00 lines, same symptom as wrong cwd). `unittest` must run from
  `build-p7c/unittest/` (fixture CSVs sit next to the binary; elsewhere = 16 file-not-found).
- NEW: PGO cross-build-dir reuse requires the rename+`-fprofile-dir` dance (see AUDIT above);
  "no warning" ≠ "profile used" — verify with `-Wmissing-profile` probe before believing any
  A2 number.
- set_row REBIND on any column shape change; m_rng ctor-draw rule (draw point frozen);
  TwoRunsSuite+RngGoldens tripwires.
- -ffp-contract=off always in gated binaries; NEVER -ffast-math/-march=native in GATE builds
  (a2n experiment is USER-class, T1-proven, not adopted); NO FP reassociation (phase 8).
  Vectorization = elementwise maps only; every += accumulator stays scalar.
- Checkpoints may break freely (re-bless from A0 ref in same commit + round-trip smoke).
- Untracked PROMPT.md/email = user files, leave alone.
- ≤5 fix iterations per divergence, then revert+log.
