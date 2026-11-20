# Phase 7b HANDOFF — A2 yardstick rebuild + T1-safe perf levers (branch `p7b-levers` off tag `p7-green` = 324958dc)

Crash-recovery file. Update EVERY green commit. Pull-forward clause fired:
p7-green A0-class bench 12.47 s vs the FROZEN phase-0 A2 yardstick 11.34 s, but that
A2 was built from the BASELINE commit (cfc3c8e0) — meaningless against the rewritten
tree. Job A: rebuild the A2 (LTO+PGO) yardstick ON THE CURRENT TREE, prove T1-clean,
stamp the true structural ratio. Job B: T1-safe levers (arena, rng-column compaction,
forcedEIR early-out, perf-record hoists), one per commit, measured vs the NEW A2,
reverted+logged if ≤ noise floor (max(1.0%, 2×MAD)). Full p7 handoff: `git show p7-green:HANDOFF.md`.

## Environment (same machine, /tmp/opencode intact)
- Gate build: `build-p7b/` = Release `-DCMAKE_CXX_FLAGS=-ffp-contract=off`
  `-DOM_BOXTEST_BINARY=ON` (copy of build-p7 config). `cmake --build build-p7b -j16`.
  unittest MUST run from `build-p7b/unittest/` → `./unittest` (94). ctest: `cd build-p7b && ctest -j16` (52).
- A1 build: `build-p7b-a1/` = Release `-ffp-contract=off -flto=auto`, tests OFF.
- PGO pass-1: `build-p7b-pgogen/` = Release `-ffp-contract=off -flto=auto
  -fprofile-generate=/tmp/opencode/om-pgo` (tests ON). Corpus: ctest -j24 + unittest +
  example@2000 hot runs seeds 1-4. Profiles are shards in that dir; `-fprofile-use` merges.
- A2 build: `build-p7b-a2/` = Release `-ffp-contract=off -flto=auto
  -fprofile-use -fprofile-correction -Wno-missing-profile`, tests OFF. Stage-B levers are
  rebuilt INTO build-p7b-a2 (same profiles, same flags — only code differs).
- T1 reference A0 = `/tmp/opencode/om-ref/A0/openMalaria` (cfc3c8e0, recipe.json beside it;
  verified 2026-11-20 p7 session). parity.py: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria
  --new <bin> --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a -j16` (53) /
  `--tiers b --seeds 1-1` etc.
- BENCH: run from `/tmp/opencode/om-bench` (XSD symlinks live there; benches elsewhere fail).
  `/tmp/opencode/om-bench/bench.sh <bin> /home/razlock/git/fork/util/example/example_scenario.xml`
  = /usr/bin/time -f %U, taskset -c 5, 1 discard + 5 runs. Compute median + MAD yourself;
  noise floor = max(1.0%, 2×MAD); discard >0.5 s-above-cluster outlier but report raw 5.
- perf events BLOCKED (paranoid=2) — only `-F 999` software sampling via perf record.
- A0-class references: current tree 12.47 median; A0 baseline 12.33. Old frozen A2 11.34
  (baseline commit — superseded by row 7b-A once stamped).

## Status
- [x] branch created off p7-green
- [ ] A: A1 build + LTO T1 bit-check (tier-b seed 1)
- [ ] A: PGO corpus + A2 build + T1 bit-check (tier-b seed 1) + bench + PERF row 7b-A
- [ ] B1 arena; B2 rng compaction; B3 forcedEIR; B4 perf hoists (each: own commit + PERF row)
- [ ] FINAL: T1 a+b 63/63, ctest 52/52, unittest 94/94, cloc, census touch, final report

## Traps (inherited, all live)
- set_row REBIND on any column shape change (B2 tripwire; death-heavy ctest catches, tier-a doesn't).
- m_rng ctor-draw-point rule: LocalRng construction consumes master draws — B2 may change the
  LocalRng TYPE but the draw must stay at the same point; RngGoldens + TwoRunsSuite are the tripwires.
- -ffp-contract=off always in gated binaries; NEVER -ffast-math/-march=native; NO FP reassociation
  (phase 8). Every 7b lever must be value-identical; T1 is the proof.
- Checkpoints may break freely (re-bless from A0 reference in same commit + round-trip smoke).
- Untracked PROMPT.md/email = user files, leave alone.

## Gates log
- (stage setup) branch p7b-levers off p7-green 324958dc; env inventory verified
  (recipe check passes at parity startup; om-bench + om-ref intact).
