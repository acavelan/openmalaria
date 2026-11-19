---
description: Contract for the OpenMalaria model rewrite (plan v3.3). Load for ANY rewrite/refactor/perf task in this repo. KPIs: code LOC (primary), performance, readability. Defines the correctness ladder (T1 bit-exact / T2 equation fixtures / T3 calibrated statistics), RNG preservation rules, target layout, phase map, gates and the PERF.md ledger. Single-threaded, single-seed, all model variants kept.
---

# OpenMalaria rewrite contract (v3.3)

Canonical plan (read it for phase briefs + pinned specs):
`/home/razlock/.opencode/plan/openmalaria-rewrite.md` — this skill is the law;
the plan is the detail. Repo: `/home/razlock/git/fork`, baseline `cfc3c8e0`.

## KPIs — every phase reports all three in PERF.md

| KPI | Baseline | Target |
|---|---|---|
| **LOC** — `cloc --quiet model --include-lang="C++,C/C++ Header" \| grep '^SUM'`, col 5 (code lines only — comments/blanks free) | 18,813 | **Σ ≤ 14.0k incl. `model/xml/` (hard cap; xml/ ≤ 2.6k)**. Breach = user sign-off event, never silent drift. `unittest/` LOC is a second column — moving model logic into tests is forbidden. Comments cost nothing toward the KPI: annotate `// parity:` generously |
| **Performance** — per-scenario **CPU-seconds** (`/usr/bin/time -f %U`), single benchmark = `util/example/example_scenario.xml` (popSize 2000, iseed 1, its built-in values — do not edit), `taskset`-pinned, governor/turbo off, median-of-5, vs the A2 (PGO) build of the *current* baseline; warmup/main split; never ctest wall | phase 0 stamps A0/A2 | flags alone ≥1.10 (A2/A0); rewrite **1.15–1.35× pre-break, 1.5–2.0× post-break vs A2**; end state ≈1.8× vs plain `-O3` |
| **Readability** | — | one free function per named model quantity, unit-testable in isolation; doc comment = symbolic formula + paper ref + `// parity: File.cpp:NN` anchor; zero constraint checks in hot math |

An optimization that does not beat the noise floor (max(1.0%, 2×MAD)) is
**REVERTED and logged**, never kept on faith. Reverts are ledger rows.

## Correctness ladder — every edit passes its gate before the next starts

| Tier | Gate | Scope |
|---|---|---|
| **T1 bit-exact** | two tiers: (a) all 51 box scenarios (small pops — variant coverage only), old vs new binary, `output.bin` **bit-identical**, single continuous run: `python3 util/compareBin.py OLD.bin NEW.bin` (0 = identical); (b) **representative run: `util/example/example_scenario.xml` @ popSize 2000, seeds 1..10 via `--seed`** — box-test populations are too small to be representative, so no phase is green without tier (b). Reference-binary recipe (compiler, flags, lib versions, commit) hashed — hash mismatch = stop, do not debug "phantom" divergence. Gate binaries: `-O3 -DNDEBUG -ffp-contract=off`, no LTO/PGO unless they pass a T1 bit-check | every phase 2–7 commit |
| **T2 equations** | reference-value suites exact (Molineaux/Penny trajectories, PK/PD, DecayFunction, AgeGroupInterpolation, decision-tree fractions, CalcNvO-linear, percentiles, CFR, `RngGoldenSuite` state-words): `ctest -R unittest` | always |
| **T3 statistics** | per-scope tolerance table frozen in phase 0, **calibrated on `example_scenario.xml` @2000 × 10 seeds** (old binary): `python3 util/statCompare.py OLD_DIR NEW_DIR` + sampler battery (KS@50k, χ² tail-pooled, regime boundaries). CDF/sampler replacements allowed **only** where bit-parity grids pass; else keep GSL | **post-break only** |

Harness self-check: `python3 unittest/HarnessScriptsTest.py`.

**The ONE sanctioned break = Phase 8** (not 6): FP reassociation, `inv_*`
reciprocals, slot traversal, SIMD samplers, CalcNvO solver flip, value-changing
bug fixes. One deliberate commit re-blessing `test/expected/` AND `unittest/`
fixtures. Everything before it is T1-provable.

## RNG rules (violations diverge chaotically — there is no partial failure)

- `model/util/xoshiro.h`, `contrib/chacha.h` semantics: copy verbatim (sole
  allowed change: the ChaCha `was_zero` discard-guard — a UB fix).
- Pre-break: **never** change order of RNG draws or FP operation order. No FMA
  contraction (gate builds pin `-ffp-contract=off`), no `-ffast-math`,
  no `-march=native`/`-mtune=native` in gate binaries (FMA = divergence), no
  loop recombination, no reciprocal precompute (`x/b ≠ x*(1/b)` — post-break
  only), no sampler algorithm swaps.
- Master seeding order (`main.cpp:152-194`) frozen incl. the
  validation-before-`seed()` ordering; percentile table uses `MasterRng(0,0)` —
  load-bearing, do not "clean up".
- Frozen draw-order contracts **a–g** (plan §3): bernoulli-last deploy filters,
  reproducibility `stable_sort`s, history-ordered FP products, decision-tree
  child order + one `uniform_01` per Random node, per-day Poisson/bernoulli
  counts, validation-before-seed, mon measure-id = output column order.
  Annotate every draw site `// parity:`; any PR touching an annotated site
  re-justifies the order in its commit message.

## Target layout (config/state enum+switch idiom; no class hierarchy with one
live implementation; every variant kept — as equations over shared machinery)

    model/xml/          scnXml→plain-struct only; only code allowed to include schema; XSD codegen KEPT
    model/core/         RunContext: SimClock (in_update asserts), options u64+PAIRS table,
                        Parameters::try_get, registries, monitoring sink, master RNG — no statics
    model/sim/          tick loop
    model/humans/       SoA columns (pinned list, plan §7/§9); Host::update = 7 free functions
    model/infections/   slots ≤21/host behind ONE SlotAllocator (the clamp lives there)
    model/clinical/     compiled decision tree + flat node-array interpreter; 3-config→2-state collapse
    model/transmission/ dense per-host arrays; PerHost class gone; SimpleMPD = optional member + 4 hooks
    model/interventions/date-cursor + SIM_FUTURE; dense components by type-rank; HostView accessors; zero serialization
    model/mon/ checkpoint/ util/  mon: MeasureStore, enum ids = output column order; util: RNG core verbatim + inline GSL ports

## Phase map (one agent per phase; gate green before the next; git: branch per
phase off tag `p(N−1)-green`; ≤5 fix iterations on a red gate, then revert + escalate)

0 harness/benchmark/PGO baseline · 1 census + `// parity:` annotations → **user sign-off** ·
2 `xml/` · 3 RunContext (+ two-runs-in-one-process statics proof) · 4a decision-tree interpreter ·
4b remaining hierarchies (interventions untouched — phase 7 only) · 5 cycle breaks ·
6 SoA sidecar (layout-only, still T1) · 7 interventions rewrite · **8 the break** · 9 perf levers.

## Forbidden

- Editing `test/expected/` or `unittest/` fixtures outside an authorized
  re-bless commit (phase-0 clause, checkpoint-conversion re-bless blessed from
  the A0 reference, or the phase-8 break commit).
- Deleting/modelling-away any model variant, scenario, or suite (user veto).
- **Threads, multi-seed runs, in-process batch — the binary is one core, one
  seed, one run per invocation** (user veto). Parallelism lives outside the binary.
- Checkpoint reproducibility expectations: checkpoint = post-warmup snapshot,
  interventions never checkpointed, restore = event replay ignoring coverage
  filters (plan §6); no phase may serialize per-host intervention state.
- **Carrying a "checkpoint bytes frozen" constraint into a phase brief (user
  veto, reiterated: checkpoints may break freely — re-bless the affected box
  goldens from the A0 reference in the same commit + keep one round-trip smoke
  test; T1/T2/T3 unchanged).**
- New dependencies without a PERF.md row justifying them.
- Touching a `friend`ed class without updating its suites in the same commit.
- Trusting round-N mining claims over source: verify flags against source
  (two "hot GSL" and "golden-free solver" claims were already found false at source).

## Skills to load per task

Profiling/benchmarks → `cpp-performance`, `linux-perf`. Optimization pattern
hunting (phases 6–9) → `performance-patterns`. Any logic/fixture work →
`test-driven-development` (write the failing T2 fixture first). **Valgrind rule
(user): valgrind is extremely slow here — run memcheck/callgrind only on tiny
popSize scenarios (e.g. a small box scenario), never on the pop-2000
benchmark.**

**Ponytail (user requirement, 2026-11-19): always load and follow the `ponytail`
skill** for every coding task in this rewrite — lazy-first ladder (YAGNI →
reuse in-repo → stdlib → platform → installed dep → one line → minimum code).
It never overrides the correctness ladder, T1 gates, or frozen contracts — it
only constrains *how much* code a compliant solution may write. Comments are
LOC-free (cloc code column only), so `// parity:` notes stay cheap; code is
not.

