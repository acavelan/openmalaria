# OpenMalaria rewrite — full status, change report, and recovery guide

Branch: `rewrite` (pushed to this repo's `origin` fork).
Written: 2026-11-24. This file is the single entry point to understand or
rebuild everything from another machine. `PERF.md` and `HANDOFF.md` hold the
raw ledgers; this is the plain-language version.

---

## 1. What was done, in one paragraph

The OpenMalaria model code (`model/`, ~19k lines) was rewritten in place, in
small gated steps, with one goal: **same numbers out, far fewer lines in** —
plus a side quest to make it faster than both the original and a competing
optimisation branch. Result: **16,317 model code LOC (−13.3%)**, **bit-exact
output on every test scenario and 10 seeds**, and **1.146× faster than the
original at pop 2000, 1.214× at pop 10,000** — within ~2% of a branch that
spent its budget purely on speed. Every single step passed the same four
gates (below), so at no commit was the model ever wrong.

## 2. The rules that never changed (needed to trust the result)

- **T1 bit-exactness is the gate.** Output must be byte-identical to the
  reference binary A0. RNG draw points and draw order frozen; formula bodies
  copied verbatim; error-message strings frozen; compiled with
  `-ffp-contract=off` (no FMA contraction, ever).
- **Reference A0** = original code at commit `cfc3c8e0`, built with
  `-O3 -DNDEBUG -ffp-contract=off`, gcc 16.2.1, GSL 2.8 (see
  `util/a0-recipe.json`). `util/parity.py --recipe` hard-fails if the
  toolchain drifts — bit-exactness is only meaningful on the same recipe.
- **Four gates per change:** build with 0 warnings/0 errors; `ctest` 64/64;
  `unittest` prints `.OK!`; T1 parity `63/63 identical` (all `test/scenario*.xml`
  + example scenario × seeds 1–10).
- **LOC KPI:** `cloc --quiet model --include-lang="C++,C/C++ Header" | grep ^SUM`
  column 5. Baseline 18,813 → HEAD 16,317 (both re-measured, not quoted).
- All model variants kept (Descriptive/Dummy/Empirical/Molineaux/Penny
  infections; Falciparum/Common/Vivax within-host; all health-system kinds).
- Perf candidates must be measured in the **shipping build form** (LTO+PGO
  recipe), not just plain `-O3` — learned the hard way (E6 below).

## 3. Branches on the remote

| branch | what it is |
|---|---|
| `rewrite` | everything, through phase 8 + perf rounds 1–6. **This is the product.** |
| `eval/bs-all` | baseline + the 7 upstream DRAFT perf PRs (476–482), for ladder comparisons; `PERF.md` round 3 verified the PR author's claims on this branch |
| `eval/e6-pyrogen` | rejected experiment E6 (kept for reference, see §6) |
| `eval/e9-humans` | rejected experiment E9 (kept for reference, see §6) |
| `eval/instr-counters` | diagnostic-only counters; produced the measurement behind E8. Never merge |
| `bsperf/*` | the PR author's original branch-PR heads (mirror of remote `bscode`) |

Tag `rewrite-status-2026-11-24` marks the pushed HEAD.

## 4. The rewrite, phase by phase (LOC is model code LOC)

| milestone | LOC | Δ | what happened |
|---|---|---|---|
| baseline `cfc3c8e0` | 18,813 | — | original code |
| phase 0–7g | 18,248 | −565 | modernisation, small-file folds, `Clock.h` created, `Global.h` deleted, dead code removed (grep-proved), trivial accessors killed, folders flattened, `Infection` redesigned (7f) |
| p8-structure | 18,257 | +9 | `model/core/` killed; `RunExtensions` dissolved into `RunContext` (+9 = honest un-cycling tax) |
| p8-checkpoint | 17,632 | −625 | old checkpointing deleted wholesale, redone xml-style as self-contained `model/checkpoint/` module; 12 `Resume_*` tests prove resume == uninterrupted run |
| p8-sweep | 17,545 | −87 | ponytail audit sweep: dead data, single-caller wrappers, twin loops |
| p8-files | 17,192 | −353 | 16 merge groups, 140→115 files, `#pragma once` everywhere, checkpoint decls out of model headers |
| p8-poly | 17,078 | −114 | `Infection` and `WHInterface` Kind-switches → real class hierarchies, ctor layering keeps RNG order exact |
| p8-collapse | 16,759 | −319 | STL collapses (find_if/erase_if/max_element/accumulate), a real bug fixed (bogus assert in `transientEIRknown`), clang-format reflow |
| p8-q3 poly wave 2 | 16,799 | +40 | `LSTMDrug` + `ClinicalModel` switch → classes (honest shape cost, user-ordered) |
| p8-xml | 16,764 | −35 | six measured folds in the xml loaders, first-error order preserved |
| naming waves 1–2 | — | reverted | short names tried, reverted by user order — see next row |
| p8-ns namespace flatten | **16,233** | −488 | one `namespace OpenMalaria` (util/Host/WithinHost/Transmission/... namespaces deleted); naming-wave churn reverted; collision audit in the commit message |
| perf hunt E-series | **16,317** | +84 | see §6 (net code added buying the speed) |

Every row: gates green at merge, main-session re-verified. `HANDOFF.md` has
per-phase detail; `git log --oneline rewrite` is the full ledger (~200
rewrite commits, each individually gated).

## 5. Performance result

Bench: `util/example/example_scenario.xml`, `taskset -c 5`, serialized
interleaved rounds, median-of-5 user seconds. Same-window comparisons only.

| pop | original A0 | competing branch | **rewrite (shipping recipe)** | vs original | gap left |
|---|---|---|---|---|---|
| 2,000 (KPI) | 12.24 s | 10.47 s | **10.68 s** | **1.146×** | 2.0% |
| 10,000 | 55.6 s | 44.97 s | **45.82 s** | **1.214×** | 1.9% |

The "competing branch" is upstream DRAFT PRs 476–482 applied to the original
(`eval/bs-all` + the LTO/PGO recipe). The rewrite matched its algorithmic
content *and* kept −13% code.

## 6. The perf experiments (E-series), honest accounting

| id | change | verdict | effect |
|---|---|---|---|
| 480-hoist | ageFactor hoist (from PR 480) | merged `752957ca` | −1.4..−2.1% |
| 481-recipe | Release=LTO default-on + `util/pgo-build.sh` (PGO trained on the KPI scenario, `-ffp-contract=off` kept) | merged `752957ca` | ≈−8% |
| E1 | fuse vectorUpdate's two slab sweeps into one host-major pass | merged | −2.4% |
| E4 | batched stable compaction in `Population::update` | merged | −4.2% @10k |
| E5 | bind per-host intervention rows once per host | merged | −1.7% @10k |
| E6 | pyrogen zero-numerator fast path | **reverted** — measured *slower* in the shipping PGO form (+1.7% @10k) despite −0.8% plain. Kept at `eval/e6-pyrogen`. This is why candidates now ship-gate in recipe form |  |  |
| E7 | heap backed by transparent huge pages, in-code (`mallopt` + `MADV_HUGEPAGE`) | merged | −6.1% @10k |
| E8 | run-based `HumanColumns::compact`: deaths cluster at the age-ordered tail, so move only the runs after each gap (memmove×16, move-loop×9) instead of re-gathering 25×N rows every 5-day step for 140 simulated years | merged | −2.1% @10k |
| E9 | same trick for the `humans` vector | **rejected** — measured flat, not worth +21 LOC. Kept at `eval/e9-humans` | 0 |
| dual-PGO | train profiles on pop 2000 *and* 10000 | not shipped | −0.7% @10k but +0.2% @2k; KPI pop is 2000 |

Attribution method that finally worked (round 6, details in `PERF.md`):
LD_PRELOAD call-count shim (`util/callcount.c`) proved every library/RNG call
count is identical across trees; flat retired-instruction sampling
(`perf record -e cpu_core/instructions/`) named the hot lanes despite LTO
renaming; scratch-worktree counters (`eval/instr-counters`) turned "diffuse"
into the E8 finding.

The residual ~2% gap vs the competing branch is named (accessor
re-materialisation through the context pointer) and **declined**: invasive,
~1% ceiling, T1 risk on 140-year sims. Devirtualisation is vetoed by user
order. This is the accepted bottom.

## 7. Recovering this session on another machine

Prereqs that bit us: **gcc 16.2.1** and **GSL 2.8** (bit-exactness is
toolchain-pinned; `util/a0-recipe.json` records what A0 was built with).
CMake ≥ 3.20, xerces-c, zlib.

```bash
git clone <your-fork-url> om && cd om
git checkout rewrite

# 1. Default (LTO-on) build:
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -DNDEBUG -ffp-contract=off"
cmake --build build -j
# gates:
ctest --test-dir build -j          # want: 100% tests passed out of 64
(cd build/unittest && ./unittest)  # want: .OK!

# 2. T1 parity: rebuild reference A0 from the pinned commit + recipe:
git worktree add ../A0 cfc3c8e0
cmake -S ../A0 -B ../A0/build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -DNDEBUG -ffp-contract=off" \
      -DCMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE=OFF
cmake --build ../A0/build -j
python3 util/parity.py --old ../A0/build/openMalaria --new build/openMalaria \
        --recipe util/a0-recipe.json --tiers a,b --seeds 1-10 -j 6
# want: "63/63 identical; 0 problems"

# 3. Shipping (fast) binary — two-pass PGO recipe, ~15 min:
util/pgo-build.sh build-pgo          # trains on util/example @popSize 2000
# same three gates apply to build-pgo

# 4. Benchmark (KPI pop 2000; scaling: SCEN=util/bench-scen-10000.xml):
util/bench-om.sh /tmp/bench.log build-pgo/openMalaria
# want (this machine, gcc 16.2.1): ~10.7s @2k, ~45.8s @10k
```

Cross-window timing comparisons are invalid — always run reference and
candidate in the same window (the bench script interleaves them for you).

## 8. What is NOT pushed (and why it doesn't matter)

- **Binaries / build dirs** (`build-r1`, `build-pgo-ex6`, …): gitignored, all
  reproducible from the recipes above.
- **The A0 reference binary** itself: rebuild it (step 2). A pushed foreign
  binary would be worse than useless — parity hard-fails on toolchain drift
  by design.
- **Local scratch worktrees** (`fork-wt-hunt1..7`, `fork-wt-inst`, `fork-wt-bs`,
  …): throwaway clones of the same commits that are now on the remote;
  recreate with `git worktree add` if ever needed.
- **Scenario copies** for scaling pops: `util/bench-scen-10000.xml` is pushed;
  any other pop is `sed 's/popSize="2000"/popSize="N"/'` on the example.

## 9. Where everything lives

- `PERF.md` — rounds 1–6 ledger: every measurement, method, and declined lever.
- `HANDOFF.md` — per-phase working notes (agent history; §1 header is stale
  vs this file, trust this file for status).
- `util/parity.py`, `util/compareBin` — T1 runner.
- `util/pgo-build.sh` — shipping recipe (LTO + PGO, `-ffp-contract=off`).
- `util/bench-om.sh`, `util/bench-scen-10000.xml` — benchmark harness.
- `util/callcount.c` — LD_PRELOAD call-count shim (attribution tool).
- `util/a0-recipe.json` — the reference build's exact toolchain record.

## 10. Appendix: "SoA should be faster" — why it had to be earned, not assumed

Our model stores host state as Structure-of-Arrays columns (`HostColumns`: one
vector per state field; a host is a row index), so the intuition is that
scanning state should beat the old row-object style nearly for free. True in
principle; three measured corrections make it accurate.

**1. SoA is a bet on access patterns, not a free speedup.** It pays when loops
touch few fields across many rows. OpenMalaria's hottest loops (per-host
within-host update, `probTransmissionToMosquito`) walk MOST FIELDS OF ONE ROW
at a time — column-major storage, row-major access: every read becomes
column-base + row-index instead of one `this`-relative offset, and no
cross-row SIMD is possible because the T1 contract freezes RNG draw order and
the pyrogen-division recurrence is serial. The model's own semantics block
SoA's headline win in exactly these loops. Worse, flatness only holds for POD
elements: 9 of 25 columns hold per-row containers (vectors/maps), so a density
read chases THREE levels — column base → row container header → its heap
buffer — partly defeating the layout precisely where it's hottest.

**2. The original code was ALREADY columnar** (`HostColumns` is upstream;
`rewrite` kept it). So the round-3 "structural tax" was never SoA vs AoS — it
was four concrete things, three of which we fixed:

| tax item | why it hurt | fixed by |
|---|---|---|
| gather-compaction per step | re-gathered O(25 cols × N rows) EVERY 5-day step for 10,293 steps, while theirs shifts only the short tail past the clustered deaths | E4+E8 (run-based memmove of column runs — itself a pure SoA win) |
| two-pass vectorUpdate | two host-major slab sweeps where one pass does | E1 |
| per-(host,species) row walks | getters re-walked ctx→column→row→iv per call, not LTO-inlinable across 3 TUs | E5 |
| scattered heap | 25 growing allocations spread the columns over more 4K pages than their denser rows | E7 (in-code THP) |
| ctx pointer chains | state moved from file-scope statics (one rip-relative load) to `this→ctx→subobject` re-materialisation | OPEN — ~1% ceiling, declined |

**3. After the fixes, SoA visibly wins its own scorecard** (pop 10,000, same
recipe, same windows, ours vs patched baseline): dTLB misses 0.05M vs 277M
(≈5,400×), cache misses 34M vs 61M, IPC 2.69 vs 2.64 — and vs the original
A0: 6.2M mallocs/run vs 164M. What remains of the gap is instruction COUNT,
not memory behaviour: 662.8B vs 635.3B retired instructions (+4.3%) —
per-element address re-materialisation in the row-major hot lanes — which the
memory lead cashes back down to the 1.9% wall-clock residual.

**Verdict:** SoA delivered exactly what it promises (locality, low TLB/cache
pressure, whole-column memmove runs, one-pass column algorithms) and cost
exactly what it promises (addressing overhead for row-major access, extra
indirection for container columns). The tax we paid was algorithmic and
codegen-level — not the layout — and at larger populations, where the memory
lead grows and the instruction surplus doesn't, the layout is the side to be
on. That the layout also cut per-column code is why it survived every phase of
the LOC rewrite.
