# Phase-1 SIGNOFF — census + parity annotations (branch `p1-census`, base `p0-green`)

## Deliverables & counts

| Deliverable | Artefact | Counts |
|---|---|---|
| statics | `census/statics.csv` (+ raw, AST evidence) | 347 rows: tag init 253 / runtime 75 / checkpointed 17 / n-a 2; **run_context yes 56** (plan §5 "≈55"); test-reset 11; needs-review **2** (both `n/a`: dead-file `genotypes`, NDEBUG-only `in_update`) |
| hierarchies | `census/hierarchies.csv` | 17 polymorphic roots, 92 classes covered, zero orphans (parser-verified); vtable-heavy 9 roots with plan §5 dispositions |
| XML glue | `census/xml_glue.csv` | 140 referenced scnXml types + 1 fn + 33 schema includes; 95 flagged phase-2 `hoist`; `scnXml::IPTDescription` DEAD (XSD-removed) |
| draws + annotations | `census/draws.csv` + `// parity:` in `model/` | **139 draw/contract sites; 108 comments inserted at every site** (tight-loop groups share one; bernoulli-last keeps its own); 21 carry §3 letters a–g (a9 d4 e3 f2 c1 b1 g1); streams host-local 34 / local 84 / master 7 / fixed-seed 1 / gsl-shared 6 / n-a 7 |
| latent bugs | `census/latent_bugs.md` | 13 plan items re-verified + 5 new finds |
| mon measure ids | `census/mon_measure_ids.csv` | 90 defs rows + 12 obsolete rows; compare gate PASS vs A0 (18 distinct ids, order intact) |

**LOC policy:** comments are code-LOC-exempt per user ruling 2026-11-19. Verified:
model code LOC **unchanged by phase 1** — cloc at phase-1 HEAD~ vs HEAD both 18,839
(row 0's 18,813 predates phase 0's own +26 seed/UB-fix lines, recorded in row 0d).
Comment lines 9,119 → 9,227 (+108). Diff proof: 108 insertions, 0 deletions,
every added line `// parity:` — scripted PASS (`./util/census/draws.sh`).

## Gates run on the annotated tree

- comment-only diff proof: **PASS** (0 non-comment additions, 0 removals)
- build-p1 rebuild: **OK** (Release, `-ffp-contract=off`, compile-commands on)
- `ctest -j24`: **52/52 passed** (9.5 s)
- tier-b T1 spot check: `parity.py --tiers b --seeds 1-2` vs A0 reference
  (recipe hash ok): **2/2 bit-identical**
- `monids.sh compare` old-vs-new measure-id streams: **PASS** (scenario10 + example)

## 5 most surprising load-bearing facts (phase-1 findings)

1. **Object-construction order IS stream order**: `master_RNG` never draws —
   every entity `LocalRng` copy-seeds from it in its ctor (`Human.cpp:111`,
   `VectorModel.cpp:195`). Reordering any Human/VectorModel construction
   sequence changes downstream draws.
2. **The MAX_INFECTIONS clamp has two RNG-visible modes**: at cap,
   `importInfection` skips the genotype draw entirely, while `update()`
   draw-then-drops (caller Poisson already consumed). Both must be reproduced
   by the phase-6 slot allocator (`census/latent_bugs.md` #4).
3. **`MasterRng(0,0)` at `PerHost.h:89`** (availability percentiles) is an
   intentional, results-frozen fixed seed ("make sure old results will not
   change") — annotated, do not clean up.
4. **Deploy paths are deliberately asymmetric** (§3a): timed path consumes
   beta/heterogeneity draws BEFORE availability filters; continuous path keeps
   `uniform_01` strictly last. Both frozen; phase 7 must not harmonise.
5. **`AnophelesModel` is not a `VectorModel` subclass** (composition — parser
   truth, plan §9 confirmed), and the decision tree's randomness is a single
   `uniform_01` per node into a `map<double,…>` of cumulative keys — numeric
   key order (not insertion) is the contract (`CMDecisionTree.cpp:314,:764`).

## Needs-review statics (complete list, 2 rows)

- `model/Host/WithinHost/Infection/DescriptiveIPTInfection.h:95 genotypes` — n/a (dead source, see new-bugs N1)
- `model/sim.h:265 in_update` — n/a (NDEBUG-only; keep Debug asserts, not release state)

## New latent bugs (full detail in census/latent_bugs.md)

- N1 `DescriptiveIPTInfection.{h,cpp}` dead **and un-buildable** (names
  XSD-removed `scnXml::IPTDescription`) — delete at phase 2.
- N2 `sim::in_update` exists only in Debug builds (census scoped not-in-build).
- N3 shared mutable GSL integration workspaces (`gsl_intgr_wksp`,
  `gsl_intgr_conv_wksp`) — phase-9 GSL porting hazard.
- N4 timed-vs-continuous deploy draw/filter asymmetry frozen (annotated).
- N5 master-stream ctor-seeding order fact (see fact 1).

## mon measure-id map confirmation

`defs[]` index == internal Measure id (90 = MeasureCount); `outId` == the
`measure` token of writeText / OMOUTBIN v2 measure blocks; store order =
`reportedMeasures` **sorted by outId** (`mon/init.cpp`) — not defs order;
`nNMFTreatments` outId=-1 never output; 12 obsolete names rejected in
`<enable>`; ctsout report order = alphabetical registered names. Drift gate:
`./util/census/monids.sh compare OLD NEW`.

## LOC delta table (cloc, `--include-lang="C++,C/C++ Header"`)

| tree point | model code | model comment |
|---|---|---|
| row 0 baseline (cfc3c8e0) | 18,813 | 9,117 |
| phase-0 end (0d, 625bc3d6) | 18,839 (+26 seed/UB) | 9,119 |
| **phase-1 HEAD (this)** | **18,839 (Δ0)** | **9,227 (+108 parity comments, LOC-exempt)** |

## Re-generation commands (all re-runnable, idempotent)

```
./util/census/statics.sh build-p1 && python3 util/census/statics_tag.py   # census/statics*.csv
python3 util/census/hierarchies.py                                       # census/hierarchies.csv
python3 util/census/xmlglue.py                                           # census/xml_glue.csv
./util/census/draws.sh                                                   # census/draws.csv + annotations + proof
./util/census/monids.sh                                                  # census/mon_measure_ids.csv
./util/census/monids.sh compare /tmp/opencode/om-ref/A0/openMalaria build-p1/openMalaria
```

## Needs from user (explicit confirmations)

1. **XSD decision (F7) stays: keep codegen + validation**, phase 2 builds
   `model/xml/` conversion glue only (budget 1.6–2.2k). Confirm.
2. **Phase-2 go-ahead** on `p1-census`-successor branch (xml conversion split +
   CMake GLOB_RECURSE switch), T1-gated. Confirm.
