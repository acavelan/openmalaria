# T2 Fixture Manifest — phase-0 harvest from the Rust port clone

Harvested 2026-11-19 from `/tmp/opencode/om-refactor` (branch clone of the Rust port).
Provenance method: every `cpp_*` dataset was traced to its generator
(`scripts/generate_validation_suite.py`) and cross-checked against literals already
committed in this repo's `unittest/` suites. Datasets were **not** trusted by name.

JSON telemetry itself: `validation_suite_data.json` (clone root), `generated_at: 2026-11-13 02:48 UTC`.
Binary-generating script: `scripts/generate_validation_suite.py::_run_cpp_xml` (L59–75) invokes
the compiled legacy C++ binary `build/openMalaria` and parses its `output.txt` / `ctsout.txt`.

## Table A — NEW fixtures harvested into `unittest/fixtures/` (genuinely binary-generated, not previously in this repo)

| fixture file | what it pins | provenance (path + evidence) | planned consuming suite | gaps / notes |
|---|---|---|---|---|
| `fixtures/ctsout_vecfull_p120_s7.csv` (10.1 KB, 244 rows) | `ctsout.txt` spot series (EIR, human infectiousness, immunity Y, transmitting humans) — full-model deterministic cycle, VecFullTest popSize=120 iseed=7, every 30 days | clone `validation_suite_data.json :: deterministic.continuous_ctsout.cpp_*`; script L249–280 runs `build/openMalaria` live and reads back `ctsout.txt`. Port claim: max rel err 4.70e-6 vs Rust on all 16 columns | new `CtsoutSpotSuite` (T2.5 spot-check of tick loop + entomology + immunity); T3 harness input | downsampled stride-6 + per-column rounding (5/3 dp) → compare rtol ≤ 1e-4, never bitwise; full-resolution ctsout goldens remain a gap (must generate from binary at popSize-2000 benchmark scenario) |
| `fixtures/multiseed_cpp_fan.csv` (9.8 KB, 210 rows) | C++ 30-seed survey percentiles (p5/p25/p50/p75/p95) for 5 archetypes (vec_full, itn_control, molineaux_pkpd, vivax_relapse, pev_mda) | clone JSON `multiseed_ensembles.archetypes[*].cpp_fan`; script L521–682: 150 live `build/openMalaria` runs; port claims max_rel_err = 0.0 over 1 206 510 survey cells, KS D = 0 | T3 statCompare-sampler battery + post-break ensemble regression (plan §5 gate) | 3-dp rounding → rtol 1e-3; only 2 pinned measures (nPatent, nUncomp); fan axis downsampled to ≤55 surveys; other 88 measures have no C++ ensemble reference |
| `fixtures/multiseed_cpp_totals.csv` (3.7 KB, 150 rows) | C++ per-seed totals of the pinned measure per archetype (seeds = iseed 1..30) | same source; `totals_cpp` key, same 150 live runs | T3 ensemble gate; sanity pre-flight before long runs | `qq_cpp`/mean/sd in JSON are derivable from these — not duplicated; 4-dp → rtol 1e-4 |

## Table B — C++-generated reference values ALREADY in this repo (referenced in place; NOT copied into `fixtures/` — `diff -rq` proves the clone's `unittest/` is byte-identical to ours)

| existing file | what it pins | provenance evidence | consuming suite (port-side tolerances, for reference) |
|---|---|---|---|
| `unittest/MolineauxCirDens.txt` | Molineaux pairwise 86-day circulating-density trajectory, seed (1095, 721347520444481703), mass 71.43 | repo-native golden; consumed as golden by clone `python/tests/test_deterministic_components.py` (rtol 1e-5, atol 1e-8) and `crates/om-within-host/tests/cpp_golden_infections.rs` | `MolineauxInfectionSuite` (T1/T2) |
| `unittest/PennyCirDens.txt`, `PennySeqDens.txt` | Penny 137-day circulating + sequestered densities, same seed | repo-native; clone tests rtol 1e-5, atol 1e-8 | `PennyInfectionSuite` |
| `unittest/PennyInfectionSuite.h` L57–59 | Penny thresholds N/C/V = 30939.0831 / 350.53964 / 315.21947 — matches JSON `penny_golden.thresholds_cpp` verbatim | literal match fork `PennyInfectionSuite.h` ↔ clone JSON | `PennyInfectionSuite` |
| `unittest/MolineauxStatsOrig`, `MolineauxStatsPairwise`, `+RG` | 200-cohort quantile tables (c5/q1/med/q3/c95 × 9 stats) | repo-native; clone `cpp_golden_infections.rs` compares rel+abs tol 1e-4 | `MolineauxInfectionSuite` cohort tests |
| `unittest/EmpiricalInfectionSuite.h` L69–124 | Empirical auto-regression step densities (2.3427793, 0.42220637, 1.2290052, 0.4272347, 0.042220637, 0.16787395) — literal match with clone tests (rtol 1e-6) | literal match fork `.h` ↔ clone `.rs`/`.py` | `EmpiricalInfectionSuite` |
| `unittest/PkPdComplianceSuite.h`, `PkPdSingleCompartmentDrugs.json` | PK/PD daily concentrations + drug factors for MQ, PPQ, PPQ2 (Hodel), PPQ3 (Tarning), AR→DHA, 50 kg — the JSON `pkpd.profiles.cpp_*` arrays are verbatim copies (checked `0.378440101`, `0.03174581`, `0.08022449`, `0.075305088`, `0.0001825220`) | literal match fork ↔ clone script L156–191 | `PkPdComplianceSuite`; port claim max conc rel err 6.2e-6 |
| `unittest/CMDecisionTreeSuite.h` | decision-tree fractions / outcome asserts (31 assertions) | repo-native golden suite | `CMDecisionTreeSuite` |
| `unittest/DecayFunctionSuite.h/.py` | DecayFunction shapes (exp/Weibull/Hill/smooth-compact/linear/step) | repo-native; clone JSON `decay_series` adds only Rust-vs-analytical comparison — no new C++ data | `DecayFunctionSuite` |
| `unittest/AgeGroupInterpolationSuite.h/.py` | piecewise + linear age interpolation | repo-native; clone JSON `age_interp` is Rust-only — no new C++ data | `AgeGroupInterpolationSuite` |
| `unittest/MLCS_*.txt` (7 files) | mosquito life-cycle survey-out goldens | repo-native | `MosqLifeCycleSuite` |
| `unittest/ChaChaSuite.h`, `XoshiroSuite.h` | RNG state-word vectors | repo-native — **skipped per task instructions** | RNG suites |

## Table C — clone JSON datasets REJECTED (not C++-generated; verified in generator source)

- `decay_series` (+`analytical`): Rust `om.DecayFunction` vs Python-closed-form (script L110–134) — no C++ involved.
- `age_interp` (L137–144), `fourier` (L147–151), `pkpd.genotype_*_survival` (L221–231): Rust-only evaluations.
- `molineaux_golden.cpp_log10`, `penny_golden.cpp_cir/seq_log10`: lossy log10/5-dp re-encodings of Table-B files (L301–316) — duplicate, lower fidelity.
- `pkpd.profiles.cpp_*`: verbatim duplicates of Table-B files (evidence above).
- `stochastic_components.*` (variates Gaussian/Gamma/Beta/LogNormal, Poisson PMF, RDT calibration, vivax broods/latency, decay heterogeneity fan): drawn from **Rust** RNG/samplers or theory formulas (L322–518); zero C++ runs.
- `headline_metrics`, `qq_cpp/mean/sd` (derivable), `sample_spaghetti_rust`: not reference data / Rust-side.

## Table D — excluded, provenance NOT confirmable

- `crates/om-monitoring/tests/monitoring_tests.rs::cpp_g6_matches_golden_ostream_examples` (9 value→string pairs). Embedded Rust literals; no generator, no binary-run evidence. Plausibly C++ `ostream precision(6)` behaviour but not provable from files → excluded from fixtures (formatter output is pinned bitwise by T1 `output.bin` anyway).

## Gap table — plan §4 T2 areas vs available cpp_* goldens

| T2 area | C++ golden exists? | where | action |
|---|---|---|---|
| Molineaux trajectory + cohort stats | yes | Table B (repo `unittest/`) | none |
| Penny densities + thresholds | yes | Table B | none |
| PK/PD reference values | yes | Table B | none |
| DecayFunction shapes | yes (suite expectations; analytical origin) | Table B | optional: harvest binary-emitted decay curves if a post-break re-bless wants binary origin rather than formula origin |
| AgeGroupInterpolation | yes (analytical origin) | Table B | same note as DecayFunction |
| Decision-tree fractions | yes | Table B `CMDecisionTreeSuite.h` | none |
| **CalcNvO-linear** | **GAP** — no golden in repo or clone (`core_tests.rs` has structural tests only, no pinned values) | — | generate from C++ binary: solver step/reference values for the linear-Net Reproductive-output case (later phase-0 binary-harvest step) |
| Percentile table | **CLOSED 2026-11-19** — was: no golden anywhere | `PercentileCFRSuite.h` (in-process golden: `PerHostAnophParams::calcAvailabilityPercentiles()` with MasterRng(0,0), 1 species, gamma CV=0.3, factor 1.0, bit-exact 100-entry table + p=100/out-of-range behavior) | no binary harvest needed — reachable via public headers |
| CFR curve | **CLOSED 2026-11-19** — was: no golden | `PercentileCFRSuite.h`: `getCommunityCFR` is a pure 2-line formula (`ClinicalModel.cpp:115-118`) reading the file-global `OM::Clinical::oddsRatioThreshold` (external linkage, set directly in-test); 7 bit-exact literals. Age-specific CFR *interpolation* stays covered by `AgeGroupInterpolator` machinery/suite | no binary harvest needed |
| RNG primitives + GSL-backed distributions | **CLOSED 2026-11-19** — was: only state-word vectors (`XoshiroSuite`, `ChaChaSuite`), no per-draw distribution/state lock | `goldens/rng_dist_v50.csv` + `RngGoldenSuite.h` (see section below) + `SamplerBatterySuite.h` (KS/moments) | phase-9 inline sampler ports must keep both green |
| Full-run continuous spot checks | harvested (new, low-res) | Table A | full-resolution ctsout golden for the popSize-2000 benchmark scenario still to be generated |
| Per-measure ensemble statistics (88 measures) | partial (2 measures only) | Table A | extend multi-seed harvest at phase-0 freeze if T3 needs more measures |

## 2026-11-19 test-infra additions (unit-side, no model/ changes)

New files:

- `unittest/golden-gen/RngGoldenAdapter.h` — minimal adapter around the `util/random.h` engine + GSL
  plumbing mirroring `model/util/random.h` verbatim (file:line cites in comments); seeds
  Xoshiro256P(1,2,3,4) and exposes post-draw state via the model's own public
  `Xoshiro256P::binary_checkpoint`.
- `unittest/golden-gen/rng_goldens_gen.cpp` + `build-goldens.sh` — the only supported way to
  regenerate `unittest/goldens/rng_dist_v50.csv` (14,012 rows: 1000 draws × 14 series —
  uniform_01, gauss(polar), log_normal, gamma a=0.5/2.5, beta (0.5,0.7)/(2.0,3.0), poisson
  3.2/25.0, weibull, binomial 10·0.3/60·0.05, exponential, max_multi_log_normal — plus 12
  `gsl_cdf_ugaussian_Pinv` spots; per row: value %.17g + 4 hex state words).
  Gaussian note: the model path is `gsl_ran_gaussian` (polar Box-Muller, random.h:155);
  `gsl_ran_gaussian_ziggurat` is on no model code path (verified by grep), so no ziggurat
  tail seed applies — recorded in the CSV header. `gsl_ran_binomial`/`gsl_ran_exponential`
  have no random.h wrapper; pinned via direct GSL calls, recorded likewise.
- `unittest/RngGoldenSuite.h` — replays the CSV through the adapter; bitwise asserts on value
  AND all 4 state words per draw (RED-checked: corrupting one CSV digit fails with
  "gauss(…) draw 500: value bits differ"). Phase-9 regression gate for the primitive layer.
- `unittest/SamplerBatterySuite.h` — 50k-sample KS(>0.005, in-file Q_KS series impl, no
  scipy) + mean/std 5-sigma z-tests (kurtosis-aware std SE) against GSL CDFs/exact moments,
  through the REAL `util/random.h` API (binomial/exponential via the adapter); regime
  boundaries per task. Total unittest binary runtime stays ~0.15 s.
- `unittest/PercentileCFRSuite.h` — see CLOSED rows above.

Triage — `DescriptiveInfectionSuite.h` (revived):

- Why it was dead: it was a test-LESS skeleton (zero `test*` methods; old comment said
  determineDensities was "not so easy" to test) whose setUp no longer compiles against the
  current API — compile evidence (g++ 16.2.1, C++20):
  `Infection_init is not a member of UnittestUtil`,
  `initParameters is not a member of DescriptiveInfection` (now `init(const Parameters&)`),
  `gsl has not been declared` (old RNG bootstrap removed),
  `no matching function for call to DescriptiveInfection::DescriptiveInfection()`.
  So: the tested-INTENT functionality (DescriptiveInfection) still exists, but the file
  exercised REMOVED API surface — hence revive-by-rewrite, not un-comment-as-is.
- What it now pins: `infectionDuration(LocalRng&)` — pure log-normal draw + integer
  time-step rounding — 12 iterations of duration value + 4 post-draw state words (golden
  generated with the same adapter). Remaining gap: `determineDensities()` still untested
  (needs Human/InfectionImmunity scaffolding + blessed trajectory data; same class of work
  as Penny/Molineaux suites — recommend a density-trajectory golden from the binary).

## Counts

- JSON datasets classified: **27** → harvested-as-new: **3** (the three `fixtures/*.csv`), duplicates-of-repo-goldens: **5 groups**, Rust-only/analytical: **12**, excluded-unproven: **1** (g6 pairs).
- Repo reference datasets inventoried: **11** (Table B, incl. 2 RNG suites skipped by instruction).
- Gap areas after 2026-11-19: **1 hard gap** (CalcNvO-linear) + **2 soft gaps** (full-res ctsout; wider-measure ensembles) + **1 new soft gap** (DescriptiveInfection::determineDensities trajectory golden). Percentile table and CFR gaps closed in-process; RNG distribution/state locking added.

## Bit-exact per-platform CI correction (2026-11-26)

**The allowance is zero.** The initial `a86e2dc9` cross-platform rounding allowance
was rejected by the user and removed. The approved contract is exact matching
against independently generated, frozen references for each platform.

`RngGoldenSuite` now compares the **real model wrappers** with the frozen GSL
adapter AND the selected stored baseline bit-for-bit for every draw. Percentiles
likewise match both the independent same-platform gamma/order-statistic reference
and the stored native table exactly. CFR remains bit-exact. All 14,012 RNG rows
retain the original schedule and state words; uniform/discrete expectations and
the original CSV are unchanged.

[`goldens/platforms/README.md`](goldens/platforms/README.md) records the five
platform selections, original source commits, native capture run, reproduction
recipe and per-platform provenance/hashes. The RNG generator is extracted from
the pre-rewrite fixture commit; the percentile driver calls the pinned original
model. Builds never regenerate expectations from the model being tested.
The original percentile literals now live unchanged in
`goldens/availability_percentiles.inc`.

Additional regressions exercise GSL draws after forced RNG-vector relocation and
the initial clinical state. The former catches stale self-pointers even when the
engine's raw checkpoint bytes still match; the latter catches an uninitialized
enum previously read during checkpoint restore.
