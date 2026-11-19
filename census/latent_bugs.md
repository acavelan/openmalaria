# Latent-bug ledger — Phase 1 re-verification (plan v3.3 §3/§8)

Every item re-checked against this tree (branch `p1-census`, build-p1 Release)
on 2026-11-19. Statuses: FIXED (phase 0) / LIVE (behaves as described, frozen
by parity) / DEAD (no live code path) / FIX-AT-BREAK (phase 8).

## Re-verified from the plan's ledger

1. **`surveySimulatedEIR_i/_l` uninit** — **FIXED** phase 0: ctor zero-init,
   value-neutral (53/53 scenarios bit-identical, PERF.md 0d). Closed.
2. **cspline dead path** — **DEAD, confirmed**: `gsl_interp.h`/`gsl_spline.h`
   includes (`util/AgeGroupInterpolation.cpp:29-32`) *and* the whole
   `AgeGroupSplineInterpolation` class (:234-273) sit inside `/* */` block
   comments; no `cspline.cpp` exists anywhere; setter accepts only
   `linear`/`none` and the XSD enum allows only those. No live `gsl_spline`
   anywhere ⇒ no cspline port needed. Phase 4b deletes the dead block +
   includes (census/hierarchies.csv `AgeGroupInterpolation` row).
3. **12a continuous-vs-resumed divergence** — **LIVE, gated**: 22 measures
   differ ≤8 % between single and checkpoint-resumed runs. Evidence:
   `/tmp/opencode/single-run-goldens.log` (4.4 MB, phase-0 sweep) + PERF.md
   0d (single-run gate mode, 53/53 bin-identical) + PERF.md "Known" line.
   Checkpoint veto active (goldens bless single-run only); re-bless scheduled
   phase 8; `--checkpoint-stop` drops from `test/CMakeLists.txt` then.
4. **MAX_INFECTIONS clamp** (`WHInterface.h:205`, cap = 21) — **LIVE**, and it
   has *two different* RNG-visible semantics, both frozen:
   * `CommonWithinHost::importInfection` (:127): `if (numInfs < MAX_INFECTIONS)`
     guards the genotype **draw itself** — at-cap import consumes **zero**
     stream values (vs `Genotypes::sampleGenotype` + `createInfection` when
     under cap).
   * `CommonWithinHost::update` (:148-150): caller has **already** drawn the
     per-day Poisson/bernoulli counts; `nNewInfs_{i,l} = min(...)` clamps the
     *infection count* while the draws stay consumed (`nNewInfsIgnored`).
   Phase-6 SlotAllocator must reproduce both (skip-draw-at-cap on import;
   draw-then-drop in update) or T1 breaks.
5. **`MasterRng(0,0)` percentile seed** — **LIVE, load-bearing**:
   `Transmission/PerHost.h:89-93` — `MasterRng seed(0, 0); seed.seed(0, 0);
   LocalRng rng(seed);` feeding `calcAvailabilityPercentiles` ("Using a
   different seed to make sure old results will not change"). Annotated
   `// parity:`; do not "clean up" (census/draws.csv master row).
6. **ExactNv0 box test** — **VERIFIED**: `test/CMakeLists.txt:65`
   (`OM_BOXTEST_NAMES`), NC list empty, runs via run.py `-- --checkpoint-stop`;
   checked against a blessed golden (plan §3 ground truth holds).
7. **Output-format contract (g6 / C locale)** — **VERIFIED as default-stream
   contract**: no `setf`/`setprecision`/`imbue`/`setlocale` anywhere in
   `mon/Monitoring.cpp`, `mon/Continuous.cpp`, `main.cpp` — all output is
   default `std::ostream` (g6) under the C locale; column-token structure in
   `writeText` (Monitoring.cpp) and ctsout (Continuous.cpp). Comparators
   parse at 1e-6 ⇒ fewer printed digits silently eats the tolerance margin
   (Rust `format_cpp_g6` lesson). Rewrite keeps default formatting until the
   phase-8 re-bless.
8. **`removeFromCohort` stale bits (Human.cpp)** — **RE-EXAMINED**: current
   tree is map-based: `addToCohort` sets `subPopExp[id]=expiry` (:135),
   `removeFromCohort` = `subPopExp.erase(id)` (:139-142), `isInSubPop` reads
   the map (Human.h:158-160). No bitvector ⇒ the "stale bits" form found by
   the Rust port is not present in this revision; remaining subtlety is
   expiry-time semantics (membership = `t > sim::now()` evaluated at insert
   time). Final semantic call at phase 7 (interventions rewrite); census
   flags it via `removeAtIds` (statics.csv) and the cohort call sites.
9. **`print_progress(int lastPercent…)` by-value** (`main.cpp:54`) — **LIVE**,
   cosmetic only (progress printed repeatedly; `lastPercent` never updated at
   caller). Progress text never reaches goldens ⇒ free to fix any time; not a
   parity risk.
10. **`VectorModel.cpp:381` assert** — **LIVE**: `assert(simulationMode ==
    forcedEIR)` guards the `EIR_*` reset block although species data can be
    non-empty in dynamic mode — Release builds compile it out; Debug runs may
    fire spuriously. Fix at phase 4b (assert → condition) — value-neutral.
11. **ChaCha `discard(n<16)` at ctr==0** (`contrib/chacha.h:161-165`) —
    **LIVE**: `if (idx + n >= 16 && ctr % 16 != 0) generate_block();` — with
    `idx==0, n<16` and a never-generated block, `discard` reads stale bytes.
    Only reachable through explicit discard paths (none in the current model
    flow — master/local streams use `seed()` then sequential draws). Rust's
    `was_zero` guard is the adopted UB-fix **at the phase-8 break** (changes
    values only where C++ was UB).
12. **`std::sort` NaN UB in percentile/solver paths** — **LIVE risk**
    (`PerHost.h` percentile sort, solver paths): a NaN comparator result is
    UB; current goldens never produce NaN there. Fix at break with
    `total_cmp`-style comparator (plan §3). Note census: `map<double,…>`
    iteration = numeric key order must survive any replacement (decision-tree
    `Branches_t`).
13. **`logDFT` OOB read** (`util/vectors.h:113`) — **LIVE** in-tree (used by
    `InfectionIncidenceModel` periodic forcing): Rust port found the
    past-the-end read. Fix at phase-8 break with regression test (plan §3).

## NEW finds from the phase-1 census

N1. **`Host/WithinHost/Infection/DescriptiveIPTInfection.{h,cpp}` is dead
    source, worse than uncompiled**: absent from `model/CMakeLists.txt`
    `Model_CPP` (only conditional-on-flag `StreamValidator.cpp` is
    conditionally excluded — this one is not in any target, incl. unittest);
    referenced nowhere else; and it names `scnXml::IPTDescription`, a type
    removed from the schema after `scenario_30.xsd` ⇒ it cannot compile
    against the current codegen even if re-added. Its static
    `DescriptiveIPTInfection::genotypes` carries scope `not-in-build` in
    census/statics.csv. Recommend deleting the pair at phase 2 (document in
    SIGNOFF; never in any LOC baseline).

N2. **`sim::in_update` is `#ifndef NDEBUG`-only** (`sim.h:265`): the
    update-phase assert state exists only in Debug builds; Release (gate)
    builds have no such member. Census row scoped `not-in-build`. Keep the
    assertions in the rewrite (cheap in Debug), but RunContext must not make
    it a gated release-state field.

N3. **Shared mutable GSL integration workspaces**: `gsl_intgr_wksp`
    (`PkPd/Drug/LSTMDrugThreeComp.cpp:128`) and `gsl_intgr_conv_wksp`
    (`LSTMDrugConversion.cpp:134`) are process-global pointers allocated at
    static-init and mutated by every `gsl_integration_qag` call (shared
    limit/cache state across all drug evaluations). Tagged `runtime`
    (hand-verified) in census/statics.csv. Phase-9 GSL porting must either
    keep the exact shared-workspace call sequence or switch to per-call
    workspaces behind a T2 gate (results can shift in the last bits).

N4. **Deploy-path availability/draw asymmetry** (plan §9 quirk, do not fix):
    `TimedHumanDeployment::deploy` consumes its beta/heterogeneity draws
    (`Deployments.h:300` etc.) **before** the availability/subPop filters
    (:312), while `ContinuousHumanDeployment::filterAndDeploy` checks
    availability/date/subPop first and keeps `uniform_01` last ("RNG call
    should be last test", :675-678). Both orders are frozen — annotated
    `// parity: a` at both paths; phase 7 must not "harmonise" them.

N5. **Census artefact worth knowing** (not a bug): `master_RNG` performs no
    direct draws — every entity stream is copy-seeded from it at
    construction (`Human.cpp:111`, `VectorModel.cpp:195`,
    `RNG(RNG<S>& source)` in `util/random.h:100`), so **object-construction
    order is stream order**. That is why census/draws.csv carries master-stream
    rows at the ctor sites and contract-f freezes the seed point.
