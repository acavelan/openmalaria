# Phase 4b HANDOFF — Config/State de-polymorph (branch `p4b-configstate` off `p4a-green`)

Crash-recovery file. Update EVERY green commit. Baseline p4a-green: model **19,160** / unittest **3,760**; A0 12.33s.

## Environment
- Build: `build-p4b/` = Release `-ffp-contract=off` `-DOM_BOXTEST_BINARY=ON` (configured).
- unittest MUST run from `build-p4b/unittest/` (csv resources there) → `./unittest` (91 tests).
- ctest: `cd build-p4b && ctest -j16` (52 tests, ~9.3s).
- T1 tier-a: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new build-p4b/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a -j16` (53). Tier-b: `--tiers b --seeds 1-10` (10).
- bench: `sh /tmp/opencode/om-bench/bench.sh build-p4b/openMalaria p4b` (A0-class build, vs A0 12.33).

## Done (per-hierarchy LOC deltas)
- **DecayFunction** (9 subclasses + factory in xml/util.cpp → 1 enum-tagged value type): Δ **−101** model LOC. Gates: ctest 52/52, unittest 91/91, T1 tier-a 53/53. 4 virtuals + 9 vtables gone; factory keeps error order verbatim (Emax child built BETWEEN present-checks and range-checks — preserved); RTL arg-order of old `max(f1->sAOD, f2->sAOD)` mirrored explicitly (f2 draws first) — T1-invisible in corpus (no operator/emax decays) but pinned. Checkpoint bytes: bool+double unchanged. No suite changes (API kept: makeObject/hetSample/eval/sampleAgeOfDecay/checkpointSample).


## Queue (risk order, small first)
1. ~~DecayFunction~~ DONE (see above)
2. AgeGroupInterpolation (4 classes; subclasses in xml/util.cpp too) — NEXT
3. InfectionIncidenceModel (4→1 merged config — §5 canonical)
4. Sampler (gamma/lognormal enum)
5. EmergenceRateEstimator (solver-vs-fitter tag)
6. PathogenesisModel + Submodels
7. ClinicalModel rest (3-config→2-state collapse)
8. LSTMDrug PK variants (enum + free fns)
9. AnophelesModel/SimpleMPD (optional member + 4 hooks)
10. TransmissionModel (PerHost stays for phase 5)
11. WHInterface + WHFalciparum/Descriptive/Vivax/Dummy (largest; per-species infection models inside)
12. immunity components (InfectionImmunity* — find rows in census)

## Learned / traps
- DecayFunction: checkpointSample writes bool+double only (hetFactor); restore = config-clone with factor set at EVERY node (Operator/Emax children share it). sampleAgeOfDecay ignores hetFactor (all kinds).
- Not mine: HumanInterventionComponent descendants (phase 7), TimedDeployment (7), PerHostInterventionData (5→7), base_exception.
- ext lifetime only after ensure_extensions; no header-inline deref of ctx.ext; gcc ctor-arg right-to-left eval was real (4a).
- Checkpoint bytes frozen: classes serializing via &stream keep field order.
- TwoRunsSuite = statics tripwire; rerun if file-statics move. Census refresh at end: ./util/census/statics.sh build-p4b + statics_tag.py.
