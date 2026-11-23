# MAIN-SESSION STATUS (2026-11-20, authoritative — everything below is per-agent history)
`rewrite` HEAD = p8-collapse merged onto p8-poly (merge 488ec0ea + PERF rows). ALL phase-8
waves 1–2 merged & verified: structure, checkpoint-redo, files, Q2-polymorphism, collapse+reflow.
**model 16,233 code LOC / 115 files (baseline 18,813 → −2,580).** One namespace: single `namespace OpenMalaria` (p8-ns merged; naming waves 1–2 REVERTED — RunContext/xml_scenario_error/TransmissionModel/ClinicalModel names are back and current; Context/xml_error/Model spellings in rows above that line are historical). Kept nested: scnXml/xsd, Continuous, Deploy, GT, Error/Messages, UnitTest.
**Draft-PR eval DONE 2026-11-23** (PERF.md section): upstream #476/#477/#482 already in tree (superseded); #478 duplicates our tooling; #479 ported, bit-exact but +0.16% (reject); #480 ported, bit-exact, −1.38% (accept, commit `d64310c1` in `fork-wt-c`, cherry-pick ready); #481 recipe measured on current tree: LTO −1.70%, LTO+PGO −7.93% bit-exact (script needs `-ffp-contract=off` deviation); combined 480+LTO+PGO = 11.25 s vs B0 12.36 s = **1.100×**. Bench B0 stamped **12.36 s** (MAD 0.01). A0 copy in /tmp/opencode was tmp-GC'd again — rebuild outside /tmp before the next T1-vs-A0 gate; candidate parities this round anchored on B0 (itself T1 a+b 63/63 vs A0 at merge). Gates on merged tree: build 0 err,
ctest 64/64 (incl Resume_*), unittest .OK!, T1 a+b seeds 1-10 = 63/63 vs A0. Verify build =
`build-mrg/`; repo-root `build/schema/*.xsd` must exist (parity harness resolves XSDs there —
restored from build-mrg; a missing one = harness red, not code red). No p8 worktrees/branches remain.
NEXT (user call): ~~8-Q1 de-GSL~~ **REJECTED 2026-11-20, do not re-attempt** (measured: verbatim GSL
port = +617 LOC for −2.6 %; user keeps GSL — PERF row p8-q1), 8-Q3 LSTMDrug/ClinicalModel poly wave
= DONE (merged; pattern in the P8-poly HANDOFF section at end of file / PERF p8-Q2 row), cap-vs-purity decision
(−2,759 to the 14.0k model cap; only [phase8] FP work or a [string-audit] xml rework can close it).
Bench numbers are stale vs this tree — re-measure A0-class (+ optionally the PGO recipe,
PERF `7d-A-amended`) before any perf claim.

# P8-collapse HANDOFF — collapse-sweep: bug fix, punch-list, xml hard-read, clang-format reflow (branch `p8-collapse` off `rewrite` HEAD)

SOURCE OF WORK: `/tmp/opencode/ponytail_collapse.md` (ranked punch-list #1–#23).
User order: fix the audited assert bug, execute the punch-list (T1-gated),
dedicated hard-read pass on `model/xml/interventions.cpp`, then ONE
whitespace-only clang-format reflow commit.

## Environment
- Worktree `/tmp/opencode/p8collapse`, build `build-cl/` = Release `-O3 -DNDEBUG
  -ffp-contract=off -DOM_BOXTEST_BINARY=ON`; ignored symlink `build -> build-cl`
  REQUIRED (parity.py resolves generated XSD via REPO/build/schema).
- Gates per commit: build 0 errors; `ctest` in build-cl (64);
  `cd build-cl/unittest && ./unittest` → ".OK!" (MUST run from that dir —
  suites open data files relative to CWD); T1 `python3 util/parity.py --old
  /tmp/opencode/om-ref/A0/openMalaria --new build-cl/openMalaria --recipe
  /tmp/opencode/om-ref/A0/recipe.json --tiers a -j12` (53/53); + `--tiers b
  --seeds 1-4` on deployment-parse commits. Commits unsigned.
- Baseline at branch start: all gates green, model 17,192 code LOC / 115 files.
- Trap re-confirmed: gate outputs are stale-binary traps — check the build
  compiled clean BEFORE trusting ctest/parity.

## Result: model 17,192 -> 16,939 code LOC (−253), 115 files unchanged.
Final gates at HEAD 8738f0c1: build 0 err, ctest 64/64, unittest .OK!,
T1 tier-a 53/53 bit-exact vs A0.

| Commit | item | LOC |
|---|---|---|
| 5e22eb3f | BUGFIX punch #23: bogus `assert(simulationMode==forcedEIR)` inside the transientEIRknown branch (TransmissionModel.cpp:307) — assert builds abort on any scenario with transient EIR + non-forced mode; deleted (release builds never saw it; assert-build-only = T1-safe) | −1 |
| a382cc38 | #2: 9 VecStat cts callbacks -> one `ctsCbVecStat(VecStat, stream)` | −48 |
| 820b9904 | #3: 5 identical nhh reduce loops -> applyReduces member-pointer lambda (+ loop-invariant hoist inside) | −25 |
| 2f4f17c1 | #4: dead XML_ASSERT/OM_ASSERT_EQ*/do_assertion* family deleted from errors.h (grep: zero uses) | −19 |
| c61178a8 | #5: ctsCbAlpha/P_B/P_CD -> shared meanOver left-fold helper | −9 |
| e5d3125d | #11: cts column-header block -> ctsHeads(prefix) lambda (dropped never-registered ctsIRSEffects) | −20 |
| a3cad600 | #7: nhhInstances expiry purge -> std::erase_if | −8 |
| 831c5489 | #9/#10/#18/#20: range-fors, find_if, tautological asserts in TransmissionModel.cpp | −13 |
| 275a5141 | #12: NeonatalMortality running-max -> max_element | −4 |
| c0efd522 | #13: updateSurveyNumbers scan -> find_if + ternary | −2 |
| b9fcb66c / 5e38a49a / a7f8eba4 / 98d7b2bd | #14 vaccinate() scan -> ranges::find_if; #15 isValidCondition -> any_of; #17 vectors::sum -> accumulate + scale range-for; #22 SpeciesIndexChecker insert-then-check (batched gates) | −11 |
| 99d0de06 | #16: shared consts for byte-frozen duplicated error messages in xml/interventions.cpp (multiset-diff proof of identical emitted strings; +2 LOC, punch est −3: literals were inline throw args) | +2 |
| 57fb3558 | xml hard-read: ITN init_kill ladder -> lambda; read_avail_limits ladder dup fold | −17 |
| 226b5c9c | #21 + extras: range-fors in xml/interventions.cpp (incl. redundant `size()>0` guards dropped) | −8 |
| 8738f0c1 | clang-format reflow (see below) | −70 |

## Skipped / rejected (with proof)
- #1 CalcNvO, #6 WHVivax cout traces, WithinHost/LSTMDrug/Clinical bodies:
  user veto / other agent owns. Untouched by design.
- #8 Deployments.h copula try/catch unwrap: **UNSAFE** — `set_gsl_handler`
  (run.cpp:226) converts GSL errors to traced_exception; `rng().beta`
  (ux≥1 branch) and `gsl_cdf_beta_Pinv` can throw outside the alpha<0||beta<0
  guard (guard misses alpha==0 and lives only in the ux<1 branch); unwrapping
  would change message + exit code on those inputs. Proof in commit msg
  a3cad600.
- #19 writeRows/countRows: filters diverge (ValueType-enum vs isDouble);
  shared lambda would be LONGER than the two bodies.
- Percentile hoist (getEntoAvailabilityPercentile out of deploy loops):
  REJECTED — its throw is conditioned on reaching a human in the age window;
  hoisting fires on inputs that today deploy fine (validation-order law).
- Line-count-neutral/unreadable folds rejected: nhh_mods_setup 5-map resize
  ladder, AnophelesModel currentNhh loop, xml :312/:373 loops.
- Other model/xml/*.cpp duplicate messages: only 1–2-site inline-arg dups —
  dedup would cost LOC; noted, not folded.

## Reflow commit 8738f0c1 (formatting-only)
- `.clang-format` at repo root: LLVM base, IndentWidth 4, ColumnLimit 120,
  UseTab Never, PointerAlignment Left, SortIncludes Never, ReflowComments
  false, **FixNamespaceComments false** (else clang-format 23 injects
  `} // namespace X` = non-whitespace), **BreakStringLiterals false**.
- 134 of 144 checked-in files formatted; 10 left byte-identical:
  9 where clang-format must inject non-whitespace to re-wrap (macro
  continuation backslashes: errors.h, ExtraAsserts.h, MolineauxInfectionSuite.h,
  PkPdComplianceSuite.h, ClinicalModel.h, LSTMDrugType.h, LSTMMedicate.h,
  IRS.h, ITN.h) + unittest/TestPaths.h (configure_file template: clang-format
  reads `@CMAKE_VAR@` as an Obj-C attribute, inserts a space, breaks @ONLY
  substitution — the ONLY place whitespace is semantic here).
- PROOF (run before commit): for every touched file
  `git show HEAD~:<f> | tr -d '[:space:]' | md5sum` == working-file md5 → 0
  mismatches. Generated schema/configured headers untouched.
- Model LOC after reflow: **16,939 / 115**; unittest 3,590.
- NOTE for future formats: re-exclude those 10 files (or wrap in
  `// clang-format off`, which is a non-whitespace change — separate commit).

---

# P8-files HANDOFF — merge files that belong (branch `p8-files` off `rewrite` HEAD = p8-sweep green ad640d22)

User order (verbatim): "don't touch CalcNvO, review the structure again and try
to merge files that belong".

## Environment
- Worktree `/tmp/opencode/p8files`, build `build-fl/` = Release `-O3 -DNDEBUG
  -ffp-contract=off -DOM_BOXTEST_BINARY=ON`; ignored symlink `build -> build-fl`
  present (parity.py resolves generated XSD via REPO/build/schema).
- Gates per merge commit: `cmake --build build-fl -j14` (0 errors);
  `cd build-fl && ctest -j12` (64); `cd build-fl/unittest && ./unittest`
  (".OK!"); T1 `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria
  --new build-fl/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json
  --tiers a -j12` (53); + `--tiers b --seeds 1-4` on commits moving
  WithinHost/Clinical/PkPd BODIES. Final: tiers a + b seeds 1-10 = 63/63.
  Commits unsigned (`--no-gpg-sign`). No bench (structural, include-level only).
- Baseline at HEAD before first merge: build ✓, ctest 64/64, unittest .OK!,
  tier-a 53/53, model 17,545 code LOC / 140 files.
- TU-count tool (demography-trap check): transitive .cpp-include closure,
  RunContext.h cone = 58/65 TUs; anything inside that cone stays a cpp.

## Merge ledger (one commit per group)

| # | merge | files Δ | gates |
|---|---|---|---|
| G12b | interventions: IRS/ITN/Vaccine/InterventionManager.cpp -> their headers (class A x4, cones 3-6 TU: each file was the header's non-XML behaviour — XML ctors stay xml/interventions.cpp (schema boundary); all moved defs inline; Vaccine.h +Host/Human.h, using-namespace OM::util dropped (was TU-local); bodies + error strings verbatim; one components-mega-file REJECTED: distinct component concepts, Components.h is their TOC) | -4 | build, ctest 64/64, unittest .OK!, T1a 53/53 |
| G12 | mon/Continuous.cpp -> Continuous.h (class A: the 4 live functions — 3 registerCallback overloads + update — inlined into the namespace's header, 14-TU cone; init() already lives xml-side; dead xml_scenario_error using dropped; Monitoring.h merge REJECTED: Continuous.h includes Population.h, pushing it into Monitoring.h's 44-TU cone) | -1 | build, ctest 64/64 (Resume_* compares ctsout), unittest .OK!, T1a 53/53 |
| G11 | Host/WithinHost/State.h -> PathogenesisModel.h (class B: the Pathogenesis::State enum + StatePair live in the Pathogenesis namespace the model file declares — pre-7g they were one directory; 4 includers repointed; no new includes reach the 37-TU cone beyond Parameters.h/random.h which it already had via WHInterface/Infection.h) | -1 | build, ctest 64/64, unittest .OK!, T1a 53/53 + T1b 1-4 |
| G10 | run.h -> RunContext.h (class B: run_simulation is the run-core entry declared beside the context it takes; decl-only, +10 LOC to the lean header, no new includes) | -1 | build, ctest 64/64 (TwoRunsSuite), unittest .OK!, T1a 53/53 |
| G9 | PkPd: LSTMTreatments.h -> LSTMMedicate.h (class B: the loader declaration for exactly the registry LSTMMedicate.h defines — Schedule/DosageTable/TreatmentRegistry + LSTMTreatments::init read as one 'drug course config' file; decl-only, scnXml fwd-decl kept at global scope) | -1 | build, ctest 64/64, unittest .OK!, T1a 53/53 |
| G8 | PkPd: LSTMDrugOneComp/ThreeComp/Conversion.cpp -> LSTMDrug.cpp (class A cpp->cpp: all four files are the body of the ONE LSTMDrug final class, split by Kind when the old hierarchy collapsed; every variant equation + both GSL workspaces verbatim, per-file include sets unioned; LSTMDrug.h cone is 59 TU, so bodies stay out of headers) | -3 | build, ctest 64/64, unittest .OK! (LSTMPkPd/PkPdCompliance suites), T1a 53/53 + T1b 1-4 |
| G7 | Host/WithinHost/Infection/DummyInfection.cpp -> Infection.cpp (class A cpp->cpp: Dummy is the dispatcher's own simplest Kind; Infection.h is a 43-TU header, so bodies stay in the cpp — Dummy variant + infection_new/infection_update dispatch now one file; stale empty 'checkpointing' section comment replaced; bodies verbatim) | -1 | build, ctest 64/64, unittest .OK!, T1a 53/53 + T1b 1-4 |
| G6 | Host/NeonatalMortality.cpp -> .h (class A: model math + update loop inline; run-state globals become C++20 inline vars, external linkage kept for the extern refs in xml/healthSystem.cpp:init and checkpoint/statics.cpp:staticCheckpoint — both files untouched; 6-TU cone) | -1 | build, ctest 64/64 (Resume_* incl. neonatal wire), unittest .OK!, T1a 53/53 + T1b 1-4 |
| G5 | Clinical/CMDecisionTree.cpp -> CMDecisionTree.h (class A: the flat node-array interpreter is the header's own API — exec_decision_tree inline; anon-namespace helpers keep internal linkage; compile_decision_tree stays in xml/healthSystem.cpp (schema boundary); ClinicalModel.h include dropped (unreferenced); all consumers already under ClinicalModel/WHInterface cone) | -1 | build, ctest 64/64, unittest .OK!, T1a 53/53 + T1b 1-4 |
| G4 | Clinical/ESCaseManagement.h -> EventScheduler.h (only consumer; same EventScheduler chain — ES trees + execute() sit after esClinicalUpdate). healthSystem.cpp/HumanInterventionComponents.h/unittest suite re-pointed (include line only; no fixture edits) | -1 | build, ctest 64/64, unittest .OK!, T1a 53/53 + T1b 1-4 |
| G1 | util/SpeciesIndexChecker.cpp -> .h (both methods inline; std:: ostringstream qualification only) | -1 | build, ctest 64/64, unittest .OK!, T1a 53/53 |
| G3 | 5-day chain: DecisionTree5Day.{h,cpp}+ImmediateOutcomes.{h,cpp}+CM5DayCommon.{h,cpp} -> CM5DayCommon.{h,cpp} (bodies verbatim; local enum CaseType -> CM5DayCommon::CaseType, value-identical; trees statics stay extern for xml/healthSystem.cpp) | -4 | build, ctest 64/64, unittest .OK!, T1a 53/53 + T1b 1-4 |
| G2 | util/DecayFunction.cpp -> .h (formula switch + sampler inlined; bodies verbatim except std:: qualifications) | -1 | build, ctest 64/64 (incl. T2 DecayFunctionSuite), unittest .OK!, T1a 53/53 |

## Rejected (belonging/TU tests)
- errors.cpp(118 LOC/67 TU), Genotypes.cpp(68/59), Diagnostic.cpp(54/58),
  LSTMDrugType.cpp(78/59), LSTMModel.cpp(83/58), PopulationAgeStructure.cpp(126/58),
  UnitParse.cpp(127/39), PerHost.cpp(106/33), InfectionIncidenceModel.cpp(156/33):
  cpp→header = demography-class (body lands in the wide RunContext/PerHost cone).
- IvState.cpp→IvState.h: 199 LOC + documented leaf header (cycle break); fold
  would drag ITN/IRS/GVI registries into the leaf. Keep cpp.
- Continuous.h+Monitoring.h: Continuous.h includes Population.h (16 TU) —
  merging pushes it into Monitoring.h's 44-TU cone. Keep separate.
- ImmediateOutcomes.h+ESCaseManagement.h: different health-system families
  (5-day vs EventScheduler). ESCaseManagement.h merged into EventScheduler.h
  instead (its only consumer + same ES concept).
- transmission.h→TransmissionModel.h: it is the xml/entomology.cpp factory TOC
  incl. AnophelesModel+EmergenceRateEstimator includes; folding drags vector
  machinery into the model header. Keep.
- NeonatalMortality.h+InfectionIncidenceModel.h: different concepts (neonatal
  death event vs FOI config model); InfectionIncidenceModel.h cone = 33 TU.
- xml/util.h incl. UnitParse.h: DocumentLoader/XMLChecker pull <schema>;
  UnitParse.h cone = 39 TU → schema would hit all of them (generated-schema
  bloat rule). XMLChecker+latentP+DocumentLoader merge together instead
  (all schema-reading, ≤5 TU); name kept DocumentLoader.h (xml/util.h would
  mislabel xml/util.cpp, whose API deliberately lives in consumer headers).
- ModelNameProvider.h→ModelOptions.h: different concepts (name registry vs
  option bits); fold pushes ModelOptions.h (340 lines) into Parameters.h cone.
- One interventions "components" mega-file: IRS/ITN/GVI/Vaccine are distinct
  component concepts; Components.h is already their TOC. Per-file folds only.
- ClinicalModel.cpp/Episode.cpp→headers: cycle (Episode.h←ClinicalModel.h,
  EventScheduler.h→ClinicalModel.h). Keep cpp.
- Exclusions honoured: CalcNvO.{h,cpp}, EmergenceRateEstimator.h,
  model/checkpoint/ internals, demography.cpp, mod.h/checkpoint usings,
  ComponentId.h leaf, RunContext.h lean, all variants alive, formula bodies +
  error strings verbatim.

---

# P8-sweep HANDOFF — ponytail audit sweep (branch `p8-sweep` off `rewrite` HEAD incl. p8-CK merge; green-stop at HEAD)

SOURCE OF WORK: `/tmp/opencode/ponytail_audit.md` (ranked ledger, auditor-verified).
All T1-safe ledger items executed except the exclusions below.

## Environment
- Worktree `/tmp/opencode/p8sweep`, build `build-sw/` = Release `-O3 -DNDEBUG
  -ffp-contract=off -DOM_BOXTEST_BINARY=ON`. Ignored symlink `build -> build-sw`
  REQUIRED (parity.py resolves generated XSD via REPO/build/schema).
- Gates per commit: `cmake --build build-sw -j14`; `cd build-sw && ctest -j12`
  (64); `cd build-sw/unittest && ./unittest` (".OK!"); T1
  `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new
  build-sw/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a -j12`
  (53); + `--tiers b --seeds 1-4` on within-host/clinical commits.
  Commits unsigned (`-c commit.gpgsign=false`, headless pinentry).

## Result: model 17,632 -> 17,545 code LOC (-87). FINAL GATES AT HEAD:
build 0 err, ctest 64/64 (incl. 12 Resume_*), unittest .OK!,
T1 tiers a+b seeds 1-10 = 63/63 bit-identical vs A0.

| Commit | item | LOC | gates |
|---|---|---|---|
| HostColumns erase_at fold (audit #3) | 25×col.erase -> lambda + packed member list | -18 | a |
| WithinHost twin folds (#5,#7,#8,#13,#14) | i/l new-infection twins fold over 2-elem Side array (draw-order proof in msg: each side draws fully before the next starts, nothing draws between); tautological `else if (opt_vaccine_genotype == false)` x4 -> else; y_lag blocks -> `WHInterface::cache_y_lag()`; InfGenotypeSorter -> lambda | -31 | a+b1-4 |
| EmpiricalInfection CSV (#6) | 7-field unroll -> loop; field/value order + error strings verbatim | -9 | a+b1-4 |
| EventScheduler (#9,#10) | const-true inHospital branches dropped (getCommunityCFR still live at subsequent-days site + CM5DayCommon); pTreatment/treatmentEffectMult (=1.0) killed, **uniform_01 DRAW KEPT** as bare call | -7 | a+b1-4 |
| Wrapper inlines (#11,#12,#15,#17,#19-part) | getLatestReport/getLatestState -> public members; getDOB -> `cols->dob[row]` at run.cpp; getBaseModelOptions -> `options.reset()`; LSTMDrugType::getIndex -> public `index` at 10 sites (LSTMDrug::getIndex is a DIFFERENT API — kept) | -18 | a |
| Deletions (#16,#18,#19) | WHInterface move ctor (grep proof: std::move only on unique_ptr holders); totalInfectionrateVariance (comment-only use); nYearsWarmupData inlined at single use | -4 | a |

## Exclusions / honest skips
- CalcNvO.cpp dead solver variants (`solve_Nv0_linear`,
  `CalcXP_pre_performance_patch`, stale `maybe_unused` on `solve_Nv0_rootfind`):
  **AWAITING USER DECISION — untouched by design, not skipped by oversight.**
- Audit #2 create*/Infection factory merge: the checkpoint redo already killed
  `checkpointedInfection` + its wrappers; the 4 surviving `create*Infection`
  wrappers are live glue feeding the `RunContext::createInfection` fn-pointer
  (5 CommonWithinHost sites + unittest) — grep proves nothing dead, removal
  would be a factory redesign (out of scope). Stale comment Infection.h:103
  still names the gone `infection_new_checkpointed` (comment-only, left).
- [phase8] / [string-audit] rows and everything vetoed: not touched. Formula
  bodies, error strings, checkpoint wire untouched (no touched field
  serialization: Resume_* green throughout).

---

# Phase 8-CK HANDOFF — kill checkpointing, redo it xml-style (branch `p8-checkpoint` off `rewrite` HEAD = `p7g-green` e34997f1)

User order (verbatim): "kill checkpointing entirely (one commit) then redo the
checkpointing similar to xml/".

## Merge state: rewrite p8-structure merged at HEAD
core/ gone, RunExtensions dissolved, SlotAllocator -> Host/WithinHost/. Resolutions: old-module
## checkpoint functions rewrite had sed'd (checkpoint_common, Monitoring/
## Continuous::checkpoint, loadFromCheckpoint, InterventionManager::checkpoint)
## stay DELETED (kill); factory field `checkpointedInfection` dropped from
## RunContext + the 4 init sites (new module uses createInfection / direct
## ctors); module ctx.ext-> sites rewritten to ctx. (importedInfections and
## monitoring are now pointers: -> / *). Gates at merge: build, ctest 64/64,
## unittest .OK!, tier-a 53/53. (The "merge traps" doc from rewrite is in the
## second half of this file, below the ---.)

## Environment
- Worktree `/tmp/opencode/p8ckpt`, build `build-ck/` = Release `-O3 -DNDEBUG
  -ffp-contract=off -DOM_BOXTEST_BINARY=ON` (~60 s clean, `cmake --build
  build-ck -j14`). Untracked symlink `build -> build-ck` exists because
  util/parity.py resolves the generated `scenario_current.xsd` via
  `REPO/build/schema` (hardcoded); without it parity tier-a fails with rc=67
  for BOTH binaries — environment trap, not a code bug.
- unittest: `cd build-ck/unittest && ./unittest` → ".OK!" (93 cases; was 94 —
  CheckpointSuite's single test was killed with the feature).
- ctest: `cd build-ck && ctest -j12` → 57 (52 box tests + 5 new Resume_*:
  test/resume.py split-run equivalence, see C1). The old box-test
  `--checkpoint-stop` loop was killed with the feature; resume.py is the new
  round-trip gate (`util/parity.py` never contained checkpoint code; its
  tier-a=53 scenarios, tier-b=example×seeds both stay).
- T1: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new
  build-ck/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a
  -j12` → 53/53. Every commit.
- Commits: `git commit --no-gpg-sign` (repo sets gpgsign but existing commits
  are unsigned; pinentry dies headless).
- Pre-kill reference copy of the old code: `/tmp/opencode/ck-ref/` (+ git
  history: parent of the KILL commit has every deleted wire function).

## KILL (done — commit `9b1d46f5`)
Deleted the whole v2 checkpoint machine: model/checkpoint.{h,cpp},
model/util/checkpoint{,_containers}.{h,cpp} (712 code LOC), ~116 operator&
sites + checkpoint()/staticCheckpoint chains across ~40 files, CLI
--checkpoint/-c/--checkpoint-file/--checkpoint-stop, run.cpp save/restore
branches, util::checkpoint_error + Error::Checkpoint (later exit codes shift
−1), unittest/CheckpointSuite.h, Continuous::init(isCheckpoint), StreamValidator
checkpoint hooks (validator itself kept, re-pointed at draw-stream use),
box-test `--checkpoint-stop` + run.py resume loop.
LOC: model 18,248 → 16,719 code (−1,529); unittest 3,791 → 3,641.
Kept: `Xoshiro256P::binary_checkpoint` + `RandomBase::checkpoint` state-word
API (RngGolden/SlotAllocator/DescriptiveInfection suites use it; module reuses).
Old wire format is DEAD: checkpoint bytes may break freely (plan clause);
re-bless/round-trip goldens = the new resume tests below.
Gates on kill commit: build ✓, ctest 52/52, unittest .OK!, T1 tier-a 53/53.

## REDO design (landed C1–C4; C5 = this wrap-up)
New module `model/checkpoint/` (mirror of `model/xml/` conventions): free
functions, per-subsystem files, ONE orchestrator listing the section order.
- `Stream.h` — header magic+version u32, templated rd/wr primitives for
  trivially-copyable scalars, length-checked list helpers. No framework.
- one .cpp per area: clock+rngs / transmission(+vector cohorts) / humans
  (Population+SoA+slots+WH+Infection+clinical+PK/PD per human) / statics /
  mon / interventions; orchestrator = `checkpoint.{h,cpp}`
  `write_checkpoint(RunContext&, ..., std::ostream&)` / `read_checkpoint(...)`.
- CLI: `--save-to FILE [--save-at DAY]` (default: start of main phase) and
  `--resume FILE`. Single file, uncompressed, one version int.
- Config registries (drugs/treatments/scenario tables) rebuilt from the
  scenario on load; only true run state goes on the wire (verified per case).
- State classes lose `private:` where the module needs them (repo already went
  all-public for Population/Infection/Anopheles; same idiom elsewhere).
- NO per-host intervention state on the wire (skill law 2026-11-20:
  interventions never checkpointed; restore replays deployment events). Proven
  days are therefore ≤ last deploy time; during-intervention days are PROVEN
  ONLY if event replay reproduces per-host records bit-exactly — if it doesn't
  (rng draw order!), HANDOFF states that honestly as the remaining gap.

## Steps (TDD ladder; every commit: build + ctest + unittest .OK! + tier-a 53)
- [x] P0 KILL (one commit).
- [x] C1 module + CLI + resume tests. Landed: model/checkpoint/ (Stream.h
      Archive rd/wr-on-one-class, sections.h, checkpoint.{h,cpp} TOC,
      statics.cpp, transmission.cpp, humans.cpp), Error::Checkpoint,
      CLI --save-to/--save-at/--resume, Continuous::init(isCheckpoint)
      resume-open + seek, run.cpp branch (skip createInitialHumans/init2/
      warmup/calibration/reset per saved phase). test/resume.py per scenario:
      full run vs save@main → resume+save@(main+50d) → resume, plus truncated
      file must fail; compares output(.bin) AND ctsout bytes.
      Wire = magic OMCk/v1, phase, clock(t0,t1,interv,end,estEnd), statics,
      monitoring(RuntimeState+cts offset), transmission(v2-superset field
      list; estimator/MPD/NonVector), humans(SoA cols+infIncidence+WH all
      kinds incl. Vivax broods+Molineaux variant counts+Pathogenesis+
      LSTMDrug/LSTM queues+Clinical/Episode), master RNG last. Restores reuse
      real ctors (createInfection / infection_new / emplace) whose draws are
      then overwritten from the wire. Budget-checked reads: every corrupt
      case = clean rc73 error, never UB. Also fixed: kill leftovers
      Genotypes.cpp assert `state.current_mode` (only breaks assert builds)
      and Infection-descriptive restore (factory unset → direct infection_new).
- [x] C4 corrupt-checkpoint suite, proven per scenario in resume.py (21
      checks/scenario): truncated, junk, bad magic, bad format version
      (byte bumped), wrong scenario (store count mismatch) — all clean
      rc73, no output file, no UB (budget-checked reads). Old v2
      boost/gzstream files hit the magic check.
- [x] C2 scripted save days per scenario in resume.py: main-start,
      mid-main (mainStart+50d, chained), late-main (end-25d, chained),
      mid-warmup (mainStart/2). Chains exercise resume→save→resume.
- [x] C3 (no-intervention corpus): Resume_* now 12 — + scenario1 (NonVector
      EIR branch), VecTest, Genotypes (SAMPLE_TRACKING mode), ExactNv0
      (exact solver — NOTE: seed-sensitive; resume.py uses the scenario's
      own iseed, NOT --seed 112), SimpleMPDTest (MPD state),
      HetVecDailyEIR, MolPairwise (multi-clone variants).
      Interventions-heavy = still the documented gap (skill law).
- [x] C5 final wrap-up: suite green (ctest 64/64, unittest .OK!, tier-a
      53/53), LOC report: raw model LOC 29,275 (kill commit) → 30,438 →
      30,426 at merge with rewrite (checkpoint module +1,163 incl. licence
      headers); PERF.md NOT touched (main-session ledger).
- Remaining (documented gap, not scheduled): interventions checkpoint/replay
  — see Known limits.

## Proven round-trips (update every green commit)
C2/C3 (all byte-identical vs uninterrupted, output.bin AND ctsout.txt,
scenario default seed): 4 save-day chains each (main-start, mid-main,
late-main end-25d, mid-warmup) + truncated-file rejection for:
- NoInterv (Descriptive WH + cts ~1.5k lines), Vivax (broods), Empirical,
  Molineaux (growing variant counts), Penny
- 1 (NonVector EIR), VecTest, Genotypes (tracking mode), ExactNv0 (exact
  emergence estimator), SimpleMPDTest (SimpleMPD state), HetVecDailyEIR,
  MolPairwise
- corrupt (per scenario, in resume.py): truncated, junk, bad magic, bad
  version byte, wrong scenario → all rc73, no output file, no UB.
Gates at C1: build ✓, ctest 57/57, unittest .OK!, tier-a 53/53.
Gates at C2/C3: ctest 64/64, unittest .OK! (model binary byte-identical to
the C1 parity-verified one).
Final gates (HEAD): build ✓, ctest 64/64 (12 Resume_* + 52 legacy),
unittest .OK!, tier-a 53/53 bit-identical vs A0.

## Known limits (honest partial)
- NO interventions section (skill law 2026-11-20): --resume on a scenario
  that deploys anything is UNSUPPORTED. Deployment schedules
  (Interventions::Deployments / InterventionManager) and global IvState are
  not on the wire and deployments are not replayed after resume, so state
  created by a deployment before the save point is lost. Exception:
  per-infection drug dose records (TreatInfectionsWithDrug / LSTMDrug
  queues) ARE on the wire with the infection. Proven days are all
  intervention-free scenarios.
- Intervention-derived side state not on the wire (same gap): human cohort
  intervention columns, AnophelesModel SimpleDecayingValue sugar-feeding
  values, vaccine/vector control registries. All config-rebuilt = empty.
- cts resume overwrites file content past the saved offset (rows the later
  part of a previous longer run wrote are replaced; fine for save+exit use).

---

# P8-structure HANDOFF — killed model/core/, relocated SlotAllocator, dissolved RunExtensions (branch `p8-structure` off `rewrite` HEAD bac09107)

User order (verbatim): "kill the core folder, move outside or in other files if
that makes sens. The createslot file is also weird, where does it really
belong? same for runextension... delate"

## Environment
- Gate build: `build-st/` = Release `-O3 -DNDEBUG -ffp-contract=off`
  `-DOM_BOXTEST_BINARY=ON`. `cmake --build build-st -j16`.
- Trap (cost a detour): parity.py hardcodes `REPO/build/schema` for the XSD —
  the worktree needs the ignored symlink `build -> build-st` (`build*/` is
  gitignored; already created here).
- unittest: run from `build-st/unittest/` (`./unittest` → ".OK!", 94 cases).
- ctest: `cd build-st && ctest -j12` (52).
- T1: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new
  build-st/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a
  -j12` (53); `--tiers b --seeds 1-10` (10) → final 63/63.
- NO bench this phase (structural, perf-neutral; nothing smelled slow —
  ctor now allocates 9 small heap objects up-front, free at this scale).
- HARD boundary honoured: never touched /tmp/opencode/p8ckpt or
  /home/razlock/git/fork's checkout. Commits unsigned (`-c
  commit.gpgsign=false`): gpg pinentry can't run headless.

## Done — all gates green per step (build 0 err, ctest 52/52, unittest .OK!
94, T1 tier-a 53/53; steps 2–3 also tier-b seeds 1–4; final a+b seeds 1–10
= 63/63 on the final binary)

| Step | model code LOC / files | Notes |
|---|---|---|
| baseline (rewrite HEAD) | 18,248 / 139 | |
| S1 kill model/core/ | 18,248 / 139 | RunContext.h, run.{h,cpp}, RunExtensions.h → model/ top level (siblings of Clock.h — where run state belongs); sed include/comment `core/` prefixes; guards Hmod_core_* → Hmod_*; GLOB build, zero CMakeLists edits |
| S2 SlotAllocator | 18,248 / 139 | model/Host/ → model/Host/WithinHost/ — real owner is WHInterface (one allocator per host; Common/Descriptive drive it, InfectionIncidenceModel uses clamp_count) |
| S3 dissolve RunExtensions | 18,257 / 138 | struct + ext indirection + ensure_extensions() deleted; every field a plain RunContext member; `ctx.ext->X` → `ctx.X` (173 sites, 32 files); RunExtensions.h file deleted |

Final: **18,257 / 138 (+9 LOC, −1 file)**. Honest note: net is +9, not ≤0.
Breakdown: RunContext.h 33→99 code LOC (field decls + 9 pointer members + 8
std includes can't shrink — the data had to move somewhere), run.cpp
181→215 (ctor/dtor/Treatments::select, +34), RunExtensions.h −70, −29 stray
include/boilerplate elsewhere. The +34 in run.cpp is the gcc-16.2
deleter-on-incomplete-type tax; if a later phase un-cycles one of the heavy
headers, its field becomes a value member and the tax shrinks.

## Where everything ended up and why
- RunContext.h, run.h, run.cpp → `model/` top level: they ARE the run core,
  siblings of Clock.h; "core/" was a folder of one-off files.
- RunExtensions.h → deleted; fields dissolved into RunContext.
- SlotAllocator.h → `model/Host/WithinHost/`: WHInterface owns one per host.
- `WithinHost::Treatments::select()` → defined in run.cpp (needed complete
  RunContext + Treatments; declaration stays in Treatments.h).

## S3 shape (read before merging p8-checkpoint)
- RunContext.h gains only 2 non-leaf-risky includes: Diagnostic.h,
  ComponentId.h (both proven cycle-free leaves). Everything else of the old
  ext stays behind forward declarations.
- 9 raw-pointer members (new[]ed in the out-of-line ctor **in old
  RunExtensions declaration order** = old ensure_extensions order = law):
  decision_trees, whTreatments, humanComponents, continuous, timed,
  importedInfections, removeAtIds (`new vector<ComponentId>[SubPopRemove::NUM]`,
  index syntax unchanged), anophParams, monitoring. ~RunContext deletes them
  reverse-order (mirrors `delete ext`). All other ex-fields are value members
  (diagnostic_set, identifierMap, function pointers, whMock, itn/irs/gvi/
  vaccineComponents, vaccineReportComponent, entoAvailabilityPercentiles,
  cts* output state, nextTimed).
- WHMockHooks is now `RunContext::WHMockHooks` (nested); unittest installers
  (UnittestUtil.h initTime, WHMock.cpp) write ctx() fields directly —
  ensure_extensions() call removed.
- Behaviour-neutral fix: `mon::Continuous::update(const RunContext&)` wrote
  `ctsStreamOffset`/the stream — the ext raw pointer made that legal through a
  const ref (a const-lie). Signature is now `RunContext&`. T1 unaffected.

## Merge traps vs p8-checkpoint (branch in /tmp/opencode/p8ckpt; resolves in main session)
p8-structure deletes files p8ckpt edits → delete/modify + contextual conflicts:
1. `model/core/RunExtensions.h` (their base): DELETED here. Their removal of
   the factory pointers `createInfection`/`checkpointedInfection` must be
   re-applied to `model/RunContext.h` (member block, S3-tree lines ~119–122)
   AND to the RunContext ctor/dtor in `model/run.cpp` (~156–186) AND to the
   4 Infection/*.cpp init sites (`ctx.createInfection = …;
   ctx.checkpointedInfection = …;` — DummyInfection.cpp:41-42,
   EmpiricalInfection.cpp:80-81, MolineauxInfection.cpp:164-165,
   PennyInfection.cpp:157-158) + the restore read site
   CommonWithinHost.cpp:369 (`ctx.checkpointedInfection(ctx, stream)`).
   If restore dies, `checkpointedInfection` (and maybe the whole field)
   goes; `createInfection` is also read by CommonWithinHost.cpp:118/147/152/
   168/173 + 2 unittest suites — keep it unless their diff says otherwise.
2. `model/run.cpp`: both branches rewrite its top half. Their target: the
   checkpoint file logic in run_simulation (base ~lines 220–300:
   checkpointFileName / readCheckpoint / writeCheckpoint calls +
   `#include "checkpoint.h"` line 45). My block: 7 includes replacing
   RunExtensions.h at line 46 (immediately adjacent — expect a conflict at
   the include list) and ctor/dtor/select inserted right before
   run_simulation. Their edits inside run_simulation body should apply
   cleanly after the sed'd `ctx.ext->` lines are accounted for.
3. sed'd checkpoint-touch lines p8ckpt patches will miss contextually:
   InterventionManager.cpp:41/45 `(*ctx.importedInfections) & stream` (was
   `ctx.ext->importedInfections & stream` — the intervention checkpoint
   wire; if their branch deletes intervention save/restore, drop these lines
   entirely), InterventionManager.cpp:67, xml/interventions.cpp:1322,
   PerHost.cpp:131/135/139, monitoring runtime sites (Monitoring.cpp ~20×
   `auto& runtime = *ctx.monitoring;` was `ctx.ext->monitoring`).
4. `unittest/UnittestUtil.h`: ensure_extensions(ctx()) line removed here;
   any of their fixtures still calling it (their tree) must drop it at merge.
5. If p8ckpt deletes checkpoint save/restore helpers from checkpoint.cpp,
   nothing here touches that file — safe.

## Left open / notes
- `.opencode/skills/openmalaria-rewrite/SKILL.md` target-layout still lists
  `model/core/` — contract-doc update is the main session's call.
- PERF.md untouched (main-session ledger). Report the +9 here.

**8-Q1 dead-weight note (from G16, NOT touched here):** in the engine block now in util/random.h — commented-out `Sseq` ctor/seed + `operator<<`/`>>` blocks, `generate_block()`/`chacha_core()` declarations with no definitions anywhere, and `rotl`'s file-static placement; plus random.h's `<chacha.h>` include + `MasterRng` typedef are only referenced by checkpoint RNG export — prune candidates with proof there.

**Class-design-pass anchors (queued after p8-files; definitions were moved between files here, designs untouched — these are the cohabited sections to split/rework):**
- `model/Host/WithinHost/Infection/Infection.cpp` — Kind-dispatch ctor/update + the Dummy kind section (was DummyInfection.cpp, G7) in one file.
- `model/PkPd/LSTMDrug.cpp` — all four Kind bodies (OneComp/ThreeComp/Conversion sections marked `// ———`) + shared ctor (G8); `updateCached`/GSL workspace file-scope state lives here.
- `model/Clinical/CM5DayCommon.cpp` — DecisionTree5Day + ImmediateOutcomes + shared CM5Day helpers (G3 5-day chain merge); `i`/`l` twin code blocks per kind live in the WH files.
- `model/Host/WithinHost/WHInterface.h` — VivaxBrood + SlotAllocator (G15) + WHInterface now share one header; brood/allocator design review happens there.
- `model/Transmission/TransmissionModel.h` — flat Kind-switched class; vector/non-vector halves are comment-marked.

**Rulings on Human/HostColumns/Population (user-approved, final — do not re-litigate):**
- Human→Population merge: REJECTED. Human is the SoA row-handle facade (typed accessors over HumanColumns + row-object identity, doB/rng/infIncidence/withinHostModel/clinicalModel members), not a data blob. Include graph: Human.h is a widely-used leaf read by ~26 TUs; Population.h is an 11-TU consumer that includes Human.h — merging Human into Population would bloat every Human reader with population/demography internals for zero LOC gain. The SoA question doesn't remove the handle: callers name per-human behaviour (ageYears, rng, infect…), not column indices.
- HostColumns.h → owning container's header (Population.h holds the only `HumanColumns cols` member): REJECTED by include-cycle, verified: HostColumns.h is included by Human.h, WHInterface.h, PerHost.h, HumanComponents.h, IvState.h, Vaccine.h + Population.h; Population.h includes Human.h, so the merged header would require Human.h -> Population.h -> Human.h. The vocabulary a row-handle and column-readers below Population consume cannot live in the container header above them. HostColumns.h stays the deliberate leaf it is (3 records, ~116 LOC, trivial includes) — matching the ComponentId.h leaf rule.

---

# P8-poly HANDOFF — Infection + WHInterface Kind-dispatch -> real class hierarchies (branch `p8-poly` off `rewrite` HEAD = p8-files green 11c0f585)

User order: convert the two model-option variant families to true classes.
Laws honoured: formula bodies VERBATIM (user veto on reassociation; only
`this->`/`inf.` re-qualifiers and payload-struct absorption); RNG draw points
AND order unchanged (draws moved into derived ctors at the exact old stream
slots; base ctors draw nothing); checkpoint module stays free-function — kind
tag + `static_cast`, no serialization decls in classes; resume ctor order kept
base -> WHFalciparum -> concrete; ALL variants kept (deletion vetoed);
CalcNvO untouched; mock hooks re-homed via installer pattern.

## Environment
- Worktree `/tmp/opencode/p8poly`, build `build-poly/` = Release `-O3 -DNDEBUG
  -ffp-contract=off -DOM_BOXTEST_BINARY=ON`; ignored symlink `build ->
  build-poly` REQUIRED (parity.py XSD trap). Commits `--no-gpg-sign`.
- Gates per commit: build 0 err; `cd build-poly && ctest -j12` (64);
  `cd build-poly/unittest && ./unittest` (".OK!"); T1
  `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new
  build-poly/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a`
  (53) + tier-b seeds 1-4 on dispatch commits; final a + b seeds 1-10 = 63/63.

## Commits
| commit | scope | model code LOC / files |
|---|---|---|
| baseline 11c0f585 | | 17,192 / 115 |
| `12461d2c` | Infection: base (update/expired/clear virtuals, `const Kind kind`, factory) + Dummy/Descriptive/Empirical/Molineaux/Penny classes; payload structs absorbed by value (per-infection mallocs die); per-kind TUs; checkpoint static_cast switch; unittest suites re-homed | 17,192 / 115 |
| `c9732f64` | WHInterface: abstract base (common wire block + virtual treatment()) -> WHFalciparum (infections list/pkpd) -> WHCommon/WHDescriptive; WHVivax : base; all Kind-switches deleted; checkpoint per-kind static_casts; unittest WHMock became a class (hooks unchanged: installer still re-homes ctx.whMock + the hook table; the mock's non-hook methods keep the old `unimplemented_exception("not needed in unit test")` branches) | 17,078 / 115 |

## Design notes (read before extending)
- Base classes carry ONLY what the checkpoint wire writes for every kind;
  kind-specific state lives in the derived class that reads it. Checkpoint
  `static_cast` is justified by the stored kind tag (same rule that made the
  old union/`inf.X` accesses safe).
- WH ctor layering = wire order: base ctor (no draws) -> WHFalciparum ctor
  (falciparum immunity draws) -> concrete ctor (own draws). Do not reorder.
- The model side has ZERO mock code: `WHInterface::init` (scenario.cpp)
  selects production variants only. The unittest mock is `WHMock final :
  WHInterface` in unittest/WHMock.cpp; makeWHMock installs the hook table on
  ctx.whMock and returns the mock — treatment()/treatSimple()/treatPkPd()/
  getTotalDensity()/diagnosticResult() go through the hooks, every other
  virtual throws unimplemented_exception (old Mock branches verbatim).
  `treatment()` is virtual on the base for exactly this (prod body shared by
  Falciparum+Vivax; Mock overrides with the hook).
- InfectionImmunitySuite uses `WHDescriptive*` directly (setParams ->
  WHFalciparum::setParams); SlotAllocatorSuite constructs WHDescriptive.

## PERF (measured, closed)
BENCH (taskset -c 5, %U, example_scenario@popSize2000, warmup discarded,
10 INTERLEAVED reps, same window both arms): p7g median **12.21 s** (MAD
0.03) vs p8-poly median **12.36 s** (MAD 0.02) = **+1.3%** — inside the
≤ +1.5% expectation. An interim +27% reading was a measurement artifact
(arms measured in different windows during parallel -j12 parity load); the
static-dispatch attribution experiment was void for the same reason. Full
table + retraction: BENCH_p8poly.md.

## Second wave (LSTMDrug, ClinicalModel): NOT started — deliberate stop
Both are still flat Kind classes (headers unchanged since the pragma sweep).
Given the open +27% perf question above, converting two more families before
the orchestrator rules on the perf finding would compound the same TU-split
effect for unknown gain. Pattern is proven twice (12461d2c, c9732f64): base =
wire/common state, derived = kind payloads, ctor layering = old stream order,
checkpoint static_casts on the kind tag. Anchors for whoever picks it up:
model/PkPd/LSTMDrug.{h,cpp} (freeLSTMDrug* = the dtor bodies),
model/Clinical/ClinicalModel.{h,cpp} (Kind ctor + treatment-episode guard).

---

# P8-q3 HANDOFF — LSTMDrug + ClinicalModel Kind-dispatch -> class hierarchies (branch `p8-q3` off `rewrite` HEAD 1b8661a0)

User order behind the wave: model-option variants as clearly independent
classes, per-approach functions; shared data structure on the base where
genuine. p8-poly pattern applied; laws honoured: formula bodies VERBATIM (FP
veto), RNG draw points+order unchanged (base ctors draw nothing; derived ctor
bodies = old per-Kind case bodies verbatim), checkpoint STAYS free-function
self-contained (kind tag + `static_cast`, zero serialisation decls in the
classes), ALL variants kept, error strings verbatim, mock hooks untouched
(nothing converted is mocked — WHMock covers WHInterface only).

## Commits
| commit | scope | model code LOC / files |
|---|---|---|
| baseline 1b8661a0 | | 16,759 / 115 |
| `f9f7bd44` | LSTMDrug -> base + `LSTMDrugOneComp`/`LSTMDrugThreeComp`/`LSTMDrugConversion` (same TU; ~59-TU cone keeps bodies in .cpp) | 16,759 / 115 (delta 0) |
| `dcd89ea5` | ClinicalModel -> base + `CMDecisionTree5Day`/`CMImmediateOutcomes`/`CMEventScheduler` | 16,805 / 115 (+46) |

## LSTMDrug shape
- Base holds ONLY what every variant shares: `const Kind kind`, `typeData`,
  `doses` (+ `medicate`), `vol_dist`; non-virtual `getIndex`; virtual
  getConcentration/calculateDrugFactor/updateConcentration + virtual dtor.
  The old `freeLSTMDrug*` dtor dispatch was already retired in p8-files
  (`vector<unique_ptr<LSTMDrug>>` owns destruction); the virtual dtor is its
  class-shaped replacement.
- Old switch-body -> method map: ctor case X = derived X ctor body;
  getConcentration case = override body (base dead `return 0.0` fell away);
  `drugFactor_X`/`updateConcentration_X` = `calculateDrugFactor`/
  `updateConcentration` overrides verbatim; per-kind helpers (updateCached,
  calculateFactor_fC / parent+metabolite getters, calculateFactor_conv,
  set*Parameters) private in the owning class; `metaboliteData` moved to
  LSTMDrugConversion; GSL workspaces + func_* integrands unchanged file-scope.
- Ctor/stream proof: each derived ctor body is the old case body line-for-
  line (OneComp: Vd→elim; ThreeComp: Vd→elim→k12→k21→k13→k31→ka; Conversion:
  Vd→ka→elim→conv→metab elim→metab Vd), base draws nothing ⇒ stream
  bit-identical (PkPdCompliance/LSTMPkPd suites + T1 63/63 confirm).
- Factory: `DrugRegistry::createInstance` — one create-line per variant,
  selection condition unchanged (conversion_rate / k12).
- Checkpoint (`checkpoint/humans.cpp`): kind tag written from a local copy
  (`kind` is const now); read tag discarded — the factory rebuilds the same
  kind from the same type config; per-kind field order/names identical.

## ClinicalModel inventory (model/Clinical/)
Converted: the one instance family — Kind { DecisionTree5Day,
ImmediateOutcomes, EventScheduler }, selected by `opt_event_scheduler`/
`opt_imm_outcomes` (healthSystem XML). All four dispatch points became
virtuals: update() tail -> `doClinicalUpdate()`; the two
`kind==ImmediateOutcomes` branches inside `clinicalUpdate5Day` ->
`uncomplicatedEvent()` (public: shared free function dispatches through it);
`isExistingCase()` ES-clock pick -> `timeOfLastTreatment()` (ES override).
ES payload (esPgState + es_* column accessors) -> CMEventScheduler; ES ctor
guard (verbatim error) -> CMEventScheduler ctor; `esClinicalUpdate` now takes
`CMEventScheduler&` (body byte-identical). Override bodies are one-liners
onto the untouched free-function bodies, in ClinicalModel.cpp (the TU
resolving the include cycle; same-file rule).
Ruled out (already clearly-identifiable / not instance variants):
- `DTKind` switch in CMDecisionTree.h — compiled decision-tree node
  interpreter (phase-4a frozen design), node kinds ≠ model-option variants.
- CM5DayCommon / DecisionTree5Day / ImmediateOutcomes / ClinicalEventScheduler
  / ESCaseManagement — per-approach STATIC config classes, each with its own
  init/setHealthSystem already; no per-instance kind state.
- ClinicalModel::init/setHS static selectors — scenario config loading (which
  XML section feeds which family), the factory input itself; not behaviour
  dispatch.

## Wire notes (checkpoint bytes unchanged both commits; wire may break freely
per contract, kept identical anyway)
- LSTMDrug: tag+doses+vol_dist+per-kind fields byte-identical.
- ClinicalModel: esPgState slot symmetric — ES writes/reads its payload via
  static_cast; 5-day variants write/read the old always-NONE int via a local.

## Gates (per commit + final)
build 0 err; ctest 64/64 incl Resume_* 12/12; unittest .OK!; T1 tier-a 53/53;
tier-b seeds 1-4 after each (both commits touch factory/ctor); FINAL tier-a +
tier-b seeds 1-10 = 63/63 identical. No bench (structural; storage and
per-day dispatch counts unchanged by construction — ClinicalModel update()
now takes ONE virtual call where the free functions took the same work;
LSTMDrug per-drug virtuals replace switch jumps, same call count).

## Left open / notes
- PERF.md untouched (main-session ledger): honest LOC delta +46 total.
- `unittest/` needed zero re-homing (no suite constructed LSTMDrug or
  ClinicalModel instances directly; only statics/factories, unchanged).
- Parallel-agent landmine respected: random.h untouched; no RNG call site
  moved except the verbatim ctor-case relocation inside these two families.

---

# p8-xml: model/xml/ compaction wave (branch p8-xml, base f7385b04)

Domain: model/xml/** only. Model LOC 16,799 → 16,764 (−35), 7 commits, all
gates green per commit (build 0 err, ctest 64/64 incl Resume_*, unittest .OK!,
T1 tier-a 53/53; commits touching interventions/pharmacology/healthSystem
parse paths also tier-b seeds 1-4; FINAL tier-a + tier-b seeds 1-10 = 63/63).

## Collapses (per commit)
- pharmacology.cpp: mirror locus-restriction error throws → one `mismatch`
  lambda; message composed from same fragments, name-pair swapped per site;
  erase()-before-throw order and loop positions unchanged (−6).
- interventions.cpp: two literal-only `ostringstream msg;` blocks → direct
  throw with adjacent-literal split (bytes identical by concatenation
  identity) (−2).
- interventions.cpp: 9× `UnitParse::readDate(..., STEPS /*backwards compat*/)`
  → `readDate_steps` wrapper; no error strings touched, args/order unchanged
  (−6).
- HumanInterventionComponents.h: print_component 12-deep ternary →
  Component::Type-indexed 13-entry name table. Enum is dense (RECRUIT_ONLY…
  CLEAR_IMMUNITY); IPT entry carries the former default fall-through bytes
  "clear immunity"; kind is only ever assigned enum members (Components.h
  set<T> + CLEAR_IMMUNITY site), so indexing is total (−10).
- interventions.cpp nhh_mods_setup: 5× resize-if ladder → loop over the five
  same-type effect maps; loop order == textual order, operator[] insertion
  per map unchanged (−4).
- healthSystem.cpp setHS: 3× presence-check+dispatch blocks → one
  `need(optional, name)` lambda; message composed byte-identically
  ("Expected " + name + " section in healthSystem data (…)"); branch order
  and check positions unchanged (−7).

## Rejected with reasons
- entomology.cpp fourier/monthly expIDFT-rescale tail (2×3 stmts): measured —
  helper costs +2 net LOC. Not done.
- util.cpp OptionCodeMap 40 assignments → braced init: same line count, pure
  cosmetics.
- healthSystem getHealthSystemACRByName/getHealthSystemTreatmentOptional drug
  ladders → member-fn-ptr table: table+lookup is as long as the ladder, and
  `&Class::getX` is ambiguous (const/non-const overloads).
- ITN/IRS initialInsecticide+qnorm3 3-statement dup → template helper: net
  ≈0.
- util.cpp DecayFunction kind ladder: if-chain is already minimal.
- interventions.cpp VaccineParams phenotype block: messages all distinct,
  loop already minimal.
- decision-tree compile_* two-branch node funcs (diagnostic/uncomplicated/
  severe/cohort): each has a distinct middle step; shared helper would be
  longer and the zone is order-critical (gcc arg-eval parity comment).
- scenario.cpp PRINT_GENOTYPES header loop (carried from recon): ~4 LOC,
  output-shape risk.

## Dead-code sweep
grep-zero found none: averageEIR (4 users), BetaSampler/WeibullSampler
setParams (entomology 500-502, vivax 104), all xml helpers referenced.
DocumentLoader.h carries no stale diagnostics (XMLChecker+latentP only).

## [string-audit] big-fish observations (NOT touched, main-session material)
- interventions.cpp factors::SurvivalFactor::init: the isDeterrent vs else
  arms hold ~15 near-identical comment lines each plus two bounds-appender
  msg blocks (6 checks / 3 checks sharing literal fragments). A shared
  table would re-derive message order → out of scope here.
- healthSystem drug readers are a trio (ACR / TreatmentOptional /
  TreatmentByName) over the same 6-drug name space; one registry would
  replace all three but is the mini-parser-table class explicitly excluded.

## Landmines respected
random.h, model/Host/WithinHost, model/PkPd, model/Clinical, XSDs, generated
scnXml untouched; build symlink uncommitted (gitignored); PERF.md/plan
untouched; no RNG draw site moved (all folds are pure reads/message plumbing).
