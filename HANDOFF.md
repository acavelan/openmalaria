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
