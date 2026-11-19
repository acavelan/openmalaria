# Phase 5 HANDOFF — cycle breaks + Infection perf fix (branch `p5-cycles` off `p4b-green`)

Crash-recovery file. Update EVERY green commit. Baseline p4b-green: model **18,599** / unittest 3,760; A0 bench 12.33s; p4b bench 12.80s (session re-measure `p5-before` = **12.77s** mad 0.000). Full phase-4b handoff: `git show p4b-green:HANDOFF.md`.

## Environment
- Build: `build-p5/` = Release `-ffp-contract=off` `-DOM_BOXTEST_BINARY=ON`.
- unittest MUST run from `build-p5/unittest/` → `./unittest` (91 tests). ctest: `cd build-p5 && ctest -j16` (52).
- T1 tier-a: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new build-p5/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a -j16` (53). Tier-b: `--tiers b --seeds 1-10` (10).
- bench: `sh /tmp/opencode/om-bench/bench.sh /home/razlock/git/fork/build-p5/openMalaria TAG` (ABSOLUTE path!).
- perf: hardware events BLOCKED (perf_event_paranoid=2) — `perf record -e task-clock -F 3997` only (no -g). p4b profile: /tmp/opencode/p5-perf/p5-sw.data. Hot: dgemm 18% (CalcNvO GSL — phase-8 flip kills), getPEpisode 7.6%, probTransmissionToMosquito_falciparum 6.5%, update_descriptive+determineDensities+decayInfection ~9% (Infection loop).

## Steps
1. [x] Infection payload — commit a3e1352a.
2. [x] Cycle break + dynamic_cast removal — commit 8d517521.
   `model/Transmission/VectorIntervention.h` (TransmissionModel + VectorControlIntervention
   bases moved verbatim + `asVectorIntervention`/`asITNComponent`/`asIRSComponent` virtuals),
   `model/interventions/ComponentId.h` (`Component::Type`+`ComponentId`, leaf: Global.h+limits).
   PerHost.h pulls ONLY the leaf (verified `g++ -H`). `DeployedVectorIntervention(const void*)`
   → `asVectorIntervention()`. Also fixed: `~PerHost()` + move ops declared OUT-OF-LINE in
   PerHost.cpp (user-declared dtor suppressed implicit moves — move was silently copying all
   per-host vectors! caught by -Wall on `std::move(humans.begin(), ...)` in Population.cpp);
   unused DecayFunction include (cycle member) in PerHost.h → sampler.h; dup include
   TransmissionModel.h in Deployments.h dropped.
3. [x] Clinical seam: only clinical→interventions call = `t.deploys[aux].deploy(...)` in
   Clinical/CMDecisionTree.cpp on a concrete owned value (`FlatTree::deploys`) — legal leaf
   use, documented with seam comment. Removed unused `Clinical/ClinicalModel.h` include from
   interventions/Deployments.h. interventions→clinical runtime edges (DecisionTreeComponent
   running `Clinical::exec_decision_tree`; changeHS → `ClinicalModel::setHS` statics) live in
   model/xml wiring layer, NOT interventions/ — phase-7 `ctx.ext` cleanup, documented here.
   Box tests: T2 conflict never materialized (checkpoint wire kept byte-compatible), so the
   sanctioned 49-scenario single-run re-bless was NOT needed; 12/12a keep real checkpoint coverage.
4. [ ] FINAL GATE:: T1 a+b 63/63, ctest, unittest, bench p5-final, cloc, census statics
   refresh (nothing moved so far — no static changes), PERF.md row 5.

## Learned / traps (inherited + new)
- m_rng(ctx.master_RNG) CONSTRUCTS-BY-DRAWING — only where old code did (4b h10 trap).
- ext lifetime only after ensure_extensions; no header-inline deref of ctx.ext; TwoRunsSuite = ext tripwire, rerun whenever RunExtensions changes.
- INHERITED-SCOPE lookup: qualify statics per body after merges; convert from `git show HEAD:file` anchors.
- Error strings byte-identical; `// parity:` comments travel with code.
- `grep -v LSTMDrug` hides `inf->Kn` uses in LSTMDrugType.cpp/LSTMDrugConversion.cpp — check ALL of model/ before deleting "dead" members (Kn kept, lazy unique_ptr).
- Checkpoint payload classes behind unique_ptr: keep read_X/write_X byte ORDER identical — wire stays compatible, no golden re-bless (this phase's alternative to box-test conversion).
- New per-host payload classes: if the class has a user-declared destructor, declare move ctor/assignment explicitly (= default in .cpp) or std::move silently deep-copies.
- ESCaseManagementSuite.h + AnophelesModelSuite.h remain excluded+stale (pre-existing); deferred: setHS-overwritten clinical statics → ctx.ext (phase 7); Human.h still includes interventions/HumanComponents.h (phase-7 territory, not the PerHost cycle).
- Untracked PROMPT.md/email = user files, leave.

## Model-only census statics (refresh if statics moved)
293 rows: ConstValue 23, TimeSpec 21, TrialVariant 20, ScenarioReader 16, DecisionTreeNode 16,
DecisionTree 16, GeneticMutation 15, GeneticMutationType 12, ClinicalDeploy 11, Time_evolution 11,
EnumMapping 10, GenericTime 10, PkPdD 10, SimpleTime 10, GammaVal 9, HumanMarking 9,
AntimalarialModel 8, ClinicalMethod 8, Vaccine 8, GenericMarking 7 (+ long tail). Top-4:
ScenarioReader 143, ClinicalMethod 108, DecisionTreeNode 98, DecisionTree 83. Unchanged since
4b-final (no statics moved in 5 yet).

## Gates log
- (step1 a3e1352a) ctest 52/52; unittest 91/91; T1 a 53/53; bench p5-after-infection-fix mean **12.48s mad .010** vs 12.77 → −2.3% KEEP (floor 1.0%); RSS 39.4→35.4MB. model LOC 18,614 (+15 payload decls).
- (step2+3 8d517521) build clean; unittest 91/91; ctest 52/52; T1 a 53/53. model LOC 18,658.
- (step3-4 44947048) build clean
- (FINAL GATE) T1 63/63; ctest 52/52; unittest 91/91; bench p5-final 12.51s mad .010 KEEP; LOC 18,657. (also fixed stray-quote warning CM5DayCommon.cpp:25 from 4b); unittest 91/91; ctest 52/52; T1 a 53/53; T1 b seeds 1-4 4/4.
