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
