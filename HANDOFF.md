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
