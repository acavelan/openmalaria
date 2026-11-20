# P8-structure HANDOFF — kill model/core/, relocate SlotAllocator, dissolve RunExtensions (branch `p8-structure` off `rewrite` HEAD bac09107)

User order (verbatim): "kill the core folder, move outside or in other files if
that makes sens. The createslot file is also weird, where does it really
belong? same for runextension... delate"

## Environment
- Gate build: `build-st/` = Release `-O3 -DNDEBUG -ffp-contract=off`
  `-DOM_BOXTEST_BINARY=ON`. `cmake --build build-st -j16`.
- unittest: run from `build-st/unittest/` (`./unittest` → ".OK!", 94 cases).
- ctest: `cd build-st && ctest -j12` (52).
- T1: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new
  build-st/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a
  -j12` (53); steps 2–3 touch WithinHost → ALSO `--tiers b --seeds 1-4` (57);
  final `--tiers a` + `--tiers b --seeds 1-10` = 63.
- NO bench this phase (structural, perf-neutral-required; machine busy).
- LOC baseline: model = 18,248 code LOC / 139 files. Do NOT touch PERF.md or
  the plan (main-session ledger). Prior p7g handoff: commit e34997f1 history.
- HARD boundary: never touch /tmp/opencode/p8ckpt or /home/razlock/git/fork.

## Steps
- [ ] S1 kill model/core/: RunContext.h, run.h, run.cpp, RunExtensions.h →
      model/ top level (siblings of Clock.h — that is where they belong).
      sed `#include "core/X.h"` → `"X.h"` across model+unittest. GLOB build →
      no CMakeLists edits.
- [ ] S2 SlotAllocator.h: model/Host/ → model/Host/WithinHost/ (WHInterface
      owns one per host; Common/Descriptive drive it).
- [ ] S3 DELETE RunExtensions.h: every field becomes a plain RunContext
      member; `ctx.ext->X` → `ctx.X` (~170 sites, 34 files); kill
      ensure_extensions(), the ext pointer, ~RunContext's delete; delete file.
      TRAPS: RunContext.h must stay lean — forward-declare member types;
      where gcc 16.2 instantiates deleters at class scope for incomplete
      types, use raw-pointer members deleted individually in the out-of-line
      ~RunContext in run.cpp. Construction ORDER of ensure_extensions() is
      law (registries + unittest backend installers; RNG streams — RngGoldens
      /TwoRunsSuite are tripwires). unittest fixtures (WHMock.cpp,
      UnittestUtil.h) install backends into these fields — update call sites.

## Merge traps vs p8-checkpoint (parallel branch, resolves in main session)
(TO BE FILLED IN STEP 3 with the exact field/line overlap.)
