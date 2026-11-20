# Phase 7g HANDOFF — whole-tree sweep: folds, dead code, accessors (branch `p7g-sweep` off `rewrite` HEAD 1da48376)

User order (verbatim): "Remove empty files, reduce number of files where
reasonable, kill dead code. Inline more functions. Don't hesitate to make small
cpp header-only. Clean the folders mess, why do we have both humans and Host
for example? Plenty of opportunity to reduce code. Make the code easy to
understand for humans. Aim at having simple structs and simple free functions.
You are allowed simple classes. getters and setters are not needed most of the
time. Check the full code and remove code."

## Environment
- Gate build: `build-p7g/` = Release `-O3 -DNDEBUG -ffp-contract=off`
  `-DOM_BOXTEST_BINARY=ON`. `cmake --build build-p7g -j16` (~60 s clean).
- Warnings build (report-only, throwaway): `build-warn/` = same + `-O1
  -Wall -Wextra -Wunused-function`; log at /tmp/opencode/warn.log.
- unittest: run from `build-p7g/unittest/` (`./unittest` → ".OK!", 94 cases).
- ctest: `cd build-p7g && ctest -j12` (52).
- T1: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new
  build-p7g/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a
  -j12` (53); after Infection/Transmission/clinical-touching commits also
  `--tiers b --seeds 1-4` (4); final `--tiers a` + `--tiers b --seeds 1-10` = 63.
- Bench (final only): from `/tmp/opencode/om-bench` (XSD symlink dir),
  taskset -c 5, median-of-3 min, absolute paths, example_scenario.xml@2000.
  p7f-era ≈12.2 s; >1 % off = investigate; ~13.9 s = external load, discard.

## Done already (pre-7g, do NOT redo)
No empty files; small cpps collapsed except xml/demography.cpp (KEEP — schema
drag); zero friends; no Global.h (checkpoint.h/mod.h carry OM-scope usings,
order-sensitive); Clock.h owns SimTime+sim+SimClock; Infection = struct+free
functions (7f).

## Steps
- [x] C1 folder folds: humans/HostColumns.h + infections/SlotAllocator.h →
      model/Host/ (dirs deleted); mon/init.h → mon/Monitoring.h. 18332→18320.
- [x] C2 dead-code round 1 (warnings-driven + grep): cts shells ×3, #if 0
      Vaccine block, ClinicalEventScheduler 5 dead instance members, dead enum
      labels ×2, vectors::addTo, dead ctor params ×5, dead set-but-unused ×2,
      -Wreorder ×2, demography comment block. →18282.
- [x] C3 accessor kill ×11 (getGenotypes, liver/bloodEffect, neg_m_exponent,
      molecular_weight_ratio, getPartialEIR ×2, getSpeciesIndexMap,
      GetModelName, setHadEvent/setHadRelapse). →18268.
- [x] C4 folds round 2 (PkPd/Drug, WithinHost/Pathogenesis) + scan round 2
      (numDrugTypes, NvOguess, dead fwd decls ×2, dead test struct). →18261.
- [x] C5 Population accessors (getSize/getRecentBirths/resetRecentBirths →
      public size/recentBirths, 14 sites). Read-pass remainder + wrapper sweep:
      all candidates kept with reasons (see keep-proofs). →18248.
- [x] C6 dangling doxygen block (AnophelesModel.h, comment-only) + final gates:
      ctest 52/52, unittest .OK!, T1 a+b seeds 1-10 = 63/63 on final binary,
      bench median 12.17 s (no regression), PERF.md row 7g. DONE.

## Method / proof tooling
- Trap (cost an hour): pre-rewrite paths (Environment.h, DailyDate.*,
  Scenario.cpp, Clinical/CohortTreatment.*, Clinical/Morbidity.*, Vivax/*,
  VectorBio/*, mon/NaturalScreens.*) do NOT exist in this repo — that code was
  folded into model/xml/*.cpp, model/core/, Host/, Transmission/ in p7b-p7f.
  sed/grep on them silently no-op. Resolve paths with `find model` first.
  (Also: the old DailyDate type itself is gone from the new tree — renamed in
  an earlier phase; its helpers survive.)
- /tmp/opencode/orphan.sh + two inline python scans (decl-identifier vs
  whole-tree occurrence count, comments stripped): function defs (.cpp),
  header decls, header inline defs, member fields, type names.
- model is a static lib → linker catches nothing; grep occurrence count is the
  proof, quoted in commit messages.

## Do-not-touch / keep-proofs
- xoshiro.h + contrib/chacha.h verbatim (contract) — chacha_core/generate_block
  look dead but are copy-mandated.
- CxxTest suites + test methods (runner codegen finds them); 6 `dont_testMol*`
  disabled fixtures (suite-deletion veto).
- slope()/max_killing_rate()/genotype()/origin()/N()/nowOrTs0/1 — readable
  formula names over raw fields, 7f policy; keep.
- Human age()/rng()/monitoringAgeGroup() encode SoA column lookups — keep.
- Kind-dispatch one-liners in WHInterface/ClinicalModel = variant dispatch,
  variants are user-veto-protected; keep.
- timed_deploy_* uniform signature (family dispatch) — unused params stay.
- CalcNvO solver untouched (no-solver-touch rule).
- Checkpoint wire frozen: VivaxBrood hadEvent/hadRelapse are write-only
  in-model but wire-streamed — kept as public fields, setters killed.

## Traps
- -ffp-contract=off always; no reassociation ever (phase 8).
- mon/Monitoring.h now carries the 3 scenario-load decls (readSurveyDates /
  initReporting/initCohorts) — defs live in xml/monitoring.cpp; scnXml fwd
  decls added to Monitoring.h (order-sensitive section: keep near top).
- Episode/TransmissionModel ctor init-lists reordered to match declaration
  order (initializers independent) — re-verify when touching either.
- Untracked PROMPT.md/email/pgodata-p7d/ = user files, leave alone.
