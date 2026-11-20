# Phase 8-CK HANDOFF — kill checkpointing, redo it xml-style (branch `p8-checkpoint` off `rewrite` HEAD = `p7g-green` e34997f1)

User order (verbatim): "kill checkpointing entirely (one commit) then redo the
checkpointing similar to xml/".

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

## REDO design (in progress — see checkboxes)
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
- [x] C4 (partial): magic/version/budget/count checks → clean errors, proven
      in resume.py (truncated) + manual junk/bad-magic/wrong-scenario runs.
      v2 (boost/gzstream) file → magic error (message could name it nicer).
- [x] C2 scripted save days per scenario in resume.py: main-start,
      mid-main (mainStart+50d, chained), late-main (end-25d, chained),
      mid-warmup (mainStart/2). Chains exercise resume→save→resume.
- [x] C3 (no-intervention corpus): Resume_* now 12 — + scenario1 (NonVector
      EIR branch), VecTest, Genotypes (SAMPLE_TRACKING mode), ExactNv0
      (exact solver — NOTE: seed-sensitive; resume.py uses the scenario's
      own iseed, NOT --seed 112), SimpleMPDTest (MPD state),
      HetVecDailyEIR, MolPairwise (multi-clone variants).
      Interventions-heavy = still the documented gap (skill law).
- [ ] C3 corpus: example, DecisionTree(5Day), Vivax, HetVecDailyEIR /
      vector-full (AnophelesModel + ExactNv0), Molineaux/Penny, Empirical,
      interventions-heavy (IvState global + replay), MSAT/Genotypes.
- [ ] C4 bad/old/truncated checkpoint file → clean model error (never UB);
      version mismatch message.
- [ ] C5 final: full resume-equivalence suite green, LOC before/after report,
      tier-a clean. PERF.md NOT touched (main-session ledger).

## Proven round-trips (update every green commit)
C2/C3 (all byte-identical vs uninterrupted, output.bin AND ctsout.txt,
scenario default seed): 4 save-day chains each (main-start, mid-main,
late-main end-25d, mid-warmup) + truncated-file rejection for:
- NoInterv (Descriptive WH + cts ~1.5k lines), Vivax (broods), Empirical,
  Molineaux (growing variant counts), Penny
- 1 (NonVector EIR), VecTest, Genotypes (tracking mode), ExactNv0 (exact
  emergence estimator), SimpleMPDTest (SimpleMPD state), HetVecDailyEIR,
  MolPairwise
- corrupt: truncated → rc73; junk/bad magic → rc73; wrong scenario → rc73
  (measure store count mismatch) — no UB, no output file.
Gates at C1: build ✓, ctest 57/57, unittest .OK!, tier-a 53/53.
Gates at C2/C3: ctest 64/64, unittest .OK! (model binary byte-identical to
the C1 parity-verified one — no model/ change, so tier-a not re-run).

## Known limits (honest partial — extend in C3+)
- NO interventions section yet (skill law): --resume on a scenario that
  deploys anything is UNSUPPORTED (nextTimed/IvState/continuous deployments
  not on the wire; deployments never replayed either). Proven days are all
  intervention-free scenarios.
- SimpleDecayingValue: writes deploy_t + het-presence flag; RESUMING a
  scenario with a SAMPLED decay deployment throws a clean
  "sampled decay deployment not yet supported" error.
- Cohorts intervention-column state is config-rebuilt (empty), not wire.
- cts resume overwrites file content past the saved offset (rows the later
  part of a previous longer run wrote are replaced; fine for save+exit use).
