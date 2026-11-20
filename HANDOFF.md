# Phase 7 HANDOFF — interventions rewrite (branch `p7-interventions` off `p6-green` merge 9c077bc8)

Crash-recovery file. Update EVERY green commit. Baseline p6: model **18,895** /
unittest 3,794 LOC; `interventions/` dir 2,837 LOC (cloc); p6 bench **12.58s** mad
.010 (A0 12.33s). Full phase-6 handoff: `git show 4ce2157d:HANDOFF.md` (p6 6.5).

## Environment (REBUILT 2026-11-20 — /tmp was wiped once; recreate same way)
- Build: `build-p7/` = Release `-DCMAKE_CXX_FLAGS=-ffp-contract=off`
  `-DOM_BOXTEST_BINARY=ON`. `cmake --build build-p7 -j16`.
- unittest MUST run from `build-p7/unittest/` → `./unittest` (94). ctest:
  `cd build-p7 && ctest -j16` (52).
- T1 reference A0 = baseline commit `cfc3c8e0` built Release+`-ffp-contract=off`:
  `git worktree add /tmp/opencode/om-ref-src cfc3c8e0` + cmake-build → copy binary
  to `/tmp/opencode/om-ref/A0/openMalaria`; regenerate `/tmp/opencode/om-ref/A0/recipe.json`
  with the 4 keys parity.py::check_recipe probes (compiler, gsl, gsl_so_sha256,
  machine; extra keys are informational). Rebuilt + verified: tier-a 53/53 old==new
  at p6 code, so the A0 rebuild is trustworthy.
- T1 tier-a: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria
  --new build-p7/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a -j16` (53).
  Tier-b: `--tiers b --seeds 1-4` per cluster, 1-10 at final.
- bench.sh at /tmp/opencode/om-bench/bench.sh was LOST with /tmp — reconstruct per
  PERF.md §"CPU-s" before the final gate: `/usr/bin/time -f %U`, taskset-pinned,
  example_scenario popSize 2000 iseed 1, median-of-5, warmup/main split. p6 = 12.58s.
- Checkpoint smoke idiom (run.py mechanism, verified): in a dir with
  `scenario_current.xsd -> schema/scenario_50.xsd`: run
  `BIN --resource-path test --scenario SC --output-format bin --output output.bin --checkpoint-stop`
  TWICE — run 1 writes checkpoint(+marker `checkpoint`)+stops, run 2 auto-resumes and
  writes output.bin. compareBin old vs new. scenario12/12a lack ITNs, so 2ITNs smoke above
  is the per-host-wire canary until 7.5 deletes it.
- perf events BLOCKED (paranoid=2). Untracked PROMPT.md/email = user files, leave.

## Plan (order = plan §7 row 7 / §9; each step ends green+committed)
1. [x] 7.1 flat timed deployments — see Gates log.
2. [~] 7.2 dense components by type-rank:
   [x] 7.2a-i per-host VECTOR-intervention state dissolved: HumanITN/HumanIRS/
   HumanGVI/PerHostInterventionData classes DEAD; `iv` column
   (`vector<vector<IvRecord>>` in HostColumns, history order = contract c);
   free fns `iv_deploy/iv_update/iv_*/iv_checkpoint` in interventions/IvState.{h,cpp}
   (component headers friend them for private params); DecayFunction: `compute` →
   `compute_with(age,factor)` + public `eval_with(age,factor)` +
   `sample_het_factor(rng|NormalSample)` (het-clone replacement, bit-identical);
   HumanVectorInterventionComponent now carries NO per-host API (makeHumanPart/
   asITN/asIRS virtuals deleted); hasActiveInterv = kind+deploy scan, no registry.
   Wire = NEW uniform record layout (wire free; resume bit-verified vs A0, see gates).
   [ ] 7.2a-ii vaccine: PerHumanVaccine/PerEffectPerHumanVaccine → record column
   (Human::vaccine member + 12 getFactor sites + its wire).
   [ ] 7.2b config-side HumanInterventionComponent tree flatten (deploy switch).
3. [ ] 7.3 ITN 3-state columns NEVER/ACTIVE/DISCARDED + §9 scalars; draw order EXACT
   (1 shared NormalSample → hole/rip/insecticide-decay (0 extra draws), attrition age
   (0–1 draws), initial insecticide; update = Poisson(holes) then Poisson(rips) only if
   n_holes*rip_rate>0; IRS order = insecticide THEN het_sample).
4. [ ] 7.4 copula g per (host, subPop-id — NOT component) → column (replaces
   `copulaGaussianSamples` std::map; iteration/wire order was ComponentId-ascending).
5. [ ] 7.5 zero intervention checkpoint bytes: delete `checkpointIntervs` from
   PerHost::operator& + PerHostInterventionData wire fns (wire NOT frozen — user ruling;
   scenario 12/12a goldens may re-bless from A0 in same commit if they move — expected
   NO change: 12/12a have no per-host intervention state and the wire removal is
   symmetric). HostView gains `iv` ref (+mon §9; deviation OK with rationale).
   Delete HumanInterventionComponent class tree + per-timestep per-entity scans.
6. [ ] 7.6 FINAL GATE: T1 a+b 63/63 seeds 1-10, ctest 52/52, unittest 94/94, bench +
   ITN-heavy secondary (`/usr/bin/time -f %U` VecFull or LifeNet old vs new), cloc,
   census refresh, PERF row 7, report (kind tables, LOC by substep, phase-8 handoff).

## Facts found this phase (archaeology, saves a re-read)
- Timed deploys (7.1, done): flat `TimedDeploymentKind` + variant payload in
  Deployments.h; deploy bodies verbatim (header inline for human/vector kinds,
  `timed_deploy_*` in xml/interventions.cpp for HS/EIR/NHH); `check_timed_date` at
  load keeps old base-ctor throw order; `replayed_on_restore(kind)` seam;
  builders `timed_deployment(date, kind[, payload])` (struct is move-only, not an
  aggregate). ChangeHS/EIR payload raw ptrs deleted on deploy; never-deployed = leak
  at exit (old behaviour preserved deliberately). Sentinel = `{sim::future(), DUMMY}`.
- PerHost (p7 target, model/Transmission/PerHost.{h,cpp}): members NOW = cols/row ptrs
  (transmission cols, p6), `copulaGaussianSamples` map, `activeComponents`
  vector<unique_ptr<PerHostInterventionData>> + static relAvailAge. `update` = virtual
  loop over activeComponents, called ONLY from `TransmissionModel::calculateEIR_vector`
  (TransmissionModel.cpp:340). Effect accessors (entoAvailabilityHetVecItv/probMosqBiting/
  probMosqResting/relMosqFecundity) = ordered FP products over activeComponents —
  `// parity: c` on the first; order = deploy history (redeploy mutates in place,
  NEVER reorders). hasActiveInterv(type) loops + getComponent(...).componentType().
  Wire: PerHost::operator& streams transmission cols + copula map + checkpointIntervs
  (count + per-component id/wire). initialise = 4 draws/species (ento sample, prob
  biting, 2× resting) — the het_sample §9 4-draw contract.
- PerHostInterventionData subclasses live in interventions/{ITN,IRS,GVI}.h(.cpp)
  (HumanITN/HumanIRS/HumanGVI) + HumanVaccine in Vaccine.h is a Human member, not
  activeComponents. Deploy chain: HumanIntervention::deploy → component deploy →
  `PerHost::deployComponent` (matches by component id: redeploy in place, else
  makeHumanPart push). 2ITNs scenario = two ITN components on one host — history
  order + multi-instance-per-type must keep working.
- `human.vaccine.getFactor(TBV/PEV/BSV, genotype)` consumed in TransmissionModel.cpp:447
  + within-host/incidence sites (~10 total; grep getFactor).
- Population.cpp `hasActiveInterv` survey callbacks (cts*Coverage) must migrate with
  the state columns.

## Traps (inherited from p6 — all still live)
- set_row REBIND: every submodel binding (cols,row) must rebind in Human::set_row;
  the new iv column follows the same rule (death-heavy ctest = canary; tier-a misses it).
- Compaction read-by-position inside Population::update erase loop (Human::row stale).
- m_rng/LocalRng construction consumes master draws — column appends happen exactly at
  old ctor points. TwoRunsSuite = statics tripwire (rerun on RunExtensions change).
- Error strings + cerr/cout texts byte-identical (PRINT_INTERVENTIONS, VERBOSE lines,
  "cum deployment: ..." — all re-proven in 7.1). `// parity:` comments travel.
- vector erase(end()) UB trap for map→vector migrations (p6 6.2b).
- Draw-order landmines: copula g cached by subPop id; ITN first-creation uses
  clock.now(), redeploy uses nowOrTs1(); Screen/DecisionTree components call
  within-host/clinical code at deploy time (recursion keeps order); TimedHumanDeployment
  uses anoph_ento_AVAIL sums, ContinuousHumanDeployment uses anoph_ento_RAW (quirk —
  do NOT "fix").

## Gates log
- (7.2a-i, this commit) per-host vector-intervention state → `iv` record column.
  Deleted: HumanITN/HumanIRS/HumanGVI/PerHostInterventionData + makeHumanPart/
  asITN/asIRS virtuals + PerHost::deployComponent/checkpointIntervs (−589 lines,
  model 18,895→18,854). New: HostColumns IvRecord+iv (append_row/erase_at mirrored),
  interventions/IvState.{h,cpp} (deploy/update/4 effects/checkpoint free fns,
  component classes friend them), DecayFunction eval_with/sample_het_factor.
  Draw orders moved verbatim (ITN: shared NormalSample→hole/rip/insecticide-het→
  attrition→clamped insecticide, ctor deployTime=now(), redeploy=nowOrTs1(); IRS:
  insecticide FIRST then het, per old ctor member-init; GVI: single het).
  build clean; unittest 94/94; ctest 52/52 (incl. checkpoint two-run suites);
  T1 a 53/53; T1 b 1-4 4/4; 2ITNs new-binary two-run resume == A0 single run
  (new uniform record wire round-trips: id,kind,deploy,disposal,n_holes,hole_index,
  insecticide,hole_rate,rip_rate,factor; decay ptr recovered from registry on load).
- (7.1, earlier commit) flat TimedDeploymentKind + variant payloads; class tree deleted
  (TimedDeployment/Derivatives + 5 xml-local classes → read/builders + free deploy
  fns); RunExtensions timed = value vector; deploy loop + loadFromCheckpoint flat.
  A0 reference REBUILT from cfc3c8e0 (recipe regenerated; p6-code tier-a 53/53 proves
  equivalence). build clean; unittest 94/94; ctest 52/52; T1 a 53/53; T1 b seeds1-4
  4/4; NEW checkpoint+ITN 2-run smoke (2ITNs --checkpoint-stop×2, old vs new)
  bit-identical. RngGoldenAdapter anchor re-pointed to Deployments.h:195.
