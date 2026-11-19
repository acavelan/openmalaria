# Phase 6 HANDOFF — SoA sidecar + SlotAllocator + Host::update 7 fns (branch `p6-soa` off `p5-green`)

Crash-recovery file. Update EVERY green commit. Baseline p5-green: model **18,657** / unittest 3,760 LOC; A0 bench 12.33s, p5 bench **12.51s** mad .010 (residual +1.5%; target ≤12.33). Full phase-5 handoff: `git show p5-green:HANDOFF.md`.

## Environment
- Build: `build-p6/` = Release `-ffp-contract=off` `-DOM_BOXTEST_BINARY=ON` (configured).
- unittest MUST run from `build-p6/unittest/` → `./unittest` (94 tests since 6.1: +SlotAllocatorSuite×3). ctest: `cd build-p6 && ctest -j16` (52).
- T1 tier-a: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new build-p6/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a -j16` (53). Tier-b: `--tiers b --seeds 1-10` (per-cluster duty: 1-4).
- bench: `sh /tmp/opencode/om-bench/bench.sh /home/razlock/git/fork/build-p6/openMalaria TAG` (ABSOLUTE path!).
- perf: events BLOCKED (paranoid=2) — `perf record -e task-clock -F 3997` only.

## Steps (order = plan)
1. [x] 6.1 SlotAllocator — single clamp owner `model/infections/SlotAllocator.h`
   (constexpr MAX_INFECTIONS=21 + try_import/clamp_new/free_slots/clamp_count);
   replaced sites: CommonWithinHost import+update, DescriptiveWithinHost import+
   update, InfectionIncidenceModel.cpp numNewInfections clamp, WHInterface.cpp
   checkpoint validation; WHInterface.h member constant DELETED. T2 suite
   `unittest/SlotAllocatorSuite.h` (mode-A golden words + at-cap 0-draw rel. pin
   + clamp table). NOTE goldens are FULL-SUITE-CONTEXT fixed (registration order
   frozen; --abort-on-fail) — ctor consumes exactly 1 draw, below-cap import +1.
   UnittestUtil gained setNumInfs/getNumInfs/resetModelOptions (friend).
2. [ ] 6.2 columns — `model/humans/HostColumns.h` struct, ONE vector per
   pinned column; owned by `Population` as member `cols`; zero globals.
   [x] cluster a (rng): `HumanColumns{ vector<LocalRng> rng; erase_at(r) }`;
   Human gained `cols`/`row` + accessor `rng()` (old `human.rng` syntax kept
   as `human.rng()`); ctor appends its row and copy-seeds the rng IN THE CTOR
   BODY at the old member-init slot (census fact 1 — T1-proven). Population
   mirrors erase (cols.erase_at) + renumber pass after the erase loop;
   push_back rows are correct by construction. UnittestUtil gained testCols().
   `Human` gains `HumanColumns* cols; size_t row;` (ctor appends its row EXACTLY
   where it constructs today — rng column appends `LocalRng(ctx.master_RNG)`
   IN PLACE of the old member init: same ctor point = same master consumption,
   census fact 1/N5; RNG move-ctor self-fixes its gsl self-pointers, verified).
   Submodel classes that host migrated members (PerHost, WHInterface,
   ClinicalModel) store `HumanColumns* cols; size_t row;` bound by the Human ctor;
   Population is the ONLY place that changes vector size (createInitialHumans /
   update erase+push_back / checkpoint read) — after the erase loop it MUST
   renumber (`humans[i].set_row(i)` incl. submodels) and mirror every
   erase/push into cols. Tests build standalone Humans via UnittestUtil::createHuman
   — give it an overload taking cols& (or keep static mini-columns there).
   Cluster order: rng → dob/death/monAgeGroup/nextCtsDist/cohort (bitmask u32 +
   sparse (id,expiry) pairs; replace std::map with per-row sorted vector<pair> —
   iteration order ComponentId-ascending is preserved, mon order safe) →
   [x] cluster b (dob, death, cohort_set u32 + cohorts sparse pairs,
   mon_age_group u16, next_cts_dist u32): old Human members DELETED;
   accessors keep syntax (monitoringAgeGroup()/nextCtsDist() new;
   mon::updateAgeGroup takes uint16_t&). Wire replicated byte-shape for
   cohorts: count(size_t) then (id as size_t, expiry int) ascending — read via
   ComponentId(stream). death column never()⇔alive, NOT on wire (old `dead`
   wasn't either). mon field size_t→u16 shrinks the wire (legal: write/read
   symmetric; cp files not byte-frozen vs p5).
   immunity (falc: m_cumulative_h/_Y/_Y_lag; vivax: cumPrimInf) → clinical
   (m_tLastTreatment + es 3×SimTime + esPreviousDensity f64) → transmission
   (anophEntoAvailabilityRaw/anophEntoAvailability/anophProbMosqBiting/
   anophProbMosqResting per species + relativeAvailabilityHet + outsideTransmission).
   Each cluster: DELETE old member, fix call sites (compiler = checker), tier-a,
   tier-b 1-4, HANDOFF, commit. Intervention-effect columns stay in PerHost
   (phase-7 landing zone — do NOT touch polymorphism).
3. [ ] 6.3 scratch (consumer-owned, per-sweep): EIR_per_genotype_i/l stack vectors
   in Host::update loop + probTransGenotype_i/l in TransmissionModel::updateKappa
   are per-host reallocations inside clean sweeps (run.cpp loop, updateKappa loop)
   → hoist into a scratch struct owned by the sweep caller. age_years/availability
   are recomputed per use, no per-call alloc — leave unless bench shows cost.
4. [ ] 6.4 HostView + 7 free fns (§9 order: host_is_dead, host_update_mon_cohorts,
   host_calculate_eir, host_sample_new_infections, host_within_host_update,
   host_clinical_update, host_update_infant_deaths) — bodies MOVED verbatim from
   Host::update (Human.cpp), // parity comments travel. HostView{cols ref, row,
   rng ref, scratch, ...} per §9. het_sample = Human.cpp hetSample (contract:
   max 4 draws, 3 leading ifs + else-if chain — verified vs §9, annotate).
5. [ ] 6.5 FINAL GATE: T1 a+b 63/63; ctest+unittest; bench p6-final vs 12.51
   (recover toward ≤12.33; revert offender if worse > noise floor max(1%,2×MAD));
   cloc; census statics/draws refresh if moved; PERF.md row 6; column table +
   7 handoff notes for phase 7.

## Traps (inherited + p6-specific)
- Compaction read-by-position: inside Population::update's erase loop,
  Human::row is STALE after any erase (renumber runs after the loop) — read
  cols[i]/cols.death[i] by POSITION there, never via Human accessors.
- vector erase() invalidates the end() sentinel (std::map loops didn't): the
  old `it = v.erase(it)` + cached-end pattern reads past the NEW end →
  erase(end()) = UB memmove that tramples the heap (root cause of the 6.2b
  2ITNs/IRS30/Rach5IC SIGSEGVs; index iteration is the fix). Audit any other
  map→vector loop migrations for this before touching them.
- m_rng(ctx.master_RNG) CONSTRUCTS-BY-DRAWING — only where old code did. In 6.2
  the append of the rng row IS the construction moment; never default-then-assign.
- TwoRunsSuite = ext/statics tripwire — rerun on any RunExtensions change.
- Error strings byte-identical; `// parity:` comments travel with moved code.
- Checkpoint: columns keep SEMANTIC round-trip (byte order now free, but
  scenario 12/12a two-run goldens are the canary — they passed 6.1 unchanged).
- Human move ops: adding raw ptr/size members keeps move=default fine; if a member
  gains a user-declared dtor re-declare moves (p5 lesson).
- PerHost move ops are out-of-line on purpose (p5) — keep.
- `grep -v LSTMDrug` blind spot still applies.
- Untracked PROMPT.md/email = user files, leave.

## Gates log
- (6.1, this commit) build clean; unittest 94/94; ctest 52/52; T1 a 53/53; T1 b seeds1-4 4/4. Model LOC 18,673 (+16 vs p5: SlotAllocator struct; WHInterface constant removed). unittest LOC 3,787 (+27: new suite + 3 helpers).
- (6.2a) rng column; same gates green (see git log e2c87de0).
- (6.2b, this commit) dob/death/cohort_set+cohorts/mon_age_group/next_cts_dist columns; old members deleted; compaction read-by-position fix; updateCohortSet index-loop fix (see Traps). build clean; unittest 94/94 (run from build-p6/unittest/); ctest 52/52; T1 a 53/53; T1 b seeds1-4 4/4. bench 12.59s mad .010 (+0.08s vs p5 — under 1% gate, recover in 6.3/6.4). Model LOC 18,797. Debug tooling built during the hunt (throwaway): /tmp/opencode/p5dbg (p5 twin w/ prints), /tmp/opencode/build-p6-rwdi (RelWithDebInfo twin), /tmp/opencode/build-p6-asan — p5dbg worktree REMOVED after use; if a hunt recurs: gdb-disable-ASLR + python watchpoint scripts in /tmp/opencode/irs/watch*.py worked where ASan's quarantine lied.
