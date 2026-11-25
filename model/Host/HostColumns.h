/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2025-2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#pragma once

#include "Clock.h"
#include "interventions/ComponentId.h"
#include "util/random.h"
#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace OpenMalaria {

class DecayFunction;

/** Per-host vector-intervention record (ITN/IRS/GVI), stored in DEPLOY-HISTORY
 * order (parity: c — the effect FP products multiply in this exact order;
 * redeploy updates the matching record in place, never reorders).
 * §9 3-state: NEVER = no record for the component; ACTIVE = deploy != never();
 * DISCARDED = deploy == never() after ITN attrition. decay/decay_factor
 * replace the old per-host DecayFunction hetSample clone (eval_with applies
 * the factor; pointer targets the component's config). */
struct IvRecord {
    /// ComponentId has no default ctor; max() = never-a-valid-id (records are
    /// filled field-by-field at deploy/checkpoint-read before use).
    IvRecord() : id(std::numeric_limits<size_t>::max()) {}
    ComponentId id;
    Component::Type kind;                 // type-rank (ITN/IRS/GVI)
    SimTime deploy = never();             // time of last (re)deployment
    SimTime disposal = never();           // ITN attrition time
    int n_holes = 0;                      // ITN: total holes
    double hole_index = 0.0;              // ITN: number+size measure
    double initial_insecticide = 0.0;     // ITN/IRS mg/m² (resampled per deploy)
    double hole_rate = 0.0;               // ITN holes/step (sampled once per human)
    double rip_rate = 0.0;                // ITN rips/hole/step (ditto)
    const DecayFunction* decay = nullptr; // component config (never owned)
    double decay_factor = 0.0;            // per-human het sample of decay
};

/** Per-human vaccine dose record (old PerEffectPerHumanVaccine): one per
 * vaccine COMPONENT deployed on this human, insertion order = the order
 * vaccine_factor multiplies effects in (parity: c; redeploy updates the
 * matching record in place, never reorders). */
struct VacRecord {
    VacRecord() : id(std::numeric_limits<size_t>::max()) {}
    ComponentId id;
    uint32_t num_doses = 0;               // doses administered (EPI-schedule guard)
    SimTime time_last = never();          // time of last vaccination
    std::vector<double> per_genotype_eff; // efficacy at last deployment (undecayed)
    const DecayFunction* decay = nullptr; // component config (never owned)
    double decay_factor = 0.0;            // per-human het sample (sampled once)
};

/** Persistent per-host SoA columns (plan §7 row 6 pinned list; owner module,
 * passed by reference, never global).
 *
 * Row i of every column belongs to Population::humans[i]. Only Population
 * changes the row count (create/checkpoint-read append; update erases) and
 * keeps Human::row in sync after every erase (renumber pass).
 *
 * RNG RULE (census fact 1): the rng row is constructed — copy-seeded from the
 * master — exactly at the owning Human's old member-init moment. Moving the
 * storage is legal; moving the construction moment is fatal (stream order). */
struct HumanColumns {
    /** Per-host RNG; 4×u64 Xoshiro256P state each (copy-seed ctor runs at the
     * Human's construction point — Human ctor body, old member-init slot). */
    std::vector<LocalRng> rng;

    /// Date of birth (old Human::dateOfBirth).
    std::vector<SimTime> dob;
    /// Time of death; never() while alive (old Human::dead bool — the
    /// wire never carried `dead`, so this column is deliberately not on it).
    std::vector<SimTime> death;
    /// Monitoring cohort-set bitmask (old Human::cohortSet).
    std::vector<uint32_t> cohort_set;
    /// Sub-population memberships (old Human::subPopExp): per row, (id,expiry)
    /// pairs sorted ascending by id — same iteration order as the old map.
    std::vector<std::vector<std::pair<ComponentId, SimTime>>> cohorts;
    /// Monitoring age-group index (old Human::monitoringAgeGroup, size_t).
    std::vector<uint16_t> mon_age_group;
    /// Next continuous-distribution index (old Human::nextCtsDist).
    std::vector<uint32_t> next_cts_dist;
    /// Immunity accumulators (old WHInterface m_cumulative_h/_Y/_Y_lag, f64).
    std::vector<double> imm_h, imm_Y, imm_Y_lag;
    /// Primary blood-stage infection count (old WHInterface cumPrimInf, u32).
    std::vector<uint32_t> cum_prim_inf;
    /// Clinical timestamps (old ClinicalModel m_tLastTreatment + ES-kind
    /// esCaseStartTime/esTimeOfRecovery/esTimeLastTreatment, SimTime).
    std::vector<SimTime> t_last_treatment, es_case_start, es_time_recovery, es_time_treatment;
    /// ES total density last time step (old ClinicalModel esPreviousDensity).
    std::vector<double> es_prev_density;
    /// Per-species transmission params (old PerHost anoph* vectors; nested:
    /// one vector<double> per row, resized by PerHost::initialise).
    std::vector<std::vector<double>> anoph_ento_raw, anoph_ento_avail, anoph_prob_biting, anoph_prob_resting;
    /// Availability heterogeneity factor (old relativeAvailabilityHet;
    /// signaling_NaN until PerHost::initialise, like the old member).
    std::vector<double> rel_avail_het;
    /// Outside-transmission flag (old bool PerHost::outsideTransmission;
    /// u8 keeps the 1-byte wire shape).
    std::vector<uint8_t> outside_trans;
    /// Vector-intervention deployment records (old PerHost::activeComponents
    /// polymorphic list; history order — see IvRecord).
    std::vector<std::vector<IvRecord>> iv;
    /// Vaccine dose records (old Human::vaccine PerHumanVaccine::effects list;
    /// insertion order — see VacRecord).
    std::vector<std::vector<VacRecord>> vac;
    /// Gaussian-copula deployment samples (old PerHost::copulaGaussianSamples
    /// map, 7.4): one gaussian per (human, subPop component) the first time a
    /// correlated deployment attempts that sub-population; value never changes.
    std::vector<std::map<ComponentId, double>> copula_g;

    /** Add a row for a newly constructed Human, at the exact moment of its old
     * member inits: rng is copy-seeded from `master` HERE (RNG RULE); dob set;
     * everything else at construction defaults (alive, no cohorts). */
    size_t append_row(SimTime dateOfBirth, MasterRng& master) {
        rng.emplace_back(master);
        dob.push_back(dateOfBirth);
        death.push_back(never());
        cohort_set.push_back(0u);
        cohorts.emplace_back();
        mon_age_group.push_back(0u);
        next_cts_dist.push_back(0u);
        imm_h.push_back(0.0);
        imm_Y.push_back(0.0);
        imm_Y_lag.push_back(0.0);
        cum_prim_inf.push_back(0u);
        t_last_treatment.push_back(never());
        es_case_start.push_back(never());
        es_time_recovery.push_back(never());
        es_time_treatment.push_back(never());
        es_prev_density.push_back(std::numeric_limits<double>::quiet_NaN());
        anoph_ento_raw.emplace_back();
        anoph_ento_avail.emplace_back();
        anoph_prob_biting.emplace_back();
        anoph_prob_resting.emplace_back();
        rel_avail_het.push_back(std::numeric_limits<double>::signaling_NaN());
        outside_trans.push_back(0u);
        iv.emplace_back();
        vac.emplace_back();
        copula_g.emplace_back();
        return rng.size() - 1;
    }

    /// Drop row `r` from every column (Population-only; mirrors humans.erase).
    /// Adding a column: append to one of the e(...) calls below.
    void erase_at(size_t r) {
        const auto e = [r](auto& col) { col.erase(col.begin() + r); };
        e(rng);
        e(dob);
        e(death);
        e(cohort_set);
        e(cohorts);
        e(mon_age_group);
        e(next_cts_dist);
        e(imm_h);
        e(imm_Y);
        e(imm_Y_lag);
        e(cum_prim_inf);
        e(t_last_treatment);
        e(es_case_start);
        e(es_time_recovery);
        e(es_time_treatment);
        e(es_prev_density);
        e(anoph_ento_raw);
        e(anoph_ento_avail);
        e(anoph_prob_biting);
        e(anoph_prob_resting);
        e(rel_avail_het);
        e(outside_trans);
        e(iv);
        e(vac);
        e(copula_g);
    }
};

}
