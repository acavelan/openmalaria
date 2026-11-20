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
#include "Context.h"
#include "Host/HostColumns.h"
#include "Transmission/PerHost.h"
#include "IncidenceModel.h"
#include "interventions/HumanComponents.h"
#include <map>

namespace OM {
namespace Clinical {
class Model;
}
namespace WithinHost {
class WHInterface;
}
namespace Transmission {
class Model;
}

namespace Host {

/** Interface to all sub-models storing data per-human individual.
 *
 * Still contains some data, but most is now contained in sub-models. */
class Human {
  public:
    /// The run this human belongs to (clock: ctx->clock). Set at
    /// construction; a pointer so Human stays move-assignable.
    Context* ctx;

    /// This human's SoA columns (owned by Population; same lifetime as ctx).
    HumanColumns* cols;
    /// This human's row index in every column (kept in sync by Population;
    /// renumbered after every structural erase — see Population::update).
    size_t row;

    /** Initialise all variables of a human datatype.
    * 
    * @param dateOfBirth date of birth (usually start of next time step)
    * @param columns persistent columns receiving this human's row (the rng
    * row is copy-seeded from the master HERE — the old member-init moment) */
    Human(Context& ctx, SimTime dateOfBirth, HumanColumns& columns);

    /// Allow move construction
    Human(Human&&) = default;
    Human& operator=(Human&&) = default;

    /// Disable copying
    Human(const Human&) = delete;
    Human& operator=(const Human&) = delete;

    /** Get human's age with respect to some time. */
    inline SimTime age(SimTime time) const;

    /** Return true if human is a member of the sub-population. */
    inline bool isInSubPop(interventions::ComponentId id) const;

    inline uint32_t getCohortSet() const;

    /** Return true if dead, false otherwise */
    bool isDead() const;

    /** After this call, isDead() will be true */
    void kill();

    /** Add the human to an intervention component's sub-population for the given duration. 
     * A duration of zero implies no effect. */
    void addToCohort(interventions::ComponentId id, SimTime duration);

    /** Remove the human from an intervention component's sub-population */
    void removeFromCohort(interventions::ComponentId id);

    /** Act on "remove from sub-population on first ..." events.
     * This is only for use during a human update. */
    void removeFirstEvent(interventions::SubPopRemove::RemoveAtCode code);

    /** Remove host from expired cohorts */
    void updateCohortSet();

    /** Contains per-species vector data (VectorModel only). */
    Transmission::PerHost perHostTransmission;

    /** The WithinHostModel models parasite density and immunity */
    unique_ptr<WithinHost::WHInterface> withinHostModel;

    /** The IncidenceModel translates per-host EIR into new infections */
    unique_ptr<IncidenceModel> infIncidence;

    /** The Model encapsulates pathogenesis (sickness status),
    * case management (medicating drugs)
    * and clinical outcomes (morbidity, reporting). */
    unique_ptr<Clinical::Model> clinicalModel;

    /** Per-host RNG: SoA column row, not a member (accessor below keeps the
     * old `human.rng()` call syntax; 4×u64 Xoshiro state per row). */
    inline util::LocalRng& rng();

    /** Set the column row index (Population-only; see Population::update). */
    void set_row(size_t r);

    /** Vaccine factor: one minus the (decayed) efficacy, multiplied over this
     * human's dose records of the given vaccine type in insertion order
     * (1 for no effect, 0 for full effect). Interventions/Vaccine.cpp. */
    double vaccine_factor(interventions::Vaccine::Types type, uint32_t genotype = 0) const;

    /** True if this human has any vaccine dose record. vaccine_factor() returns
     * exactly 1.0 when false (empty product), so hot loops can hoist the call. */
    inline bool has_vaccines() const { return !cols->vac[row].empty(); }

    /** Monitoring age group (column row; accessor keeps read/write syntax). */
    inline uint16_t& monitoringAgeGroup();
    /// ditto (read-only use, e.g. episode recording)
    inline uint16_t& monitoringAgeGroup() const;

    /** The next continuous distribution in the series (column row). */
    inline uint32_t& nextCtsDist();
};

void summarize(Context& ctx, Human& human, bool surveyOnlyNewEp);

/** Per-sweep scratch reused across humans within a time step (owned by the
 * sweep caller — the run loop). Holds the EIR-per-genotype out-vectors that
 * getEIR overwrites via assign() each host, so the buffer allocates once per
 * run instead of once per host per tick, plus the scalars the §9 host-update
 * stages hand to each other (all written before read within one host's pass). */
struct HostScratch {
    std::vector<double> EIR_per_genotype_i, EIR_per_genotype_l;
    SimTime age0;
    double age1, eir_i, eir_l, eir;
    int nNewInfs_i, nNewInfs_l;
};

/** Per-host view for the 7 §9 update stages (pinned phases-6-7 shape).
 * Phase-6 form carries Human& because the submodels still hang off it;
 * phase 7 dissolves them into columns and the view shrinks toward
 * {cols,row,rng,scratch,iv,mon}. */
struct HostView {
    Context& ctx;
    Human& human;
    Transmission::Model& transmission;
    HostScratch& scratch;
};

/// The 7 §9 stages, in call order; Host::update is the driver. Bodies moved
/// verbatim from the old Host::update (draw/FP order sacred).
bool host_is_dead(HostView& v);
void host_update_mon_cohorts(HostView& v);
void host_calculate_eir(HostView& v);
void host_sample_new_infections(HostView& v);
void host_within_host_update(HostView& v);
void host_clinical_update(HostView& v);
void host_update_infant_deaths(HostView& v);

void update(Context& ctx, Human& human, Transmission::Model& transmission, HostScratch& scratch);

inline SimTime Human::age(SimTime time) const { return time - cols->dob[row]; }

inline util::LocalRng& Human::rng() { return cols->rng[row]; }

inline uint16_t& Human::monitoringAgeGroup() { return cols->mon_age_group[row]; }

inline uint16_t& Human::monitoringAgeGroup() const { return cols->mon_age_group[row]; }

inline uint32_t& Human::nextCtsDist() { return cols->next_cts_dist[row]; }

inline bool Human::isInSubPop(interventions::ComponentId id) const {
    // Definition (unchanged from the old std::map): member if (id,expiry) is
    // present with expiry >= nowOrTs0(). Few spans per host; linear is fine.
    for (auto const& span : cols->cohorts[row]) {
        if (span.first == id)
            return span.second >= ctx->clock.nowOrTs0();
    }
    return false; // no history of membership
}

inline uint32_t Human::getCohortSet() const { return cols->cohort_set[row]; }

inline bool Human::isDead() const { return cols->death[row] != sim::never(); }

inline void Human::kill() {
    // A human cannot be revived: only Population::update erases the row.
    cols->death[row] = ctx->clock.nowOrTs1();
}

}
}
