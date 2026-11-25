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

#ifndef Hmod_humans_HostColumns
#define Hmod_humans_HostColumns

#include "Global.h"
#include "interventions/ComponentId.h"
#include "util/random.h"
#include <utility>
#include <vector>

namespace OM {

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
    std::vector<util::LocalRng> rng;

    /// Date of birth (old Human::dateOfBirth).
    std::vector<SimTime> dob;
    /// Time of death; sim::never() while alive (old Human::dead bool — the
    /// wire never carried `dead`, so this column is deliberately not on it).
    std::vector<SimTime> death;
    /// Monitoring cohort-set bitmask (old Human::cohortSet).
    std::vector<uint32_t> cohort_set;
    /// Sub-population memberships (old Human::subPopExp): per row, (id,expiry)
    /// pairs sorted ascending by id — same iteration order as the old map.
    std::vector<std::vector<std::pair<interventions::ComponentId, SimTime>>> cohorts;
    /// Monitoring age-group index (old Human::monitoringAgeGroup, size_t).
    std::vector<uint16_t> mon_age_group;
    /// Next continuous-distribution index (old Human::nextCtsDist).
    std::vector<uint32_t> next_cts_dist;

    /** Add a row for a newly constructed Human, at the exact moment of its old
     * member inits: rng is copy-seeded from `master` HERE (RNG RULE); dob set;
     * everything else at construction defaults (alive, no cohorts). */
    size_t append_row( SimTime dateOfBirth, util::MasterRng& master ) {
        rng.emplace_back( master );
        dob.push_back( dateOfBirth );
        death.push_back( sim::never() );
        cohort_set.push_back( 0u );
        cohorts.emplace_back();
        mon_age_group.push_back( 0u );
        next_cts_dist.push_back( 0u );
        return rng.size() - 1;
    }

    /// Drop row `r` from every column (Population-only; mirrors humans.erase).
    void erase_at( size_t r ) {
        rng.erase( rng.begin() + r );
        dob.erase( dob.begin() + r );
        death.erase( death.begin() + r );
        cohort_set.erase( cohort_set.begin() + r );
        cohorts.erase( cohorts.begin() + r );
        mon_age_group.erase( mon_age_group.begin() + r );
        next_cts_dist.erase( next_cts_dist.begin() + r );
    }
};

}
#endif
