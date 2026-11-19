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
#include "util/random.h"#include <vector>

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

    /// Drop row `r` from every column (Population-only; mirrors humans.erase).
    void erase_at( size_t r ) { rng.erase( rng.begin() + r ); }
};

}
#endif
