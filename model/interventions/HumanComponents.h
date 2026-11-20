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
#include "interventions/Interfaces.h"
#include "util/sampler.h"
#include "util/DecayFunction.h"

namespace OM {
namespace Host {
    class Human;
}
namespace interventions {
    using util::LocalRng;

// ———  sub-population removal  ———

namespace SubPopRemove{
    /** Identifiers for types of sub-population removal option, followed by NUM
     * in last place (as a counter).
     * 
     * ON_FIRST_BOUT: remove humans from some sub-populations at the start of
     * their first clinical events following inclusion in the sub-population.
     * 
     * ON_FIRST_INFECTION: remove patent humans from some sub-population(s)
     * during each survey.
     * 
     * ON_FIRST_TREATMENT: remove humans from some sub-populations when
     * assigning each course of treatment.
     */
    enum RemoveAtCode {
        ON_FIRST_BOUT,
        ON_FIRST_INFECTION,
        ON_FIRST_TREATMENT,
        NUM };
}

/** For each RemoveAtCode (excluding NUM), this is a list of all sub-population
 * identifiers for which the option is enabled. */

// ———  vaccines  ———

namespace Vaccine{
    enum Types { PEV, BSV, TBV, NumVaccineTypes };
}

/** Per vaccine effect (type), per human details.
 * Phase 7.2a-ii: PerEffectPerHumanVaccine and PerHumanVaccine are gone —
 * per-human state is the VacRecord column (Host/HostColumns.h), read by
 * Human::vaccine_factor and written by interventions::vaccinate
 * (interventions/Vaccine.cpp). */

}
}
