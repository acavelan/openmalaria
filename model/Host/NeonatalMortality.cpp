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

#include "Host/NeonatalMortality.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Diagnostic.h"

namespace OM::Host {
namespace {
constexpr double gEst = 0.011, critPrev2025 = 0.25, critPrevPrim = 0.19, pBirthPrim = 0.3;
constexpr double y = pBirthPrim * gEst, z = -1.0 / critPrev2025;
const SimTime ageLb = sim::fromYearsI(20), ageUb = sim::fromYearsI(25);
}

NeonatalMortality makeNeonatalMortality(SimTime interval) {
    if (interval <= 0 || interval > 150)
        throw util::xml_scenario_error("neonatal mortality interval must be in [1, 150] days");
    return {0.0, vector<double>(150 / interval, 0.0)};
}

void updateNeonatalRisk(NeonatalMortality& state, double prev2025, size_t index) {
    double maxPrev = prev2025;
    //update the vector containing the prevalence by gestational age
    state.prevByGestationalAge[index] = prev2025;
    for(size_t i = 0; i < state.prevByGestationalAge.size(); ++i) {
        if (state.prevByGestationalAge[i] > maxPrev) {
            maxPrev = state.prevByGestationalAge[i];
        }
    }

    // equation (2) p 75 AJTMH 75 suppl 2
    double prevPG= maxPrev / (critPrevPrim + maxPrev);
    // equation (1) p 75 AJTMH 75 suppl 2 including 30% multiplier
    state.riskFromMaternalInfection = y * (1.0-exp(prevPG * z));
}

void update(NeonatalMortality& state, vector<Human>& population,
        const WithinHost::Diagnostic& diagnostic, SimTime stepStart, SimTime interval) {
    // ———  find potential mothers and their prevalence  ———
    // For individuals in the age range 20-25, we sum:
    int nCounter=0;	// total number
    int pCounter=0;	// number with patent infections, needed for prev in 20-25y

    for (Human& human : population){
        // diagnosticDefault() gives patency after the last time step's
        // update, so it's appropriate to use age at the beginning of this step.
        SimTime age = human.age(stepStart);

        if( age >= ageUb ) continue;
        if( age < ageLb ) break;	// Not interested in younger individuals.

        nCounter ++;
        if( human.withinHostModel->diagnosticResult(human.rng, diagnostic) ){
            pCounter ++;
        }
    }

    // ———  calculate risk of neonatal mortality  ———
    //default value for prev2025, for use when there are no 20-25 year olds
    double prev2025 = 0.25;
    if( nCounter > 0 )
        prev2025 = double(pCounter) / nCounter;

    size_t index = util::mod_nn(stepStart / interval, state.prevByGestationalAge.size());
    updateNeonatalRisk(state, prev2025, index);
}
}
