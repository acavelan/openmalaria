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

#ifndef Hmod_InfantMortality
#define Hmod_InfantMortality

#include "Global.h"
#include "Parameters.h"

namespace OM {
namespace Clinical {

/// Infant mortality reporting (state owned by Model, checkpointed).
struct InfantMortality {
    /** Initialise parameters. */
    void init( const OM::Parameters& parameters ) {
        infantDeaths.resize(sim::stepsPerYear());
        infantIntervalsAtRisk.resize(sim::stepsPerYear());
        nonMalariaMortality=parameters[Parameter::NON_MALARIA_INFANT_MORTALITY];
    }
    
    void preMainSimInit() {
        std::fill(infantDeaths.begin(), infantDeaths.end(), 0);
        std::fill(infantIntervalsAtRisk.begin(), infantIntervalsAtRisk.end(), 0);
    }
    
    void reportRisk(size_t index, bool isDoomed) {
        infantIntervalsAtRisk[index] += 1;     // baseline
        if (isDoomed)
            infantDeaths[index] += 1;  // deaths
    }
    
    /// Calculate infant mortality as deaths/1000 livebirths for the whole main-
    /// simulation period (not as deaths/1000 years-at-risk per survey).
    /// 
    /// This mimicks field data on all-cause mortality in infants.
    /// Uses the kaplan-meier method because the demography was set up to provide
    /// a stable age-distribution but unfortunately does not accurately describe
    /// death rates. The kaplan-meier estimate is the product of the proportion of
    /// infants survivng at each interval.
    double allCause() const {
        double infantPropSurviving=1.0;       // use to calculate proportion surviving
        for( size_t i = 0; i < sim::stepsPerYear(); i += 1 ){
            // multiply by proportion of infants surviving at each interval
            infantPropSurviving *= double(infantIntervalsAtRisk[i] - infantDeaths[i])
                / double(infantIntervalsAtRisk[i]);
        }
        // Child deaths due to malaria (per 1000), plus non-malaria child deaths. Deaths per 1000 births is the return unit.
        return (1.0 - infantPropSurviving) * 1000.0 + nonMalariaMortality;
    }
    
    /// Checkpointing (order unchanged from the old static checkpoint)
    template<class S>
    void operator& (S& stream) {
        infantDeaths & stream;
        infantIntervalsAtRisk & stream;
    }
    
    // Infant death summaries (checkpointed).
    vector<int> infantDeaths;
    vector<int> infantIntervalsAtRisk;
    /// Non-malaria mortality in under 1year olds. Set by init().
    double nonMalariaMortality = 0.0;
};

} }
#endif
