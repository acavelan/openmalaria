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
#include "util/random.h"
#include "Host/Human.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/WHInterface.h"

#include <cmath>
#include <vector>

namespace scnXml {
    class Clinical;
}
namespace OM {
    class Population;
namespace Host {

using util::LocalRng;

class NeonatalMortality {
public:
  /// Initialisation
  static void init( RunContext& ctx, const scnXml::Clinical& clinical );
  
  /** Called for each birth; returns true if infant dies due to mother's
   * infection. */
  static bool eventNeonatalMortality(LocalRng& rng);
  
  /** Calculate risk of a neonatal mortality based on humans 20-25 years old. */
  static void update (RunContext& ctx, vector<Human>& population);
  
};

// ———  run state (init in xml/healthSystem.cpp; wire in checkpoint/statics.cpp;
// both reference these by extern — inline keeps external linkage)  ———

//Goodman estimated for neonatal mortality due to malaria in pregnancy
const double gEst = 0.011;
//Critical value of Prev20-25 for neonatal mortality
const double critPrev2025 = 0.25;
//Critical value for estimating prevalence in primigravidae
const double critPrevPrim = 0.19;
//Proportion of births with primigravid mothers
const double pBirthPrim = 0.3;

// optimised constants:
const double y = pBirthPrim * gEst;
const double z = -1.0 / critPrev2025;

/// Probability for a newborn to die (indirect death) because the mother is infected.
/// Depends on the prevalence of parasitaemia in mother at some previous t.
inline double riskFromMaternalInfection = 0.0;
/// Array of stored prevalences of mothers over last 5 months
inline std::vector<double> prevByGestationalAge;

/// Lower and upper bounds for potential mothers (as in model description)
inline SimTime ageLb = sim::fromYearsI(20), ageUb = sim::fromYearsI(25);

// The model is parameterised based on patency levels; the diagnostic
// used for this may be important.
inline const WithinHost::Diagnostic* neonatalDiagnostic = 0;

inline bool NeonatalMortality::eventNeonatalMortality(LocalRng& rng) {
  // parity: order — local RNG stream: draw call order fixes stream state
  return rng.uniform_01() <= riskFromMaternalInfection;
}

inline void NeonatalMortality::update (RunContext& ctx, vector<Host::Human> &population) {
    // ———  find potential mothers and their prevalence  ———
    // For individuals in the age range 20-25, we sum:
    int nCounter=0;	// total number
    int pCounter=0;	// number with patent infections, needed for prev in 20-25y
    
    for (Human& human : population){
        // diagnosticDefault() gives patency after the last time step's
        // update, so it's appropriate to use age at the beginning of this step.
        SimTime age = human.age(ctx.clock.ts0());
        
        // Note: since we're using a linked list, we have to iterate until we reach
        // the individuals we're interested in. Due to population structure, it's
        // probably quickest to start iterating from the oldest.
        if( age >= ageUb ) continue;
        if( age < ageLb ) break;	// Not interested in younger individuals.
        
        nCounter ++;
        if( human.withinHostModel->diagnosticResult(human.rng(), *neonatalDiagnostic) ){
            pCounter ++;
        }
    }
    
    // ———  calculate risk of neonatal mortality  ———
    //default value of prev2025, for use when there are no 20-25 year olds
    double prev2025 = 0.25;
    if( nCounter > 0 )
        prev2025 = double(pCounter) / nCounter;
    
    double maxPrev = prev2025;
    //update the vector containing the prevalence by gestational age
    size_t index = ctx.clock.moduloSteps(ctx.clock.ts0(), prevByGestationalAge.size());
    prevByGestationalAge[index] = prev2025;
    // max is order-free (plain > comparisons, no FP combine) — punch-list #12
    maxPrev = std::max(maxPrev, *std::max_element(prevByGestationalAge.begin(), prevByGestationalAge.end()));
    
    // equation (2) p 75 AJTMH 75 suppl 2
    double prevPG= maxPrev / (critPrevPrim + maxPrev);
    // equation (1) p 75 AJTMH 75 suppl 2 including 30% multiplier
    riskFromMaternalInfection = y * (1.0-exp(prevPG * z));
}

} }
