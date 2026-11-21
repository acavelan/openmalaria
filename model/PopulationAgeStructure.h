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

#include <vector>
#include "Clock.h"

namespace scnXml {
class Demography;
}
namespace OpenMalaria {
class RunContext;

/** Per-run demography state (was: AgeStructure statics; lives in
     * RunContext::age). */
struct AgeStructureState {
    /*! Field-data age groups: ageGroupBounds[i] is the lower bound of
        group i (bounds has ngroups+1 entries), ageGroupPercent[i] the
        percentage of the population in group i. Set by
        estimateRemovalRates(), read by the fit objective (file-local in
        PopulationAgeStructure.cpp). */
    std::vector<double> ageGroupBounds, ageGroupPercent;
    /// Parameters defining the smooth target age-distribution curve.
    double mu0 = 0.0, mu1 = 0.0, alpha0 = 0.0, alpha1 = 0.0;
    /// rho is growth rate of human population.
    static constexpr double rho = 0.0;
    /** Target cumulative percentage of population by age, from oldest
        age to youngest: cumAgeProp[size()+1-i] is the proportion of people
        aged i time steps or older. Set by init()/calcCumAgeProp(). */
    std::vector<double> cumAgeProp;
};

/** Encapsulates code just setting up the age structure (i.e.  cumAgeProp). */
class AgeStructure {
  public:
    /** Set up cumAgeProp from XML data. */
    static void init(RunContext& ctx, const scnXml::Demography& demography);

    /** Return maximum individual lifetime in intervals that AgeStructure can handle. */
    static size_t getMaxTStepsPerLife(const RunContext& ctx);

    /** Return the expected population size of individuals aged ageTSteps or
        * older, based on a total population size of targetPop. */
    static int targetCumPop(RunContext& ctx, size_t ageTSteps, int targetPop);

  private:
    static void estimateRemovalRates(RunContext& ctx, const scnXml::Demography& demography);

    /** Takes the best-fitting demography parameters estimated by
        * estimateRemovalRates and calculates the age structure (cumAgeProp). */
    static void calcCumAgeProp(RunContext& ctx);
};
}
