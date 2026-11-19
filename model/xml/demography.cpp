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

// model/xml/ conversion layer — demography module: AgeStructure removal-rate
// estimation from scnXml::Demography. Moved verbatim from
// PopulationAgeStructure.cpp; declaration stays in PopulationAgeStructure.h.

#include "Global.h"
#include "PopulationAgeStructure.h"
#include "util/errors.h"
#include <schema/demography.h>

#include <cmath>

namespace OM {
    using namespace OM::util;

    // Free function defined in PopulationAgeStructure.cpp (external linkage):
    void minimizeCalc_rss (double (*func) (double,double), double param1, double param2);

void AgeStructure::estimateRemovalRates( RunContext& ctx, const scnXml::Demography& demography ){
    // mu1, alpha1: These are estimated here
    // mu0, alpha0: These are fixed alpha0=4.0, mu0 calculated from the other parameters
    // rho - population growth rate (input)
    
    //Get lower and upper age bounds for age groups and cumulative precentage of population from field data
    double sumperc = 0.0;
    const scnXml::DemogAgeGroup::GroupSequence& group = demography.getAgeGroup().getGroup();
    
    size_t ngroups = group.size() + 1;
    ageGroupBounds.resize(ngroups+1,0);
    ageGroupPercent.resize(ngroups,0);
    
    //Add age group for first month of life
    ageGroupBounds[0] = 0.0;
    ageGroupBounds[1] = 1.0 / 12.0;
    ageGroupPercent[0] = 0.0;
    for(size_t i = 1;i < ngroups; i++) {
	ageGroupBounds[i+1] = group[i-1].getUpperbound();
	ageGroupPercent[i] = group[i-1].getPoppercent();
	sumperc += ageGroupPercent[i];
    }
    sumperc = 100.0 / sumperc; // multiplier to get percentages
    for(size_t i = 0;i < ngroups; i++) {
	ageGroupPercent[i]  = ageGroupPercent[i] * sumperc;
    }

    /* NOTE: unused --- why?
    double tol = 0.00000000001;
    int maxcal = 500000;
    int npar = 2;
    int iw = 3; */
    double p1 = 0.371626412;
    double p2 = 0.841209593;
    // returns "double rss":
    minimizeCalc_rss (&setDemoParameters, p1, p2);
}
} // namespace OM
