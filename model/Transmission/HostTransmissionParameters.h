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

#include "util/AgeGroupInterpolation.h"
#include "util/sampler.h"

namespace OM::Transmission {

struct HostSpeciesConfiguration {
    util::DistributionParameters availability;
    util::BetaParameters biting;
    util::BetaParameters findingRestSite;
    util::BetaParameters survivingRest;
};

struct HostSpeciesParameters {
    /** @brief Probabilities of finding a host and surviving a feeding cycle
     * 
     * These parameters describe the mean and heterogeneity of α_i, P_B_i,
     * P_C_i and P_D_i across the human population. */
    //@{
    /** Availability rate (α_i) */
    double entoAvailabilityFactor = 1.0;

    /** Availability heterogeneity factor
     * This multiplies the entoAvailabilityFactor. */
    unique_ptr<util::Sampler> entoAvailability;

    /** Probability of mosquito successfully biting host (P_B_i) */
    util::BetaSampler probMosqBiting;

    /** Probability of mosquito escaping human and finding a resting site without
     * dying, after biting the human (P_C_i). */
    util::BetaSampler probMosqFindRestSite;

    /** Probability of mosquito successfully resting after finding a resting site
     * (P_D_i). */
    util::BetaSampler probMosqSurvivalResting;
    //@}
    
};

struct HostTransmissionParameters {
    util::AgeGroupInterpolator availabilityByAge;
    std::vector<HostSpeciesParameters> species;
    std::vector<double> availabilityPercentiles;

    double availabilityPercentile(int p) const {
        if (p >= 0 && p < 100) return availabilityPercentiles[p];
        if (p == 100) return std::numeric_limits<double>::max();
        throw util::xml_scenario_error("Availability percentiles in intervention specifications can only be in [0 and 100]");
    }
};

HostSpeciesParameters makeHostSpeciesParameters(const HostSpeciesConfiguration&);
void calculateAvailabilityPercentiles(HostTransmissionParameters&);

}
