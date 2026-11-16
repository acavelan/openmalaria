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

#include "Transmission/HostTransmissionParameters.h"

namespace OM::Transmission {
HostSpeciesParameters makeHostSpeciesParameters(const HostSpeciesConfiguration& config) {
    HostSpeciesParameters result;
    result.entoAvailability = util::createSampler(1.0, config.availability);
    result.probMosqBiting.setParamsMV(config.biting.mean, config.biting.variance);
    result.probMosqFindRestSite.setParamsMV(config.findingRestSite.mean, config.findingRestSite.variance);
    result.probMosqSurvivalResting.setParamsMV(config.survivingRest.mean, config.survivingRest.variance);
    return result;
}

void calculateAvailabilityPercentiles(HostTransmissionParameters& parameters) {
    int nSamples = 100000;
    vector<double> samples(nSamples, 0.0);
    // Independent stream: percentile initialization must not affect simulation RNGs.
    util::MasterRng seed(0, 0);
    seed.seed(0, 0);
    util::LocalRng rng(seed);
    for (int i = 0; i < nSamples; ++i) {
        double avail = 0.0;
        for (const auto& species : parameters.species)
            avail += species.entoAvailabilityFactor * species.entoAvailability->sample(rng);
        samples[i] = avail;
    }
    sort(samples.begin(), samples.end());
    parameters.availabilityPercentiles.resize(100);
    for (int i = 0; i < 100; ++i)
        parameters.availabilityPercentiles[i] = samples[int(i*nSamples/100)];
    parameters.availabilityPercentiles[0] = 0.0;
}

}
