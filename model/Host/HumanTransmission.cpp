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

#include "Host/HumanTransmission.h"
#include "interventions/HumanComponents.h"

namespace OM::Host {
void initialize(HumanTransmission& host, util::LocalRng& rng, double availabilityFactor,
        const Transmission::HostTransmissionParameters& parameters) {
    host.relativeAvailability = availabilityFactor;
    host.species.resize(parameters.species.size());
    for (size_t i = 0; i < host.species.size(); ++i) {
        const auto& base = parameters.species[i];
        auto& species = host.species[i];
        species.sampledAvailability = base.entoAvailability->sample(rng);
        species.availability = base.entoAvailabilityFactor * species.sampledAvailability * availabilityFactor;
        species.biteProbability = base.probMosqBiting.sample(rng);
        auto pRest1 = base.probMosqFindRestSite.sample(rng);
        auto pRest2 = base.probMosqSurvivalResting.sample(rng);
        species.restingProbability = pRest1 * pRest2;
    }
}

// Checkpoints retain the four original vectors in field-major order.
static constexpr auto fields = {&MosquitoInteraction::sampledAvailability,
    &MosquitoInteraction::availability, &MosquitoInteraction::biteProbability,
    &MosquitoInteraction::restingProbability};
void checkpoint(HumanTransmission& host, ostream& stream) {
    for (auto field : fields) {
        host.species.size() & stream;
        for (auto& species : host.species) species.*field & stream;
    }
    host.relativeAvailability & stream;
    host.outsideTransmission & stream;
}
void checkpoint(HumanTransmission& host, istream& stream) {
    for (auto field : fields) {
        size_t count;
        count & stream;
        util::checkpoint::validateListSize(count);
        if (field == *fields.begin()) host.species.resize(count);
        else if (count != host.species.size()) throw util::checkpoint_error("inconsistent mosquito species counts");
        for (auto& species : host.species) species.*field & stream;
    }
    host.relativeAvailability & stream;
    host.outsideTransmission & stream;
}

double relativeAvailabilityAge(const HumanTransmission& host, double ageYears,
        const Transmission::HostTransmissionParameters& parameters) {
    return host.outsideTransmission ? 0.0 : parameters.availabilityByAge.eval(ageYears);
}
double relativeAvailabilityHetAge(const HumanTransmission& host, double ageYears,
        const Transmission::HostTransmissionParameters& parameters) {
    return host.relativeAvailability * relativeAvailabilityAge(host, ageYears, parameters);
}
double availability(const HumanTransmission& host, const interventions::HumanInterventions& state, size_t species, SimTime effectTime) {
    double value = host.species[species].availability;
    for (const auto& component : state.vectorInterventions)
        value *= component->relativeAttractiveness(species, effectTime);
    return value;
}
double bitingProbability(const HumanTransmission& host, const interventions::HumanInterventions& state, size_t species, SimTime effectTime) {
    double value = host.species[species].biteProbability;
    for (const auto& component : state.vectorInterventions)
        value *= component->preprandialSurvivalFactor(species, effectTime);
    return value;
}
double restingProbability(const HumanTransmission& host, const interventions::HumanInterventions& state, size_t species, SimTime effectTime) {
    double value = host.species[species].restingProbability;
    for (const auto& component : state.vectorInterventions)
        value *= component->postprandialSurvivalFactor(species, effectTime);
    return value;
}
double fecundity(const interventions::HumanInterventions& state, size_t species, SimTime effectTime) {
    double value = 1.0;
    for (const auto& component : state.vectorInterventions)
        value *= component->relFecundity(species, effectTime);
    return value;
}
double availabilityToBite(const HumanTransmission& host, const interventions::HumanInterventions& state,
        size_t species, SimTime effectTime) {
    return availability(host, state, species, effectTime) * bitingProbability(host, state, species, effectTime);
}
double availabilityWithAge(const HumanTransmission& host, const interventions::HumanInterventions& state,
        size_t species, double ageYears, const Transmission::HostTransmissionParameters& parameters, SimTime effectTime) {
    return availability(host, state, species, effectTime) * relativeAvailabilityAge(host, ageYears, parameters);
}
}
