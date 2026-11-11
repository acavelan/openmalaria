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
#include "Global.h"
#include "Transmission/HostTransmissionParameters.h"

namespace OM::interventions { struct HumanInterventions; }
namespace OM::Host {
struct MosquitoInteraction {
    double sampledAvailability;
    double availability;
    double biteProbability;
    double restingProbability;
};

struct HumanTransmission {
    bool outsideTransmission = false;
    double relativeAvailability = numeric_limits<double>::signaling_NaN();
    vector<MosquitoInteraction> species;
};

void initialize(HumanTransmission&, util::LocalRng&, double availabilityFactor,
    const Transmission::HostTransmissionParameters&);
void checkpoint(HumanTransmission&, ostream&);
void checkpoint(HumanTransmission&, istream&);

double relativeAvailabilityAge(const HumanTransmission&, double ageYears,
    const Transmission::HostTransmissionParameters&);
double relativeAvailabilityHetAge(const HumanTransmission&, double ageYears,
    const Transmission::HostTransmissionParameters&);
double availability(const HumanTransmission&, const interventions::HumanInterventions&, size_t species, SimTime effectTime);
double bitingProbability(const HumanTransmission&, const interventions::HumanInterventions&, size_t species, SimTime effectTime);
double restingProbability(const HumanTransmission&, const interventions::HumanInterventions&, size_t species, SimTime effectTime);
double fecundity(const interventions::HumanInterventions&, size_t species, SimTime effectTime);
double availabilityToBite(const HumanTransmission&, const interventions::HumanInterventions&, size_t species, SimTime effectTime);
double availabilityWithAge(const HumanTransmission&, const interventions::HumanInterventions&, size_t species,
    double ageYears, const Transmission::HostTransmissionParameters&, SimTime effectTime);
}
