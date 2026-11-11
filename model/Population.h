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
#include "Host/Human.h"

namespace OM {
struct Population {
    size_t targetSize = 0;
    int recentBirths = 0;
    vector<Host::Human> humans;
};

void createInitialHumans(Population&, const ModelParameters&, util::MasterRng&);
void update(Population&, const ModelParameters&, util::MasterRng&);
void checkpoint(Population&, istream&, const ModelParameters&, util::MasterRng&,
    const vector<unique_ptr<interventions::HumanInterventionComponent>>& definitions);
void checkpoint(Population&, ostream&);

void registerContinousPopulationCallbacks(const Transmission::HostTransmissionParameters& parameters);

/// Delegate to print the number of hosts
void ctsHosts (Population &population, ostream& stream);
/// Delegate to print cumulative numbers of hosts under various age limits
void ctsHostDemography (Population &population, ostream& stream);
/// Delegate to print the number of births since last count
void ctsRecentBirths (Population &population, ostream& stream);
/// Delegate to print the number of patent hosts
void ctsPatentHosts (Population &population, ostream& stream);
/// Delegate to print immunity's cumulativeh parameter
void ctsImmunityh (Population &population, ostream& stream);
/// Delegate to print immunity's cumulativeY parameter (mean across population)
void ctsImmunityY (Population &population, ostream& stream);
/// Delegate to print immunity's cumulativeY parameter (median across population)
void ctsMedianImmunityY (Population &population, ostream& stream);
/// Delegate to print the mean age-based availability reduction of each human relative to an adult
void ctsMeanAgeAvailEffect (Population &population, ostream& stream, const Transmission::HostTransmissionParameters& parameters);
void ctsITNCoverage (Population &population, ostream& stream);
void ctsIRSCoverage (Population &population, ostream& stream);
void ctsGVICoverage (Population &population, ostream& stream);

}
