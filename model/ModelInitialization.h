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
#include "ModelParameters.h"

namespace scnXml { class Scenario; class Mosq; }
namespace OM {
struct Model;
class Parameters;
void initialize(Model& model, scnXml::Scenario& scenario);
Host::IncidenceParameters readIncidenceParameters(const Parameters& parameters);
Transmission::HostSpeciesParameters readHostSpeciesParameters(const scnXml::Mosq& mosq);
void calculateAvailabilityPercentiles(Transmission::HostTransmissionParameters& parameters);
}
