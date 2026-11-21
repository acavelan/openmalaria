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

#include "Transmission/TransmissionModel.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/Anopheles/EmergenceRateEstimator.h"
#include "util/ModelOptions.h"

#include <memory>
#include <vector>

namespace scnXml {
class Entomology;
class AnophelesParams;
}

namespace OM {
namespace Transmission {
using namespace OM::util;

inline SimulationMode readMode(const string& str) {
    if (str == "forced")
        return forcedEIR;
    else if (str == "dynamic")
        return dynamicEIR;
    else
        // Note: originally 3 (transientEIRknown) could be specified; now it's
        // set automatically.
        throw util::xml_scenario_error(string("mode attribute invalid: ").append(str));
}

// Factory functions are defined in model/xml/entomology.cpp.
TransmissionModel* createNonVectorModel(RunContext& ctx, const scnXml::Entomology& entoData);
bool anophelesCompare(const scnXml::AnophelesParams& a1, const scnXml::AnophelesParams& a2);
Anopheles::AnophelesModel* createAnophelesModel(RunContext& ctx, size_t i, const scnXml::AnophelesParams& anoph,
                                                vector<double>& initialisationEIR, int populationSize,
                                                int interventionMode);
TransmissionModel* createVectorModel(RunContext& ctx, const scnXml::Entomology& entoData, int populationSize);
TransmissionModel* createTransmissionModel(RunContext& ctx, const scnXml::Entomology& entoData, int populationSize);
} // namespace Transmission
} // namespace OM
