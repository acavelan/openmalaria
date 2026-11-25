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
#include <exception>
#include <string>
#include <memory>

namespace OpenMalaria {
class RunContext;
class Population;
class TransmissionModel;
inline constexpr int SCHEMA_VERSION = 50;

/// The engine receives runtime objects, never an XSD tree. Members are ordered
/// to destroy transmission before population, as in the original driver.
struct ScenarioRuntime {
    std::unique_ptr<Population> population;
    std::unique_ptr<TransmissionModel> transmission;
    bool surveyOnlyNewEp;
    SimTime humanWarmupLength, estEndTime;
};

/// Parse, validate and initialise in the frozen dependency/RNG order.
ScenarioRuntime initialiseScenario(RunContext& ctx, const std::string& path, bool resuming);

/// Report an XSD exception, if this is one. Keeps library types and their
/// diagnostics inside xml/; main retains its existing exit-code handling.
bool reportXmlException(std::exception_ptr error);

} // namespace OpenMalaria
