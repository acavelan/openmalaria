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
#include "Population.h"

namespace OM {
namespace Transmission { class TransmissionModel; }

struct Reporting {
    int newInfections = 0;
};

// Owns the migrated subsystems. Time, interventions and other reporting state
// still use their existing APIs until their separate migrations are complete.
struct Model {
    explicit Model(size_t populationSize);
    ~Model();
    Model(const Model&) = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) = delete;
    Model& operator=(Model&&) = delete;

    ModelParameters parameters;
    util::MasterRng masterRng{0, 0};
    Population population;
    std::unique_ptr<Transmission::TransmissionModel> transmission;
    Reporting reporting;

    SimTime humanWarmupLength = 0;
    bool surveyOnlyNewEp = false;
};

void update(Model& model);
}
