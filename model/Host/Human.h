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
#include "Host/HumanTransmission.h"
#include "InfectionIncidenceModel.h"
#include "interventions/HumanComponents.h"

namespace OM {
namespace Clinical { class ClinicalModel; }
namespace WithinHost { class WHInterface; }
struct ModelParameters;

namespace Host {
struct Human {
    Human(SimTime birthDate, const ModelParameters&, util::MasterRng&);

    HumanTransmission transmission;
    unique_ptr<WithinHost::WHInterface> withinHostModel;
    IncidenceState infIncidence;
    unique_ptr<Clinical::ClinicalModel> clinicalModel;
    util::LocalRng rng;
    interventions::HumanInterventions interventions;
    size_t monitoringAgeGroup = 0;
    SimTime birthDate = sim::never();
    bool dead = false;

    SimTime age(SimTime date) const { return date - birthDate; }
    bool isInSubPop(interventions::ComponentId id) const {
        auto it = interventions.cohortExpiry.find(id);
        return it != interventions.cohortExpiry.end() && it->second >= sim::nowOrTs0();
    }
    void addToCohort(interventions::ComponentId id, SimTime duration);
    void removeFirstEvent(const vector<interventions::ComponentId>& rules, SimTime date);
    void updateCohortSet();
    void checkpoint(istream&, const vector<unique_ptr<interventions::HumanInterventionComponent>>&);
    void checkpoint(ostream&);
};
} }
