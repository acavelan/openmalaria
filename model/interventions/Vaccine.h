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

#include "interventions/HumanComponents.h"

namespace scnXml { class VaccineDescription; }
namespace OM::interventions {
struct VaccineComponent : HumanInterventionComponent {
    VaccineComponent(ComponentId id, Vaccine::Types type, bool report = false)
        : HumanInterventionComponent(id), type(type), reportDeployments(report) {}
    VaccineComponent(ComponentId, const scnXml::VaccineDescription&, Vaccine::Types, bool report);

    void deploy(Host::Human&, mon::Deploy::Method, VaccineLimits) const override;
    Component::Type componentType() const override;
    void print_details(ostream&) const override;
    double getInitialEfficacy(LocalRng&, size_t numPrevDoses, uint32_t genotype) const;

    Vaccine::Types type;
    util::DecayFunction decayFunc;
    vector<vector<double>> perGenotypeInitialMeanEfficacy;
    double efficacyB = 1.0;
    // Preserve reporting of only the first vaccine configured in this model.
    bool reportDeployments;
};
}
