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

// Phase 7.2b: flat human-intervention component config. The
// HumanInterventionComponent class tree is replaced by ComponentConfig: an id, a
// Component::Type kind and a variant payload (plan §5 tag+switch idiom).
// Deployment and detail-printing are the deploy_component/print_component
// switches (model/xml/HumanInterventionComponents.h, one-TU pattern).

#ifndef Hmod_interventions_Components
#define Hmod_interventions_Components

#include "interventions/Interfaces.h"
#include "interventions/Vaccine.h"
#include "interventions/ITN.h"
#include "interventions/IRS.h"
#include "interventions/GVI.h"
#include <variant>

namespace scnXml {
    class Screen;
    class DTTreatSimple;
    class DTTreatPKPD;
    class DecisionTree;
}

namespace OM {
namespace WithinHost { class Diagnostic; }
namespace Clinical { class FlatTree; }
namespace interventions {

/// No state: reports deployment, nothing else. Also the payload for
/// CLEAR_IMMUNITY (its deploy reads no parameters either).
struct RecruitOnlyParams {};

/// Screening: deploys one of two nested intervention lists.
struct ScreenParams {
    const WithinHost::Diagnostic& diagnostic;
    HumanIntervention positive, negative;
    ScreenParams( RunContext&, const scnXml::Screen& );
};

/// Simple treatment: no PK/PD, just remove parasites.
struct TreatSimpleParams {
    SimTime timeLiver, timeBlood;
    TreatSimpleParams( RunContext&, const scnXml::DTTreatSimple& );
};

/// Treatment via the PK/PD model (single "treatment", cf. CMDTTreatPKPD).
struct TreatPKPDParams {
    size_t schedule;    ///< index of the schedule
    size_t dosage;      ///< index of the dosage table
    double delay_d;     ///< delay in days
    TreatPKPDParams( RunContext&, const scnXml::DTTreatPKPD& );
};

/// Execute a "case management" decision tree.
struct DecisionTreeParams {
    const Clinical::FlatTree& tree;
    DecisionTreeParams( RunContext&, const scnXml::DecisionTree& );
};

/** One human intervention component (plan: dense, by type-rank).
 * Owned by RunExtensions::humanComponents, indexed by id.id; never removed,
 * so payload pointers (the per-type registries) stay valid. */
struct ComponentConfig {
    ComponentId id;
    Component::Type kind;
    SimTime duration = sim::never();    ///< sub-pop removal expireAfter
    std::variant<RecruitOnlyParams, ScreenParams, TreatSimpleParams,
        TreatPKPDParams, DecisionTreeParams, VaccineParams, ITNParams,
        IRSParams, GVIParams> data;};

/// Deploy component c to human (replaces the old virtual deploy()).
/// Draws/mon reports happen in the same order as the deleted class bodies.
void deploy_component( RunContext& ctx, const ComponentConfig& c, Host::Human& human,
    mon::Deploy::Method method, VaccineLimits vaccLimits );

/// One-line component description (PRINT_INTERVENTIONS; old print_details()).
void print_component( std::ostream& out, const ComponentConfig& c );

} }
#endif
