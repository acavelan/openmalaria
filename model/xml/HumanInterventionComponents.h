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

// XML parsing and payload construction only. Runtime deployment and printing
// live in interventions/Components.cpp. Included by xml/interventions.cpp.

// The includes here are more for documentation than required.
#include "RunContext.h"
#include "interventions/Components.h"
#include "interventions/InterventionManager.h"
#include "interventions/IvState.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "Clinical/EventScheduler.h"
#include "Clinical/CM5DayCommon.h"
#include "Host/Human.h"
#include "Transmission/TransmissionModel.h"
#include "Host/WithinHost/Diagnostic.h"
#include "PkPd/LSTMMedicate.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include <schema/healthSystem.h>
#include <schema/interventions.h>

namespace OpenMalaria {

// ———  HumanIntervention  ———

// Components sort into standard type-rank order (contract b).
static void sort_components(RunContext& ctx, std::vector<ComponentId>& components) {
    auto kind = [&ctx](ComponentId id) { return (*ctx.humanComponents)[id.id].kind; };
    /** Sort components according to a standard order.
     *
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     *
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort(components.begin(), components.end(),
                     [&kind](ComponentId a, ComponentId b) { return kind(a) < kind(b); });
}
static ComponentId component_id(RunContext& ctx, const std::string& textId) {
    ComponentId id = InterventionManager::getComponentId(ctx, textId);
    if (id.id >= ctx.humanComponents->size())
        throw base_exception("invalid component id"); // used before defined
    return id;
}
HumanIntervention::HumanIntervention(RunContext& ctx,
                                     const xsd::cxx::tree::sequence<scnXml::Component>& componentList) {
    components.reserve(componentList.size());
    for (auto it = componentList.begin(), end = componentList.end(); it != end; ++it)
        components.push_back(component_id(ctx, it->getId()));
    sort_components(ctx, components);
}
HumanIntervention::HumanIntervention(RunContext& ctx, const xsd::cxx::tree::sequence<scnXml::Component>& componentList,
                                     const xsd::cxx::tree::sequence<scnXml::Condition>& conditionList)
    : HumanIntervention(ctx, componentList) {
    for (auto it = conditionList.begin(), end = conditionList.end(); it != end; ++it) {
        double min = it->getMinValue().present() ? it->getMinValue().get() : -numeric_limits<double>::infinity();
        double max = it->getMaxValue().present() ? it->getMaxValue().get() : numeric_limits<double>::infinity();
        conditions.push_back(setupCondition(ctx, it->getMeasure(), min, max, it->getInitialState()));
    }
}
HumanIntervention::HumanIntervention(RunContext& ctx, const xsd::cxx::tree::sequence<scnXml::DTDeploy>& componentList) {
    components.reserve(componentList.size());
    for (auto it = componentList.begin(), end = componentList.end(); it != end; ++it)
        components.push_back(component_id(ctx, it->getComponent()));
    sort_components(ctx, components);
}

// ———  Utilities  ———

TriggeredDeployments::TriggeredDeployments(RunContext& ctx, const scnXml::TriggeredDeployments& elt) {
    lists.reserve(elt.getDeploy().size());
    for (auto it = elt.getDeploy().begin(), end = elt.getDeploy().end(); it != end; ++it) {
        lists.push_back(SubList(ctx, *it));
    }
}
TriggeredDeployments::SubList::SubList(RunContext& ctx, const scnXml::Deploy& elt)
    : HumanIntervention(ctx, elt.getComponent()), minAge(ctx.clock.from_years_n(elt.getMinAge())), maxAge(future()),
      coverage(elt.getP()) {
    if (elt.getMaxAge().present()) maxAge = ctx.clock.from_years_n(elt.getMaxAge().get());

    if (minAge < zero() || maxAge < minAge) {
        throw xml_scenario_error("triggered intervention must have 0 <= minAge <= maxAge");
    }

    if (coverage < 0.0 || coverage > 1.0) {
        throw xml_scenario_error("triggered intervention must have 0 <= coverage <= 1");
    }

    // Zero coverage: optimise
    if (coverage <= 0.0 || minAge >= maxAge) components.clear();
}
// ———  Component payload constructors  ———

ScreenParams::ScreenParams(RunContext& ctx, const scnXml::Screen& elt)
    : diagnostic(diagnostics::get(ctx, elt.getDiagnostic())), positive(ctx, elt.getPositive()),
      negative(ctx, elt.getNegative()) {}

/// Simple treatment: no PK/PD, just remove parasites
TreatSimpleParams::TreatSimpleParams(RunContext& ctx, const scnXml::DTTreatSimple& elt) {
    //NOTE: this code is currently identical to that in CMDTTreatSimple
    try {
        SimTime durL = UnitParse::readShortDuration(ctx.clock, elt.getDurationLiver(), UnitParse::NONE),
                durB = UnitParse::readShortDuration(ctx.clock, elt.getDurationBlood(), UnitParse::NONE);
        SimTime neg1 = -ctx.clock.one_TS();
        if (durL < neg1 || durB < neg1) {
            throw xml_scenario_error("treatSimple: cannot have durationBlood or durationLiver less than -1");
        }
        if (ModelOptions::option(VIVAX_SIMPLE_MODEL)) {
            if (durL != zero() || durB != neg1)
                throw unimplemented_exception("vivax model only supports timestepsLiver=0, timestepsBlood=-1");
            // Actually, the model ignores these parameters; we just don't want somebody thinking it doesn't.
        }
        timeLiver = durL;
        timeBlood = durB;
    } catch (const format_error& e) {
        throw xml_scenario_error(string("treatSimple: ").append(e.message()));
    }
}

//NOTE: this is similar to CMDTTreatPKPD, but supporting only a single "treatment"
TreatPKPDParams::TreatPKPDParams(RunContext& ctx, const scnXml::DTTreatPKPD& elt)
    : schedule(ctx.treatments.findSchedule(elt.getSchedule())), dosage(ctx.treatments.findDosages(elt.getDosage())),
      delay_d(elt.getDelay_h() / 24.0 /* convert to days */) {}

DecisionTreeParams::DecisionTreeParams(RunContext& ctx, const scnXml::DecisionTree& elt)
    : tree(compile_decision_tree(ctx, elt, false)) {}

}
