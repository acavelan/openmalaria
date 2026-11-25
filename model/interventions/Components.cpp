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
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "interventions/Components.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/CMDecisionTree.h"
#include "Host/Human.h"
#include "RunContext.h"
#include "interventions/IvState.h"
#include "mon/Monitoring.h"
#include "util/random.h"

namespace OpenMalaria {

void HumanIntervention::deploy(RunContext& ctx, Human& human, Deploy::Method method, VaccineLimits vaccLimits) const {
    for (size_t cond : conditions) {
        // Abort deployment if any condition is false.
        if (!checkCondition(ctx, cond))
            return;
    }

    for (ComponentId id : components) {
        const ComponentConfig& component = (*ctx.humanComponents)[id.id];
        // we must report first, since it can change cohort and sub-population
        // which may affect what deployment does (at least in the case of reporting deployments)
        human.addToCohort(id, component.duration);
        deploy_component(ctx, component, human, method, vaccLimits);
    }
}

void HumanIntervention::print_details(std::ostream& out) const {
    out << "human:";
    for (ComponentId id : components) {
        out << '\t' << id.id;
    }
}

void TriggeredDeployments::deploy(RunContext& ctx, Human& human, Deploy::Method method,
                                  VaccineLimits vaccLimits) const {
    for (auto it = lists.begin(), end = lists.end(); it != end; ++it) {
        it->deploy(ctx, human, method, vaccLimits);
    }
}

void TriggeredDeployments::SubList::deploy(RunContext& ctx, Human& human, Deploy::Method method,
                                           VaccineLimits vaccLimits) const {
    // This may be used within a time step; in that case we use human age at
    // the beginning of the step since this gives the expected result with age
    // limits of 0 to 1 time step.
    //TODO: why does the above contradict what we do? Is this due to a date being shifted later?
    SimTime age = human.age(ctx.clock.nowOrTs1());
    if (age >= minAge && age < maxAge) {
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if (coverage >= 1.0 || human.rng().bernoulli(coverage)) {
            HumanIntervention::deploy(ctx, human, method, vaccLimits);
        }
    }
}

void deploy_component(RunContext& ctx, const ComponentConfig& c, Human& human, Deploy::Method method,
                      VaccineLimits vaccLimits) {
    switch (c.kind) {
    case Component::RECRUIT_ONLY: // reports to monitoring, nothing else
        recordDeploy(measure("nMassRecruitOnly"), measure("nCtsRecruitOnly"), human, method);
        break;
    case Component::SCREEN: {
        const ScreenParams& p = std::get<ScreenParams>(c.data);
        recordDeploy(measure("nMassScreenings"), measure("nCtsScreenings"), human, method);
        if (human.withinHostModel->diagnosticResult(human.rng(), p.diagnostic))
            p.positive.deploy(ctx, human, method, vaccLimits);
        else
            p.negative.deploy(ctx, human, method, vaccLimits);
        break;
    }
    case Component::TREAT_SIMPLE: {
        const TreatSimpleParams& p = std::get<TreatSimpleParams>(c.data);
        recordDeploy(measure("nMDAs"), measure("nCtsMDA"), human, method);
        human.withinHostModel->treatSimple(human, p.timeLiver, p.timeBlood);
        break;
    }
    case Component::TREAT_PKPD: {
        const TreatPKPDParams& p = std::get<TreatPKPDParams>(c.data);
        recordDeploy(measure("nMDAs"), measure("nCtsMDA"), human, method);
        double age = inYears(human.age(ctx.clock.nowOrTs1()));
        human.withinHostModel->treatPkPd(p.schedule, p.dosage, age, p.delay_d);
        break;
    }
    case Component::CM_DT: {
        const DecisionTreeParams& p = std::get<DecisionTreeParams>(c.data);
        DTExecCtx ec(ctx, human, inYears(human.age(ctx.clock.nowOrTs1())), human.clinicalModel->latestState);
        CMDTOut out = exec_decision_tree(p.tree, ec);
        if (out.treated) {
            recordDeploy(measure("nMDAs"), measure("nCtsMDA"), human, method);
        }
        if (out.screened) {
            recordDeploy(measure("nMassScreenings"), measure("nCtsScreenings"), human, method);
        }
        break;
    }
    case Component::PEV:
    case Component::BSV:
    case Component::TBV:
        if (vaccinate(ctx, human, c.id, vaccLimits) && ctx.vaccineReportComponent == c.id) {
            recordDeploy(measure("nMassVaccinations"), measure("nEPIVaccinations"), human, method);
        }
        if (method == Deploy::TREAT)
            recordEvent(measure("nTreatDeployments"), human);
        break;
    case Component::ITN:
        iv_deploy(ctx, human, c);
        recordDeploy(measure("nMassITNs"), measure("nEPI_ITNs"), human, method);
        break;
    case Component::IRS:
        iv_deploy(ctx, human, c);
        recordDeploy(measure("nMassIRS"), measure("nCtsIRS"), human, method);
        break;
    case Component::GVI:
        iv_deploy(ctx, human, c);
        recordDeploy(measure("nMassGVI"), measure("nCtsGVI"), human, method);
        break;
    case Component::CLEAR_IMMUNITY:
        human.withinHostModel->clearImmunity();
        //NOTE: not reported
        break;
    default: // Component::IPT: no implementation exists
        throw SWITCH_DEFAULT_EXCEPTION;
    }
}

void print_component(std::ostream& out, const ComponentConfig& c) {
    // Indexed by Component::Type (every kind assignable to c.kind is a member
    // of that enum). IPT entry reproduces the former ternary's default:
    // "clear immunity" (Component::IPT has no implementation; the old chain
    // fell through to its last branch for it too).
    static const char* const NAMES[] = {
        "Recruit only", "Screen", "treatSimple", "treatPKPD",      "Decision tree",
        "PEV",          "BSV",    "TBV",         "clear immunity", // IPT: former default fall-through
        "ITN",          "IRS",    "GVI",         "clear immunity"};
    out << c.id.id << '\t' << NAMES[c.kind];
}

}
