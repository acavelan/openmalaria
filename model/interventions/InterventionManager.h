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
#include "util/errors.h"
#include "util/CommandLine.h"
#include "Population.h"
#include "Host/ImportedInfections.h"
#include "interventions/Deployments.h"
#include "interventions/Interfaces.h"
#include "Transmission/TransmissionModel.h"
namespace scnXml {
class Interventions;
}

namespace OpenMalaria {
class Population;
class Human;

/** Management of interventions deployed on a per-time-step basis. */
class InterventionManager {
  public:
    /** Read XML descriptions. */
    static void init(RunContext& ctx, const scnXml::Interventions& intervElt, TransmissionModel& transmission);

    /** Call after loading a checkpoint, passing the intervention-period time.
     * 
     * Serves to replace health-system and EIR where changeHS/changeEIR
     * interventions have been used. */

    /** @brief Deploy interventions
     *
     * Timed interventions are deployed for this time step.
     * 
     * Continuous interventions are deployed as humans reach the target ages.
     * Unlike with vaccines, missing one schedule doesn't preclude the next. */
    static void deploy(RunContext& ctx, vector<Human>& population, TransmissionModel& transmission);

    /** Get a numeric ComponentId from the textual identifier used in the XML.
     * 
     * If textId is unknown, an xml_scenario_error is thrown. */
    static ComponentId getComponentId(RunContext& ctx, const std::string textId);
};

// ———  implementation (was InterventionManager.cpp)  ———

// ———  InterventionManager  ———

// static functions:

inline ComponentId InterventionManager::getComponentId(RunContext& ctx, const string textId) {
    auto it = ctx.identifierMap.find(textId);
    if (it == ctx.identifierMap.end()) {
        ostringstream msg;
        msg << "unable to find an intervention component with id \"" << textId
            << "\" (wrong name, no definition or used before definition?)";
        throw xml_scenario_error(msg.str());
    }
    return it->second;
}

inline void InterventionManager::deploy(RunContext& ctx, vector<Human>& population, TransmissionModel& transmission) {
    if (ctx.clock.interv < zero())
        return;

    // deploy imported infections (not strictly speaking an intervention)
    ctx.importedInfections->import(ctx, population);

    // deploy timed interventions
    SimTime now = ctx.clock.interv_date();
    while ((*ctx.timed)[ctx.nextTimed].date <= now) {
        if (CommandLine::option(CommandLine::VERBOSE)) {
            cout << "deploy::" << now << "::Timed";
            print_timed_details(cout, (*ctx.timed)[ctx.nextTimed]);
            cout << endl;
        }

        (*ctx.timed)[ctx.nextTimed].deploy(ctx, population, transmission);
        ctx.nextTimed += 1;
    }

    // deploy continuous interventions
    for (Human& human : population) {
        uint32_t nextCtsDist = human.nextCtsDist();
        // deploy continuous interventions
        while (nextCtsDist < ctx.continuous->size()) {
            if (CommandLine::option(CommandLine::VERBOSE)) {
                cout << "deploy::" << now << "::Continuous";
                (*ctx.continuous)[nextCtsDist].print_details(cout);
                cout << endl;
            }
            if (!(*ctx.continuous)[nextCtsDist].filterAndDeploy(ctx, human))
                break; // deployment (and all remaining) happens in the future
            nextCtsDist = ++human.nextCtsDist();
        }
    }
}

}
