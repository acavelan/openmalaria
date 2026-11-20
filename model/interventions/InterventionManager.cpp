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

#include "core/RunContext.h"
#include "core/RunExtensions.h"
#include "interventions/InterventionManager.h"
#include "Population.h"
#include "util/CommandLine.h"
#include "interventions/Deployments.h"

namespace OM
{
namespace interventions
{

// ———  InterventionManager  ———

// static functions:

void InterventionManager::checkpoint (RunContext& ctx, std::istream& stream) {
    using namespace OM::util::checkpoint;
    ctx.ext->importedInfections & stream;
}
void InterventionManager::checkpoint (RunContext& ctx, std::ostream& stream) {
    using namespace OM::util::checkpoint;
    ctx.ext->importedInfections & stream;
}

ComponentId InterventionManager::getComponentId(RunContext& ctx, const string textId)
{
    auto it = ctx.ext->identifierMap.find(textId);
    if (it == ctx.ext->identifierMap.end())
    {
        ostringstream msg;
        msg << "unable to find an intervention component with id \"" << textId
            << "\" (wrong name, no definition or used before definition?)";
        throw util::xml_scenario_error(msg.str());
    }
    return it->second;
}


void InterventionManager::deploy(RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel &transmission)
{
    if (ctx.clock.interv < sim::zero()) return;

    // deploy imported infections (not strictly speaking an intervention)
    ctx.ext->importedInfections.import(ctx, population);

    // deploy timed interventions
    SimTime now = ctx.clock.interv_date();
    while (ctx.ext->timed[ctx.ext->nextTimed].date <= now)
    {
        if (util::CommandLine::option(util::CommandLine::VERBOSE))
        {
            cout << "deploy::" << now << "::Timed";
            print_timed_details(cout, ctx.ext->timed[ctx.ext->nextTimed]);
            cout << endl;
        }

        ctx.ext->timed[ctx.ext->nextTimed].deploy(ctx, population, transmission);
        ctx.ext->nextTimed += 1;
    }

    // deploy continuous interventions
    for (Host::Human &human : population)
    {
        uint32_t nextCtsDist = human.nextCtsDist();
        // deploy continuous interventions
        while (nextCtsDist < ctx.ext->continuous.size())
        {
            if (util::CommandLine::option(util::CommandLine::VERBOSE))
            {
                cout << "deploy::" << now << "::Continuous";
                ctx.ext->continuous[nextCtsDist].print_details(cout);
                cout << endl;
            }
            if (!ctx.ext->continuous[nextCtsDist].filterAndDeploy(ctx, human)) break; // deployment (and all remaining) happens in the future
            nextCtsDist = ++human.nextCtsDist();
        }
    }
}

} // namespace interventions
} // namespace OM
