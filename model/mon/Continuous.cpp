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

#include "RunContext.h"
#include "mon/Continuous.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"

#include <vector>
#include <map>
#include <fstream>

namespace OM {
    namespace mon {
    namespace Continuous {
        using util::xml_scenario_error;

        void registerCallback (RunContext& ctx, string optName, string titles, function<void(ostream&)> f){
            auto& registered = ctx.ctsRegistered;
            assert(registered.count(optName) == 0); // name clash/registered twice?
            function<void(const Population&, ostream&)> _f = [f](const Population&, ostream& ostream){ f(ostream); };
            registered[optName] = {titles, _f};
        }

        void registerCallback (RunContext& ctx, string optName, string titles, function<void(const vector<Host::Human> &, ostream&)> f){
            auto& registered = ctx.ctsRegistered;
            assert(registered.count(optName) == 0); // name clash/registered twice?
            function<void(const Population&, ostream&)> _f = [f](const Population &p, ostream& ostream){ f(p.humans, ostream); };
            registered[optName] = {titles, _f};
        }

        void registerCallback (RunContext& ctx, string optName, string titles, function<void(Population &, ostream&)> f){
            auto& registered = ctx.ctsRegistered;
            assert(registered.count(optName) == 0); // name clash/registered twice?
            registered[optName] = {titles, f};
        }

        void update (RunContext& ctx, Population &population){
            auto& ctsOStream = ctx.ctsOStream;
            auto& toReport = ctx.ctsToReport;
            const SimTime ctsPeriod = ctx.ctsPeriod;
            const bool duringInit = ctx.ctsDuringInit;
            if( ctsPeriod == sim::zero() )
                return;	// output disabled
            if( !duringInit ){
                if( ctx.clock.interv < sim::zero()
                    || mod_nn(ctx.clock.interv, ctsPeriod) != sim::zero() )
                    return;
            } else {
                if( mod_nn(ctx.clock.now(), ctsPeriod) != sim::zero() )
                    return;
                ctsOStream << ctx.clock.inSteps(ctx.clock.now()) << '\t';
            }

            if( duringInit && ctx.clock.interv < sim::zero() ){
                ctsOStream << "nan";
            }else{
                // NOTE: we could switch this to output dates, but (1) it would be
                // breaking change and (2) it may be harder to use.
                ctsOStream << ctx.clock.inSteps(ctx.clock.interv);
            }
            for (const auto& report : toReport) std::get<1>(report)(population, ctsOStream);

            // We must flush often to avoid temporarily outputting partial lines
            // (resulting in incorrect real-time graphs).
            ctsOStream << '\n' << flush;

            ctx.ctsStreamOffset = ctsOStream.tellp();
        }
    } // namespace Continuous
    } // namespace mon
}
