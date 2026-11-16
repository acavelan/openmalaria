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

#include "mon/Continuous.h"
#include "util/errors.h"

#include <vector>
#include <map>
#include <fstream>

namespace OM {
    namespace mon {
    namespace Continuous {

        void init(State& state, const string& filename, bool isCheckpoint) {
            if (state.ctsPeriod == sim::zero()) return;

            if( isCheckpoint ){
                // When loading a check-point, we resume reporting to this file.
                state.ctsOStream.open(filename, ios::binary|ios::ate|ios::in|ios::out);
                if( state.ctsOStream.fail() )
                    throw util::checkpoint_error ("Continuous: resume error (no file)");
                return; // checkpoint() sets the saved output position later
            }

            state.ctsOStream.open(filename, ios::binary|ios::out);
            state.ctsOStream << "##\t##\n"; // live-graph needs a delimiter specifier when it is not a comma
            if( state.duringInit ) state.ctsOStream << "simulation time\t";
            state.ctsOStream << "timestep"; // TODO: change to days or remove or leave?
            for (const auto& report : state.toReport) state.ctsOStream << report.titles;
            state.ctsOStream << '\n' << flush;
            state.streamOffset = state.ctsOStream.tellp();
        }

        void checkpoint (State& state, ostream& stream){
            if( state.ctsPeriod == sim::zero() )
                return;	// output disabled

            state.streamOffset & stream;
        }
        void checkpoint (State& state, istream& stream){
            if( state.ctsPeriod == sim::zero() )
                return;	// output disabled

            /* We attempt to resume output correctly after a reload by recording
            * the last position, and relocating there.
            * 
            * (Keeping results in memory until end of sim would be another,
            * slightly safer, option, but loses real-time output.) */
            state.streamOffset & stream;
            // We skip back to the last write-point, so anything written after the
            // last checkpoint will be repeated:
            state.ctsOStream.seekp( state.streamOffset, ios_base::beg );

            if( state.ctsOStream.fail() )
                throw util::checkpoint_error ("Continuous: resume error (bad pos/file)");
        }

        void registerCallback(State& state, string optName, string titles, function<void(ostream&)> f) {
            registerCallback(state, std::move(optName), std::move(titles),
                [f = std::move(f)](Population&, ostream& out) { f(out); });
        }

        void registerCallback(State& state, string optName, string titles, function<void(const vector<Host::Human>&, ostream&)> f) {
            registerCallback(state, std::move(optName), std::move(titles),
                [f = std::move(f)](Population& population, ostream& out) { f(population.humans, out); });
        }

        void registerCallback(State& state, string optName, string titles, function<void(Population&, ostream&)> f) {
            assert(state.registered.count(optName) == 0); // name clash/registered twice?
            state.registered[optName] = {std::move(titles), std::move(f)};
        }

        void update (State& state, Population &population){
            if( state.ctsPeriod == sim::zero() )
                return;	// output disabled
            if( !state.duringInit ){
                if( sim::intervTime() < sim::zero()
                    || mod_nn(sim::intervTime(), state.ctsPeriod) != sim::zero() )
                    return;
            } else {
                if( mod_nn(sim::now(), state.ctsPeriod) != sim::zero() )
                    return;
                state.ctsOStream << sim::inSteps(sim::now()) << '\t';
            }

            if( state.duringInit && sim::intervTime() < sim::zero() ){
                state.ctsOStream << "nan";
            }else{
                // NOTE: we could switch this to output dates, but (1) it would be
                // breaking change and (2) it may be harder to use.
                state.ctsOStream << sim::inSteps(sim::intervTime());
            }
            for (const auto& report : state.toReport) report.write(population, state.ctsOStream);

            // We must flush often to avoid temporarily outputting partial lines
            // (resulting in incorrect real-time graphs).
            state.ctsOStream << '\n' << flush;

            state.streamOffset = state.ctsOStream.tellp();
        }
    } // namespace Continuous
    } // namespace mon
}
