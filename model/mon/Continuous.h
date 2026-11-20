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

#ifndef Hmod_Output_Continuous
#define Hmod_Output_Continuous

#include "Clock.h"
#include "RunContext.h"
#include <string>
#include "Population.h"
#include "Host/Human.h"

#include <functional>

namespace scnXml{ class Monitoring; }
namespace OM {
    class Population;
namespace mon {
namespace Continuous {
    
    /** Functions dealing with continuous output data.
     *
     * Requirements:
     *  (1) frequency of and which data is output should be controllable
     *  (2) format should be compatible with LiveGraph and (German) Excel.
     */
	/** Load XML description of options and open the continuous-output file.
	 *
	 * Callbacks should be registered before init() is called.
	 * @param isCheckpoint  Resuming: reopen the existing output file without
	 * truncating or rewriting the header (the checkpoint read seeks it). */
	void init (RunContext& ctx, const scnXml::Monitoring& monitoring,
	    bool isCheckpoint=false);

	/** Register a callback function which produces output.
	 *
	 * This function will be called to generate output, if enabled in XML.
	 * It may output more than one statistic, if for example vector output
	 * is wanted instead of a single value. It should then title these in
	 * the form "name(index)".
	 *
	 * @param optName	Name of this output, (used for XML on/off options)
	 * @param titles	Titles for the output table; each should be preceeded by a \t
	 * @param outputCb A callback function, which when called, outputs its
	 * data to the passed stream, with each entry preceeded by '\t'.
	 */
    void registerCallback (RunContext& ctx, string optName, string titles, function<void(ostream&)> f);

    void registerCallback (RunContext& ctx, string optName, string titles, function<void(const vector<Host::Human> &, ostream&)> f);

    void registerCallback (RunContext& ctx, string optName, string titles, function<void(Population &, ostream&)> f);
	
	/// Generate time-step's output. Called at beginning of time step.
        /// Passed population since some callbacks use this to generate output.
	void update (RunContext& ctx, Population &population);


        inline void registerCallback (RunContext& ctx, string optName, string titles, function<void(ostream&)> f){
            auto& registered = ctx.ctsRegistered;
            assert(registered.count(optName) == 0); // name clash/registered twice?
            function<void(const Population&, ostream&)> _f = [f](const Population&, ostream& ostream){ f(ostream); };
            registered[optName] = {titles, _f};
        }

        inline void registerCallback (RunContext& ctx, string optName, string titles, function<void(const vector<Host::Human> &, ostream&)> f){
            auto& registered = ctx.ctsRegistered;
            assert(registered.count(optName) == 0); // name clash/registered twice?
            function<void(const Population&, ostream&)> _f = [f](const Population &p, ostream& ostream){ f(p.humans, ostream); };
            registered[optName] = {titles, _f};
        }

        inline void registerCallback (RunContext& ctx, string optName, string titles, function<void(Population &, ostream&)> f){
            auto& registered = ctx.ctsRegistered;
            assert(registered.count(optName) == 0); // name clash/registered twice?
            registered[optName] = {titles, f};
        }

        inline void update (RunContext& ctx, Population &population){
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
} // namespace OM
#endif
