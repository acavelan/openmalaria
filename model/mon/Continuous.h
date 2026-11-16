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

#include "Global.h"
#include <string>
#include "Population.h"
#include "Host/Human.h"

#include <functional>
#include <fstream>
#include <map>

namespace OM {
    struct Population;
namespace mon {
namespace Continuous {
    
    struct State {
        struct Output {
            string titles;
            std::function<void(Population&, ostream&)> write;
        };
        map<string, Output> registered;
        vector<Output> toReport;
        fstream ctsOStream;
        streamoff streamOffset = 0;
        SimTime ctsPeriod = sim::zero();
        bool duringInit = false;
    };

    /** Functions dealing with continuous output data.
     *
     * Requirements:
     *  (1) frequency of and which data is output should be controllable
     *  (2) format should be compatible with LiveGraph and (German) Excel.
     */
    // Open configured output, or resume it before checkpoint() restores the position.
    void init(State& state, const string& filename, bool isCheckpoint);

    void checkpoint(State& state, ostream& stream);
    void checkpoint(State& state, istream& stream);
        
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
    void registerCallback (State& state, string optName, string titles, function<void(ostream&)> f);

    void registerCallback (State& state, string optName, string titles, function<void(const vector<Host::Human> &, ostream&)> f);

    void registerCallback (State& state, string optName, string titles, function<void(Population &, ostream&)> f);
	
	/// Generate time-step's output. Called at beginning of time step.
        /// Passed population since some callbacks use this to generate output.
	void update (State& state, Population &population);

} // namespace Continuous
} // namespace mon
} // namespace OM
#endif
