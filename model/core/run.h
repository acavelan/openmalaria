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

#ifndef Hmod_core_Run
#define Hmod_core_Run

#include <string>

namespace OM {
class RunContext;

/** Load a scenario and run the complete simulation (init, warmup, EIR
 * calibration, main phase, survey output) — the whole of the old main()
 * body after command-line parsing, in the frozen order (plan §3).
 *
 * All state lives in ctx, owned by the caller, so one process can run the
 * model any number of times independently. Exceptions are not caught:
 * main() (and TwoRunsSuite) own the error reporting. Assumes
 * util::CommandLine::parse() has already been called. */
void run_simulation(OM::RunContext& ctx, const std::string& scenarioFile);

}
#endif
