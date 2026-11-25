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

#include "Clock.h"

#include "RunContext.h"
#include "Host/HostColumns.h" // advise_heap_huge() only

#include "util/errors.h"
#include "util/CommandLine.h"
#include "xml/DocumentLoader.h"

#include <limits>
#ifdef __linux__
#include <malloc.h>
#endif

using namespace OpenMalaria;

/// main() — parses the command line, runs the simulation (run.cpp),
/// reports errors (error texts and order unchanged from the old main()).
int main(int argc, char* argv[]) {
#ifdef __linux__
    // Heap huge-page setup (best-effort perf hints only, see Host/HostColumns.h;
    // values untouched, bit-exact). Equivalent of
    // GLIBC_TUNABLES=glibc.malloc.hugetlb=1 (measured -6.1% on the pop-10000
    // benchmark). This kernel gives fault-time THP to a heap VMA >= 2 MB whose
    // pages are first touched while it is flagged; brk grows it in ~128 KB
    // steps, and each step lands in a fresh sub-2 MB VMA that never sees a
    // huge page. So: single-step 64 MB extensions (TOP_PAD), no trim shrinking
    // (a re-grow after trim also creates a fresh unflagged VMA), then flag the
    // still-small heap — every page of the extension faults huge from then on.
    // mallopt() on non-glibc is a no-op returning 0; ignore the result.
    mallopt(M_TOP_PAD, 64 << 20);
    mallopt(M_TRIM_THRESHOLD, std::numeric_limits<int>::max() / 8);
#endif
    advise_heap_huge(/*minSize=*/0);

    int exitStatus = EXIT_SUCCESS;

    string scenarioFile;
    RunContext ctx; // everything this run owns (plan §5)

    try {
        scenarioFile = CommandLine::parse(argc, argv);
        run_simulation(ctx, scenarioFile);
    } catch (const cmd_exception& e) {
        if (e.getCode() == 0) {
            // this is not an error, but exiting due to command line
            std::cerr << e.what() << "; exiting..." << endl;
        } else {
            std::cerr << "Command-line error: " << e.what();
            exitStatus = e.getCode();
        }
    } catch (const traced_exception& e) {
        std::cerr << "Code error: " << e.what() << endl;
        std::cerr << e << flush;
        std::cerr << "This is likely an error in the C++ code. Please report!" << endl;
        exitStatus = e.getCode();
    } catch (const xml_scenario_error& e) {
        std::cerr << "Error: " << e.what() << endl;
        std::cerr << "In: " << scenarioFile << endl;
        exitStatus = e.getCode();
    } catch (const base_exception& e) {
        std::cerr << "Error: " << e.message() << endl;
        exitStatus = e.getCode();
    } catch (const exception& e) {
        if (reportXmlException(std::current_exception())) {
            exitStatus = Error::XSD;
        } else {
            std::cerr << "Error: " << e.what() << endl;
            exitStatus = EXIT_FAILURE;
        }
    } catch (...) {
        std::cerr << "Unknown error" << endl;
        exitStatus = EXIT_FAILURE;
    }

    // If we get to here, we already know an error occurred.
    if (errno != 0) std::perror("OpenMalaria");

    return exitStatus;
}
