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

#ifndef Hmod_core_RunExtensions
#define Hmod_core_RunExtensions

// Per-run registries that cannot be value members of RunContext.h because
// their headers need a complete RunContext (include cycle). Owned by
// RunContext::ext; created by ensure_extensions() at the start of
// run_simulation() (and by UnittestUtil for tests).

#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/Infection/CommonInfection.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace OM {

struct RunExtensions {
    /// Decision-tree library (frees memory at end of run; dedup on add).
    std::vector<std::unique_ptr<Clinical::CMDecisionTree>> decision_library;
    /// WithinHost::Treatments option registry (built while reading XML).
    std::vector<WithinHost::Treatments> whTreatments;
    /// Within-host infection factory, set by the active model's init.
    WithinHost::CommonInfection* (*createInfection)
        (RunContext&, util::LocalRng&, uint32_t, int) = nullptr;
    WithinHost::CommonInfection* (*checkpointedInfection)
        (RunContext&, std::istream&) = nullptr;
    /// Diagnostics registry (built while reading XML / from defaults).
    std::map<std::string, WithinHost::Diagnostic> diagnostic_set;
    const WithinHost::Diagnostic* monitoring_diagnostic = nullptr;
};

/// Create RunContext::ext if not yet present (idempotent).
inline void ensure_extensions( RunContext& ctx ){
    if( !ctx.ext ) ctx.ext = new RunExtensions();
}

}
#endif
