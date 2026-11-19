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
#include "mon/Monitoring.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/Infection/CommonInfection.h"
#include "Host/ImportedInfections.h"
#include "interventions/Deployments.h"
#include "interventions/GVI.h"
#include "interventions/HumanComponents.h"
#include "interventions/IRS.h"
#include "interventions/ITN.h"
#include "interventions/Vaccine.h"

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

    // ——— intervention registries (was: InterventionManager/ITN/IRS/GVI/
    //     Vaccine/removeAtIds statics) ———
    /// Map of textual identifiers to numeric identifiers for components.
    std::map<std::string, interventions::ComponentId> identifierMap;
    /// All human intervention components, indexed; owns the components.
    std::vector<std::unique_ptr<interventions::HumanInterventionComponent>> humanComponents;
    /// Continuous interventions, sorted by deployment age (weakly increasing).
    std::vector<interventions::ContinuousHumanDeployment> continuous;
    /// All timed interventions, sorted (time weakly increasing), sentinel last.
    std::vector<std::unique_ptr<interventions::TimedDeployment>> timed;
    uint32_t nextTimed = 0;
    /// Imported infections (not really an intervention; schema association).
    OM::Host::ImportedInfections importedInfections;
    /// Sub-population removal component lists, per remove-at trigger.
    std::vector<interventions::ComponentId> removeAtIds[interventions::SubPopRemove::NUM];
    /// Per-type component parameter look-ups by id (not owning: components
    /// are owned by humanComponents).
    std::vector<interventions::ITNComponent*> itnComponents;
    std::vector<interventions::IRSComponent*> irsComponents;
    std::vector<interventions::GVIComponent*> gviComponents;
    std::vector<interventions::VaccineComponent*> vaccineComponents;
    /// Which vaccine component deployment monitoring reports on.
    interventions::ComponentId vaccineReportComponent =
        interventions::ComponentId::wholePop();

    // ——— entomology parameters (was: PerHostAnophParams statics) ———
    /// Per-species host-availability/biting parameters (XML read).
    std::vector<Transmission::PerHostAnophParams> anophParams;
    std::vector<double> entoAvailabilityPercentiles;

    // ——— continuous output (was: mon/Continuous.cpp statics) ———
    /// Continuous-output file; its resume position; callback registry.
    std::fstream ctsOStream;
    std::streamoff ctsStreamOffset = 0;
    std::map<std::string, std::tuple<std::string,
        std::function<void(Population&, std::ostream&)>>> ctsRegistered;
    std::vector<std::tuple<std::string,
        std::function<void(Population&, std::ostream&)>>> ctsToReport;
    SimTime ctsPeriod = sim::zero();
    bool ctsDuringInit = false;

    /// Monitoring runtime state (survey numbers, stores, conditions; mon/).
    /// (In ext, not RunContext: the state's map keys need ComponentId
    /// complete, which would cycle RunContext.h through Interfaces.h.)
    mon::internal::RuntimeState monitoring;
};

/// Create RunContext::ext if not yet present (idempotent).
inline void ensure_extensions( RunContext& ctx ){
    if( !ctx.ext ) ctx.ext = new RunExtensions();
}

}
#endif
