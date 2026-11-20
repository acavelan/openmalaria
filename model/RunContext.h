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

#ifndef Hmod_RunContext
#define Hmod_RunContext

#include "Clock.h"
#include "util/mod.h"
#include "util/random.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Diagnostic.h"      // leaf (Clock, random)
#include "PkPd/LSTMDrugType.h"
#include "PkPd/LSTMMedicate.h"
#include "interventions/ComponentId.h"       // leaf (Clock, checkpoint)
#include "PopulationAgeStructure.h"
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace scnXml {
class Scenario;
}

namespace OM {

// Forward declarations for the registry members below: their defining headers
// need a complete RunContext (they would cycle through RunContext.h), so this
// header only ever sees pointers to them. gcc instantiates member deleters at
// class scope, so any member whose type is incomplete here is a raw pointer
// created by RunContext's out-of-line ctor and deleted by its out-of-line
// dtor (model/run.cpp).
namespace Clinical { struct FlatTree; }
namespace Host { class ImportedInfections; }
namespace Transmission { class PerHostAnophParams; }
namespace interventions {
    struct ComponentConfig;
    struct ContinuousHumanDeployment;
    struct TimedDeployment;
    struct ITNParams;
    struct IRSParams;
    struct GVIParams;
    class VaccineParams;
}
namespace mon { namespace internal { struct RuntimeState; } }
class Population;
namespace util {
class ModelNameProvider;
}
namespace WithinHost {
class Treatments;
struct Infection;
}

/** The simulation clock and scenario time constants: everything the old
 * sim class stored in globals. Constructed fresh with each RunContext, so
 * independent runs never share clock state. Values and construction order
 * match the old sim:: statics exactly. */

/** Everything one simulation run owns: clock, RNGs, registries, intervention
 * state, monitoring, continuous output. Owned by the caller of
 * run_simulation(); never a global. */
class RunContext {
public:
    RunContext();     // out-of-line: news the pointer members (run.cpp)
    ~RunContext();    // out-of-line: deletes them (incomplete types here)

    SimClock clock;
    /** Master RNG. // parity: order — constructed here (0, 0), matching the
     * old global's construction; the percentile-stream seed depends on it;
     * run_simulation reseeds from --seed/Iseed after validation; checkpointed
     * in model/checkpoint.cpp. */
    util::MasterRng master_RNG{0, 0};
    /// Genotype registry (built by Genotypes::init during scenario load).
    WithinHost::Genotypes genotypes;
    /// Drug-type registry (PkPd; built by LSTMDrugType::init).
    PkPd::DrugRegistry drugs;
    /// PK/PD treatment-course registry (LSTMTreatments::init).
    PkPd::TreatmentRegistry treatments;
    /// Demography / target age-structure state (PopulationAgeStructure.h).
    AgeStructureState age;

    // ——— Per-run registries (was: struct RunExtensions, reached via an ext
    // indirection created by ensure_extensions(); dissolved into plain
    // members). // parity: declaration order below == the old RunExtensions
    // order == creation order in RunContext's ctor; registries and unittest
    // backend installers depend on it. Pointer members: see the forward
    // declarations above; created/destroyed out-of-line in run.cpp. ———

    /// Compiled decision trees (frees memory at end of run; was:
    /// decision_library of polymorphic trees with equality dedup — the dedup
    /// only saved memory, identical trees behave identically).
    std::vector<std::unique_ptr<Clinical::FlatTree>>* decision_trees;
    /// WithinHost::Treatments option registry (built while reading XML).
    std::vector<WithinHost::Treatments>* whTreatments;
    /// Within-host infection factory, set by the active model's init.
    WithinHost::Infection* (*createInfection)
        (RunContext&, util::LocalRng&, uint32_t, int) = nullptr;
    WithinHost::Infection* (*checkpointedInfection)
        (RunContext&, std::istream&) = nullptr;
    /// Unittest backend for WithinHost::WHInterface Kind::Mock (installed by
    /// unittest/WHMock.cpp; null in production, where Kind::Mock is never
    /// created). Same fn-pointer-table idiom as createInfection above.
    struct WHMockHooks {
        void* state = nullptr;
        void (*treatment)(void* state);
        bool (*treatSimple)(void* state, SimTime timeLiver, SimTime timeBlood);
        void (*treatPkPd)(void* state, RunContext& ctx, size_t schedule,
            size_t dosages, double age, double delay_d);
        double (*getTotalDensity)(const void* state);
        bool (*diagnosticResult)(const void* state, util::LocalRng& rng,
            const WithinHost::Diagnostic& diagnostic);
    };
    const WHMockHooks* whMock = nullptr;
    /// Diagnostics registry (built while reading XML / from defaults).
    std::map<std::string, WithinHost::Diagnostic> diagnostic_set;
    const WithinHost::Diagnostic* monitoring_diagnostic = nullptr;

    // ——— intervention registries (was: InterventionManager/ITN/IRS/GVI/
    //     Vaccine/removeAtIds statics) ———
    /// Map of textual identifiers to numeric identifiers for components.
    std::map<std::string, interventions::ComponentId> identifierMap;
    /// All human intervention components, indexed by ComponentId; owns them
    /// (append-only: registries below hold stable payload pointers).
    std::vector<interventions::ComponentConfig>* humanComponents;
    /// Continuous interventions, sorted by deployment age (weakly increasing).
    std::vector<interventions::ContinuousHumanDeployment>* continuous;
    /// All timed interventions, sorted (time weakly increasing), sentinel last.
    std::vector<interventions::TimedDeployment>* timed;
    uint32_t nextTimed = 0;
    /// Imported infections (not really an intervention; schema association).
    OM::Host::ImportedInfections* importedInfections;
    /// Sub-population removal component lists: array of
    /// interventions::SubPopRemove::NUM lists (new[]ed in run.cpp).
    std::vector<interventions::ComponentId>* removeAtIds;
    /// Per-type component parameter look-ups by id (not owning: components
    /// are owned by humanComponents).
    std::vector<const interventions::ITNParams*> itnComponents;
    std::vector<const interventions::IRSParams*> irsComponents;
    std::vector<const interventions::GVIParams*> gviComponents;
    std::vector<const interventions::VaccineParams*> vaccineComponents;
    /// Which vaccine component deployment monitoring reports on.
    interventions::ComponentId vaccineReportComponent =
        interventions::ComponentId::wholePop();

    // ——— entomology parameters (was: PerHostAnophParams statics) ———
    /// Per-species host-availability/biting parameters (XML read).
    std::vector<Transmission::PerHostAnophParams>* anophParams;
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
    /// (Pointer member: its map keys need ComponentId from the heavy
    /// intervention headers, which need a complete RunContext.)
    mon::internal::RuntimeState* monitoring;
};

}
#endif
