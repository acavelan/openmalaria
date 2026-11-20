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

#ifndef Hmod_core_RunContext
#define Hmod_core_RunContext

#include "Clock.h"
#include "util/mod.h"
#include "util/random.h"
#include "Host/WithinHost/Genotypes.h"
#include "PkPd/LSTMDrugType.h"
#include "PkPd/LSTMMedicate.h"
#include "PopulationAgeStructure.h"
#include <memory>
#include <string>
#include <vector>

namespace scnXml {
class Scenario;
}

namespace OM {

struct RunExtensions;   // core/RunExtensions.h: registries needing complete RunContext

namespace util {
class ModelNameProvider;
}

/** The simulation clock and scenario time constants: everything the old
 * sim class stored in globals. Constructed fresh with each RunContext, so
 * independent runs never share clock state. Values and construction order
 * match the old sim:: statics exactly. */

/** Everything one simulation run owns: the clock now; further run state
 * (parameters, interventions, monitoring, RNGs, ...) moves in here through
 * the rest of phase 3. Owned by the caller of run_simulation(); never a
 * global. */
class RunContext {
public:
    ~RunContext();      // out-of-line: deletes ext (incomplete type here)
    
    SimClock clock;
    /** Master RNG. // parity: order — constructed here (0, 0), matching the
     * old global's construction; the percentile-stream seed depends on it;
     * run_simulation reseeds from --seed/Iseed after validation. */
    util::MasterRng master_RNG{0, 0};
    /// Genotype registry (built by Genotypes::init during scenario load).
    WithinHost::Genotypes genotypes;
    /// Drug-type registry (PkPd; built by LSTMDrugType::init).
    PkPd::DrugRegistry drugs;
    /// PK/PD treatment-course registry (LSTMTreatments::init).
    PkPd::TreatmentRegistry treatments;
    /// Registries whose types need a complete RunContext; created by
    /// ensure_extensions() (model/core/RunExtensions.h) at run start,
    /// owned and deleted by ~RunContext (raw pointer: gcc 16 instantiates
    /// unique_ptr's deleter on incomplete types at class scope).
    RunExtensions* ext = nullptr;
    /// Demography / target age-structure state (PopulationAgeStructure.h).
    AgeStructureState age;
};

}
#endif
