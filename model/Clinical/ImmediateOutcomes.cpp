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

#include "core/RunContext.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/random.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;

// These parameters are set by setHealthSystem() and do not need checkpointing.
double ImmediateOutcomes::cureRateUCOfficial[ImmediateOutcomes::NumCaseTypes];
double ImmediateOutcomes::cureRateUCSelfTreat[NumCaseTypes];
WithinHost::TreatmentId ImmediateOutcomes::treatmentUC[NumCaseTypes];
bool ImmediateOutcomes::useDiagnosticUC = false;




// ———  per-human, update  ———

void ImmediateOutcomes::uncomplicatedEvent (RunContext& ctx,
    Human& human,
    Episode::State pgState
) {
    // If last treatment prescribed was in recent memory, consider second line.
    CaseType regimen = FirstLine;
    if (m_tLastTreatment + healthSystemMemory > ctx.clock.ts0()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = SecondLine;
    }

    latestReport.update (ctx, human, Episode::State( pgState ) );

    // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
    double x = human.rng.uniform_01();
    if( x < accessUCAny[regimen] * m_treatmentSeekingFactor ){
        // UC1: official care OR self treatment
        // UC2: official care only
        
        if( useDiagnosticUC ){
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
            auto diag = WithinHost::diagnostics::monitoringDiagnostic();
            if( !human.withinHostModel->diagnosticResult(human.rng, diag) )
                return; // negative outcome: no treatment
        }
        
        m_tLastTreatment = ctx.clock.ts0();
        mon::recordEvent(measures[regimen], human);
        
        double p = ( x < accessUCSelfTreat[regimen] * m_treatmentSeekingFactor ) ?
            cureRateUCSelfTreat[regimen] : cureRateUCOfficial[regimen];
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if( human.rng.bernoulli(p) ){
            // Could report Episode::RECOVERY to latestReport,
            // but we don't report out-of-hospital recoveries anyway.
            human.withinHostModel->treatment( human, treatmentUC[regimen] );
        } else {
            // No change in parasitological status: treated outside of hospital
        }
        
        human.withinHostModel->optionalPqTreatment(human);
    } else {
        // No change in parasitological status: non-treated
    }
}

}
}
