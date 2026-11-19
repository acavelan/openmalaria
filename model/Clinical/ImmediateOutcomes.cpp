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

#include "core/RunExtensions.h"
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
double ImmediateOutcomes::cureRateUCOfficial[CM5DayCommon::NumCaseTypes];
double ImmediateOutcomes::cureRateUCSelfTreat[CM5DayCommon::NumCaseTypes];
WithinHost::TreatmentId ImmediateOutcomes::treatmentUC[CM5DayCommon::NumCaseTypes];
bool ImmediateOutcomes::useDiagnosticUC = false;




// ———  per-human, update  ———

void immOutUncomplicatedEvent (ClinicalModel& cm, RunContext& ctx,
    Human& human,
    Episode::State pgState
) {
    // If last treatment prescribed was in recent memory, consider second line.
    CM5DayCommon::CaseType regimen = CM5DayCommon::FirstLine;
    if (cm.m_tLastTreatment + ClinicalModel::healthSystemMemory > ctx.clock.ts0()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = CM5DayCommon::SecondLine;
    }

    cm.latestReport.update (ctx, human, Episode::State( pgState ) );

    // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
    double x = human.rng.uniform_01();
    if( x < CM5DayCommon::accessUCAny[regimen] * cm.m_treatmentSeekingFactor ){
        // UC1: official care OR self treatment
        // UC2: official care only
        
        if( ImmediateOutcomes::useDiagnosticUC ){
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
            auto diag = WithinHost::diagnostics::monitoringDiagnostic( ctx );
            if( !human.withinHostModel->diagnosticResult(human.rng, diag) )
                return; // negative outcome: no treatment
        }
        
        cm.m_tLastTreatment = ctx.clock.ts0();
        mon::recordEvent(CM5DayCommon::measures[regimen], human);
        
        double p = ( x < CM5DayCommon::accessUCSelfTreat[regimen] * cm.m_treatmentSeekingFactor ) ?
            ImmediateOutcomes::cureRateUCSelfTreat[regimen] : ImmediateOutcomes::cureRateUCOfficial[regimen];
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if( human.rng.bernoulli(p) ){
            // Could report Episode::RECOVERY to cm.latestReport,
            // but we don't report out-of-hospital recoveries anyway.
            human.withinHostModel->treatment( human, ImmediateOutcomes::treatmentUC[regimen] );
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
