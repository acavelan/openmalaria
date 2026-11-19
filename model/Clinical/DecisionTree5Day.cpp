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
#include "Clinical/DecisionTree5Day.h"
#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/WHInterface.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/errors.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;

enum CaseType { FirstLine, SecondLine, NumCaseTypes };

// These parameters are set by setHealthSystem() and do not need checkpointing.
const FlatTree *treeUCOfficial = 0, *treeUCSelfTreat = 0;

// ———  static, init  ———



// ———  per-human, update  ———

void treeUncomplicatedEvent (ClinicalModel& cm, RunContext& ctx, Human& human, Episode::State pgState ){    
    // If last treatment prescribed was in recent memory, consider second line.
    CaseType regimen = FirstLine;
    if (cm.m_tLastTreatment + ClinicalModel::healthSystemMemory > ctx.clock.ts0()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = SecondLine;
    }

    cm.latestReport.update (ctx, human, Episode::State( pgState ) );

    // parity: d — clinical decision fraction draw order (plan §3d)
    double x = human.rng().uniform_01();
    if( x < CM5DayCommon::accessUCAny[regimen] * cm.m_treatmentSeekingFactor ){
        // Run tree (which may deploy treatment)
        DTExecCtx ec( ctx, human, sim::inYears( human.age( ctx.clock.ts0() ) ),
                      pgState );
        CMDTOut output = exec_decision_tree(
            x < CM5DayCommon::accessUCSelfTreat[regimen] * cm.m_treatmentSeekingFactor ?
                *treeUCSelfTreat : *treeUCOfficial, ec );
        
        if( output.treated ){   // if any treatment or intervention deployed
            cm.m_tLastTreatment = ctx.clock.ts0();
            mon::recordEvent(CM5DayCommon::measures[regimen], human);
        }
        if( output.screened ){
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
        }
        
        human.withinHostModel->optionalPqTreatment(human);
    } else {
        // No care sought
    }
}

}
}
