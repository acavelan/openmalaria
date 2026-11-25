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

#include "Clinical/DecisionTree5Day.h"
#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/errors.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;

enum CaseType { FirstLine, SecondLine, NumCaseTypes };

// These parameters are set by setHealthSystem() and do not need checkpointing.
const CMDecisionTree *treeUCOfficial = 0, *treeUCSelfTreat = 0;

// ———  static, init  ———



// ———  per-human, update  ———

void DecisionTree5Day::uncomplicatedEvent ( Human& human, Episode::State pgState ){    
    // If last treatment prescribed was in recent memory, consider second line.
    CaseType regimen = FirstLine;
    if (m_tLastTreatment + healthSystemMemory > sim::ts0()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = SecondLine;
    }

    latestReport.update (human, Episode::State( pgState ) );

    // parity: d — clinical decision fraction draw order (plan §3d)
    double x = human.rng.uniform_01();
    if( x < accessUCAny[regimen] * m_treatmentSeekingFactor ){
        CMHostData hostData( human, sim::inYears(human.age(sim::ts0())), pgState );
        
        // Run tree (which may deploy treatment)
        CMDTOut output = ( x < accessUCSelfTreat[regimen] * m_treatmentSeekingFactor ) ?
            treeUCSelfTreat->exec( hostData ) :
            treeUCOfficial->exec( hostData );
        
        if( output.treated ){   // if any treatment or intervention deployed
            m_tLastTreatment = sim::ts0();
            mon::recordEvent(measures[regimen], human);
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
