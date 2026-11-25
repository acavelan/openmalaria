

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

// The 5-day case-management chain: shared config (CM5DayCommon) plus the two
// health-system kinds, DecisionTree5Day and ImmediateOutcomes (one file each
// before this merge; bodies verbatim).

#include "RunContext.h"
#include "Clinical/CM5DayCommon.h"
#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/WHInterface.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/random.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;

// ———  static pararmeters and methods  ———

// These parameters are set by setHealthSystem() and do not need checkpointing.
mon::Measure CM5DayCommon::measures[NumCaseTypes] = {
    mon::measure("nTreatments1"),  // first line official
    mon::measure("nTreatments2")   // second line official
};
double CM5DayCommon::accessUCAny[NumCaseTypes];
double CM5DayCommon::accessUCSelfTreat[NumCaseTypes];
double CM5DayCommon::accessSevere;
double CM5DayCommon::cureRateSevere;
WithinHost::TreatmentId CM5DayCommon::treatmentSevere;

bool cfr_pf_use_hospital = false;

void CM5DayCommon::init(){
    cfr_pf_use_hospital = util::ModelOptions::option( util::CFR_PF_USE_HOSPITAL );
}


// ———  per-human, construction and destruction  ———

void clinicalUpdate5Day (ClinicalModel& cm, RunContext& ctx, Human& human, double ageYears, WithinHost::Pathogenesis::StatePair &pg) {
    Episode::State pgState = static_cast<Episode::State>( pg.state );

    if (pgState & Episode::MALARIA) {
        if (pgState & Episode::COMPLICATED){
            severeMalaria (cm, ctx, human, pgState, ageYears, cm.doomed);
        }else if (ClinicalModel::indirectMortBugfix || !pg.indirectMortality) {
            // NOTE: the "not indirect mortality" bit is a historical accident.
            // Validity is debatable, but there's no point changing now.
            // (This does affect tests.)
            if (cm.kind == ClinicalModel::Kind::ImmediateOutcomes)
                immOutUncomplicatedEvent (cm, ctx, human, pgState);
            else
                treeUncomplicatedEvent (cm, ctx, human, pgState);
        }
    } else if (pgState & Episode::SICK) { // sick but not from malaria
        if (cm.kind == ClinicalModel::Kind::ImmediateOutcomes)
                immOutUncomplicatedEvent (cm, ctx, human, pgState);
            else
                treeUncomplicatedEvent (cm, ctx, human, pgState);
    }
    
    if (pg.indirectMortality && cm.doomed == NOT_DOOMED)
        cm.doomed = -ctx.clock.one_TS();
    
    if( cm.t_last_treatment() == ctx.clock.ts0() )
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_TREATMENT);
    
    if( pgState & Episode::SICK )
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_BOUT);
}

void severeMalaria (ClinicalModel& cm, RunContext& ctx, 
    Human &human,
    Episode::State pgState,
    double ageYears,
    int& doomed
) {
    // Probability of getting treatment (only part which is case managment):
    double p2 = CM5DayCommon::accessSevere * cm.m_treatmentSeekingFactor;
    // Probability of getting cured after getting treatment:
    double p3 = CM5DayCommon::cureRateSevere;
    // p4 is the hospital case-fatality rate from Tanzania
    double p4 = ClinicalModel::caseFatalityRate.eval (ageYears);
    // p5a here is the community threshold case-fatality rate
    double p5a = ClinicalModel::getCommunityCFR (p4);
    // p5b here is the in-hospital treatment-failure case-fatality rate
    // In our code it was originally the community CFR, but in the published model description it should be the hospital CFR
    double p5b = (cfr_pf_use_hospital) ? p4 : p5a;
    // p6 is P(seq) for treated patients
    double p6 = ClinicalModel::pSequelaeInpatient.eval (ageYears);
    // p7 is P(seq) when parasites aren't cleared
    // note: if p7 were not equal p6, we might have to do something similar to p5b
    double p7 = p6;

    double q[9];
    // Community deaths
    q[0] = (1 - p2) * p5a;
    // Community sequelae
    q[1] = q[0] + (1 - p2) * (1 - p5a) * p7;
    // Community survival
    q[2] = q[1] + (1 - p2) * (1 - p5a) * (1 - p7);
    // In-hospital parasitological failure deaths
    q[3] = q[2] + p2 * (1 - p3) * p5b;
    // In-hospital parasitological failure sequelae
    q[4] = q[3] + p2 * (1 - p3) * (1 - p5b) * p7;
    // In-hospital parasitological failure survivors
    q[5] = q[4] + p2 * (1 - p3) * (1 - p5b) * (1 - p7);
    // In-hospital parasitological success deaths
    q[6] = q[5] + p2 * p3 * p4;
    // In-hospital parasitological success sequelae
    q[7] = q[6] + p2 * p3 * (1 - p4) * p6;
    // In-hospital parasitological success survival
    q[8] = q[7] + p2 * p3 * (1 - p4) * (1 - p6);
    
    // Expectation of death is:
    // double exDeath = (q[6] - q[5]) + (q[3] - q[2]) + q[0];
    // expanded, separating in-hospital and without (q[0]):
    const double exHospitalDeath = p2 * (p3 * p4 + (1 - p3) * p5b);
    const double exDeath = exHospitalDeath + (1 - p2) * p5a;
    mon::recordStat(mon::measure("expectedHospitalDeaths"), human, exHospitalDeath);
    mon::recordStat(mon::measure("expectedDirectDeaths"), human, exDeath);
    
    // Expectation of sequelae is:
    // double exSeq = (q[7] - q[6]) + (q[4] - q[3]) + (q[1] - q[0]);
    // expanded and simplified, noting that p7 == p6 :
    const double exSeq = (p2 * (p3 * (1 - p4) + (1 - p3) * (1 - p5b)) + (1 - p2) * (1 - p5a)) * p6;
    mon::recordStat(mon::measure("expectedSequelae"), human, exSeq);

    // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
    double prandom = human.rng().uniform_01();
    
    //NOTE: we do not model diagnostics in this case
    if( prandom >= q[2] ){      // treated in hospital
        cm.t_last_treatment() = ctx.clock.ts0();
        mon::recordEvent(mon::measure("nTreatments3"), human);
        Episode::State stateTreated = Episode::State (pgState | Episode::EVENT_IN_HOSPITAL);
        
        if( prandom >= q[5] ){  // treatment successful at clearing parasites
            // this actually means successful treatment in this case:
            human.withinHostModel->treatment( human, CM5DayCommon::treatmentSevere );
            
            if( prandom < q[6] ){       // death (despite success in clearing parasites)
                cm.latestReport.update (ctx, human, Episode::State (stateTreated | Episode::DIRECT_DEATH));
                doomed  = DOOMED_COMPLICATED;
            }else if( prandom < q[7] ){ // Patient recovers, but with sequelae (don't report full recovery)
                cm.latestReport.update (ctx, human, Episode::State (stateTreated | Episode::SEQUELAE));
            }else /*prandom < q[8]*/{   // patient recovers fully
                cm.latestReport.update (ctx, human, Episode::State (stateTreated | Episode::RECOVERY));
            }
        }else /*prandom < q[5]*/{  // Treated but parasites not cleared (in hospital)
            // No change in parasitological status among in-hospital patients
            
            if( prandom < q[3] ){       // death
                cm.latestReport.update (ctx, human, Episode::State (stateTreated | Episode::DIRECT_DEATH));
                doomed  = DOOMED_COMPLICATED;
            }else if( prandom < q[4] ){ // sequelae without parasite clearance
                cm.latestReport.update (ctx, human, Episode::State (stateTreated | Episode::SEQUELAE));
            }else /*prandom < q[5]*/{   // full recovery from episode
                cm.latestReport.update (ctx, human, pgState);
            }
        }
    }else /*prandom < q[2]*/{
        // Not treated, thus no change in parasitological status
        
        if( prandom < q[0] ){   // outcome 0: death in community
            cm.latestReport.update (ctx, human, Episode::State (pgState | Episode::DIRECT_DEATH));
            doomed  = DOOMED_COMPLICATED;
        }else if( prandom < q[1] ){     // outcome 1: sequelae in community
            cm.latestReport.update (ctx, human, Episode::State (pgState | Episode::SEQUELAE));
        }else /*prandom < q[2]*/{       // outcome 2: full recovery from episode in community
            cm.latestReport.update (ctx, human, pgState);
        }
    }
}

// ———  DecisionTree5Day kind (old DecisionTree5Day.cpp)  ———

// These parameters are set by setHealthSystem() and do not need checkpointing.
const FlatTree *treeUCOfficial = 0, *treeUCSelfTreat = 0;

// ———  static, init  ———



// ———  per-human, update  ———

void treeUncomplicatedEvent (ClinicalModel& cm, RunContext& ctx, Human& human, Episode::State pgState ){    
    // If last treatment prescribed was in recent memory, consider second line.
    CM5DayCommon::CaseType regimen = CM5DayCommon::FirstLine;
    if (cm.t_last_treatment() + ClinicalModel::healthSystemMemory > ctx.clock.ts0()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = CM5DayCommon::SecondLine;
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
            cm.t_last_treatment() = ctx.clock.ts0();
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

// ———  ImmediateOutcomes kind (old ImmediateOutcomes.cpp)  ———

// ImmediateOutcomes parameters are set by setHealthSystem() and do not need checkpointing.
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
    if (cm.t_last_treatment() + ClinicalModel::healthSystemMemory > ctx.clock.ts0()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = CM5DayCommon::SecondLine;
    }

    cm.latestReport.update (ctx, human, Episode::State( pgState ) );

    // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
    double x = human.rng().uniform_01();
    if( x < CM5DayCommon::accessUCAny[regimen] * cm.m_treatmentSeekingFactor ){
        // UC1: official care OR self treatment
        // UC2: official care only
        
        if( ImmediateOutcomes::useDiagnosticUC ){
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
            auto diag = WithinHost::diagnostics::monitoringDiagnostic( ctx );
            if( !human.withinHostModel->diagnosticResult(human.rng(), diag) )
                return; // negative outcome: no treatment
        }
        
        cm.t_last_treatment() = ctx.clock.ts0();
        mon::recordEvent(CM5DayCommon::measures[regimen], human);
        
        double p = ( x < CM5DayCommon::accessUCSelfTreat[regimen] * cm.m_treatmentSeekingFactor ) ?
            ImmediateOutcomes::cureRateUCSelfTreat[regimen] : ImmediateOutcomes::cureRateUCOfficial[regimen];
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if( human.rng().bernoulli(p) ){
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
