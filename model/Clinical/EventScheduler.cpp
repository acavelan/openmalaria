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

#include "schema/healthSystem.h"
#include "Clinical/EventScheduler.h"
#include "util/random.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include <schema/scenario.h>

#include <limits>

namespace OM { namespace Clinical {
    using namespace OM::util;


// -----  static init  -----

void esInit( ClinicalParams& params, const OM::Parameters& parameters, const scnXml::Clinical& clinical )
{
    if( sim::oneTS() != sim::oneDay() )
        throw util::xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day time step.");
    
    params.nonMalariaFevers = util::ModelOptions::option( util::NON_MALARIA_FEVERS );
    
    params.alpha = exp( -parameters[Parameter::CFR_NEG_LOG_ALPHA] );
    if( !(0.0<=params.alpha && params.alpha<=1.0) ){
        throw util::xml_scenario_error(
            "Clinical outcomes: propDeathsFirstDay should be within range [0,1]");
    }
    params.neg_v = -parameters[Parameter::CFR_SCALE_FACTOR];
    
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc2 =
            clinical.getNonMalariaFevers().get();
        if( !nmfDesc2.getPrNeedTreatmentMF().present() ||
            !nmfDesc2.getPrNeedTreatmentNMF().present() ){
            throw util::xml_scenario_error( "prNeedTreatmentMF and "
                "prNeedTreatmentNMF elements required in model->clinical->"
                "NonMalariaFevers" );
        }
        params.NMF_need_antibiotic.set( nmfDesc2.getPrNeedTreatmentNMF().get(), "prNeedTreatmentNMF" );
        params.MF_need_antibiotic.set( nmfDesc2.getPrNeedTreatmentMF().get(), "prNeedTreatmentMF" );
    }
}

void esSetParameters( ClinicalParams& params, const scnXml::HSEventScheduler& esData) {
    const scnXml::ClinicalOutcomes& coData = esData.getClinicalOutcomes();
    
    params.maxUCSeekingMemory = sim::fromDays(coData.getMaxUCSeekingMemory());
    params.uncomplicatedCaseDuration = sim::fromDays(coData.getUncomplicatedCaseDuration());
    params.complicatedCaseDuration = sim::fromDays(coData.getComplicatedCaseDuration());
    params.extraDaysAtRisk = sim::fromDays(coData.getComplicatedRiskDuration()) - params.complicatedCaseDuration;
    if( params.uncomplicatedCaseDuration < sim::fromDays(1)
	|| params.complicatedCaseDuration < sim::fromDays(1)
	|| params.maxUCSeekingMemory < sim::zero()
	|| params.extraDaysAtRisk + params.complicatedCaseDuration < sim::fromDays(1) // at risk at least 1 day
	|| params.extraDaysAtRisk > sim::zero()        // at risk longer than case duration
    ){
	throw util::xml_scenario_error(
            "Clinical outcomes: constraints on case/risk/memory duration not met (see documentation)");
    }
    
    params.cumDailyPrImmUCTS.reserve( coData.getDailyPrImmUCTS().size() );
    double cumP = 0.0;
    for( auto it = coData.getDailyPrImmUCTS().begin(); it != coData.getDailyPrImmUCTS().end(); ++it ){
        cumP += *it;
        params.cumDailyPrImmUCTS.push_back( cumP );
    }
    if( cumP < 0.999 || cumP > 1.001 ){
        throw util::xml_scenario_error( "Event scheduler: dailyPrImmUCTS seq must add up to 1" );
    }
    params.cumDailyPrImmUCTS.back() = 1.0;
    
    params.caseFatalityRate.scale( params.alpha );
    
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !esData.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of healthSystem->EventScheduler required");
        }
        const scnXml::HSESNMF& nmfDesc = esData.getNonMalariaFevers().get();
        
        double pT = nmfDesc.getPrTreatment();
        params.logOddsAbBase = log( pT / (1.0 - pT) );
        params.logOddsAbNegTest = log(nmfDesc.getEffectNegativeTest());
        params.logOddsAbPosTest = log(nmfDesc.getEffectPositiveTest());
        params.logOddsAbNeed = log(nmfDesc.getEffectNeed());
        params.logOddsAbInformal = log(nmfDesc.getEffectInformal());
        params.oneMinusEfficacyAb = 1.0 - nmfDesc.getTreatmentEfficacy();
        params.severeNmfMortality.set( nmfDesc.getCFR(), "CFR" );
    }
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double tSF) :
        pgState (Episode::NONE),
        caseStartTime (sim::never()),
        timeOfRecovery (sim::never()),
        timeLastTreatment (sim::never()),
        previousDensity (numeric_limits<double>::quiet_NaN())
{
    if( tSF != 1.0 ){
	// p(treatment seeking) is part of tree & handled generically, so we
	// don't have a way of modifying it.
	throw xml_scenario_error("treatment seeking heterogeneity not supported");
    }
}


// -----  other methods  -----

bool ClinicalEventScheduler::isExistingCase(const ClinicalParams& params) {
        // If treated in the recent past:
        return sim::now() > timeLastTreatment && sim::now() <= timeLastTreatment + params.healthSystemMemory;
}

void ClinicalEventScheduler::doClinicalUpdate (Human& human, double ageYears, WithinHost::Pathogenesis::StatePair &pg,
                                               const ClinicalParams& params){
    WHInterface& withinHostModel = *human.withinHostModel;
    Episode::State newState = static_cast<Episode::State>( pg.state );
    util::streamValidate( (newState << 16) & pgState );
    
    if ( sim::ts0() == timeOfRecovery ) {
	if( pgState & Episode::DIRECT_DEATH ){
	    // Human dies this time step (last day of risk of death)
	    doomed = DOOMED_COMPLICATED;
	    
	    latestReport.update (human, pgState, params.healthSystemMemory);
        } else {
	    if ( pgState & Episode::COMPLICATED ) {
                const double pSequelae = params.pSequelaeInpatient.eval( ageYears );
                mon::recordStat(mon::measure("expectedSequelae"), human, pSequelae);
		if( human.rng.uniform_01() < pSequelae ){
		    pgState = Episode::State (pgState | Episode::SEQUELAE);
                }else{
		    pgState = Episode::State (pgState | Episode::RECOVERY);
                }
	    } else {
		pgState = Episode::State (pgState | Episode::RECOVERY);
            }
	    // report bout, at conclusion of episode:
	    latestReport.update (human, pgState, params.healthSystemMemory);
	    
	    // Individual recovers (and is immediately susceptible to new cases)
	    pgState = Episode::NONE;	// recovery (reset to healthy state)
	    
	    // And returns to transmission (if was removed)
	    human.perHostTransmission.outsideTransmission = false;
	}
    }
    
    // The bug fixed here is that indirectMortality was not set if no event occurred
    bool indirectMortality = false;
    if( params.indirectMortBugfix ) indirectMortality = pg.indirectMortality;
    
    if ( newState & Episode::SICK ){
        // we have some new case: is it severe/complicated?
        if ( newState & Episode::COMPLICATED ){
            if ( pgState & Episode::COMPLICATED ) {
                // previously severe: no events happen for course of medication
            } else {
                // previously healthy or UC: progress to severe
                pgState = Episode::State (pgState | newState | Episode::RUN_CM_TREE);
                indirectMortality = pg.indirectMortality;
                caseStartTime = sim::ts0();
            }
        } else {
            // uncomplicated case (UC/UC2/NMF): is it new?
            if (pgState & Episode::SICK) {
                // previously UC; nothing to do
            }else{
                if( caseStartTime < sim::ts0() ) {
                    // new UC case
                    pgState = Episode::State (pgState | newState | Episode::RUN_CM_TREE);
                    indirectMortality = pg.indirectMortality;
                    
                    double uVariate = human.rng.uniform_01();
                    size_t i = 0;       // units: days
                    for(; i < params.cumDailyPrImmUCTS.size(); ++i){
                        if( uVariate < params.cumDailyPrImmUCTS[i] ){
                            goto gotDelay;
                        }
                    }
                    assert(false);      // should have uVariate < 1 = params.cumDailyPrImmUCTS[len-1]
                    gotDelay:
                    // set start time: current time plus length of delay (days)
                    caseStartTime = sim::ts0() + sim::fromDays(i);
                }
            }
        }
        
        if (indirectMortality && doomed == NOT_DOOMED)
            doomed = -sim::oneTS(); // start indirect mortality countdown
    }
    
    if( caseStartTime == sim::ts0() && (pgState & Episode::RUN_CM_TREE) ){
        // OK, we're about to run the CM tree
        pgState = Episode::State (pgState & ~Episode::RUN_CM_TREE);
        
        // If last treatment prescribed was in recent memory, consider second line.
        if( timeLastTreatment + params.healthSystemMemory > sim::ts0() ){
            pgState = Episode::State (pgState | Episode::SECOND_CASE);
        }
	
	CMDTOut auxOut = esExecute(
	    CMHostData( human, ageYears, pgState ), params );
	
        if( auxOut.treated ){	// I.E. some treatment was given
            timeLastTreatment = sim::ts0();
            if( pgState & Episode::COMPLICATED ){
                mon::recordEvent(mon::measure("nTreatments3"), human);
            }else{
                if( pgState & Episode::SECOND_CASE ){
                    mon::recordEvent(mon::measure("nTreatments2"), human);
                }else{
                    mon::recordEvent(mon::measure("nTreatments1"), human);
                }
            }
        }
        if( auxOut.screened ){
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
        }
	
	if ( true /*FIXME auxOut.hospitalisation != CMAuxOutput::NONE*/ ) {	// in hospital
	    pgState = Episode::State (pgState | Episode::EVENT_IN_HOSPITAL);
	    
	    /*FIXME if( auxOut.hospitalisation == CMAuxOutput::DELAYED )
		++caseStartTime;*/
	}
	
	// Case fatality rate (first day of illness)
	// P(death) is some fixed input scaled by age-specific CFR.
	if( (pgState & Episode::COMPLICATED)
	    && !(pgState & Episode::DIRECT_DEATH)
	) {
            const bool inHospital = true /*FIXME auxOut.hospitalisation == CMAuxOutput::IMMEDIATE */;
	    double pDeath = params.caseFatalityRate.eval( ageYears );
	    // community fatality rate when not in hospital or delayed hospital entry
            if( !inHospital )
                pDeath = params.getCommunityCFR( pDeath );
            mon::recordStat(mon::measure("expectedDirectDeaths"), human, pDeath);
            if( inHospital )
                mon::recordStat(mon::measure("expectedHospitalDeaths"), human, pDeath);
	    if (human.rng.uniform_01() < pDeath) {
		pgState = Episode::State (pgState | Episode::DIRECT_DEATH | Episode::EVENT_FIRST_DAY);
		// Human is killed at end of time at risk
		//timeOfRecovery += params.extraDaysAtRisk;	(no point updating; will be set later: ATORWD)
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	}
	
	if( params.nonMalariaFevers ){
            if( (pgState & Episode::SICK) && !(pgState & Episode::COMPLICATED) ){
                // Have a NMF or UC malaria case
                
                /** Given a non-malaria fever, return the probability of it requiring
                 * treatment. */
                bool isMalarial = pgState & Episode::MALARIA;
                double pNeedTreat = isMalarial ?
                        params.MF_need_antibiotic.eval( ageYears ) :
                        params.NMF_need_antibiotic.eval( ageYears );
                bool needTreat = human.rng.bernoulli(pNeedTreat);
                
                // Calculate chance of antibiotic administration:
                double pTreatment = 0.0;
                                double treatmentEffectMult = 1.0;
                
                if( human.rng.uniform_01() < pTreatment ){
                    /*FIXME: impossible due to above; NMF output removed
                    Survey::current().addInt( Report::MI_NMF_TREATMENTS, human, 1 );
                    treatmentEffectMult = params.oneMinusEfficacyAb;
                    */
                }
                
                // In a severe NMF case (only when not malarial), there is a
                // chance of death:
                if( needTreat ){
                    double pDeath = params.severeNmfMortality.eval( ageYears ) * treatmentEffectMult;
                    if( human.rng.uniform_01() < pDeath ){
                        pgState = Episode::State (pgState | Episode::DIRECT_DEATH);
                    }
                }
            }
        }
    } else {
	// No new event (haven't changed state this time step).
	
	// Case fatality rate (subsequent days)
	// Complicated case & at risk of death (note: params.extraDaysAtRisk <= 0)
	if( (pgState & Episode::COMPLICATED)
	    && !(pgState & Episode::DIRECT_DEATH)
	    && (sim::ts0() < timeOfRecovery + params.extraDaysAtRisk)
	) {
	    // In complicated episodes, S(t), the probability of survival on
	    // subsequent days t, is described by log(S(t)) = -v(Y(t)/Y(t-1)),
	    // for parasite density Y(t). v_neg below is -v.
            // TODO: this model should be revisited, and if possible placed
            // within the WHFalciparum class (or a subclass).
	    if( withinHostModel.getTotalDensity() > 0.0 ){	// avoid division by zero
		double parasiteReductionEffect = withinHostModel.getTotalDensity() / previousDensity;
		double pDeath = 1.0 - exp( params.neg_v * parasiteReductionEffect );
		// community fatality rate when not in hospital
		if( !(pgState & Episode::EVENT_IN_HOSPITAL) )
		    pDeath = params.getCommunityCFR( pDeath );
                mon::recordStat(mon::measure("expectedDirectDeaths"), human, pDeath);
                if( pgState & Episode::EVENT_IN_HOSPITAL )
                    mon::recordStat(mon::measure("expectedHospitalDeaths"), human, pDeath);
		if (human.rng.uniform_01() < pDeath) {
		    pgState = Episode::State (pgState | Episode::DIRECT_DEATH);
		    // Human is killed at end of time at risk
		    timeOfRecovery = timeOfRecovery + params.extraDaysAtRisk;	// may be re-set later (see ATORWD)
		}
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	}
    }
    
    // Start of case. Not necessarily start of sickness due to treatment-seeking
    // delays and travel time.
    if( caseStartTime == sim::ts0() ){
	// Patients in hospital are removed from the transmission cycle.
	// This should have an effect from the start of the next time step.
	// NOTE: This is not very accurate, but considered of little importance.
	if (pgState & Episode::EVENT_IN_HOSPITAL)
	    human.perHostTransmission.outsideTransmission = true;
	
	if (pgState & Episode::COMPLICATED) {
	    // params.complicatedCaseDuration should to some respects be associated
	    // with medication duration, however ongoing medications after
	    // exiting hospital are OK and medications terminating before the
	    // end of hospitalisation shouldn't matter too much if the person
	    // can't recieve new infections due to zero transmission in hospital.
	    timeOfRecovery = sim::ts0() + params.complicatedCaseDuration;
	    // Time should be adjusted to end of at-risk period when patient dies:
	    if( pgState & Episode::DIRECT_DEATH )	// death may already have been determined
		timeOfRecovery = timeOfRecovery + params.extraDaysAtRisk;	// ATORWD (search keyword)
	} else {
	    timeOfRecovery = sim::ts0() + params.uncomplicatedCaseDuration;
	}
    }
    
    // Remove on first models...
    if( timeLastTreatment == sim::ts0() )
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_TREATMENT);
    
    if( pgState & Episode::SICK )
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_BOUT);
}


void ClinicalEventScheduler::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    int s;
    s & stream;
    pgState = Episode::State(s);
    caseStartTime & stream;
    timeOfRecovery & stream;
    timeLastTreatment & stream;
    previousDensity & stream;
}
void ClinicalEventScheduler::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    pgState & stream;
    caseStartTime & stream;
    timeOfRecovery & stream;
    timeLastTreatment & stream;
    previousDensity & stream;
}

} }
