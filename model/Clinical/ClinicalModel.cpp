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

#include "Clinical/ClinicalModel.h"

#include "Clinical/EventScheduler.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/DecisionTree5Day.h"
#include "Host/NeonatalMortality.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"

namespace OM { namespace Clinical {

// -----  ClinicalParams  -----

ClinicalParams::ClinicalParams( const Parameters& parameters, const scnXml::Scenario& scenario ) {
    const scnXml::Clinical& clinical = scenario.getModel().getClinical();
    try{
        indirectMortBugfix = util::ModelOptions::option (util::INDIRECT_MORTALITY_FIX);
        //NOTE: if changing XSD, this should not have a default unit:
        healthSystemMemory = UnitParse::readShortDuration(clinical.getHealthSystemMemory(), UnitParse::STEPS);
        oddsRatioThreshold = exp( parameters[Parameter::LOG_ODDS_RATIO_CF_COMMUNITY] );
        InfantMortality::init( parameters );
        cfrPfUseHospital = util::ModelOptions::option( util::CFR_PF_USE_HOSPITAL );
    }catch( const util::format_error& e ){
        throw util::xml_scenario_error( string("model/clinical/healthSystemMemory: ").append(e.message()) );
    }
    
    if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER)){
        eventScheduler = true;
        esInit( *this, parameters, clinical );
    }else{
        if( scenario.getHealthSystem().getImmediateOutcomes().present() ){
            immOutcomes = true;
            
            if( util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS) ){
                cerr << "Deprecation warning: healthSystem: use of ImmediateOutcomes can be replaced by the more flexible DecisionTree5Day (optional)" << endl;
            }
        }
        // else: decision tree 5 day
    }
}

void ClinicalParams::loadHealthSystem( const scnXml::HealthSystem& healthSystem, const interventions::InterventionRegistry& reg ){
    caseFatalityRate.set( healthSystem.getCFR(), "CFR" );
    pSequelaeInpatient.set( healthSystem.getPSequelaeInpatient(), "pSequelaeInpatient" );
    if( eventScheduler ){
        if( !healthSystem.getEventScheduler().present() ){
            throw util::xml_scenario_error ("Expected EventScheduler "
                "section in healthSystem data (initial or intervention)");
        }
        loadHSEventScheduler( *this, healthSystem.getEventScheduler().get(), reg );
    }else if( immOutcomes ){
        if ( !healthSystem.getImmediateOutcomes().present() ){
            throw util::xml_scenario_error ("Expected ImmediateOutcomes "
                "section in healthSystem data (initial or intervention)");
        }
        loadHSImmediateOutcomes( *this, healthSystem.getImmediateOutcomes().get(), reg );
    }else{
        if( !healthSystem.getDecisionTree5Day().present() ){
            throw util::xml_scenario_error ("Expected DecisionTree5Day "
                "section in healthSystem data (initial or intervention)");
        }
        loadHSDT5Day( *this, healthSystem.getDecisionTree5Day().get(), reg );
    }
}

double ClinicalParams::getCommunityCFR (double caseFatalityRatio) const{
    double x = caseFatalityRatio * oddsRatioThreshold;
    return x / (1 - caseFatalityRatio + x);
}

// -----  static methods  -----

unique_ptr<ClinicalModel> ClinicalModel::createClinicalModel (double tSF, const ClinicalParams& params) {
    if (params.eventScheduler){
        return make_unique<ClinicalEventScheduler>(tSF);
    }else if( params.immOutcomes ){
        return make_unique<ImmediateOutcomes>(tSF);
    }else{
        return make_unique<DecisionTree5Day>( tSF );
    }
}



// -----  non-static construction, destruction and checkpointing  -----

ClinicalModel::ClinicalModel () :
    doomed(NOT_DOOMED)
{}
ClinicalModel::~ClinicalModel () {
  // latestReport is reported, if any, by destructor
}


// -----  other non-static methods  -----

bool ClinicalModel::isDead( SimTime age ){
    if( age >= sim::maxHumanAge())       // too old (reached age limit)
        doomed = DOOMED_TOO_OLD;
    return doomed > NOT_DOOMED;	// killed by some means
}

void ClinicalModel::update (Human& human, double ageYears, bool newBorn, Host::NeonatalMortality& neonatal,
                            const ClinicalParams& params) {
    if (doomed < NOT_DOOMED)	// Countdown to indirect mortality
        doomed -= sim::oneTS();
    
    //indirect death: if this human's about to die, don't worry about further episodes:
    if (doomed <= DOOMED_EXPIRED) {	//clinical bout 6 intervals before
        mon::recordEvent(mon::measure("nIndDeaths"), human);
        doomed = DOOMED_INDIRECT;
        return;
    }
    if(newBorn /* i.e. first update since birth */) {
        // Chance of neonatal mortality:
        if (neonatal.eventNeonatalMortality(human.rng)) {
            mon::recordEvent(mon::measure("nIndDeaths"), human);
            doomed = DOOMED_NEONATAL;
            return;
        }
    }
    
    const bool isDoomed = doomed != NOT_DOOMED;
    WithinHost::Pathogenesis::StatePair pg = human.withinHostModel->determineMorbidity( human, ageYears, isDoomed );
    latestState = static_cast<Episode::State>(pg.state);

    doClinicalUpdate (human, ageYears, pg, params);
}

void ClinicalModel::updateInfantDeaths( SimTime age ){
    // update array for the infant death rates
    if (age < sim::oneYear()){
        size_t index = age / sim::oneTS();
        
        // Testing doomed == DOOMED_NEXT_TS gives very slightly different results than
        // testing doomed == DOOMED_INDIRECT (due to above if(..))
        bool isDoomed = doomed == DOOMED_COMPLICATED
                || doomed == DOOMED_NEXT_TS
                || doomed == DOOMED_NEONATAL;
        InfantMortality::reportRisk( index, isDoomed );
    }
}


void ClinicalModel::checkpoint (istream& stream) {
    latestReport & stream;
    doomed & stream;
}
void ClinicalModel::checkpoint (ostream& stream) {
    latestReport & stream;
    doomed & stream;
}



// ——— infant mortality ———

// Infant death summaries (checkpointed).
vector<int> infantDeaths;
vector<int> infantIntervalsAtRisk;

/// Non-malaria mortality in under 1year olds.
/// Set by init ()
double nonMalariaMortality;

void InfantMortality::init( const OM::Parameters& parameters ){
    infantDeaths.resize(sim::stepsPerYear());
    infantIntervalsAtRisk.resize(sim::stepsPerYear());
    nonMalariaMortality=parameters[Parameter::NON_MALARIA_INFANT_MORTALITY];
}

void InfantMortality::preMainSimInit() {
    std::fill(infantDeaths.begin(), infantDeaths.end(), 0);
    std::fill(infantIntervalsAtRisk.begin(), infantIntervalsAtRisk.end(), 0);
}

void InfantMortality::staticCheckpoint(istream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
}
void InfantMortality::staticCheckpoint (ostream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
}

void InfantMortality::reportRisk(size_t index, bool isDoomed) {
    infantIntervalsAtRisk[index] += 1;     // baseline
    if (isDoomed)
        infantDeaths[index] += 1;  // deaths
}

double InfantMortality::allCause(){
    double infantPropSurviving=1.0;       // use to calculate proportion surviving
    for( size_t i = 0; i < sim::stepsPerYear(); i += 1 ){
        // multiply by proportion of infants surviving at each interval
        infantPropSurviving *= double(infantIntervalsAtRisk[i] - infantDeaths[i])
            / double(infantIntervalsAtRisk[i]);
    }
    // Child deaths due to malaria (per 1000), plus non-malaria child deaths. Deaths per 1000 births is the return unit.
    return (1.0 - infantPropSurviving) * 1000.0 + nonMalariaMortality;
}


} }
