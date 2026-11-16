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
#include "mon/Monitoring.h"
#include "util/errors.h"

namespace OM { namespace Clinical {

// -----  static methods  -----

unique_ptr<ClinicalModel> ClinicalModel::createClinicalModel (double tSF, const ClinicalParameters& config) {
    if (config.kind == ClinicalParameters::Kind::EventScheduler){
        return unique_ptr<ClinicalModel>(new ClinicalEventScheduler (tSF, config));
    }else if( config.kind == ClinicalParameters::Kind::ImmediateOutcomes ){
        return unique_ptr<ClinicalModel>(new ImmediateOutcomes (tSF, config));
    }else{
        return unique_ptr<ClinicalModel>(new DecisionTree5Day( tSF, config ));
    }
}

double ClinicalModel::getCommunityCFR (double caseFatalityRatio) const{
    double x = caseFatalityRatio * config.oddsRatioThreshold;
    return x / (1 - caseFatalityRatio + x);
}



// -----  non-static construction, destruction and checkpointing  -----

ClinicalModel::ClinicalModel (const ClinicalParameters& config) :
    doomed(NOT_DOOMED), config(config)
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

void ClinicalModel::update (Human& human, double ageYears, bool newBorn, double neonatalRisk, const interventions::CohortRemovalRules& removalRules) {
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
        if (human.rng.uniform_01() <= neonatalRisk) {
            mon::recordEvent(mon::measure("nIndDeaths"), human);
            doomed = DOOMED_NEONATAL;
            return;
        }
    }
    
    const bool isDoomed = doomed != NOT_DOOMED;
    WithinHost::Pathogenesis::StatePair pg = human.withinHostModel->determineMorbidity( human, ageYears, isDoomed );
    latestState = static_cast<Episode::State>(pg.state);

    doClinicalUpdate (human, ageYears, pg, removalRules);
}

void ClinicalModel::updateInfantDeaths(SimTime age, InfantMortality& reporting) {
    // update array for the infant death rates
    if (age < sim::oneYear()){
        size_t index = age / sim::oneTS();
        
        // Testing doomed == DOOMED_NEXT_TS gives very slightly different results than
        // testing doomed == DOOMED_INDIRECT (due to above if(..))
        bool isDoomed = doomed == DOOMED_COMPLICATED
                || doomed == DOOMED_NEXT_TS
                || doomed == DOOMED_NEONATAL;
        reporting.reportRisk(index, isDoomed);
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

void InfantMortality::reset() {
    std::fill(deaths.begin(), deaths.end(), 0);
    std::fill(intervalsAtRisk.begin(), intervalsAtRisk.end(), 0);
}

void InfantMortality::reportRisk(size_t index, bool isDoomed) {
    intervalsAtRisk[index] += 1;     // baseline
    if (isDoomed)
        deaths[index] += 1;  // deaths
}

double InfantMortality::allCause() const {
    double infantPropSurviving=1.0;       // use to calculate proportion surviving
    for( size_t i = 0; i < intervalsAtRisk.size(); i += 1 ){
        // multiply by proportion of infants surviving at each interval
        infantPropSurviving *= double(intervalsAtRisk[i] - deaths[i])
            / double(intervalsAtRisk[i]);
    }
    // Child deaths due to malaria (per 1000), plus non-malaria child deaths. Deaths per 1000 births is the return unit.
    return (1.0 - infantPropSurviving) * 1000.0 + nonMalariaMortality;
}


} }
