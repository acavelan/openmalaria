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

// model/xml/ conversion layer — health system / clinical side of the host
// modules: non-malaria fevers (PathogenesisModel), neonatal mortality
// diagnostic selection and within-host treatment options from scnXml.
// Moved verbatim from Host/WithinHost/Pathogenesis/PathogenesisModel.cpp,
// Host/NeonatalMortality.cpp and Host/WithinHost/Treatments.cpp. State stays
// defined in the original files (extern re-declared below).

#include "Global.h"
#include "Parameters.h"
#include "Host/NeonatalMortality.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "Host/WithinHost/Pathogenesis/Submodels.h"
#include "Host/WithinHost/Treatments.h"
#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "schema/healthSystem.h"

#include <cmath>
#include <string>

using namespace std;

namespace OM {
namespace WithinHost {
using namespace OM::util;
namespace Pathogenesis {

// State defined in PathogenesisModel.cpp (external linkage):
extern double pg_indirRiskCoFactor, pg_severeMalThreshold, pg_comorbIntercept,
    pg_inv_critAgeComorb;
extern AgeGroupInterpolator pg_NMF_incidence;
extern bool opt_predetermined_episodes, opt_mueller_pres_model;

void PathogenesisModel::init( const Parameters& parameters, const scnXml::Clinical& clinical, bool nmfOnly ){
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc = clinical.getNonMalariaFevers().get();
        pg_NMF_incidence.set( nmfDesc.getIncidence(), "incidence" );
    }
    if( nmfOnly ) return;
    
    pg_indirRiskCoFactor = 1 - exp(-parameters[Parameter::INDIRECT_RISK_COFACTOR]);
    pg_severeMalThreshold = parameters[Parameter::SEVERE_MALARIA_THRESHHOLD] + 1;
    pg_comorbIntercept = 1 - exp(-parameters[Parameter::COMORBIDITY_INTERCEPT]);
    pg_inv_critAgeComorb = 1 / parameters[Parameter::CRITICAL_AGE_FOR_COMORBIDITY];

    if (util::ModelOptions::option (util::PREDETERMINED_EPISODES)) {
        opt_predetermined_episodes = true;
        //no separate init:
        PyrogenPathogenesis::init( parameters );
    } else {
        if (util::ModelOptions::option (util::MUELLER_PRESENTATION_MODEL)){
            opt_mueller_pres_model = true;
            MuellerPathogenesis::init( parameters );
        }else{
            PyrogenPathogenesis::init( parameters );
        }
    }
}

}

TreatmentId Treatments::addTreatment( const scnXml::TreatmentOption& desc ){
    Treatments treatment( desc );
    TreatmentId id( treatments.size() );
    treatments.push_back( std::move(treatment) );
    return id;
}

Treatments::Stages stageFromString( const std::string& str ){
    if( str == "liver" ) return Treatments::LIVER;
    if( str == "blood" ) return Treatments::BLOOD;
    if( str == "both" ) return Treatments::BOTH;
    throw util::format_error( std::string("stage must be liver, blood or both, not ").append(str) );
}
Treatments::Treatments( const scnXml::TreatmentOption& elt ) :
    TriggeredDeployments(elt), timeLiver(sim::zero()), timeBlood(sim::zero())
{
    for( auto it = elt.getClearInfections().begin(), end = elt.getClearInfections().end();
        it != end; ++it )
    {
        try{
            //NOTE: if changing XSD, this should not be called "timesteps" or have a default unit
            SimTime len = UnitParse::readShortDuration( it->getTimesteps(), UnitParse::STEPS );
            if( len < -sim::oneTS() || len == sim::zero() ){ 
                throw util::format_error( "timesteps must be ≥ 1 or have special value -1" );
            }
            Stages stage = stageFromString( it->getStage() );
            if( stage & LIVER ){
                if( timeLiver != sim::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of liver stage effect" );
                timeLiver = len;
            }
            if( stage & BLOOD ){
                if( timeBlood != sim::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of blood stage effect" );
                timeBlood = len;
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string(".../clearInfections/timesteps: ").append(e.message()) );
        }
    }
}

}

namespace Host {
using namespace OM::util;

// State defined in NeonatalMortality.cpp (external linkage):
extern std::vector<double> prevByGestationalAge;
extern const WithinHost::Diagnostic* neonatalDiagnostic;

void NeonatalMortality::init( const scnXml::Clinical& clinical ){
    SimTime fiveMonths = sim::fromDays( 5 * 30 );
    prevByGestationalAge.assign( sim::inSteps(fiveMonths), 0.0 );
    
    if( clinical.getNeonatalMortality().present() ){
        neonatalDiagnostic = &WithinHost::diagnostics::get(
            clinical.getNeonatalMortality().get().getDiagnostic() );
    }else{
        //NOTE: this is a compatibility option for older scenarios
        if( CommandLine::option( CommandLine::DEPRECATION_WARNINGS ) ){
            cerr << "Deprecation warning: specification of the diagnostic "
                "used by the Neonatal Mortality model is recommended "
                "(model/clinical/neonatalMortality)" << endl;
        }
        neonatalDiagnostic = &WithinHost::diagnostics::monitoringDiagnostic();
    }
}
}
}
