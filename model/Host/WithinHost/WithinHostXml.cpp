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

#include "Host/WithinHost/WithinHostXml.h"
#include "Host/WithinHost/WithinHostParameters.h"
#include "Host/WithinHost/WHVivax.h"
#include "Host/WithinHost/Infection/DescriptiveInfection.h"
#include "Host/WithinHost/Infection/EmpiricalInfection.h"
#include "Host/WithinHost/Infection/MolineauxInfection.h"
#include "Host/WithinHost/Infection/PennyInfection.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "PkPd/LSTMModel.h"
#include "ScenarioXml.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"

namespace OM::xml {
using namespace WithinHost;
using namespace util;
namespace {
void initMueller(const Parameters& parameters, Pathogenesis::PathogenesisParameters& config){
    config.rateMultiplier_31 = parameters[Parameter::MUELLER_RATE_MULTIPLIER]
        * sim::yearsPerStep();
    config.densityExponent_32 = parameters[Parameter::MUELLER_DENSITY_EXPONENT];
}

void initPyrogen(const Parameters& parameters, Pathogenesis::PathogenesisParameters& config){
    config.initPyroThres = parameters[Parameter::Y_STAR_0];

    double delt = 1.0 / 11;
    double smuY = -log(0.5) /
        (sim::stepsPerYear() * parameters[Parameter::Y_STAR_HALF_LIFE]);
    config.b = -smuY * delt;

    config.Ystar2_13 = parameters[Parameter::Y_STAR_SQ];

    //alpha: factor determining increase in pyrogenic threshold
    double alpha14 = parameters[Parameter::ALPHA];
    config.a = alpha14 * sim::oneTS() * delt;

    //Ystar1: critical value of parasite density in determing increase in pyrog t
    config.Ystar1_26 = parameters[Parameter::Y_STAR_1];
}

void initPathogenesis(const Parameters& parameters, const scnXml::Clinical& clinical, bool nmfOnly, Pathogenesis::PathogenesisParameters& config){
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc = clinical.getNonMalariaFevers().get();
        xml::setAgeGroups(config.pg_NMF_incidence, nmfDesc.getIncidence(), "incidence" );
    }
    if( nmfOnly ) return;

    config.pg_indirRiskCoFactor = 1 - exp(-parameters[Parameter::INDIRECT_RISK_COFACTOR]);
    config.pg_severeMalThreshold = parameters[Parameter::SEVERE_MALARIA_THRESHHOLD] + 1;
    config.pg_comorbIntercept = 1 - exp(-parameters[Parameter::COMORBIDITY_INTERCEPT]);
    config.pg_inv_critAgeComorb = 1 / parameters[Parameter::CRITICAL_AGE_FOR_COMORBIDITY];

    if (util::ModelOptions::option (util::PREDETERMINED_EPISODES)) {
        config.kind = Pathogenesis::PathogenesisParameters::Kind::Predetermined;
        //no separate init:
        initPyrogen( parameters, config );
    } else {
        if (util::ModelOptions::option (util::MUELLER_PRESENTATION_MODEL)){
            config.kind = Pathogenesis::PathogenesisParameters::Kind::Mueller;
            initMueller( parameters, config );
        }else{
            initPyrogen( parameters, config );
        }
    }
}

void initFalciparum( const OM::Parameters& parameters, const scnXml::Model& model, WithinHostParameters& config ) {
    config.sigma_i=sqrt(parameters[Parameter::SIGMA_I_SQ]);
    config.immPenalty_22=1-exp(parameters[Parameter::IMMUNITY_PENALTY]);
    config.immEffectorRemain=exp(-parameters[Parameter::IMMUNE_EFFECTOR_DECAY]);
    config.asexImmRemain=exp(-parameters[Parameter::ASEXUAL_IMMUNITY_DECAY]);

    // calculate inverses here, so we can use multiplication later (faster):
    config.invCumulativeYstar = 1.0 / parameters[Parameter::CUMULATIVE_Y_STAR];
    config.invCumulativeHstar = 1.0 / parameters[Parameter::CUMULATIVE_H_STAR];
    config.alpha_m = 1.0 - exp(-parameters[Parameter::NEG_LOG_ONE_MINUS_ALPHA_M]);
    config.decayM = parameters[Parameter::DECAY_M];

    config.y_lag_len = sim::inSteps(sim::fromDays(20)) + 1;

    //NOTE: should also call cleanup() on the PathogenesisModel, but it only frees memory which the OS does anyway
    initPathogenesis( parameters, model.getClinical(), false, config.pathogenesis );

    // Default value to use if no <parameters> element is present in XML.
    SimTime latentP = UnitParse::readShortDuration("15d", UnitParse::DAYS);
    try
    {
        if (model.getParameters().present())
        {
            latentP = UnitParse::readShortDuration(
                model.getParameters().get().getLatentp(), UnitParse::NONE);
        }
    }
    catch( const util::format_error& e )
    {
            throw util::xml_scenario_error( string("model/parameters/latentP: ").append(e.message()) );
    }

    Infection::init( latentP );
}

void initCommon( const scnXml::Scenario& scenario, WithinHostParameters& config ){
    const scnXml::Human human = scenario.getModel().getHuman();
    if( !human.getWeight().present() ){
        // Technically this is needed by the PK/PD and Molineaux models
        throw util::xml_scenario_error( "model->human->weight element required by certain models" );
    }
    xml::setAgeGroups(config.massByAge, human.getWeight().get(), "weight" );
    config.hetMassMultStdDev = human.getWeight().get().getMultStdDev();
    // hetWeightMult must be large enough that birth weight is at least 0.5 kg:
    config.minHetMassMult = 0.5 / config.massByAge.eval( 0.0 );

    config.reportInfectedOrPatentInfected = mon::isUsed(mon::measure("totalInfs")) ||
        mon::isUsed(mon::measure("totalPatentInf"));

    PkPd::LSTMModel::init( scenario );
}

void initDescriptive(WithinHostParameters& config){
    config.reportPatentInfected = (
                                mon::isUsed(mon::measure("totalPatentInf"))
                                || mon::isUsed(mon::measure("totalPatentInf_Imported"))
                                || mon::isUsed(mon::measure("totalPatentInf_Introduced"))
                                || mon::isUsed(mon::measure("totalPatentInf_Indigenous"))
                            );
}

}

void initWithinHost( const OM::Parameters& parameters, const scnXml::Scenario& scenario, WithinHostParameters& config ) {
    config.opt_vaccine_genotype = ModelOptions::option(VACCINE_GENOTYPE);
    config.reportInfectionsByGenotype = mon::isUsed(mon::measure("nInfectByGenotype")) ||
        mon::isUsed(mon::measure("nPatentByGenotype")) ||
        mon::isUsed(mon::measure("logDensByGenotype"));

    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        config.kind = WithinHostParameters::Kind::Vivax;
        WHVivax::init( parameters, scenario.getModel() );
        initPathogenesis(parameters, scenario.getModel().getClinical(), true, config.pathogenesis);
    }else{
        initFalciparum( parameters, scenario.getModel(), config );

        if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL)) {
            config.kind = WithinHostParameters::Kind::Dummy;
        } else if (util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
            config.kind = WithinHostParameters::Kind::Empirical;
            EmpiricalInfection::init();    // 1-day time step check
        } else if (util::ModelOptions::option (util::MOLINEAUX_WITHIN_HOST_MODEL)) {
            config.kind = WithinHostParameters::Kind::Molineaux;
            MolineauxInfection::init( parameters );
        } else if (util::ModelOptions::option (util::PENNY_WITHIN_HOST_MODEL)) {
            config.kind = WithinHostParameters::Kind::Penny;
            PennyInfection::init();
        } else {
            DescriptiveInfection::init( parameters );      // 5-day time step check
        }

        if( config.kind != WithinHostParameters::Kind::Descriptive ){
            initCommon( scenario, config );
        }else{
            initDescriptive(config);
        }
    }
}

}
