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

#include "ModelInitialization.h"
#include "Model.h"
#include "PopulationAgeStructure.h"
#include "Parameters.h"
#include "Transmission/transmission.h"
#include "Transmission/TransmissionModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ClinicalModel.h"
#include "interventions/InterventionManager.h"
#include "mon/Continuous.h"
#include "mon/init.h"
#include "util/ModelOptions.h"
#include "schema/scenario.h"

namespace OM {
using namespace util;
using Transmission::TransmissionModel;

Host::IncidenceParameters readIncidenceParameters(const Parameters& parameters) {
    Host::IncidenceParameters result;
    result.baseline_avail_shape_param = parameters[Parameter::BASELINE_AVAILABILITY_SHAPE];
    result.baseline_avail_offset = 0.5 * pow(result.baseline_avail_shape_param, 2);
    
    result.gamma_p=parameters[Parameter::GAMMA_P];
    result.Sinf=1-exp(-parameters[Parameter::NEG_LOG_ONE_MINUS_SINF]);
    result.Simm=parameters[Parameter::SIMM];
    result.EstarInv = 1.0/parameters[Parameter::E_STAR];
    result.Xstar_pInv = 1.0/parameters[Parameter::X_STAR_P];
    
    // Constraint on Gamma availability variance; log(1 + 0.649) is approximately 0.5.
    const double r_square_Gamma = 0.649;

    result.noPreErythrocytic = util::ModelOptions::option (util::NO_PRE_ERYTHROCYTIC);
    const bool negativeBinomial = util::ModelOptions::option (util::NEGATIVE_BINOMIAL_MASS_ACTION);
    result.vaccineGenotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    if (negativeBinomial) {
        result.kind = Host::IncidenceKind::NegativeBinomial;
        result.inf_rate_shape_param = (result.baseline_avail_shape_param+1.0) / (r_square_Gamma*result.baseline_avail_shape_param - 1.0);
        result.inf_rate_shape_param=std::max(result.inf_rate_shape_param, 0.0);
    } else if (util::ModelOptions::option (util::LOGNORMAL_MASS_ACTION)) {
        result.kind = Host::IncidenceKind::Lognormal;
        
        //! constant defining the constraint for the log Normal variance
        /// Used for the case where availability is assumed log Normally distributed
        double r_square_LogNormal = log(1.0+r_square_Gamma);
        
        result.inf_rate_shape_param = sqrt(r_square_LogNormal - 1.86 * pow(result.baseline_avail_shape_param, 2));
        result.inf_rate_shape_param=std::max(result.inf_rate_shape_param, 0.0);
        result.inf_rate_offset = 0.5 * pow(result.inf_rate_shape_param, 2);
        if( (std::isnan)(result.inf_rate_shape_param) ){
            throw util::xml_scenario_error( "bad parameter 16 (BASELINE_AVAILABILITY_SHAPE)" );
        }
    }else{
        if( util::ModelOptions::option( util::TRANS_HET ) ||
            util::ModelOptions::option( util::COMORB_TRANS_HET ) ||
            util::ModelOptions::option( util::TRANS_TREAT_HET ) ||
            util::ModelOptions::option( util::TRIPLE_HET ) )
        {
            result.kind = Host::IncidenceKind::Heterogeneity;
            cerr << "Warning: will use heterogeneity workaround." << endl;
        }
    }
    
    return result;
}
Transmission::HostSpeciesParameters readHostSpeciesParameters(const scnXml::Mosq& mosq) {
    Transmission::HostSpeciesParameters result;
    if(mosq.getAvailability().getDistr() != "const" && (util::ModelOptions::option(NEGATIVE_BINOMIAL_MASS_ACTION) || util::ModelOptions::option(LOGNORMAL_MASS_ACTION) || util::ModelOptions::option(COMORB_TRANS_HET)))
        throw util::xml_scenario_error( "HostSpeciesParameters::HostSpeciesParameters(): ModelOptions NEGATIVE_BINOMIAL_MASS_ACTION, LOGNORMAL_MASS_ACTION, COMORB_TRANS_HET are not compatible with ento-availability heterogeneity distributions other than \"const\"");

    result.entoAvailability = util::createSampler(1.0, mosq.getAvailability());
    result.probMosqBiting.setParams( mosq.getMosqProbBiting() );
    result.probMosqFindRestSite.setParams( mosq.getMosqProbFindRestSite() );
    result.probMosqSurvivalResting.setParams( mosq.getMosqProbResting() );
    return result;
}

void calculateAvailabilityPercentiles(Transmission::HostTransmissionParameters& parameters) {
    int nSamples = 100000;
    vector<double> samples(nSamples, 0.0);
    // Independent stream: percentile initialization must not affect simulation RNGs.
    util::MasterRng seed(0, 0);
    seed.seed(0, 0);
    util::LocalRng rng(seed);
    for (int i = 0; i < nSamples; ++i) {
        double avail = 0.0;
        for (const auto& species : parameters.species)
            avail += species.entoAvailabilityFactor * species.entoAvailability->sample(rng);
        samples[i] = avail;
    }
    sort(samples.begin(), samples.end());
    parameters.availabilityPercentiles.resize(100);
    for (int i = 0; i < 100; ++i)
        parameters.availabilityPercentiles[i] = samples[int(i*nSamples/100)];
    parameters.availabilityPercentiles[0] = 0.0;
}

void initialize(Model& model, scnXml::Scenario& scenario) {
    // 1) elements with no dependencies on other elements initialised here:
    WithinHost::Genotypes::init( scenario );
    model.masterRng.seed( scenario.getModel().getComputationParameters().getIseed(), 0 ); // Init RNG with Iseed
    util::ModelNameProvider modelNameProvider(scenario.getModel());

    // 2) elements depending on only elements initialised in (1):
    sim::init(scenario, modelNameProvider); // Also reads survey dates.

    // 3) elements depending on only elements initialised in (2).
    Parameters parameters( scenario.getModel().getParameters(), modelNameProvider ); // Depends on ModelNameProvider.
    util::ModelOptions::init( scenario.getModel().getModelOptions(), modelNameProvider ); // Depends on ModelNameProvider.

    // 4) elements depending on only elements initialised in (3).
    WithinHost::diagnostics::init( parameters, scenario ); // Depends on Parameters.

    // 5) elements depending on only elements initialised in (4).
    mon::initReporting( scenario ); // Reporting init depends on diagnostics and monitoring

    // Init models used by humans
    model.parameters.hostTransmission.availabilityByAge.set( scenario.getModel().getHuman().getAvailabilityToMosquitoes(), "availabilityToMosquitoes" );
    model.parameters.incidence = readIncidenceParameters(parameters);
    mon::Continuous::registerCallback("new infections", "\tnew infections", [&model](ostream& stream) {
        stream << '\t' << model.reporting.newInfections;
        model.reporting.newInfections = 0;
    });
    WithinHost::WHInterface::init( parameters, scenario );
    Clinical::ClinicalModel::init( parameters, scenario );
    Host::NeonatalMortality::init( scenario.getModel().getClinical() );
    AgeStructure::init( scenario.getDemography() );

    // Species parameters and vector RNG precede creation of the human population.
    size_t popSize = scenario.getDemography().getPopSize();

    auto& population = model.population;
    model.transmission.reset(Transmission::createTransmissionModel(scenario.getEntomology(), popSize, model.parameters.hostTransmission, model.masterRng));
    auto& transmission = *model.transmission;

    registerContinousPopulationCallbacks(model.parameters.hostTransmission);

    // Depends on transmission model (for species indexes):
    // MDA1D may depend on health system (too complex to verify)
    interventions::InterventionManager::init(scenario.getInterventions(), population, transmission );
    Clinical::ClinicalModel::setHS( scenario.getHealthSystem() ); // Depends on interventions, PK/PD (from humanPop)
    mon::initCohorts( scenario.getMonitoring() ); // Depends on interventions

    model.surveyOnlyNewEp = scenario.getMonitoring().getSurveyOptions().getOnlyNewEpisode();

}
}
