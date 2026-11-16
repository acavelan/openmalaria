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
#include "Host/WithinHost/WithinHostXml.h"
#include "ScenarioXml.h"
#include "Model.h"
#include "PopulationAgeStructure.h"
#include "Parameters.h"
#include "Transmission/TransmissionModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/ClinicalXml.h"
#include "interventions/Interventions.h"
#include "mon/Continuous.h"
#include "mon/init.h"
#include "util/ModelOptions.h"
#include "schema/scenario.h"

namespace OM {
using namespace util;
using Transmission::TransmissionModel;

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
    xml::setAgeGroups(model.parameters.hostTransmission.availabilityByAge, scenario.getModel().getHuman().getAvailabilityToMosquitoes(), "availabilityToMosquitoes");
    model.parameters.incidence = Host::makeIncidenceParameters(xml::readIncidenceConfiguration(parameters));
    mon::Continuous::registerCallback(model.reporting.continuous, "new infections", "\tnew infections", [&model](ostream& stream) {
        stream << '\t' << model.reporting.newInfections;
        model.reporting.newInfections = 0;
    });
    xml::initWithinHost( parameters, scenario, model.parameters.withinHost );
    xml::initClinical( parameters, scenario, model.parameters.clinical );
    model.reporting.infantMortality = {vector<int>(sim::stepsPerYear()),
        vector<int>(sim::stepsPerYear()), parameters[Parameter::NON_MALARIA_INFANT_MORTALITY]};
    model.population.neonatalMortality = Host::makeNeonatalMortality(sim::oneTS());
    model.parameters.neonatalDiagnostic = xml::readNeonatalDiagnostic(scenario.getModel().getClinical());
    model.parameters.cumulativeAgeProportions = makeAgeStructure(
        xml::readDemography(scenario.getDemography()), sim::maxHumanAge(), sim::oneTS());

    // Species parameters and vector RNG precede creation of the human population.
    size_t popSize = scenario.getDemography().getPopSize();

    auto& population = model.population;
    model.transmission = xml::createTransmission(scenario.getEntomology(), popSize, model.parameters.hostTransmission, model.masterRng);
    auto& transmission = *model.transmission;

    transmission.registerContinuousCallbacks(model.reporting.continuous);
    registerContinousPopulationCallbacks(model.parameters.hostTransmission, model.reporting.continuous);

    // Depends on transmission model (for species indexes):
    // MDA1D may depend on health system (too complex to verify)
    model.interventions.init(scenario.getInterventions(), population, transmission, model.parameters.clinical );
    xml::setHealthSystem( scenario.getHealthSystem(), model.interventions, model.parameters.clinical ); // Depends on interventions, PK/PD (from humanPop)
    mon::initCohorts( scenario.getMonitoring(), model.interventions ); // Depends on interventions

    model.surveyOnlyNewEp = scenario.getMonitoring().getSurveyOptions().getOnlyNewEpisode();

}
}
