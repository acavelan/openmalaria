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

#include "Model.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "util/vectors.h"
#include "util/StreamValidator.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ClinicalModel.h"
#include "Transmission/TransmissionModel.h"
#include "interventions/Interventions.h"
#include "mon/Continuous.h"

namespace OM {
Model::Model(size_t populationSize) : population{populationSize} {}
Model::~Model() = default;

static void summarizeHuman(Host::Human &human, bool surveyOnlyNewEp, const interventions::CohortRemovalRules& removalRules) {
    if( surveyOnlyNewEp && human.clinicalModel->isExistingCase() ){
        // This modifies the denominator to treat the health-system-memory
        // period immediately after a bout as 'not at risk'.
        return;
    }
    
    mon::recordStat(mon::measure("nHost"), human);
    mon::recordStat(mon::measure("sumAge"), human, sim::inYears(human.age(sim::now())) );
    bool patent = human.withinHostModel->summarize (human);
    mon::recordStat(mon::measure("nExpectd"), human, human.infIncidence.probabilityInfected);
    
    if( patent && mon::isReported() ){
        // this should happen after all other reporting!
        human.removeFirstEvent(removalRules[interventions::SubPopRemove::ON_FIRST_INFECTION], sim::nowOrTs0());
    }
}

static void updateHuman(Host::Human &human, Transmission::TransmissionModel& transmission,
        const Host::IncidenceParameters& incidence, Reporting& reporting, double neonatalRisk, const interventions::CohortRemovalRules& removalRules)
{
    // For integer age checks we use age0 to e.g. get 73 steps comparing less than 1 year old
    SimTime age0 = human.age(sim::ts0());
    if (human.clinicalModel->isDead(age0)) {
        human.dead = true;
        return;
    }
    
    util::streamValidate( age0 );

    // monitoringAgeGroup is the group for the start of the time step.
    mon::updateAgeGroup(human.monitoringAgeGroup, age0);
    
    human.updateCohortSet();

    // Age at  the end of the update period. In most cases
    // the difference between this and age at the start is not especially
    // important in the model design, but since we parameterised with
    // age1 we should stick with it.
    double age1 = sim::inYears(human.age(sim::ts1()));

    // age1 used only in relativeAvailabilityAge(); difference to age0 should be minor
    vector<double> EIR_per_genotype_i, EIR_per_genotype_l;
    if (transmission.isVector)
        interventions::updateVectorInterventions(human.interventions, human.rng, sim::ts0());
    transmission.getEIR(human, age0, age1, EIR_per_genotype_i, EIR_per_genotype_l);

    double EIR_i = util::vectors::sum(EIR_per_genotype_i);
    double EIR_l = util::vectors::sum(EIR_per_genotype_l);
    double EIR = EIR_i + EIR_l;

    double x = 0.0;
    
    if(EIR > 0)
        x = EIR_i / EIR;

    double expectedNNewInfs = Host::expectedInfections(human.infIncidence, incidence, human.rng, EIR,
        human.transmission.relativeAvailability,
        incidence.vaccineGenotype ? 1.0 : interventions::vaccineFactor(human.interventions.vaccines, interventions::Vaccine::PEV, sim::ts1()), sim::oneTS());
    double expectedNNewInfs_i = x * expectedNNewInfs;
    double expectedNNewInfs_l = (1.0-x) * expectedNNewInfs;

    // One more rng call if imported infection
    int nNewInfs_l = Host::sampleNewInfections( human.rng, expectedNNewInfs_l);
    int nNewInfs_i = Host::sampleNewInfections( human.rng, expectedNNewInfs_i);

    // age1 used when medicating drugs (small effect) and in immunity model (which was parameterised for it)
    human.withinHostModel->update(human, human.rng, nNewInfs_i, nNewInfs_l, EIR_per_genotype_i, EIR_per_genotype_l, age1);
    mon::recordEvent(mon::measure("nNewInfections"), human, nNewInfs_i+nNewInfs_l);
    mon::recordEvent(mon::measure("nNewInfections_Introduced"), human, nNewInfs_i);
    mon::recordEvent(mon::measure("nNewInfections_Indigenous"), human, nNewInfs_l);
    reporting.newInfections += nNewInfs_i+nNewInfs_l;

    // age1 used to get case fatality and sequelae probabilities, determine pathogenesis
    human.clinicalModel->update( human, age1, age0 == sim::zero(), neonatalRisk, removalRules );
    human.clinicalModel->updateInfantDeaths(age0, reporting.infantMortality);
}

void update(Model& model) {
    auto& population = model.population;
    auto& transmission = *model.transmission;
    // Monitoring. sim::now() gives time of end of last step,
    // and is when reporting happens in our time-series.
    mon::Continuous::update(model.reporting.continuous, population);
    if( sim::intervDate() == mon::nextSurveyDate() ){
        for(Host::Human &human : population.humans)
            summarizeHuman(human, model.surveyOnlyNewEp, model.interventions.removalRules);
        transmission.summarize();
        mon::concludeSurvey();
    }
    
    // Deploy interventions, at time sim::now().
    model.interventions.deploy( population.humans, transmission );
    
    // Time step updates. Time steps are mid-day to mid-day.
    // sim::ts0() gives the date at the start of the step, sim::ts1() the date at the end.
    sim::start_update();

    // This should be called before humans contract new infections in the simulation step.
    // This needs the whole population (it is an approximation before all humans are updated).
    transmission.vectorUpdate(population.humans);
    
    // NOTE: no neonatal mortalities will occur in the first 20 years of warmup
    // (until humans old enough to be pregnate get updated and can be infected).
    Host::update(population.neonatalMortality, population.humans,
        model.parameters.neonatalDiagnostic.value(), sim::ts0(), sim::oneTS());
    
    for (Host::Human& human : population.humans)
    {
        if (human.birthDate + sim::maxHumanAge() >= model.humanWarmupLength) // this is last time of possible update
            updateHuman(human, transmission, model.parameters.incidence, model.reporting, population.neonatalMortality.riskFromMaternalInfection, model.interventions.removalRules);
    }
   
    update(population, model.parameters, model.masterRng);
    
    // Doesn't matter whether non-updated humans are included (value isn't used
    // before all humans are updated).
    transmission.updateKappa(population.humans);
    transmission.surveyEIR();

    sim::end_update();

}
}
