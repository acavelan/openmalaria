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

// Engine driver — body moved verbatim from main.cpp (phase 3, plan §5: one
// init/run sequence shared by main() and unittest/TwoRunsSuite.h).

#include "Global.h"

#include "core/run.h"

#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/ModelNameProvider.h"
#include "util/StreamValidator.h"
#include "xml/DocumentLoader.h"
#include "xml/XMLChecker.h"

#include "mon/Continuous.h"
#include "mon/init.h"
#include "mon/Monitoring.h"

#include "interventions/InterventionManager.h"
#include "Clinical/ClinicalModel.h"

#include "Host/NeonatalMortality.h"
#include "checkpoint.h"
#include "core/RunExtensions.h"

#include <cerrno>

namespace OM {
    using interventions::InterventionManager;
    using Transmission::TransmissionModel;
}

using namespace OM;

namespace OM {
namespace {

void print_progress(const RunContext& ctx, int &lastPercent, SimTime &estEndTime)
// parity: was a by-value param + static lastPercent (latent bug: static never
// resets in-process); now the caller owns the counter — single-run output unchanged
{
    int percent = (ctx.clock.now() * 100) / estEndTime;
    if( percent != lastPercent ) {   // avoid huge amounts of output for performance/log-file size reasons
        lastPercent = percent;
        cerr << "\r" << percent << "%\t" << flush;
    }
}

void print_errno(const RunContext& ctx)
{
    if( errno != 0 )
    {
       char err[256];
       sprintf(err, "t = %d Please report! Error: ", int(ctx.clock.now()));
       std::perror(err);
       errno = 0;
    }
}

// Internal simulation loop
void run(RunContext& ctx, Population &population, TransmissionModel &transmission, SimTime humanWarmupLength, SimTime &endTime, SimTime &estEndTime, bool surveyOnlyNewEp, string phase)
{
    int lastPercent = -1;

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Starting " << phase << "..." << endl;

    while (ctx.clock.now() < endTime)
    {
        if (util::CommandLine::option(util::CommandLine::VERBOSE) && ctx.clock.interv_date() > 0)
            cout << "Time step: " << ctx.clock.now() / ctx.clock.one_TS() << ", internal days: " << ctx.clock.now() << " | " << estEndTime << ", Intervention Date: " << ctx.clock.interv_date() << endl;

        // Monitoring. ctx.clock.now() gives time of end of last step,
        // and is when reporting happens in our time-series.
        mon::Continuous::update( ctx, population );
        if( ctx.clock.interv_date() == mon::nextSurveyDate(ctx) ){
            for(Host::Human &human : population.humans)
                Host::summarize(ctx, human, surveyOnlyNewEp);
            transmission.summarize(ctx);
            mon::concludeSurvey(ctx);
        }
        
        // Deploy interventions, at time ctx.clock.now().
        InterventionManager::deploy( ctx, population.humans, transmission );
        
        // Time step updates. Time steps are mid-day to mid-day.
        // ctx.clock.ts0() gives the date at the start of the step, ctx.clock.ts1() the date at the end.
        ctx.clock.start_update();

        // This should be called before humans contract new infections in the simulation step.
        // This needs the whole population (it is an approximation before all humans are updated).
        transmission.vectorUpdate(population.humans);
        
        // NOTE: no neonatal mortalities will occur in the first 20 years of warmup
        // (until humans old enough to be pregnate get updated and can be infected).
        Host::NeonatalMortality::update (ctx, population.humans);
        
        Host::HostScratch scratch;   // EIR out-vectors, reused across humans/ticks
        for (Host::Human& human : population.humans)
        {
            if (human.getDOB() + ctx.clock.max_human_age >= humanWarmupLength) // this is last time of possible update
                Host::update(ctx, human, transmission, scratch);
        }
       
        population.update();
        
        // Doesn't matter whether non-updated humans are included (value isn't used
        // before all humans are updated).
        transmission.updateKappa(population.humans);
        transmission.surveyEIR();

        ctx.clock.end_update();

        if (util::CommandLine::option(util::CommandLine::PROGRESS))
            print_progress(ctx, lastPercent, estEndTime);
        print_errno(ctx);
    }

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Finishing " << phase << "..." << endl;
}

} // namespace

/// run_simulation() — loads scenario XML and runs the simulation (was main())
RunContext::~RunContext() { delete ext; }

void run_simulation(RunContext& ctx, const std::string& scenarioFile)
{
    ensure_extensions( ctx );

    bool startedFromCheckpoint;

    SimTime estEndTime, endTime;

    util::set_gsl_handler();

    unique_ptr<scnXml::Scenario> scenario = util::loadScenario(scenarioFile);

    util::XMLChecker().PerformPostValidationChecks(*scenario);

    // 1) elements with no dependencies on other elements initialised here:
    WithinHost::Genotypes::init( ctx.genotypes, *scenario );
    {
        const int seedOverride = util::CommandLine::getSeedOverride();
        const unsigned iseed = seedOverride >= 0 ? unsigned (seedOverride) :
            scenario->getModel().getComputationParameters().getIseed();
        // parity: f — validation completes BEFORE ctx.master_RNG.seed(): load->check->seed (plan §3f); never seed before checks
        ctx.master_RNG.seed( iseed, 0 ); // Init RNG with Iseed (--seed overrides)
    }
    util::ModelNameProvider modelNameProvider(scenario->getModel());

    // 2) elements depending on only elements initialised in (1):
    ctx.clock.init(ctx, *scenario, modelNameProvider); // Also reads survey dates.

    // 3) elements depending on only elements initialised in (2).
    Parameters parameters( scenario->getModel().getParameters(), modelNameProvider ); // Depends on ModelNameProvider.
    util::ModelOptions::init( ctx.clock, scenario->getModel().getModelOptions(), modelNameProvider ); // Depends on ModelNameProvider.

    // 4) elements depending on only elements initialised in (3).
    WithinHost::diagnostics::init( ctx, parameters, *scenario ); // Depends on Parameters.

    // 5) elements depending on only elements initialised in (4).
    mon::initReporting( ctx, *scenario ); // Reporting init depends on diagnostics and monitoring

    // Init models used by humans
    Transmission::PerHost::init( ctx, scenario->getModel().getHuman().getAvailabilityToMosquitoes() );
    Host::InfectionIncidenceModel::init( ctx, parameters );
    WithinHost::WHInterface::init( ctx, parameters, *scenario );
    Clinical::ClinicalModel::init( ctx, parameters, *scenario );
    Host::NeonatalMortality::init( ctx, scenario->getModel().getClinical() );
    AgeStructure::init( ctx, scenario->getDemography() );

    // 3) elements depending on other elements; dependencies on (1) are not mentioned:
    // Transmission model initialisation depends on Transmission::PerHost and
    // genotypes (both from Human, from Population::init()) and
    // monitoring age groups (from Surveys.init()):
    // Note: PerHost dependency can be postponed; it is only used to set adultAge
    size_t popSize = scenario->getDemography().getPopSize();

    std::unique_ptr<Population> population = std::unique_ptr<Population>(new Population(ctx, popSize));
    std::unique_ptr<TransmissionModel> transmission = std::unique_ptr<TransmissionModel>(Transmission::createTransmissionModel(ctx, scenario->getEntomology(), popSize));

    registerContinousPopulationCallbacks(ctx);

    // Depends on transmission model (for species indexes):
    // MDA1D may depend on health system (too complex to verify)
    interventions::InterventionManager::init(ctx, scenario->getInterventions(), *population, *transmission );
    Clinical::ClinicalModel::setHS( ctx, scenario->getHealthSystem() ); // Depends on interventions, PK/PD (from humanPop)
    mon::initCohorts(ctx,  scenario->getMonitoring() ); // Depends on interventions

    bool surveyOnlyNewEp = scenario->getMonitoring().getSurveyOptions().getOnlyNewEpisode();

    ctx.clock.t0 = sim::zero();
    ctx.clock.t1 = sim::zero();
    
    // Make sure warmup period is at least as long as a human lifespan, as the
    // length required by vector warmup, and is a whole number of years.
    SimTime humanWarmupLength = ctx.clock.max_human_age;
    if(transmission->interventionMode != Transmission::forcedEIR)
        humanWarmupLength = max(humanWarmupLength, sim::fromYearsI(55)); // Data is summed over 5 years; add an extra 50 for stabilization.

    humanWarmupLength = sim::fromYearsI( static_cast<int>(ceil(sim::inYears(humanWarmupLength))) );

    // ———  End of static data initialisation  ———
    string checkpointFileName = util::CommandLine::getCheckpointName();

    if(checkpointFileName == "")
        checkpointFileName = "checkpoint";

    if(util::CommandLine::option(util::CommandLine::CHECKPOINT))
    {
        ifstream checkpointFile(checkpointFileName,ios::in);
        startedFromCheckpoint = checkpointFile.is_open();
        if(startedFromCheckpoint == false)
            errno = 0; // Cleanup errno if file doesn't exist
    }
    else
        startedFromCheckpoint = false;
    
    estEndTime = humanWarmupLength + (ctx.clock.end_date - ctx.clock.start_date) + ctx.clock.one_TS();
    assert( estEndTime + sim::never() < sim::zero() );

    if (startedFromCheckpoint)
    {
        mon::Continuous::init(ctx, scenario->getMonitoring(), true);
        readCheckpoint(ctx, checkpointFileName, endTime, estEndTime, *population, *transmission);

        /** Calculate ento availability percentiles **/
        Transmission::PerHostAnophParams::calcAvailabilityPercentiles(ctx);
    }
    else
    {
        mon::Continuous::init(ctx, scenario->getMonitoring(), false);
        population->createInitialHumans();
        transmission->init2(population->humans);
        
        /** Calculate ento availability percentiles **/
        Transmission::PerHostAnophParams::calcAvailabilityPercentiles(ctx);

        /** Warm-up phase: 
         * Run the simulation using the equilibrium inoculation rates over one
         * complete lifespan (ctx.clock.max_human_age) to reach immunological
         * equilibrium in all age classes. Don't report any events. */
        endTime = humanWarmupLength;
        run(ctx, *population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp, "Warmup");

        /** Transmission init phase:
         * Fit the emergence rate to the input EIR */
        SimTime iterate = transmission->initIterate();
        while(iterate > sim::zero())
        {
            endTime = endTime + iterate;
            // adjust estimation of final time step: end of current period + length of main phase
            estEndTime = endTime + (ctx.clock.end_date - ctx.clock.start_date) + ctx.clock.one_TS();
            run(ctx, *population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp, "EIR Calibration");
            iterate = transmission->initIterate();
        }

        /** Main phase:
         * This procedure starts with the current state of the simulation 
         * It continues updating assuming:
         * (i)         the default (exponential) demographic model
         * (ii)        the entomological input defined by the EIRs in intEIR()
         * (iii)       the intervention packages defined in Intervention()
         * (iv)        the survey times defined in Survey() */
        // reset endTime and estEndTime to their exact value after initIterate()
        estEndTime = endTime = endTime + (ctx.clock.end_date - ctx.clock.start_date) + ctx.clock.one_TS();
        ctx.clock.interv = sim::zero();
        Host::InfectionIncidenceModel::preMainSimInit();
        Clinical::InfantMortality::preMainSimInit();
        ctx.genotypes.preMainSimInit();
        population->resetRecentBirths();
        transmission->summarize(ctx); // Only to reset TransmissionModel::inoculationsPerAgeGroup
        mon::initMainSim(ctx);

        if(util::CommandLine::option (util::CommandLine::CHECKPOINT))
        {
            writeCheckpoint(ctx, startedFromCheckpoint, checkpointFileName, endTime, estEndTime, *population, *transmission);
            if( util::CommandLine::option (util::CommandLine::CHECKPOINT_STOP) )
                throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
        }
    }

    // Main phase loop
    run(ctx, *population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp, "Intervention period");
   
    cerr << '\r' << flush;  // clean last line of progress-output
    
    for(Host::Human &human : population->humans)
        human.clinicalModel->flushReports();

    mon::writeSurveyData(ctx);
    
#ifdef OM_STREAM_VALIDATOR
    util::StreamValidator.saveStream();
#endif
    
    // simulation's destructor runs
}

}
