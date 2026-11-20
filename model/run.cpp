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

#include "Clock.h"

#include "RunContext.h"

#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/ModelNameProvider.h"
#include "util/StreamValidator.h"
#include "xml/DocumentLoader.h"
#include "xml/XMLChecker.h"

#include "checkpoint/checkpoint.h"
#include "mon/Continuous.h"
#include "mon/Monitoring.h"

#include "interventions/InterventionManager.h"
#include "Clinical/ClinicalModel.h"

#include "Host/NeonatalMortality.h"
#include "Transmission/transmission.h"
// RunContext's pointer members + Treatments::select need these complete
// (they cannot be pulled in through RunContext.h — include cycle):
#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/ImportedInfections.h"
#include "interventions/Components.h"
#include "interventions/Deployments.h"
#include "Transmission/PerHost.h"
#include <cassert>

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

// --save-to FILE: write the checkpoint and exit successfully (main maps
// cmd_exception code 0 to a clean exit; destructors close the output files).
void save_and_exit(RunContext& ctx, Population &population, TransmissionModel &transmission,
        checkpoint::Phase phase, SimTime &endTime, SimTime &estEndTime)
{
    checkpoint::writeCheckpoint (ctx, population, transmission, phase,
        endTime, estEndTime, util::CommandLine::getSaveName());
    throw util::cmd_exception ("checkpoint written to " +
        util::CommandLine::getSaveName(), util::Error::None);
}

// --save-at DAY (not "main", handled at the main-phase reset): fire once at
// the first loop boundary at or after the given simulation day.
void maybe_save(RunContext& ctx, Population &population, TransmissionModel &transmission,
        bool mainPhase, SimTime &endTime, SimTime &estEndTime)
{
    const string spec = util::CommandLine::getSaveName().empty() ?
        string() : util::CommandLine::getSaveAt();
    if (spec.empty() || spec == "main") return;
    if (ctx.clock.now() >= sim::fromDays (std::stoi (spec)))
        save_and_exit (ctx, population, transmission,
            mainPhase ? checkpoint::Phase::Main : checkpoint::Phase::Warmup,
            endTime, estEndTime);
}

// Internal simulation loop
void run(RunContext& ctx, Population &population, TransmissionModel &transmission, SimTime humanWarmupLength, SimTime &endTime, SimTime &estEndTime, bool surveyOnlyNewEp, string phase)
{
    int lastPercent = -1;
    const bool mainPhase = phase == "Intervention period";

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Starting " << phase << "..." << endl;

    while (ctx.clock.now() < endTime)
    {
        maybe_save (ctx, population, transmission, mainPhase, endTime, estEndTime);

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
            if (human.cols->dob[human.row] + ctx.clock.max_human_age >= humanWarmupLength) // this is last time of possible update
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

/// RunContext ctor/dtor: create/destroy the registry pointer members.
/// // parity: order — identical to the old ensure_extensions()/RunExtensions
/// member-init order (declaration order of RunContext.h); registries and
/// unittest backend installers depend on it.
RunContext::RunContext() :
    decision_trees(new std::vector<std::unique_ptr<Clinical::FlatTree>>),
    whTreatments(new std::vector<WithinHost::Treatments>),
    humanComponents(new std::vector<interventions::ComponentConfig>),
    continuous(new std::vector<interventions::ContinuousHumanDeployment>),
    timed(new std::vector<interventions::TimedDeployment>),
    importedInfections(new OM::Host::ImportedInfections),
    removeAtIds(new std::vector<interventions::ComponentId>
        [interventions::SubPopRemove::NUM]),
    anophParams(new std::vector<Transmission::PerHostAnophParams>),
    monitoring(new mon::internal::RuntimeState)
{}

// (was: delete ext; destroys the old RunExtensions members in the same
// reverse-declaration order)
RunContext::~RunContext()
{
    delete monitoring;
    delete anophParams;
    delete importedInfections;
    delete timed;
    delete continuous;
    delete humanComponents;
    delete[] removeAtIds;
    delete whTreatments;
    delete decision_trees;
}

// The one accessor for a TreatmentId (plain struct now; needs whTreatments).
const WithinHost::Treatments& WithinHost::Treatments::select(
    RunContext& ctx, TreatmentId treatId )
{
    assert( treatId.id < ctx.whTreatments->size() );
    return (*ctx.whTreatments)[treatId.id];
}

/// run_simulation() — loads scenario XML and runs the simulation (was main())
void run_simulation(RunContext& ctx, const std::string& scenarioFile)
{
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
    interventions::InterventionManager::init(ctx, scenario->getInterventions(), *transmission );
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

    estEndTime = humanWarmupLength + (ctx.clock.end_date - ctx.clock.start_date) + ctx.clock.one_TS();
    assert( estEndTime + sim::never() < sim::zero() );

    const bool resuming = !util::CommandLine::getResumeName().empty();
    checkpoint::Phase phase = checkpoint::Phase::Warmup;

    mon::Continuous::init(ctx, scenario->getMonitoring(), resuming);
    if (resuming)
    {
        // Everything config-derived was just rebuilt from the scenario; only
        // run state comes off the wire (clock, statics, monitoring,
        // transmission, humans, all RNG streams).
        checkpoint::readCheckpoint (ctx, *population, *transmission, phase,
            endTime, estEndTime, util::CommandLine::getResumeName());
        Transmission::PerHostAnophParams::calcAvailabilityPercentiles(ctx);
        if (util::CommandLine::option(util::CommandLine::VERBOSE))
            cout << "Resumed checkpoint at day " << ctx.clock.now() / ctx.clock.one_TS()
                 << " (" << (phase == checkpoint::Phase::Main ? "main" : "warmup")
                 << " phase)..." << endl;
    }
    else
    {
        population->createInitialHumans();
        transmission->init2(population->humans);

        /** Calculate ento availability percentiles **/
        Transmission::PerHostAnophParams::calcAvailabilityPercentiles(ctx);

        endTime = humanWarmupLength;
    }

    if (phase == checkpoint::Phase::Warmup)
    {   // normal start, or a checkpoint saved during warmup/calibration
        /** Warm-up phase:
         * Run the simulation using the equilibrium inoculation rates over one
         * complete lifespan (ctx.clock.max_human_age) to reach immunological
         * equilibrium in all age classes. Don't report any events. */
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
         * (ii)        the entomological input defined in intEIR()
         * (iii)       the intervention packages defined in Intervention()
         * (iv)        the survey times defined in Survey() */
        // reset endTime and estEndTime to their exact value after initIterate()
        estEndTime = endTime = endTime + (ctx.clock.end_date - ctx.clock.start_date) + ctx.clock.one_TS();
        ctx.clock.interv = sim::zero();
        Host::InfectionIncidenceModel::preMainSimInit();
        Clinical::InfantMortality::preMainSimInit();
        ctx.genotypes.preMainSimInit();
        population->recentBirths = 0;
        transmission->summarize(ctx); // Only to reset TransmissionModel::inoculationsPerAgeGroup
        mon::initMainSim(ctx);

        // --save-at main: checkpoint exactly at the start of the main phase.
        if (!util::CommandLine::getSaveName().empty() &&
            util::CommandLine::getSaveAt() == "main")
            save_and_exit (ctx, *population, *transmission,
                checkpoint::Phase::Main, endTime, estEndTime);
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
