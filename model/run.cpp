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
#include "util/StreamValidator.h"
#include "xml/DocumentLoader.h"

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

namespace OpenMalaria {}

using namespace OpenMalaria;

namespace OpenMalaria {
namespace {

void print_progress(const RunContext& ctx, int& lastPercent, SimTime& estEndTime)
// parity: was a by-value param + static lastPercent (latent bug: static never
// resets in-process); now the caller owns the counter — single-run output unchanged
{
    int percent = (ctx.clock.now() * 100) / estEndTime;
    if (percent != lastPercent) { // avoid huge amounts of output for performance/log-file size reasons
        lastPercent = percent;
        cerr << "\r" << percent << "%\t" << flush;
    }
}

void print_errno(const RunContext& ctx) {
    if (errno != 0) {
        char err[256];
        sprintf(err, "t = %d Please report! Error: ", int(ctx.clock.now()));
        std::perror(err);
        errno = 0;
    }
}

// --save-to FILE: write the checkpoint and exit successfully (main maps
// cmd_exception code 0 to a clean exit; destructors close the output files).
void save_and_exit(RunContext& ctx, Population& population, TransmissionModel& transmission, Phase phase,
                   SimTime& endTime, SimTime& estEndTime) {
    writeCheckpoint(ctx, population, transmission, phase, endTime, estEndTime, CommandLine::getSaveName());
    throw cmd_exception("checkpoint written to " + CommandLine::getSaveName(), Error::None);
}

// --save-at DAY (not "main", handled at the main-phase reset): fire once at
// the first loop boundary at or after the given simulation day.
void maybe_save(RunContext& ctx, Population& population, TransmissionModel& transmission, bool mainPhase,
                SimTime& endTime, SimTime& estEndTime) {
    const string spec = CommandLine::getSaveName().empty() ? string() : CommandLine::getSaveAt();
    if (spec.empty() || spec == "main") return;
    if (ctx.clock.now() >= fromDays(std::stoi(spec)))
        save_and_exit(ctx, population, transmission, mainPhase ? Phase::Main : Phase::Warmup, endTime, estEndTime);
}

// Internal simulation loop
void run(RunContext& ctx, Population& population, TransmissionModel& transmission, SimTime humanWarmupLength,
         SimTime& endTime, SimTime& estEndTime, bool surveyOnlyNewEp, string phase) {
    int lastPercent = -1;
    const bool mainPhase = phase == "Intervention period";

    if (CommandLine::option(CommandLine::VERBOSE)) cout << "Starting " << phase << "..." << endl;

    HostScratch scratch; // EIR out-vectors, reused across humans and ticks
    while (ctx.clock.now() < endTime) {
        maybe_save(ctx, population, transmission, mainPhase, endTime, estEndTime);

        if (CommandLine::option(CommandLine::VERBOSE) && ctx.clock.interv_date() > 0)
            cout << "Time step: " << ctx.clock.now() / ctx.clock.one_TS() << ", internal days: " << ctx.clock.now()
                 << " | " << estEndTime << ", Intervention Date: " << ctx.clock.interv_date() << endl;

        // Monitoring. ctx.clock.now() gives time of end of last step,
        // and is when reporting happens in our time-series.
        Continuous::update(ctx, population);
        if (ctx.clock.interv_date() == nextSurveyDate(ctx)) {
            for (Human& human : population.humans)
                summarize(ctx, human, surveyOnlyNewEp);
            transmission.summarize(ctx);
            concludeSurvey(ctx);
        }

        // Deploy interventions, at time ctx.clock.now().
        InterventionManager::deploy(ctx, population.humans, transmission);

        // Time step updates. Time steps are mid-day to mid-day.
        // ctx.clock.ts0() gives the date at the start of the step, ctx.clock.ts1() the date at the end.
        ctx.clock.start_update();

        // This should be called before humans contract new infections in the simulation step.
        // This needs the whole population (it is an approximation before all humans are updated).
        transmission.vectorUpdate(population.humans);

        // NOTE: no neonatal mortalities will occur in the first 20 years of warmup
        // (until humans old enough to be pregnate get updated and can be infected).
        NeonatalMortality::update(ctx, population.humans);

        for (Human& human : population.humans) {
            if (human.cols->dob[human.row] + ctx.clock.max_human_age >=
                humanWarmupLength) // this is last time of possible update
                update(ctx, human, transmission, scratch);
        }

        population.update();

        // Doesn't matter whether non-updated humans are included (value isn't used
        // before all humans are updated).
        transmission.updateKappa(population.humans);
        transmission.surveyEIR();

        ctx.clock.end_update();

        if (CommandLine::option(CommandLine::PROGRESS)) print_progress(ctx, lastPercent, estEndTime);
        print_errno(ctx);
    }

    if (CommandLine::option(CommandLine::VERBOSE)) cout << "Finishing " << phase << "..." << endl;
}

} // namespace

/// RunContext ctor/dtor: create/destroy the registry pointer members.
/// // parity: order — identical to the old ensure_extensions()/RunExtensions
/// member-init order (declaration order of RunContext.h); registries and
/// unittest backend installers depend on it.
RunContext::RunContext()
    : decision_trees(std::make_unique<std::vector<std::unique_ptr<FlatTree>>>()),
      whTreatments(std::make_unique<std::vector<Treatments>>()),
      humanComponents(std::make_unique<std::vector<ComponentConfig>>()),
      continuous(std::make_unique<std::vector<ContinuousHumanDeployment>>()),
      timed(std::make_unique<std::vector<TimedDeployment>>()),
      importedInfections(std::make_unique<ImportedInfections>()),
      removeAtIds(std::make_unique<std::vector<ComponentId>[]>(SubPopRemove::NUM)),
      anophParams(std::make_unique<std::vector<PerHostAnophParams>>()), monitoring(std::make_unique<RuntimeState>()) {}

RunContext::~RunContext() = default;

// The one accessor for a TreatmentId (plain struct now; needs whTreatments).
const Treatments& Treatments::select(RunContext& ctx, TreatmentId treatId) {
    assert(treatId.id < ctx.whTreatments->size());
    return (*ctx.whTreatments)[treatId.id];
}

/// run_simulation() — loads scenario XML and runs the simulation (was main())
void run_simulation(RunContext& ctx, const std::string& scenarioFile) {
    set_gsl_handler();
    const bool resuming = !CommandLine::getResumeName().empty();
    auto [population, transmission, surveyOnlyNewEp, humanWarmupLength, estEndTime] =
        initialiseScenario(ctx, scenarioFile, resuming);
    SimTime endTime;
    Phase phase = Phase::Warmup;
    if (resuming) {
        // Everything config-derived was just rebuilt from the scenario; only
        // run state comes off the wire (clock, statics, monitoring,
        // transmission, humans, all RNG streams).
        readCheckpoint(ctx, *population, *transmission, phase, endTime, estEndTime, CommandLine::getResumeName());
        PerHostAnophParams::calcAvailabilityPercentiles(ctx);
        if (CommandLine::option(CommandLine::VERBOSE))
            cout << "Resumed checkpoint at day " << ctx.clock.now() / ctx.clock.one_TS() << " ("
                 << (phase == Phase::Main ? "main" : "warmup") << " phase)..." << endl;
    } else {
        population->createInitialHumans();
        transmission->init2(population->humans);

        /** Calculate ento availability percentiles **/
        PerHostAnophParams::calcAvailabilityPercentiles(ctx);

        endTime = humanWarmupLength;
    }

    if (phase == Phase::Warmup) { // normal start, or a checkpoint saved during warmup/calibration
        /** Warm-up phase:
         * Run the simulation using the equilibrium inoculation rates over one
         * complete lifespan (ctx.clock.max_human_age) to reach immunological
         * equilibrium in all age classes. Don't report any events. */
        run(ctx, *population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp, "Warmup");

        /** Transmission init phase:
         * Fit the emergence rate to the input EIR */
        SimTime iterate = transmission->initIterate();
        while (iterate > zero()) {
            endTime = endTime + iterate;
            // adjust estimation of final time step: end of current period + length of main phase
            estEndTime = endTime + (ctx.clock.end_date - ctx.clock.start_date) + ctx.clock.one_TS();
            run(ctx, *population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp,
                "EIR Calibration");
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
        ctx.clock.interv = zero();
        InfectionIncidenceModel::preMainSimInit();
        InfantMortality::preMainSimInit();
        ctx.genotypes.preMainSimInit();
        population->recentBirths = 0;
        transmission->summarize(ctx); // Only to reset TransmissionModel::inoculationsPerAgeGroup
        initMainSim(ctx);

        // --save-at main: checkpoint exactly at the start of the main phase.
        if (!CommandLine::getSaveName().empty() && CommandLine::getSaveAt() == "main")
            save_and_exit(ctx, *population, *transmission, Phase::Main, endTime, estEndTime);
    }

    // Main phase loop
    run(ctx, *population, *transmission, humanWarmupLength, endTime, estEndTime, surveyOnlyNewEp,
        "Intervention period");

    cerr << '\r' << flush; // clean last line of progress-output

    for (Human& human : population->humans)
        human.clinicalModel->flushReports();

    writeSurveyData(ctx);

#ifdef OM_STREAM_VALIDATOR
    StreamValidator.saveStream();
#endif

    // simulation's destructor runs
}

}
