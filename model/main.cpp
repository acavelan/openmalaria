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

#include "Global.h"
#include "Model.h"
#include "Transmission/TransmissionModel.h"
#include "ModelInitialization.h"
#include "ScenarioXml.h"

#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/ModelNameProvider.h"
#include "util/StreamValidator.h"
#include "util/DocumentLoader.h"
#include "util/XMLChecker.h"

#include "mon/Continuous.h"
#include "mon/init.h"
#include "mon/Monitoring.h"

#include "interventions/Interventions.h"
#include "Clinical/ClinicalModel.h"

#include "Host/NeonatalMortality.h"
#include "checkpoint.h"

#include "schema/scenario.h"

#include <cerrno>

namespace OM {
    using interventions::Interventions;
    using Transmission::TransmissionModel;
}

using namespace OM;

void print_progress(int lastPercent, SimTime &estEndTime)
{
    int percent = (sim::now() * 100) / estEndTime;
    if( percent != lastPercent ) {   // avoid huge amounts of output for performance/log-file size reasons
        lastPercent = percent;
        cerr << "\r" << percent << "%\t" << flush;
    }
}

void print_errno()
{
    if( errno != 0 )
    {
       char err[256];
       sprintf(err, "t = %d Please report! Error: ", int(sim::now()));
       std::perror(err);
       errno = 0;
    }
}

// Internal simulation loop
void run(Model& model, SimTime &endTime, SimTime &estEndTime, const string& phase)
{
    int lastPercent = -1;

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Starting " << phase << "..." << endl;

    while (sim::now() < endTime)
    {
        if (util::CommandLine::option(util::CommandLine::VERBOSE) && sim::intervDate() > 0)
            cout << "Time step: " << sim::now() / sim::oneTS() << ", internal days: " << sim::now() << " | " << estEndTime << ", Intervention Date: " << sim::intervDate() << endl;

        update(model);

        if (util::CommandLine::option(util::CommandLine::PROGRESS))
            print_progress(lastPercent, estEndTime);
        print_errno();
    }

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Finishing " << phase << "..." << endl;
}

/// main() — loads scenario XML and runs simulation
int main(int argc, char* argv[])
{
    int exitStatus = EXIT_SUCCESS;
    bool startedFromCheckpoint;

    string scenarioFile;
    string checkpointFileName;
    SimTime estEndTime, endTime;
    
    try {
        util::set_gsl_handler();
        
        scenarioFile = util::CommandLine::parse (argc, argv);
        unique_ptr<scnXml::Scenario> scenario = util::loadScenario(scenarioFile);

        util::XMLChecker().PerformPostValidationChecks(*scenario);

        Model model(scenario->getDemography().getPopSize());
        initialize(model, *scenario);
        auto* population = &model.population;
        auto* transmission = model.transmission.get();

        sim::s_t0 = sim::zero();
        sim::s_t1 = sim::zero();
        
        // Make sure warmup period is at least as long as a human lifespan, as the
        // length required by vector warmup, and is a whole number of years.
        SimTime& humanWarmupLength = model.humanWarmupLength;
        humanWarmupLength = sim::maxHumanAge();
        if(transmission->interventionMode != Transmission::forcedEIR)
            humanWarmupLength = max(humanWarmupLength, sim::fromYearsI(55)); // Data is summed over 5 years; add an extra 50 for stabilization.

        humanWarmupLength = sim::fromYearsI( static_cast<int>(ceil(sim::inYears(humanWarmupLength))) );

        // ———  End of static data initialisation  ———
        checkpointFileName = util::CommandLine::getCheckpointName();

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
        
        estEndTime = humanWarmupLength + (sim::endDate() - sim::startDate()) + sim::oneTS();
        assert( estEndTime + sim::never() < sim::zero() );

        if (startedFromCheckpoint)
        {
            xml::initContinuous(model.reporting.continuous, scenario->getMonitoring(), util::CommandLine::getCtsoutName(), true);
            readCheckpoint(checkpointFileName, endTime, estEndTime, model);

            /** Calculate ento availability percentiles **/
            calculateAvailabilityPercentiles(model.parameters.hostTransmission);
        }
        else
        {
            xml::initContinuous(model.reporting.continuous, scenario->getMonitoring(), util::CommandLine::getCtsoutName(), false);
            createInitialHumans(*population, model.parameters, model.masterRng);
            transmission->init2(population->humans);
            
            /** Calculate ento availability percentiles **/
            calculateAvailabilityPercentiles(model.parameters.hostTransmission);

            /** Warm-up phase: 
             * Run the simulation using the equilibrium inoculation rates over one
             * complete lifespan (sim::maxHumanAge()) to reach immunological
             * equilibrium in all age classes. Don't report any events. */
            endTime = humanWarmupLength;
            run(model, endTime, estEndTime, "Warmup");

            /** Transmission init phase:
             * Fit the emergence rate to the input EIR */
            SimTime iterate = transmission->initIterate();
            while(iterate > sim::zero())
            {
                endTime = endTime + iterate;
                // adjust estimation of final time step: end of current period + length of main phase
                estEndTime = endTime + (sim::endDate() - sim::startDate()) + sim::oneTS();
                run(model, endTime, estEndTime, "EIR Calibration");
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
            estEndTime = endTime = endTime + (sim::endDate() - sim::startDate()) + sim::oneTS();
            sim::s_interv = sim::zero();
            model.reporting.newInfections = 0;
            model.reporting.infantMortality.reset();
            WithinHost::Genotypes::preMainSimInit();
            population->recentBirths = 0;
            transmission->summarize(); // Only to reset TransmissionModel::inoculationsPerAgeGroup
            mon::initMainSim();

            if(util::CommandLine::option (util::CommandLine::CHECKPOINT))
            {
                writeCheckpoint(startedFromCheckpoint, checkpointFileName, endTime, estEndTime, model);
                if( util::CommandLine::option (util::CommandLine::CHECKPOINT_STOP) )
                    throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
            }
        }

        // Main phase loop
        run(model, endTime, estEndTime, "Intervention period");
       
        cerr << '\r' << flush;  // clean last line of progress-output
        
        for(Host::Human &human : population->humans)
            human.clinicalModel->flushReports();

        mon::writeSurveyData(model.reporting.infantMortality);
        
    # ifdef OM_STREAM_VALIDATOR
        util::StreamValidator.saveStream();
    # endif
        
        // simulation's destructor runs
    } catch (const OM::util::cmd_exception& e) {
        if( e.getCode() == 0 ){
            // this is not an error, but exiting due to command line
            std::cerr << e.what() << "; exiting..." << endl;
        }else{
            std::cerr << "Command-line error: "<<e.what();
            exitStatus = e.getCode();
        }
    } catch (const ::xsd::cxx::tree::exception<char>& e) {
        if (errno == 2) // Errno value 2 corresponds to "No such file or directory" i.e. this is not a mismatch between scenario contents and schema rules.
        {
            // We don't print the content of the xsdcxx exception since the error messages are not relevant to the user in this case.
            std::cerr << "Parsing scenario file failed.\n\tLikely the XSD schema file is not present at an expected location/with expected filename." << std::endl;
        }
        else
        {
            // e.what() returns a concise error message e.g. "instance document parsing failed".
            // Whereas inserting e into a stream inserts the (possibly long and repetitive) list of error messages arising from the parsing attempt.
            std::cerr << "XSD error: " << e.what() << '\n' << e << endl;

            // We print this *after* the above since users are more likely to take note of the tail of stderr/stdout output than the head.
            std::cerr << "Parsing scenario file failed.\n\tLikely the XSD schema is named/located correctly but the scenario does not conform to the schema." << std::endl;
        }
        exitStatus = OM::util::Error::XSD;
    } catch (const OM::util::checkpoint_error& e) {
        std::cerr << "Checkpoint error: " << e.what() << endl;
        std::cerr << e << flush;
        exitStatus = e.getCode();
    } catch (const OM::util::traced_exception& e) {
        std::cerr << "Code error: " << e.what() << endl;
        std::cerr << e << flush;
        std::cerr << "This is likely an error in the C++ code. Please report!" << endl;
        exitStatus = e.getCode();
    } catch (const OM::util::xml_scenario_error& e) {
        std::cerr << "Error: " << e.what() << endl;
        std::cerr << "In: " << scenarioFile << endl;
        exitStatus = e.getCode();
    } catch (const OM::util::base_exception& e) {
        std::cerr << "Error: " << e.message() << endl;
        exitStatus = e.getCode();
    } catch (const exception& e) {
        std::cerr << "Error: " << e.what() << endl;
        exitStatus = EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << endl;
        exitStatus = EXIT_FAILURE;
    }
    
    // If we get to here, we already know an error occurred.
    if( errno != 0 )
        std::perror( "OpenMalaria" );
    
    return exitStatus;
}
