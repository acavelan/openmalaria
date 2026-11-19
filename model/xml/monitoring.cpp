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

// model/xml/ conversion layer — monitoring module: scnXml::Monitoring /
// scnXml::Scenario survey+reporting readers. Moved verbatim from
// mon/init.cpp (initAgeGroups/initReporting/readSurveyDates/initCohorts) and
// mon/Continuous.cpp (Continuous::init). Declarations stay in the mon headers.

#include "Global.h"
#include "core/RunExtensions.h"
#include "Host/WithinHost/Diagnostic.h"
#include "interventions/InterventionManager.h"
#include "mon/Continuous.h"
#include "mon/Monitoring.h"
#include "mon/init.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include <schema/monitoring.h>
#include <schema/scenario.h>

#include <algorithm>
#include <bit>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

// Continuous.cpp file-scope state (external linkage; declared here so the
// moved init can reach it; definitions remain in mon/Continuous.cpp):
namespace OM { namespace mon { namespace Continuous {
    extern std::fstream ctsOStream;
    extern std::streamoff streamOffset;
    extern std::map<std::string, std::tuple<std::string, std::function<void(Population&, std::ostream&)>>> registered;
    extern std::vector< std::tuple<std::string, std::function<void(Population&, std::ostream&)>> > toReport;
    extern SimTime ctsPeriod;
    extern bool duringInit;
} } }

namespace OM {
namespace mon {

using internal::SurveyDate;



using internal::SurveyDate;

namespace {

constexpr uint32_t maxCohortNumber = uint32_t{1} << 21;

static void initAgeGroups(RunContext& ctx, const scnXml::Monitoring& monitoring)
{
    auto& runtime = ctx.ext->monitoring;
    const scnXml::MonAgeGroup::GroupSequence& groups = monitoring.getAgeGroup().getGroup();
    if (!(monitoring.getAgeGroup().getLowerbound() <= 0.0)) {
        throw util::xml_scenario_error("Expected survey age-group lowerbound of 0");
    }

    runtime.ageGroupUpperBound.resize(groups.size() + 1);
    for (size_t i = 0; i < groups.size(); ++i) {
        runtime.ageGroupUpperBound[i] = ctx.clock.from_years_d(groups[i].getUpperbound());
    }
    runtime.ageGroupUpperBound[groups.size()] = sim::future();
}

} // namespace

void initReporting(RunContext& ctx, const scnXml::Scenario& scenario)
{
    auto& runtime = ctx.ext->monitoring;
    runtime.reportIMR = -1;
    const scnXml::MonitoringOptions& optsElt = scenario.getMonitoring().getSurveyOptions();
    std::vector<OutMeasure> reportedMeasures;
    reportedMeasures.reserve(optsElt.getOption().size());

    std::set<int> outIds;
    for (const scnXml::MonitoringOption& optElt : optsElt.getOption()) {
        if (!optElt.getValue()) continue;

        OutMeasure om = findOutMeasure(optElt.getName());
        if (om.measure == invalidMeasure) {
            if (isObsoleteMeasure(optElt.getName())) {
                throw util::xml_scenario_error("obsolete survey option: " + std::string(optElt.getName()));
            }
            throw util::xml_scenario_error("unrecognised survey option: " + std::string(optElt.getName()));
        }
        if ((om.measure == measure("sumlogDens") || om.measure == measure("logDensByGenotype")) &&
            WithinHost::diagnostics::monitoringDiagnostic( ctx ).allowsFalsePositives())
        {
            throw util::xml_scenario_error("measure " + std::string(optElt.getName()) + " may not be used when monitoring diagnostic sensitivity < 1");
        }

        auto applyCategory = [&](Dim flag, const auto& option, const char* label) {
            if (!option.present()) return;
            if (option.get() && !(om.dims & flag)) {
                throw util::xml_scenario_error("measure " + std::string(optElt.getName()) + " does not support categorisation by " + label);
            }
            if (!option.get()) om.dims &= ~flag;
        };
        applyCategory(Dim::Age, optElt.getByAge(), "age group");
        applyCategory(Dim::Cohort, optElt.getByCohort(), "cohort");
        applyCategory(Dim::Species, optElt.getBySpecies(), "species");
        applyCategory(Dim::Genotype, optElt.getByGenotype(), "genotype");
        applyCategory(Dim::Drug, optElt.getByDrugType(), "drug type");

        if (optElt.getOutputNumber().present()) om.outId = optElt.getOutputNumber().get();
        if (om.outId < 0) {
            throw util::xml_scenario_error("monitoring output number must not be negative");
        }
        if (!outIds.insert(om.outId).second) {
            throw util::xml_scenario_error("monitoring output number " + std::to_string(om.outId) + " used more than once");
        }
        if (om.measure == allCauseIMR) {
            runtime.reportIMR = om.outId;
            continue;
        }
        reportedMeasures.push_back(om);
    }
    std::sort(reportedMeasures.begin(), reportedMeasures.end(),
        [](const OutMeasure& lhs, const OutMeasure& rhs) { return lhs.outId < rhs.outId; });
    const size_t nSpecies = scenario.getEntomology().getVector().present()
        ? scenario.getEntomology().getVector().get().getAnopheles().size() : 1;
    const size_t nDrugs = scenario.getPharmacology().present()
        ? scenario.getPharmacology().get().getDrugs().getDrug().size() : 1;
    internal::initStores(ctx, reportedMeasures, nSpecies, nDrugs);
}

SimTime readSurveyDates(RunContext& ctx, const scnXml::Monitoring& monitoring)
{
    const SimClock& clock = ctx.clock;
    auto& runtime = ctx.ext->monitoring;
    const scnXml::Surveys::SurveyTimeSequence& survs = monitoring.getSurveys().getSurveyTime();
    std::map<SimTime, bool> surveys;

    auto addSurvey = [&surveys](SimTime date, bool reporting) {
        if (reporting) surveys[date] = true;
        else surveys.insert(std::make_pair(date, false));
    };

    for (size_t i = 0; i < survs.size(); ++i) {
        const scnXml::SurveyTime& surv = survs[i];
        try {
            std::string s = surv;
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), s.end());
            SimTime cur = UnitParse::readDate(clock, s, UnitParse::STEPS);
            const bool reporting = surv.getReported();
            if (surv.getRepeatStep().present() != surv.getRepeatEnd().present()) {
                throw util::xml_scenario_error("surveyTime: use of repeatStep or repeatEnd without other");
            }
            if (surv.getRepeatStep().present()) {
                const SimTime step = UnitParse::readDuration(clock, surv.getRepeatStep().get(), UnitParse::NONE);
                if (step < clock.one_TS()) {
                    throw util::xml_scenario_error("surveyTime: repeatStep must be >= 1");
                }
                const SimTime end = UnitParse::readDate(clock, surv.getRepeatEnd().get(), UnitParse::NONE);
                while (cur < end) {
                    addSurvey(cur, reporting);
                    cur = cur + step;
                }
            } else {
                addSurvey(cur, reporting);
            }
        } catch (const util::format_error& e) {
            throw util::xml_scenario_error(std::string("surveyTime: ").append(e.message()));
        }
    }

    runtime.surveyDates.clear();
    runtime.surveyDates.reserve(surveys.size());
    runtime.nSurveys = 0;
    for (const auto& [date, reporting] : surveys) {
        const size_t num = reporting ? runtime.nSurveys++ : NOT_USED;
        runtime.surveyDates.push_back({date, num});
    }
    if (runtime.surveyDates.empty()) {
        throw util::xml_scenario_error("Scenario defines no surveys; at least one is required.");
    }
    if (!runtime.surveyDates.back().isReported()) {
        std::cerr << "Warning: the last survey is unreported. Having surveys beyond the last reported survey is pointless." << std::endl;
    }

    if (util::CommandLine::option(util::CommandLine::PRINT_SURVEY_TIMES)) {
        std::cout << "Survey\tsteps\tdate\n";
        for (const SurveyDate& surveyDate : runtime.surveyDates) {
            if (!surveyDate.isReported()) continue;
            std::cout
                << (surveyDate.num + 1) << '\t'
                << clock.inSteps(surveyDate.date - clock.start_date) << '\t'
                << surveyDate.date << std::endl;
        }
    }

    runtime.nCohorts = 1;
    if (monitoring.getCohorts().present()) {
        runtime.nCohorts = static_cast<uint32_t>(1) << monitoring.getCohorts().get().getSubPop().size();
    }
    initAgeGroups(ctx, monitoring);
    return runtime.surveyDates.back().date;
}

void initCohorts(RunContext& ctx, const scnXml::Monitoring& monitoring)
{
    auto& runtime = ctx.ext->monitoring;
    runtime.cohortSubPopNumbers.clear();
    runtime.cohortSubPopIds.clear();
    if (!monitoring.getCohorts().present()) return;

    const scnXml::Cohorts monCohorts = monitoring.getCohorts().get();
    uint32_t nextId = 0;
    for (auto it = monCohorts.getSubPop().begin(), end = monCohorts.getSubPop().end(); it != end; ++it) {
        const interventions::ComponentId compId = interventions::InterventionManager::getComponentId(ctx, it->getId());
        const bool inserted = runtime.cohortSubPopIds.insert(std::make_pair(compId, nextId)).second;
        if (!inserted) {
            throw util::xml_scenario_error(
                std::string("cohort specification uses sub-population \"").append(it->getId()).append("\" more than once"));
        }
        const auto number = it->getNumber();
        if (number < 1 || number > maxCohortNumber ||
            !std::has_single_bit(static_cast<uint32_t>(number))) {
            throw util::xml_scenario_error(
                std::string("cohort specification assigns sub-population \"").append(it->getId())
                    .append("\" a number which is not a power of 2 between 1 and ")
                    .append(std::to_string(maxCohortNumber)));
        }
        runtime.cohortSubPopNumbers.push_back(number);
        nextId += 1;
    }
}

namespace Continuous {

using util::xml_scenario_error;

        /* Initialise: enable outputs registered and requested in XML.
         * Search for Continuous::registerCallback to see outputs available. */
        void init (RunContext& ctx, const scnXml::Monitoring& monitoring, bool isCheckpoint) {
            auto& ctsOStream = ctx.ext->ctsOStream;
            auto& registered = ctx.ext->ctsRegistered;
            auto& toReport = ctx.ext->ctsToReport;
            SimTime& ctsPeriod = ctx.ext->ctsPeriod;
            bool& duringInit = ctx.ext->ctsDuringInit;
            auto& streamOffset = ctx.ext->ctsStreamOffset;
            const scnXml::Monitoring::ContinuousOptional& ctsOpt = monitoring.getContinuous();
            if( !ctsOpt.present() ) {
                ctsPeriod = sim::zero();
                return;
            }
            const auto& options = ctsOpt.get();
            try{
                //NOTE: if changing XSD, this should not have a default unit:
                ctsPeriod = UnitParse::readShortDuration( ctx.clock, options.getPeriod(), UnitParse::STEPS );
                if( ctsPeriod < ctx.clock.one_TS() )
                    throw util::format_error("must be >= 1 time step");
            }catch( const util::format_error& e ){
                throw xml_scenario_error( string("monitoring/continuous/period: ").append(e.message()) );
            }

            if( options.getDuringInit().present() )
                duringInit = options.getDuringInit().get();

            for (const auto& option : options.getOption()) {
                auto registeredOutput = registered.find(option.getName());
                if (registeredOutput == registered.end()) {
                    throw xml_scenario_error("monitoring.continuous: no output " + string(option.getName()));
                }
                if (option.getValue()) toReport.push_back(registeredOutput->second);
            }

            const string filename = util::CommandLine::getCtsoutName();

            if( isCheckpoint ){
                // When loading a check-point, we resume reporting to this file.
                ctsOStream.open(filename, ios::binary|ios::ate|ios::in|ios::out);
                if( ctsOStream.fail() )
                    throw util::checkpoint_error ("Continuous: resume error (no file)");
                return; // checkpoint() sets the saved output position later
            }

            ctsOStream.open(filename, ios::binary|ios::out);
            ctsOStream << "##\t##\n"; // live-graph needs a delimiter specifier when it is not a comma
            if( duringInit ) ctsOStream << "simulation time\t";
            ctsOStream << "timestep"; // TODO: change to days or remove or leave?
            for (const auto& report : toReport) ctsOStream << std::get<0>(report);
            ctsOStream << '\n' << flush;
            streamOffset = ctsOStream.tellp();
        }

} // namespace Continuous

} // namespace mon
} // namespace OM
