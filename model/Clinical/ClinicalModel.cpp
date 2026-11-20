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

// p8-files: header fold rejected — ClinicalModel.cpp resolves the
// ClinicalModel<->EventScheduler/CM5DayCommon include cycle; bodies stay here.
#include "Context.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/CM5DayCommon.h"

#include "Clinical/EventScheduler.h"
#include "Clinical/CM5DayCommon.h"
#include "Host/NeonatalMortality.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "util/UnitParse.h"

namespace OM {
namespace Clinical {

bool opt_event_scheduler = false;
bool opt_imm_outcomes = false;

bool ClinicalModel::indirectMortBugfix;
SimTime ClinicalModel::healthSystemMemory{sim::never()};

//log odds ratio of case-fatality in community compared to hospital
double oddsRatioThreshold;

util::AgeGroupInterpolator ClinicalModel::caseFatalityRate;
util::AgeGroupInterpolator ClinicalModel::pSequelaeInpatient;

// -----  static methods  -----

unique_ptr<ClinicalModel> ClinicalModel::createClinicalModel(Context& ctx, double tSF) {
    // One create-line per variant; the kind flags are set at scenario load
    // (ClinicalModel::init, xml/healthSystem.cpp).
    if (opt_event_scheduler)
        return unique_ptr<ClinicalModel>(new CMEventScheduler(ctx, tSF));
    else if (opt_imm_outcomes)
        return unique_ptr<ClinicalModel>(new CMImmediateOutcomes(ctx, tSF));
    else
        return unique_ptr<ClinicalModel>(new CMDecisionTree5Day(ctx, tSF));
}

double ClinicalModel::getCommunityCFR(double caseFatalityRatio) {
    double x = caseFatalityRatio * oddsRatioThreshold;
    return x / (1 - caseFatalityRatio + x);
}

// -----  non-static construction, destruction and checkpointing  -----

ClinicalModel::ClinicalModel(Context& ctx, double tSF, const Kind kind)
    : kind(kind), ctx(ctx), m_treatmentSeekingFactor(tSF), latestReport(ctx), doomed(NOT_DOOMED) {}

// parity: guard = old ClinicalModel ctor's Kind::EventScheduler check, verbatim
// (no RNG draws in any ctor; ctor layering cannot shift the stream).
CMEventScheduler::CMEventScheduler(Context& ctx, double tSF) : ClinicalModel(ctx, tSF, Kind::EventScheduler) {
    if (tSF != 1.0) {
        // p(treatment seeking) is part of tree & handled generically, so we
        // don't have a way of modifying it. (parity: old ES ctor)
        throw util::xml_scenario_error("treatment seeking heterogeneity not supported");
    }
}

// Variant dispatch: each override calls its variant's free-function body
// (this TU is the one that resolves the ClinicalModel<->EventScheduler/
// CM5DayCommon include cycle, so the hooks live here). Bodies unchanged.
void CMDecisionTree5Day::doClinicalUpdate(Human& human, double ageYears, WithinHost::Pathogenesis::StatePair& pg) {
    clinicalUpdate5Day(*this, ctx, human, ageYears, pg);
}

void CMDecisionTree5Day::uncomplicatedEvent(Human& human, Episode::State pgState) {
    treeUncomplicatedEvent(*this, ctx, human, pgState);
}

void CMImmediateOutcomes::doClinicalUpdate(Human& human, double ageYears, WithinHost::Pathogenesis::StatePair& pg) {
    clinicalUpdate5Day(*this, ctx, human, ageYears, pg);
}

void CMImmediateOutcomes::uncomplicatedEvent(Human& human, Episode::State pgState) {
    immOutUncomplicatedEvent(*this, ctx, human, pgState);
}

void CMEventScheduler::doClinicalUpdate(Human& human, double ageYears, WithinHost::Pathogenesis::StatePair& pg) {
    esClinicalUpdate(*this, ctx, human, ageYears, pg);
}

SimTime CMEventScheduler::timeOfLastTreatment() const { return es_time_treatment(); }

bool ClinicalModel::isExistingCase() {
    // If treated in the recent past:
    const SimTime t = timeOfLastTreatment();
    return ctx.clock.now() > t && ctx.clock.now() <= t + healthSystemMemory;
}

// -----  other non-static methods  -----

bool ClinicalModel::isDead(SimTime age) {
    if (age >= ctx.clock.max_human_age) // too old (reached age limit)
        doomed = DOOMED_TOO_OLD;
    return doomed > NOT_DOOMED; // killed by some means
}

void ClinicalModel::update(Human& human, double ageYears, bool newBorn) {
    if (doomed < NOT_DOOMED) // Countdown to indirect mortality
        doomed -= ctx.clock.one_TS();

    //indirect death: if this human's about to die, don't worry about further episodes:
    if (doomed <= DOOMED_EXPIRED) { //clinical bout 6 intervals before
        mon::recordEvent(mon::measure("nIndDeaths"), human);
        doomed = DOOMED_INDIRECT;
        return;
    }
    if (newBorn /* i.e. first update since birth */) {
        // Chance of neonatal mortality:
        if (Host::NeonatalMortality::eventNeonatalMortality(human.rng())) {
            mon::recordEvent(mon::measure("nIndDeaths"), human);
            doomed = DOOMED_NEONATAL;
            return;
        }
    }

    const bool isDoomed = doomed != NOT_DOOMED;
    WithinHost::Pathogenesis::StatePair pg = human.withinHostModel->determineMorbidity(human, ageYears, isDoomed);
    latestState = static_cast<Episode::State>(pg.state);

    doClinicalUpdate(human, ageYears, pg);
}

void ClinicalModel::updateInfantDeaths(SimTime age) {
    // update array for the infant death rates
    if (age < sim::oneYear()) {
        size_t index = age / ctx.clock.one_TS();

        // Testing doomed == DOOMED_NEXT_TS gives very slightly different results than
        // testing doomed == DOOMED_INDIRECT (due to above if(..))
        bool isDoomed = doomed == DOOMED_COMPLICATED || doomed == DOOMED_NEXT_TS || doomed == DOOMED_NEONATAL;
        InfantMortality::reportRisk(index, isDoomed);
    }
}

// ——— infant mortality ———

// Infant death summaries (checkpointed).
vector<int> infantDeaths;
vector<int> infantIntervalsAtRisk;

/// Non-malaria mortality in under 1year olds.
/// Set by init ()
double nonMalariaMortality;

void InfantMortality::init(Context& ctx, const OM::Parameters& parameters) {
    infantDeaths.resize(ctx.clock.steps_per_year);
    infantIntervalsAtRisk.resize(ctx.clock.steps_per_year);
    nonMalariaMortality = parameters[Parameter::NON_MALARIA_INFANT_MORTALITY];
}

void InfantMortality::preMainSimInit() {
    std::fill(infantDeaths.begin(), infantDeaths.end(), 0);
    std::fill(infantIntervalsAtRisk.begin(), infantIntervalsAtRisk.end(), 0);
}

void InfantMortality::reportRisk(size_t index, bool isDoomed) {
    infantIntervalsAtRisk[index] += 1; // baseline
    if (isDoomed)
        infantDeaths[index] += 1; // deaths
}

double InfantMortality::allCause() {
    double infantPropSurviving = 1.0; // use to calculate proportion surviving
    for (size_t i = 0; i < infantIntervalsAtRisk.size(); i += 1) {
        // multiply by proportion of infants surviving at each interval
        infantPropSurviving *= double(infantIntervalsAtRisk[i] - infantDeaths[i]) / double(infantIntervalsAtRisk[i]);
    }
    // Child deaths due to malaria (per 1000), plus non-malaria child deaths. Deaths per 1000 births is the return unit.
    return (1.0 - infantPropSurviving) * 1000.0 + nonMalariaMortality;
}

}
}
