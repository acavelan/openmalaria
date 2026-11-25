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
#include "RunContext.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/CM5DayCommon.h"

#include "Clinical/EventScheduler.h"
#include "Host/NeonatalMortality.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "util/UnitParse.h"

namespace OpenMalaria {

bool opt_event_scheduler = false;
bool opt_imm_outcomes = false;

bool ClinicalModel::indirectMortBugfix;
SimTime ClinicalModel::healthSystemMemory{never()};

//log odds ratio of case-fatality in community compared to hospital
double oddsRatioThreshold;

AgeGroupInterpolator ClinicalModel::caseFatalityRate;
AgeGroupInterpolator ClinicalModel::pSequelaeInpatient;

// -----  static methods  -----

unique_ptr<ClinicalModel> ClinicalModel::createClinicalModel(RunContext& ctx, double tSF) {
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

ClinicalModel::ClinicalModel(RunContext& ctx, double tSF, const Kind kind)
    : kind(kind), ctx(ctx), m_treatmentSeekingFactor(tSF), latestReport(ctx), doomed(NOT_DOOMED) {}

// parity: guard = old ClinicalModel ctor's Kind::EventScheduler check, verbatim
// (no RNG draws in any ctor; ctor layering cannot shift the stream).
CMEventScheduler::CMEventScheduler(RunContext& ctx, double tSF) : ClinicalModel(ctx, tSF, Kind::EventScheduler) {
    if (tSF != 1.0) {
        // p(treatment seeking) is part of tree & handled generically, so we
        // don't have a way of modifying it. (parity: old ES ctor)
        throw xml_scenario_error("treatment seeking heterogeneity not supported");
    }
}

// Variant dispatch: each override calls its variant's free-function body
// (this TU is the one that resolves the ClinicalModel<->EventScheduler/
// CM5DayCommon include cycle, so the hooks live here). Bodies unchanged.
void CMDecisionTree5Day::doClinicalUpdate(Human& human, double ageYears, StatePair& pg) {
    clinicalUpdate5Day(*this, ctx, human, ageYears, pg);
}

void CMDecisionTree5Day::uncomplicatedEvent(Human& human, Episode::State pgState) {
    treeUncomplicatedEvent(*this, ctx, human, pgState);
}

void CMImmediateOutcomes::doClinicalUpdate(Human& human, double ageYears, StatePair& pg) {
    clinicalUpdate5Day(*this, ctx, human, ageYears, pg);
}

void CMImmediateOutcomes::uncomplicatedEvent(Human& human, Episode::State pgState) {
    immOutUncomplicatedEvent(*this, ctx, human, pgState);
}

void CMEventScheduler::doClinicalUpdate(Human& human, double ageYears, StatePair& pg) {
    esClinicalUpdate(*this, ctx, human, ageYears, pg);
}

SimTime CMEventScheduler::timeOfLastTreatment() const { return es_time_treatment; }

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
        recordEvent(measure("nIndDeaths"), human);
        doomed = DOOMED_INDIRECT;
        return;
    }
    if (newBorn /* i.e. first update since birth */) {
        // Chance of neonatal mortality:
        if (NeonatalMortality::eventNeonatalMortality(human.rng())) {
            recordEvent(measure("nIndDeaths"), human);
            doomed = DOOMED_NEONATAL;
            return;
        }
    }

    const bool isDoomed = doomed != NOT_DOOMED;
    StatePair pg = human.withinHostModel->determineMorbidity(human, ageYears, isDoomed);
    latestState = static_cast<Episode::State>(pg.state);

    doClinicalUpdate(human, ageYears, pg);
}

void ClinicalModel::updateInfantDeaths(SimTime age) {
    // update array for the infant death rates
    if (age < oneYear()) {
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

void InfantMortality::init(RunContext& ctx, const Parameters& parameters) {
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

// Episode reporting shares the clinical model's Human/WH include cone.
Episode::~Episode() { report(); }

void Episode::flush() {
    report();
    time = never();
}

void Episode::update(RunContext& ctx, const Human& human, Episode::State newState) {
    bool toReport = false;
    if (healthSystemMemoryFix && (time + ClinicalModel::hsMemory() <= ctx.clock.ts0()))
        toReport = true;
    else if (!healthSystemMemoryFix && (time + ClinicalModel::hsMemory() < ctx.clock.ts0()))
        toReport = true;

    if (toReport) {
        infectionType = human.withinHostModel->getInfectionOrigin();
        report();
        time = ctx.clock.ts0();
        surveyPeriod = eventSurveyNumber(ctx);
        ageGroup = human.monitoringAgeGroup();
        cohortSet = human.getCohortSet();
        state = newState;
    } else {
        state = Episode::State(state | newState);
    }
}

void Episode::report() {
    if (time < zero()) return; // Nothing to report

    // parity: Episode.cpp:67-115 — report order is the output column order.
    if (state & MALARIA) {
        if (state & COMPLICATED) {
            recordEvent(*rctx, measure("nSevere"), surveyPeriod, ageGroup, cohortSet);
            if (state & SEVERE)
                recordEvent(*rctx, measure("nSevereWithoutComorbidities"), surveyPeriod, ageGroup, cohortSet);
        } else {
            recordEvent(*rctx, measure("nUncomp"), surveyPeriod, ageGroup, cohortSet);
            if (infectionType == InfectionOrigin::Indigenous)
                recordEvent(*rctx, measure("nUncomp_Indigenous"), surveyPeriod, ageGroup, cohortSet);
            else if (infectionType == InfectionOrigin::Introduced)
                recordEvent(*rctx, measure("nUncomp_Introduced"), surveyPeriod, ageGroup, cohortSet);
            else
                recordEvent(*rctx, measure("nUncomp_Imported"), surveyPeriod, ageGroup, cohortSet);
        }
        if (state & EVENT_IN_HOSPITAL) {
            if (state & DIRECT_DEATH) {
                recordEvent(*rctx, measure("nDirDeaths"), surveyPeriod, ageGroup, cohortSet);
                recordEvent(*rctx, measure("nHospitalDeaths"), surveyPeriod, ageGroup, cohortSet);
                if (state & EVENT_FIRST_DAY) {
                    recordEvent(*rctx, measure("Clinical_FirstDayDeaths"), surveyPeriod, ageGroup, cohortSet);
                    recordEvent(*rctx, measure("Clinical_HospitalFirstDayDeaths"), surveyPeriod, ageGroup, cohortSet);
                }
            } else if (state & SEQUELAE) {
                recordEvent(*rctx, measure("nSeq"), surveyPeriod, ageGroup, cohortSet);
                recordEvent(*rctx, measure("nHospitalSeqs"), surveyPeriod, ageGroup, cohortSet);
            } else if (state & RECOVERY) {
                recordEvent(*rctx, measure("nHospitalRecovs"), surveyPeriod, ageGroup, cohortSet);
            }
        } else {
            if (state & DIRECT_DEATH) {
                recordEvent(*rctx, measure("nDirDeaths"), surveyPeriod, ageGroup, cohortSet);
                if (state & EVENT_FIRST_DAY)
                    recordEvent(*rctx, measure("Clinical_FirstDayDeaths"), surveyPeriod, ageGroup, cohortSet);
            } else if (state & SEQUELAE) {
                recordEvent(*rctx, measure("nSeq"), surveyPeriod, ageGroup, cohortSet);
            }
        }
    } else if (state & SICK) {
        recordEvent(*rctx, measure("nNMFever"), surveyPeriod, ageGroup, cohortSet);
        if (state & DIRECT_DEATH) recordEvent(*rctx, measure("nNmfDeaths"), surveyPeriod, ageGroup, cohortSet);
    }
}

}
