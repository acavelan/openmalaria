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

// p8-files: header fold rejected — this TU reads Human/WHInterface, so inlining
// into Episode.h would drag that chain into the whole clinical cone.
#include "RunContext.h"
#include "Clinical/Episode.h"
#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"

namespace OpenMalaria {

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
    if (time < zero()) // Nothing to report
        return;

    // Reports malarial/non-malarial UC fever dependent on cause, not diagnosis.
    if (state & MALARIA) {
        // Malarial fevers: report bout
        if (state & COMPLICATED) {
            recordEvent(*rctx, measure("nSevere"), surveyPeriod, ageGroup, cohortSet);
            if (state & SEVERE) {
                recordEvent(*rctx, measure("nSevereWithoutComorbidities"), surveyPeriod, ageGroup, cohortSet);
            }
        } else { // UC or UC2
            recordEvent(*rctx, measure("nUncomp"), surveyPeriod, ageGroup, cohortSet);
            if (infectionType == InfectionOrigin::Indigenous)
                recordEvent(*rctx, measure("nUncomp_Indigenous"), surveyPeriod, ageGroup, cohortSet);
            else if (infectionType == InfectionOrigin::Introduced)
                recordEvent(*rctx, measure("nUncomp_Introduced"), surveyPeriod, ageGroup, cohortSet);
            else
                recordEvent(*rctx, measure("nUncomp_Imported"), surveyPeriod, ageGroup, cohortSet);
        }

        // Report outcomes of malarial fevers
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
                if (state & EVENT_FIRST_DAY) {
                    recordEvent(*rctx, measure("Clinical_FirstDayDeaths"), surveyPeriod, ageGroup, cohortSet);
                }
            } else if (state & SEQUELAE) {
                recordEvent(*rctx, measure("nSeq"), surveyPeriod, ageGroup, cohortSet);
            }
            // Don't care about out-of-hospital recoveries
        }
    } else if (state & SICK) {
        // Report non-malarial fever and outcomes
        recordEvent(*rctx, measure("nNMFever"), surveyPeriod, ageGroup, cohortSet);

        if (state & DIRECT_DEATH) {
            recordEvent(*rctx, measure("nNmfDeaths"), surveyPeriod, ageGroup, cohortSet);
        }
    }
}

}
