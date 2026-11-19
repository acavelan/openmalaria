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

#ifndef Hmod_CM5DayCommon
#define Hmod_CM5DayCommon

#include "core/RunContext.h"
#include "Host/WithinHost/Pathogenesis/State.h"
#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"

namespace OM {
namespace Clinical {
using Host::Human;

/**
 * Case-management CONFIG shared by the 5-day case management models
 * (DecisionTree5Day and ImmediateOutcomes kinds). Non-polymorphic: static
 * config only; per-human state lives on ClinicalModel, behaviour is the
 * clinicalUpdate5Day free function below.
 */
class CM5DayCommon
{
public:
    static void init();

    enum CaseType { FirstLine, SecondLine, NumCaseTypes };
    static mon::Measure measures[NumCaseTypes];
    static double accessUCAny[NumCaseTypes];
    static double accessUCSelfTreat[NumCaseTypes];
    static double accessSevere;
    static double cureRateSevere;
    static WithinHost::TreatmentId treatmentSevere;
};

/** Common 5-day clinical update (old CM5DayCommon::doClinicalUpdate; the
 * kind of cm selects the uncomplicatedEvent implementation, parity: the old
 * virtual call). */
void clinicalUpdate5Day (ClinicalModel& cm, RunContext& ctx, Human& human,
                         double ageYears, WithinHost::Pathogenesis::StatePair &pg);

/** Called when a severe/complicated (with co-infection) malaria sickness occurs.
 *
 * Note: sets doomed = 4 if patient dies. */
void severeMalaria (ClinicalModel& cm, RunContext& ctx, Human& human,
                 Episode::State pgState, double ageYears, int& doomed);

}
}
#endif
