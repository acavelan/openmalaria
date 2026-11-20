

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

// The 5-day case-management family: shared config (CM5DayCommon) plus the
// two 5-day health-system kinds DecisionTree5Day and ImmediateOutcomes
// (their headers were one-include away from this one; merged here).

#pragma once

#include "Context.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Clinical/ClinicalModel.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"

namespace scnXml {
class HSDT5Day;
class HSImmediateOutcomes;
}
namespace OM {
namespace Clinical {
using Host::Human;

/**
 * Case-management CONFIG shared by the 5-day case management models
 * (DecisionTree5Day and ImmediateOutcomes kinds). Non-polymorphic: static
 * config only; per-human state lives on ClinicalModel, behaviour is the
 * clinicalUpdate5Day free function below.
 */
class CM5DayCommon {
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
 * variant selects the uncomplicatedEvent implementation via the virtual —
 * parity: the old virtual call). */
void clinicalUpdate5Day(ClinicalModel& cm, Context& ctx, Human& human, double ageYears,
                        WithinHost::Pathogenesis::StatePair& pg);

/** Called when a severe/complicated (with co-infection) malaria sickness occurs.
 *
 * Note: sets doomed = 4 if patient dies. */
void severeMalaria(ClinicalModel& cm, Context& ctx, Human& human, Episode::State pgState, double ageYears, int& doomed);

/**
 * This models case management with 5-day time steps, optionally including
 * PK/PD modeling for uncomplicated cases.
 * 
 * Uncomplicated cases: access, otherwise known as "seeking any type of
 * treatment", is determined by a fixed-function decision, which may be
 * modified by a treatment-seeking factor. Treatment decisions (type of
 * treatment, use of diagnostics, effectiveness) are determined by a
 * programmable decision tree.
 * 
 * Severe cases: all decisions and outcomes are calculated via a fixed-function
 * probability tree, using the same logic for handling severe cases as has long
 * been used.
 */
class DecisionTree5Day {
  public:
    /** Load health system data from initial data or an intervention's data (both from XML).
     * (Re)loads all data affected by this healthSystem element. */
    static void setHealthSystem(Context& ctx, const scnXml::HSDT5Day& hsDescription);
};

/// Uncomplicated-case handling for the DecisionTree5Day kind (old
/// DecisionTree5Day::uncomplicatedEvent).
void treeUncomplicatedEvent(ClinicalModel& cm, Context& ctx, Human& human, Episode::State pgState);

/**
 * This models case management at a 5-day time step with all-or-nothing
 * treatment.
 * 
 * Uncomplicated cases: access, otherwise known as "seeking any type of
 * treatment", is determined by a fixed-function decision, which may be
 * modified by a treatment-seeking factor. Treatment decisions are also fixed
 * function.
 * 
 * Severe cases: all decisions and outcomes are calculated via a fixed-function
 * probability tree, using the same logic for handling severe cases as has long
 * been used.
 * 
 * NOTE: this provides a strict subset of the functionality of DecisionTree5Day.
 */
class ImmediateOutcomes {
  public:
    /** Load health system data from initial data or an intervention's data (both from XML).
     * (Re)loads all data affected by this healthSystem element. */
    static void setHealthSystem(Context& ctx, const scnXml::HSImmediateOutcomes& hsioData);

    static double cureRateUCOfficial[CM5DayCommon::NumCaseTypes];
    static double cureRateUCSelfTreat[CM5DayCommon::NumCaseTypes];
    static WithinHost::TreatmentId treatmentUC[CM5DayCommon::NumCaseTypes];
    static bool useDiagnosticUC;
};

/// Uncomplicated-case handling for the ImmediateOutcomes kind (old
/// ImmediateOutcomes::uncomplicatedEvent).
void immOutUncomplicatedEvent(ClinicalModel& cm, Context& ctx, Human& human, Episode::State pgState);

}
}
