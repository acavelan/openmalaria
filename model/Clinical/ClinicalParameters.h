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

#pragma once
#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/Treatment.h"
#include "util/AgeGroupInterpolation.h"
#include "mon/OutMeasures.h"

namespace OM::Clinical {
struct ClinicalParameters {
    enum class Kind { DecisionTree5Day, ImmediateOutcomes, EventScheduler };
    Kind kind = Kind::DecisionTree5Day;
    vector<unique_ptr<CMDecisionTree>> nodes;
    const CMDecisionTree* escm_uncomplicated = nullptr;
    const CMDecisionTree* escm_complicated = nullptr;
    const CMDecisionTree* treeUCOfficial = nullptr;
    const CMDecisionTree* treeUCSelfTreat = nullptr;

    bool indirectMortBugfix = false;
    // Maximum gap between sickness bouts belonging to the same episode.
    SimTime healthSystemMemory = sim::never();
    // Community-to-hospital case-fatality odds ratio.
    double oddsRatioThreshold = 0.0;
    util::AgeGroupInterpolator caseFatalityRate; // hospital case-fatality rate by age
    util::AgeGroupInterpolator pSequelaeInpatient; // sequelae probability by age

    // Shared five-day care probabilities and immediate-outcome treatments.
    mon::Measure measures[2] = {mon::measure("nTreatments1"), mon::measure("nTreatments2")};
    double accessUCAny[2] = {};
    double accessUCSelfTreat[2] = {};
    double accessSevere = 0.0;
    double cureRateSevere = 0.0;
    WithinHost::Treatment treatmentSevere;
    bool cfr_pf_use_hospital = false;
    double cureRateUCOfficial[2] = {};
    double cureRateUCSelfTreat[2] = {};
    WithinHost::Treatment treatmentUC[2];
    bool useDiagnosticUC = false;

    // Daily event-scheduler configuration.
    bool opt_non_malaria_fevers = false;
    /// Maximum number of time steps (including first of case) an individual
    /// will remember they are sick before resetting.
    SimTime maxUCSeekingMemory = sim::never();
    /// Length of an uncomplicated case
    SimTime uncomplicatedCaseDuration = sim::never();
    /// Length of a complicated case
    SimTime complicatedCaseDuration = sim::never();
    /// Time-span for which individual is at risk of death in complicated case
    /// minus length of complicated case (must be <= 0)
    SimTime extraDaysAtRisk = sim::never();
    /// First value is probability immediate treatment, second is first +
    /// probability 1-day delay to treatment seeking, etc. Last value must be 1.
    /// Index units are days.
    vector<double> cumDailyPrImmUCTS;

    /// Parameter of S(t) for t > 0
    double neg_v = 0.0;
    /// Parameter
    double alpha = 0.0;

    /** Base log odds of treatment of non-malarial fevers in absense of a
     * malaria diagnostic and irrespective of whether treatment is needed.
     *
     * In our model, this is logit(P₀), not β₀. */
    double logOddsAbBase = numeric_limits<double>::signaling_NaN();
    /** Added to log odds treatment when a malaria diagnostic indicates no
     * parasites. Symbol in model: β₁. */
    double logOddsAbNegTest = numeric_limits<double>::signaling_NaN();
    /** Added to log odds treatment when a malaria diagnostic indicates
     * parasites. Symbol in model: β₂. */
    double logOddsAbPosTest = numeric_limits<double>::signaling_NaN();
    /** Added to log odds treatment when NMF is categorized as an illness
     * potentially leading to death (Pathogenesis::NEED_ANTIBIOTIC).
     * Symbol in model: β₃. */
    double logOddsAbNeed = numeric_limits<double>::signaling_NaN();
    /** Added to log odds treatment when given by an informal provider.
     * Symbol in model: β₄. */
    double logOddsAbInformal = numeric_limits<double>::signaling_NaN();
    /** One minus the efficacy of antibiotic/NMF treatment (i.e. a multiplier
     * for fatality-rate given that the case is treated). */
    double oneMinusEfficacyAb = numeric_limits<double>::signaling_NaN();
    /** Case fatality rate of non-malaria fevers requiring treatment given that
     * the case is not treated. */
    util::AgeGroupInterpolator severeNmfMortality;

    /// Probability that an NMF needs antibiotic treatment and could lead to death.
    util::AgeGroupInterpolator NMF_need_antibiotic;
    /** Probability that a malarial fever classified as uncomplicated requires
     * antibiotic treatment (and death could occur from non-malarial causes).
     * Unclear whether using this at the same time as comorbidity parameters
     * makes any sense. */
    util::AgeGroupInterpolator MF_need_antibiotic;
};
}
