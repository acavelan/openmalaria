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

#ifndef Hmod_ClinicalParams
#define Hmod_ClinicalParams

#include "Global.h"
#include "Parameters.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"
#include "util/AgeGroupInterpolation.h"

#include <vector>

namespace scnXml {
    class Scenario;
    class Clinical;
    class HealthSystem;
    class HSEventScheduler;
    class HSImmediateOutcomes;
    class HSDT5Day;
}
namespace OM { namespace interventions { struct InterventionRegistry; } }

namespace OM {
namespace Clinical {

class CMDecisionTree;

/// Uncomplicated-case treatment regimens (was CM5DayCommon::CaseType).
enum CaseType { FirstLine, SecondLine, NumCaseTypes };

/** All clinical / health-system scenario parameters: the former static data
 * of ClinicalModel, CM5DayCommon, ClinicalEventScheduler, ImmediateOutcomes
 * and the DecisionTree5Day.cpp / ESCaseManagement.cpp file globals.
 *
 * Owned by Model. Constructed once before the Model (before any Human is
 * created), loaded with health-system data after interventions are set up and
 * reloaded by ChangeHS intervention deployments; read-only from the per-human
 * update path. Not checkpointed: rebuilt from the scenario XML. */
struct ClinicalParams {
    ClinicalParams() = default;
    /** Read model-wide clinical parameters from the scenario (stage 1). */
    ClinicalParams( const Parameters& parameters, const scnXml::Scenario& scenario );

    /** (Re)load all health-system data affected by this healthSystem element
     * (stage 2; also done when a ChangeHS intervention is deployed).
     * Implemented per health-system kind, in the kind's .cpp file. */
    void loadHealthSystem( const scnXml::HealthSystem& healthSystem,
                           const interventions::InterventionRegistry& reg );

    /** Case fatality "rate" in the community as a function of the hospital
     * case-fatality ratio. */
    double getCommunityCFR(double caseFatalityRatio) const;

    // ———  model kind (from ModelOptions / scenario at construction)  ———
    bool eventScheduler = false;       ///< use ClinicalEventScheduler
    bool immOutcomes = false;          ///< use ImmediateOutcomes (else DecisionTree5Day)
    bool indirectMortBugfix = false;   ///< INDIRECT_MORTALITY_FIX option
    bool cfrPfUseHospital = false;     ///< CFR_PF_USE_HOSPITAL option
    bool nonMalariaFevers = false;     ///< NON_MALARIA_FEVERS option (ES only)

    // ———  shared by all clinical models (loadHealthSystem)  ———
    /// The maximum age of a sickness bout, for another bout to be considered
    /// part of the same episode.
    SimTime healthSystemMemory = sim::never();
    /// exp(LOG_ODDS_RATIO_CF_COMMUNITY): log odds ratio of case-fatality in
    /// community compared to hospital.
    double oddsRatioThreshold = 0.0;
    /// Age-specific hospital case fatality "rates"
    util::AgeGroupInterpolator caseFatalityRate;
    /// Age-specific in-hospital rates of sequelae given a severe malaria bout
    /// (out-patients currently have the same probabilities).
    util::AgeGroupInterpolator pSequelaeInpatient;

    // ———  shared 5-day case-management data (CM5DayCommon)  ———
    mon::Measure measures[NumCaseTypes] = {
        mon::measure("nTreatments1"),  // first line official
        mon::measure("nTreatments2")   // second line official
    };
    double accessUCAny[NumCaseTypes] = {};
    double accessUCSelfTreat[NumCaseTypes] = {};
    double accessSevere = 0.0;
    double cureRateSevere = 0.0;
    WithinHost::TreatmentId treatmentSevere = {};

    // ———  ImmediateOutcomes  ———
    double cureRateUCOfficial[NumCaseTypes] = {};
    double cureRateUCSelfTreat[NumCaseTypes] = {};
    WithinHost::TreatmentId treatmentUC[NumCaseTypes] = {};
    bool useDiagnosticUC = false;

    // ———  DecisionTree5Day (case-management trees, owned by the
    /// CMDecisionTree library)  ———
    const CMDecisionTree* treeUCOfficial = nullptr;
    const CMDecisionTree* treeUCSelfTreat = nullptr;

    // ———  ClinicalEventScheduler  ———
    SimTime maxUCSeekingMemory = sim::never();
    SimTime uncomplicatedCaseDuration = sim::never();
    SimTime complicatedCaseDuration = sim::never();
    SimTime extraDaysAtRisk = sim::never();
    std::vector<double> cumDailyPrImmUCTS;
    double neg_v = 0.0;
    double alpha = 0.0;
    // Non-malaria fevers (loadHealthSystem; NaN until then):
    double logOddsAbBase = std::numeric_limits<double>::quiet_NaN();
    double logOddsAbNegTest = std::numeric_limits<double>::quiet_NaN();
    double logOddsAbPosTest = std::numeric_limits<double>::quiet_NaN();
    double logOddsAbNeed = std::numeric_limits<double>::quiet_NaN();
    double logOddsAbInformal = std::numeric_limits<double>::quiet_NaN();
    double oneMinusEfficacyAb = std::numeric_limits<double>::quiet_NaN();
    util::AgeGroupInterpolator severeNmfMortality;
    util::AgeGroupInterpolator NMF_need_antibiotic;
    util::AgeGroupInterpolator MF_need_antibiotic;
    // Event-scheduler case-management trees (owned by the CMDecisionTree library):
    const CMDecisionTree* escmUncomplicated = nullptr;
    const CMDecisionTree* escmComplicated = nullptr;
};

}
}
#endif
