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

#include "Host/Human.h"
#include "Episode.h"
#include <memory>

namespace scnXml{
    class Scenario;
    class HealthSystem;
}
namespace OpenMalaria {
struct StatePair;

/// Positive values of doomed (codes). (parity: old ClinicalModel members)
enum DoomedCause {
        DOOMED_EXPIRED = -35,   // codes less than or equal to this mean "dead now"
        DOOMED_NEXT_TS = -30,   // will expire on next time step
        NOT_DOOMED = 0, // all codes greater than this mean "already dead"; codes less than this mean a count-down to death has started
        DOOMED_TOO_OLD = 1,		///< died because reached age limit
        DOOMED_COMPLICATED = 4,	///< died from severe malaria or malaria with a coinfection
        DOOMED_NEONATAL = 6,	///< died due to mother's malaria infection
        DOOMED_INDIRECT = 7,	///< died indirectly from malaria (after a delay)
};

/** The clinical model models the effects of sickness dependant on malarial
 * parasite densities and administers anti-malaria treatments via the drug
 * model (or in a simpler case, directly clearing infections).
 * 
 * So far, sickness types include uncomplicated and severe malaria cases and
 * non-malaria sickness.
 * 
 * Patient outcomes include full recovery, recovery with sequelae and death.
 * 
 * Reporting includes patient outcome and potentially drug usage and use of
 * RDTs (Rapid Diagnostic Tests) for costing purposes.
 *
 * Base of the case-management variant hierarchy (p8-poly pattern): holds the
 * state every variant shares (treatment-seeking factor, episode report/state,
 * doomed codes, column bind) and draws NO random numbers; each XML-selected
 * variant is a clearly independent derived class below with its own update
 * functions. createClinicalModel() selects the variant (one create-line
 * each). */\
class ClinicalModel
{
protected:
    /// The run this clinical model belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

public:
    /// Which case-management behaviour this instance runs (XML-selected).
    enum class Kind { DecisionTree5Day, ImmediateOutcomes, EventScheduler };

    /// @brief Static functions
    //@{
    /// First stage of initialisation
    static void init ( RunContext& ctx, const Parameters& parameters, const scnXml::Scenario& model );
  
    /** Second stage of initialisation, done after interventions are configured
     * 
     * Also done when a certain intervention is deployed. */
    static void setHS (RunContext& ctx, const scnXml::HealthSystem& healthSystem);
    
    /** Return a new ClinicalModel.
     *
     * @param tSF	treatment seeking factor, passed to CaseManagementModel */
    static unique_ptr<ClinicalModel> createClinicalModel (RunContext& ctx, double tSF);
    //@}
    
    /// Destructor
    
    /** Returns true if the human has been killed by some means.
     * 
     * Also kills the human if he/she reaches the simulation age limit.
     */
    bool isDead( SimTime age );
    
    /** Run main part of the model: determine the sickness status and any
     * treatment for the human.
     * 
     * @param ageYears Age of human.
     * @param newBorn True if human age is one time step old */
    void update (Human& human, double ageYears, bool newBorn);
    
    /** For infants, updates the infantIntervalsAtRisk and potentially
     * infantDeaths arrays. */
    void updateInfantDeaths( SimTime age );
    
    /** Special option to allow reports not to be delivered for existing cases
     * (within health-system-memory and not new cases). */
    /** Latest episode is treated (has been given treatment)?
     * (some interventions check this to prevent multiple simultaneous treatments) */
    bool isExistingCase();
    
    inline static SimTime hsMemory() {
        return healthSystemMemory;
    }
    
    /// Force all pending summaries to be reported. Should only be called when
    /// class is about to be destroyed anyway to avoid affecting output.
    inline void flushReports (){
        latestReport.flush();
    }

    /** Time of the last treatment (never() if never treated): SoA
     * column t_last_treatment (old member; accessors keep read/write syntax). */
    SimTime& t_last_treatment() const { return cols->t_last_treatment[row]; }

    /** Bind the host human's column row (Human ctor after creation; see
     * Human::set_row for renumber rebinds). */
    void bind_columns( HumanColumns& c, size_t r ) { cols = &c; row = r; }
    
    /** Construct a per-human instance (base: no RNG draws; the EventScheduler
     * variant rejects treatment-seeking heterogeneity in its own ctor). */
    ClinicalModel (RunContext& ctx, double tSF, Kind kind);

    virtual ~ClinicalModel() {}

    /// Variant behaviour (public: the shared 5-day free function
    /// clinicalUpdate5Day dispatches through uncomplicatedEvent).
    /// Base bodies never run: every variant overrides.

    /// Per-kind clinical update tail (old virtual doClinicalUpdate; the
    /// shared preamble of update() ends here).
    virtual void doClinicalUpdate (Human& human, double ageYears, StatePair& pg) {
        assert(false); // no variant updates through the bare base
    }

    /// Uncomplicated-case handling (old virtual uncomplicatedEvent); reached
    /// only by the two 5-day kinds.
    virtual void uncomplicatedEvent (Human& human, Episode::State pgState) {
        assert(false); // 5-day kinds only
    }

    /// isExistingCase() memory clock (ES keeps its own column).
    virtual SimTime timeOfLastTreatment() const { return t_last_treatment(); }

protected:

public:
    /// Per-human State (public: the variant's update free function
    /// [clinicalUpdate5Day / esClinicalUpdate] reads/writes it directly).
    const Kind kind;
    /// treatment seeking for heterogeneity (5-day kinds; ES kind: always 1)
    double m_treatmentSeekingFactor;
    /** Last episode; report to survey pending a new episode or human's death. */
    Episode latestReport;
    Episode::State latestState;

    /// The host human's SoA columns/row (bind_columns; re-bound by
    /// Human::set_row — stale mirror = wrong host after renumber).
    HumanColumns* cols = nullptr;
    size_t row = 0;

    
    /** Can indicate that the individual is dead or about to die.
     *
     * If doomed < 0, the individual is doomed to die.
     * 
     * If doomed > 0, the individual is dead, and will be removed from the
     * population at the beginning of the next time step. NOTE:
     * updateInfantDeaths() counts deaths after the fact, thus cannot remove immediately. */
    int doomed;
    
    
    /// True if bug-fix is enabled
    static bool indirectMortBugfix;
    
    /// The maximum age of a sickness bout, for another bout to be considered
    /// part of the same episode.
    /// 
    /// Used by both the clinical models in roughly the same way, but will have
    /// different values in each to match Global::interval.
    static SimTime healthSystemMemory;
    
    ///@brief Case fatality and sequelae "rate" data
    //@{
    /// Calculate the case fatality "rate" in the community as a function of
    /// that in hospitals.
    static double getCommunityCFR(double caseFatalityRatio);
    
    /// Age-specific hospital case fatality "rates"
    static AgeGroupInterpolator caseFatalityRate;
    /// Age-specific in-hospital rates of sequelae given a severe malaria bout
    /// Note: out-patients have currently have the same probabilities of sequelae
    static AgeGroupInterpolator pSequelaeInpatient;
    //@}
};

/** DecisionTree5Day health system (old ClinicalDecisionTree5Day):
 * uncomplicated-case access by fixed function, treatment by the compiled
 * decision tree. Case-management config statics live in CM5DayCommon /
 * DecisionTree5Day; bodies in CM5DayCommon.cpp. */
class CMDecisionTree5Day final : public ClinicalModel {
  public:
    CMDecisionTree5Day(RunContext& ctx, double tSF) : ClinicalModel(ctx, tSF, Kind::DecisionTree5Day) {}
    void doClinicalUpdate(Human& human, double ageYears, StatePair& pg) override;
    void uncomplicatedEvent(Human& human, Episode::State pgState) override;
};

/** ImmediateOutcomes health system (old ClinicalImmediateOutcomes): all-or-
 * nothing treatment, no decision tree. Config statics live in
 * ImmediateOutcomes; bodies in CM5DayCommon.cpp. */
class CMImmediateOutcomes final : public ClinicalModel {
  public:
    CMImmediateOutcomes(RunContext& ctx, double tSF) : ClinicalModel(ctx, tSF, Kind::ImmediateOutcomes) {}
    void doClinicalUpdate(Human& human, double ageYears, StatePair& pg) override;
    void uncomplicatedEvent(Human& human, Episode::State pgState) override;
};

/** EventScheduler health system (old ClinicalEventScheduler): continuous-time
 * bouts with scheduled recovery; case management via compiled trees. Config
 * statics live in ClinicalEventScheduler / ESCaseManagement; body in
 * EventScheduler.cpp. Owns the ES bout state (parity: old
 * ClinicalEventScheduler members; the other kinds never had these fields). */
class CMEventScheduler final : public ClinicalModel {
  public:
    /// Rejects treatment-seeking heterogeneity (old ES ctor guard).
    CMEventScheduler(RunContext& ctx, double tSF);
    void doClinicalUpdate(Human& human, double ageYears, StatePair& pg) override;
    SimTime timeOfLastTreatment() const override;

  public: // run state; wire in checkpoint/humans.cpp
    Episode::State esPgState = NONE;    ///< current sickness state
    /// ES bout clock/state (old esCaseStartTime/esTimeOfRecovery/
    /// esTimeLastTreatment SimTime members; column row refs, same syntax).
    SimTime& es_case_start() const { return cols->es_case_start[row]; }
    SimTime& es_time_recovery() const { return cols->es_time_recovery[row]; }
    SimTime& es_time_treatment() const { return cols->es_time_treatment[row]; }
    /// Total density last time step (old esPreviousDensity; column).
    double& es_prev_density() const { return cols->es_prev_density[row]; }
};


/// Infant mortality reporting
class InfantMortality {
public:
    /** Initialise parameters. */
    static void init( RunContext& ctx, const Parameters& parameters );
    static void preMainSimInit();

    static void reportRisk(size_t index, bool isDoomed);
    
    /// Calculate infant mortality as deaths/1000 livebirths for the whole main-
    /// simulation period (not as deaths/1000 years-at-risk per survey).
    /// 
    /// This mimicks field data on all-cause mortality in infants.
    /// Uses the kaplan-meier method because the demography was set up to provide
    /// a stable age-distribution but unfortunately does not accurately describe
    /// death rates. The kaplan-meier estimate is the product of the proportion of
    /// infants survivng at each interval.
    static double allCause();
};


}
