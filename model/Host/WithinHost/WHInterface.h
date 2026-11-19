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

#ifndef Hmod_WithinHost_Interface
#define Hmod_WithinHost_Interface

#include "Global.h"
#include "core/RunContext.h"
#include "humans/HostColumns.h"
#include "infections/SlotAllocator.h"
#include "util/random.h"
#include "util/sampler.h"
#include "util/errors.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/Pathogenesis/State.h"
#include "Host/WithinHost/Infection/Infection.h"

#include "Parameters.h"

#include "PkPd/LSTMModel.h"

#include <list>
#include <memory>

using namespace std;

namespace scnXml{
    class Scenario;
    class Model;
    class TreatmentOption;
    class LiverStageDrug;
}
class UnittestUtil;
class InfectionImmunitySuite;

namespace OM {
namespace Host {
    class Human;
}
namespace WithinHost {
namespace Pathogenesis {
    class PathogenesisModel;
}

class WHInterface;

using util::LocalRng;

/// True if any by-genotype reporting is enabled (read-only)
extern bool reportInfectionsByGenotype;

struct HypnozoiteReleaseDistribution {
    /// The run this sampler serves (clock accessed as ctx.clock).
    RunContext& ctx;

    HypnozoiteReleaseDistribution(RunContext& ctx, std::unique_ptr<util::Sampler> sampler, double latentRelapse) :
        ctx(ctx),
        sampler(std::move(sampler)),
        latentRelapse(latentRelapse)
        {}
    
    /// Sample the time until next release
    SimTime sampleReleaseDelay(LocalRng& rng) const {
        double liverStageMaximumDays = 16.0*30.0; // maximum of about 16 months in liver stage 
        double delay = numeric_limits<double>::quiet_NaN();       // in days
        int count = 0;
        int maxcount = 1e3;
        
        do{
            delay = sampler->sample(rng);
            count += 1;
            
            if( count >= maxcount ){
                throw util::xml_scenario_error( "<vivax><hypnozoiteRelease>  [random delay calculation causes probably an indefinite loop]:\n The hypnozoite release distribution seems off, sigma of secondRelease could be too high. We except the hypnozoite to reside a maximum of 16 months in the liver stage. Sigma choose well, dear padawan." );
            }
        } while( delay > liverStageMaximumDays || delay < 0.0  );
        
        assert( delay >= 0 && delay < liverStageMaximumDays );
        return ctx.clock.roundToTSFromDays( delay + latentRelapse );
    }
    
private:
    std::unique_ptr<util::Sampler> sampler;
    double latentRelapse;   // days
};

/**
 * A brood is the set of hypnozoites resulting from an innoculation, plus an
 * associated blood stage.
 * 
 * In this model, if a hypnozoite releases while a blood stage infection
 * initiated by another hypnozoite from the same brood is active, the newly
 * released hypnozoite does nothing, however, blood stage infections from other
 * broods have no effect.
 */
class VivaxBrood{
public:
    /** Create.
     * 
     * @param host      The host creating this. Not really needed, except to
     * prevent VivaxBrood being default-constructed by its container.
     */
    /// The run this brood belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    VivaxBrood( RunContext& ctx, LocalRng& rng, int origin, WHInterface *host );
    ~VivaxBrood();
    /** Save a checkpoint. */
    void checkpoint( ostream& stream );
    /** Create from checkpoint. */
    VivaxBrood( RunContext& ctx, istream& stream );
    
    struct UpdResult{
        UpdResult() : newPrimaryBS(false), newRelapseBS(false), newRelapsebBS(false), newBS(false) {}
        bool newPrimaryBS, newRelapseBS, newRelapsebBS, newBS, isFinished;
    };
    /**
     * Do per time step update: remove finished blood stage infections and act
     * on newly releasing hypnozoites.
     * 
     * @return pair describing whether the infection is finished (no more blood
     *  or liver stages) and whether a new primary blood-stage has started (to
     *  update cum_prim_inf_ref())
     */
    UpdResult update(LocalRng& rng);
    
    inline void setHadEvent( bool hadEvent ){ this->hadEvent = hadEvent; }
    inline bool hasHadEvent()const{ return hadEvent; }
    inline void setHadRelapse( bool hadRelapse ){ this->hadRelapse = hadRelapse; }
    inline bool hasHadRelapse()const{ return hadRelapse; }
    
    /** Equivalent to a blood stage existing. We do not model incidence of
     * gametocytes independently, thus this also tests existance of
     * gametocytes. */
    inline bool isPatent() const{
        return bloodStageClearDate > ctx.clock.latestTs0();
    }
    
    /** Fully clear blood stage parasites. */
    void treatmentBS();
    
    /** Fully clear liver stage parasites. */
    void treatmentLS();
    
private:
    
    // list of times at which the merozoite and hypnozoites release, ordered by
    // time of release, soonest last (i.e. last element is next one to release)
    vector<SimTime> releaseDates;
    
    // Either sim::never() (no blood stage) or the start of the time step on
    // which the blood stage will clear.
    SimTime bloodStageClearDate = sim::never();
    
    // Whether the primary blood stage infection has started
    bool primaryHasStarted;

    // Whether the relapse blood stage infection has started
    bool relapseHasStarted;
	
	//Whether the later category of relapses has started
	bool relapsebHasStarted;

    // Whether any clinical event has been triggered by this brood
    bool hadEvent;

    // Whether a relapse event has been triggered by this brood
    bool hadRelapse;
	
	// Infection origin
	int origin;
};

/**
 * The within-host models. These models encapsulate the infections and related
 * immunity factors of a single human, starting with infection (i.e. assuming
 * successful innoculation), including some drug action code, and outputting
 * parasite densities.
 *
 * One final class covering all variants, selected by Kind (from the scenario's
 * model choice; never serialised): Vivax (the simple vivax model), Common (the
 * shared P. falciparum model used by the Dummy/Empirical/Molineaux/Penny
 * infection models) and Descriptive (the original descriptive model). Kind
 * ::Mock is only reachable from unit tests (see unittest/WHMock.h; the
 * behaviour lives in a hook table on RunExtensions).
 *
 * Method definitions are split per kind over WHFalciparum.cpp (shared
 * P. falciparum immunity code), CommonWithinHost.cpp, DescriptiveWithinHost.cpp
 * and WHVivax.cpp.
 */
class WHInterface final {
public:
    /// The within-host model variant implemented by this instance.
    enum class Kind { Vivax, Common, Descriptive, Mock };

    /** Bind the host human's column row (Human ctor immediately after
     * creation; unit tests bind a mini-column themselves). */
    void bind_columns( HumanColumns& c, size_t r ) { cols = &c; row = r; }

    /// The run this model belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    /// @brief Static methods
    //@{
    /// Initialise static parameters (all model variants)
    static void init(RunContext& ctx, const OM::Parameters& parameters, const scnXml::Scenario& scenario);
    /// Initialise P. falciparum immunity statics (Common/Descriptive kinds).
    static void initFalciparum( RunContext& ctx, const OM::Parameters& parameters, const scnXml::Model& model );
    /// Initialise vivax statics (Vivax kind).
    static void initVivax( RunContext& ctx, const OM::Parameters& parameters, const scnXml::Model& model );
    /// Initialise Common-kind statics.
    static void initCommon(RunContext& ctx, const scnXml::Scenario& scenario);
    /// Initialise Descriptive-kind statics. Must run after monitoring is set up.
    static void initDescriptive(RunContext& ctx);
    
    /** Called when health system parameters are loaded (Vivax kind).
     * 
     * If no "LiverStageDrug" parameters are present, this is still called but with null pointer.
     */
    static void setHSParametersVivax( const scnXml::LiverStageDrug* );
    
    /** Configure a new treatment option, and return the code used to select
     * that option later. */
    static TreatmentId addTreatment( RunContext& ctx, const scnXml::TreatmentOption& desc );

    /// Create an instance using the appropriate model
    static unique_ptr<WHInterface> createWithinHostModel( RunContext& ctx, LocalRng& rng, double comorbidityFactor );
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHInterface( RunContext& ctx, LocalRng& rng, double comorbidityFactor, Kind kind );
    ~WHInterface();
    
    WHInterface(const WHInterface&) = delete;
    WHInterface& operator= (const WHInterface&) = delete;
    // Out-of-line (defaulted) so PathogenesisModel can stay incomplete here.
    WHInterface(WHInterface&&);
    WHInterface& operator=(WHInterface&&) = delete;    // ctx reference member

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }
    //@}

    /** Return the infectiousness of this human to biting mosquitoes.
     * This step is independent of parasite genetics.
     * 
     * @returns the probability of this human infecting a feeding mosquito.
     * 
     * Calculates the probability of transmitting an infection of a given
     * genotype to a mosquito and store the result in probTransGenotype[g]
     * for a given genotype g 
     * 
     * Calulates the probability for imported infections _i and
     * for local infections _l */
    double probTransmissionToMosquito(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l)const;

    /// @returns true if host has patent parasites
    bool summarize(Host::Human& human) const;

    /// Create a new infection within this human
    void importInfection(LocalRng& rng);

    /** Carry out the effects of some treatment option, optionally with intervention deployment.
     * 
     * This is equivalent to calling treatSimple, then deploying any included interventions.
     * 
     * This may be used by intervention deployment, thus should use time 
     * ctx.clock.nowOrTs1(). */
    void treatment( Host::Human& human, TreatmentId treatId );
    
    /// Conditionally gives Primaquine as a treatment. Reports as appropriate.
    void optionalPqTreatment( Host::Human& human );
    
    /** Treat a patient via the simple treatment model. Return true if any
     * blood-stage treatment is administered. Report any liver-stage treatments. */
    bool treatSimple( Host::Human& human, SimTime timeLiver, SimTime timeBlood );
    
    /** Give a patient a course of drugs, via the Pk/Pd model
     * 
     * Note: doses sizes are modified according to age via the dosage
     * table given at the time this function is called.
     *
     * @param schedule Index of a treatment schedule
     * @param dosages Index of a dosage table
     * @param age Age of human in years
     */
    void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d);

    /** Add new infections and update the parasite densities of existing
     * infections. Also update immune status.
     *
     * @param nNewInfs Number of inoculations this time-step
     * @param genotype_weights For use in selecting infection genotypes. See
     *  documentation of ctx.genotypes.sampleGenotype().
     * @param ageInYears Age of human
     * @param bsvFactor Parasite survival factor for blood-stage vaccines
     */
    void update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears);

    /** TODO: this should not need to be exposed. It is currently used by a
     * severe outcome (pDeath) model inside the EventScheduler "case
     * management" model, and case management diagnostics. */
    double getTotalDensity() const;
    
    /** Simulate use of a diagnostic test.
     *
     * Does not report for costing purposes.
     * 
     * Is used both during time step updates and during monitoring.
     * 
     * @returns true when the diagnostic is positive
     */
    bool diagnosticResult( LocalRng& rng, const Diagnostic& diagnostic ) const;
    
    /** Use the pathogenesis model to determine, based on infection status
     * and random draw, this person't morbidity.
     * 
     * This function is called after update() every time step.
     * 
     * @param human A reference to the human, used when reporting
     * @param ageYears Age of human host in years
     * @param isDoomed True if the human is already doomed to die, used for reporting
     */
    Pathogenesis::StatePair determineMorbidity( Host::Human& human, double ageYears, bool isDoomed );

    /// Special intervention: clears all immunity
    void clearImmunity();
    
    // TODO(monitoring): these shouldn't have to be exposed (perhaps use summarize to report the data):
    double getCumulative_h() const;
    double getCumulative_Y() const;

    InfectionOrigin getInfectionOrigin() const;

    /** The maximum number of infections a human can have: owned by
     * SlotAllocator (model/infections/SlotAllocator.h), the single clamp
     * owner; see its comment for the two RNG-visible modes. */

private:
    /// @brief Per-kind constructor/destructor fragments
    //@{
    /// Shared P. falciparum immunity state; draws from rng (stream order fixed).
    void ctorFalciparum( LocalRng& rng, double comorbidityFactor );
    void ctorVivax( LocalRng& rng, double comorbidityFactor );
    void ctorCommon( LocalRng& rng, double comorbidityFactor );
    void ctorDescriptive( LocalRng& rng, double comorbidityFactor );
    void dtorVivax();
    void dtorCommon();
    void dtorDescriptive();
    //@}

    /// @brief Per-kind implementations behind the public (dispatching) methods
    //@{
    double probTransmissionToMosquito_falciparum(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l)const;
    double probTransmissionToMosquito_vivax(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l)const;
    
    bool summarize_common( Host::Human& human )const;
    bool summarize_descriptive( Host::Human& human ) const;
    bool summarize_vivax(Host::Human& human) const;
    
    void importInfection_vivax(LocalRng& rng);
    void importInfection_common(LocalRng& rng);
    void importInfection_descriptive(LocalRng& rng);
    
    void treatment_falciparum( Host::Human& human, TreatmentId treatId );
    void optionalPqTreatment_vivax( Host::Human& human );
    
    bool treatSimple_falciparum( Host::Human& human, SimTime timeLiver, SimTime timeBlood );
    bool treatSimple_vivax( Host::Human& human, SimTime timeLiver, SimTime timeBlood );
    
    void treatPkPd_common(size_t schedule, size_t dosage, double age, double delay_d);
    void treatPkPd_descriptive(size_t schedule, size_t dosages, double age, double delay_d);
    void treatPkPd_vivax(size_t schedule, size_t dosages, double age, double delay_d);
    
    void update_vivax(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears);
    void update_common(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears);
    void update_descriptive(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears);
    
    double getTotalDensity_vivax() const;
    double getCumulative_h_vivax() const;
    double getCumulative_Y_vivax() const;
    
    bool diagnosticResult_falciparum( LocalRng& rng, const Diagnostic& diagnostic ) const;
    bool diagnosticResult_vivax( LocalRng& rng, const Diagnostic& diagnostic ) const;
    
    Pathogenesis::StatePair determineMorbidity_falciparum( Host::Human& human, double ageYears, bool isDoomed );
    Pathogenesis::StatePair determineMorbidity_vivax( Host::Human& human, double ageYears, bool );
    
    void clearImmunity_vivax();
    void clearImmunity_falciparum();
    //@}

    /// @brief Checkpointing chains (kind-dispatched; wire format frozen)
    //@{
    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);
    void checkpoint_falciparum (istream& stream);
    void checkpoint_falciparum (ostream& stream);
    void checkpoint_common (istream& stream);
    void checkpoint_common (ostream& stream);
    void checkpoint_descriptive (istream& stream);
    void checkpoint_descriptive (ostream& stream);
    void checkpoint_vivax (istream& stream);
    void checkpoint_vivax (ostream& stream);
    /// Descriptive kind: load one infection from a checkpoint.
    void loadInfection_descriptive(istream& stream);
    //@}

    /// @brief Shared P. falciparum immunity model
    //@{
    /** Clear infections of the appropriate stages (Common/Descriptive kinds).
     * 
     * @param stage Which stages of the infection are affected
     */
    void clearInfections( Treatments::Stages stage );
    
    /** Updates for the immunity model − assumes imm_h_ref() and imm_Y_ref()
     * have already been incremented.
     *
     * Applies decay of immunity against asexual blood stages, if present. */
    void updateImmuneStatus();

    /** @returns A multiplier describing the proportion of parasites surviving
     * immunity effects this time step.
     * 
     * Note that in the Descriptive model this multiplies log(density), but the
     * new density has no effect on future densities, wheras the Empirical model
     * multiplies the actual density (which then affects density on the following
     * time step). */
    double immunitySurvivalFactor (double ageInYears, double cumulativeExposureJ);
    
    static void setParams(double cumYStar, double cumHStar, double aM, double dM);   // for unit test only
    //@}
    
    /// Which model this instance is (never serialised).
    Kind m_kind;

    /// The host human's SoA columns/row (bind_columns; immunity state lives
    /// in the imm_* / cum_prim_inf columns).
    HumanColumns* cols = nullptr;
    size_t row = 0;
    
    /// Multiplicity of infection
    int numInfs = 0;

    /// @brief P. falciparum state (unused by the Vivax/Mock kinds)
    //@{
    //!innate ability to control parasite densities
    double _innateImmSurvFact = 0.0;

    /** Number of infections received since birth, cumulative parasite
     * density since birth and its one-step lag: SoA columns (imm_h/imm_Y/
     * imm_Y_lag; bound via bind_columns — accessor refs keep read/write
     * syntax). Old: imm_h_ref()/_Y/_Y_lag doubles. */
    double& imm_h_ref() const { return cols->imm_h[row]; }
    double& imm_Y_ref() const { return cols->imm_Y[row]; }
    double& imm_Y_lag_ref() const { return cols->imm_Y_lag[row]; }
    /// Primary blood-stage infections since birth (old cumPrimInf; column).
    uint32_t& cum_prim_inf_ref() const { return cols->cum_prim_inf[row]; }
    
    /// Total asexual blood stage density (sum of density of infections).
    double totalDensity = 0.0;
    /// Density from the point of view of diagnostics dependent on HRP2
    double hrp2Density = 0.0;
    
    /** Maximum parasite density of any infection during the previous interval.
     *
     * With 5-day time steps, this is not just the maximum density of any infection
     * at the end of the time step, but something designed to emulate the maximum
     * of 5 daily samples. */
    double timeStepMaxDensity = 0.0;
    
    /** Total asexual blood stage density over last 20 days (uses samples from
    * 10, 15 and 20 days ago).
    *
    * m_y_lag[ctx.clock.ts0().moduloSteps(y_lag_len)] corresponds to the density
    * from the previous time step (once updateInfection has been called). */
    std::vector<double> m_y_lag_i, m_y_lag_l;
    
    /// The PathogenesisModel introduces illness dependant on parasite density
    unique_ptr<Pathogenesis::PathogenesisModel> pathogenesisModel;
    
    /// End of step on which treatment expires = start of first step after expiry
    SimTime treatExpiryLiver = sim::zero(), treatExpiryBlood = sim::zero();

    /// Length of m_y_lag array. Wouldn't have to be dynamic if Global::interval was known at compile-time.
    /// set by initHumanParameters
    static int y_lag_len;
    //@}
    
    /// @brief Common-kind state (infections/opt_vaccine_genotype shared with Descriptive)
    //@{
    /// Multiplies the mean mass (for this age) as a heterogeneity factor.
    double hetMassMultiplier = 1.0;
    
    /// Encapsulates drug code for each human
    PkPd::LSTMModel pkpdModel;
    
    /** The list of all infections this human has (Common/Descriptive kinds).
     *
     * Since infection models and within host models are very much intertwined,
     * the idea is that each WithinHostModel has its own list of infections. */
    std::list<Infection*> infections;

    bool opt_vaccine_genotype = false;
    //@}
    
    /// @brief Vivax-kind state
    //@{
    list<VivaxBrood> vivaxInfections;
    
    /* Is flagged as never getting PQ: this is a heterogeneity factor. Example:
     * Set to zero if everyone can get PQ, 0.5 if females can't get PQ and
     * males aren't tested (i.e. all can get it) or (1+p)/2 where p is the
     * chance a male being tested and found to be G6PD deficient. */
    bool noPQ = false;
    
    Pathogenesis::State morbidity = Pathogenesis::NONE;
    
    // The number of primary blood stage infections to have started since this
    // human was born: SoA column, see cum_prim_inf_ref().
    
    // Either never or expiry time of treatment (Vivax kind)
    SimTime treatExpiryLiverVivax = sim::zero(), treatExpiryBloodVivax = sim::zero();

    // The probability of a clinical event from a primary infection
    double pEvent = numeric_limits<double>::quiet_NaN();
    // The probability of a clinical event from an early relapse
    double pFirstRelapseEvent = numeric_limits<double>::quiet_NaN();
	// The probability of a clinical event from a later relapse
	double pSecondRelapseEvent = numeric_limits<double>::quiet_NaN();
    
    // Used for reporting; updated by update()
    double pSevere = 0.0;
    //@}

    friend class ::UnittestUtil;
    friend class ::InfectionImmunitySuite;
};

}
}
#endif
