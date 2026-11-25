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

#include "Clock.h"
#include "RunContext.h"
#include "Host/HostColumns.h"
#include "util/random.h"
#include "util/sampler.h"
#include "util/errors.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Host/WithinHost/Infection/Infection.h"

#include "Parameters.h"

#include "PkPd/LSTMModel.h"

#include <algorithm>
#include <list>
#include <memory>

using namespace std;

namespace scnXml {
class Scenario;
class Model;
class LiverStageDrug;
}
class UnittestUtil;
class InfectionImmunitySuite;

namespace OpenMalaria {
class Human;
class PathogenesisModel;

/** Sole owner of the per-host infection-slot cap (plan §5: "MAX_INFECTIONS=21
 * lives in one place"). The only real reason for the limit is to prevent
 * bad input from making the number of infections balloon stupidly.
 *
 * The cap has TWO RNG-visible modes which must never be conflated
 * (census fact 2, latent-bugs.md #4; pinned by unittest/SlotAllocatorSuite.h):
 *  A. import mode  — try_import(): at the cap the caller skips the entire
 *     create, including the genotype draw (0 draws consumed).
 *     // parity: CommonWithinHost.cpp:111, DescriptiveWithinHost.cpp:81
 *  B. update mode  — clamp_new(): the caller's Poisson draw was already
 *     consumed, so counts are clamped draw-then-drop; _l clamps first, then
 *     _i against what remains.
 *     // parity: CommonWithinHost.cpp:131-132, DescriptiveWithinHost.cpp:104-105
 */
struct SlotAllocator {
    /// Exact constraint is: _MOI <= MAX_INFECTIONS.
    static constexpr int MAX_INFECTIONS = 21;

    /// Free slots for a host with num_infs infections (0 at the cap).
    static constexpr int free_slots(int num_infs) { return MAX_INFECTIONS - num_infs; }

    /// Mode A: true iff an import may proceed (genotype draw included).
    static constexpr bool try_import(int num_infs) { return num_infs < MAX_INFECTIONS; }

    /// Mode B: clamp to-create counts in place (_l first — order is load-bearing).
    static inline void clamp_new(int& n_i, int& n_l, int num_infs) {
        n_l = std::min(n_l, free_slots(num_infs));
        n_i = std::min(n_i, free_slots(num_infs) - n_l);
    }

    /// Clamp a stand-alone drawn count (post-draw; no rng effect).
    /// parity: InfectionIncidenceModel.cpp:257-258
    static constexpr int clamp_count(int want) { return want > MAX_INFECTIONS ? MAX_INFECTIONS : want; }
};

class WHInterface;

/// True if any by-genotype reporting is enabled (read-only)
extern bool reportInfectionsByGenotype;

struct HypnozoiteReleaseDistribution {
    /// The run this sampler serves (clock accessed as ctx.clock).
    RunContext& ctx;

    HypnozoiteReleaseDistribution(RunContext& ctx, std::unique_ptr<Sampler> sampler, double latentRelapse)
        : ctx(ctx), sampler(std::move(sampler)), latentRelapse(latentRelapse) {}

    /// Sample the time until next release
    SimTime sampleReleaseDelay(LocalRng& rng) const {
        double liverStageMaximumDays = 16.0 * 30.0;         // maximum of about 16 months in liver stage
        double delay = numeric_limits<double>::quiet_NaN(); // in days
        int count = 0;
        int maxcount = 1e3;

        do {
            delay = sampler->sample(rng);
            count += 1;

            if (count >= maxcount) {
                throw xml_scenario_error(
                    "<vivax><hypnozoiteRelease>  [random delay calculation causes probably an indefinite loop]:\n The hypnozoite release distribution seems off, sigma of secondRelease could be too high. We except the hypnozoite to reside a maximum of 16 months in the liver stage. Sigma choose well, dear padawan.");
            }
        } while (delay > liverStageMaximumDays || delay < 0.0);

        assert(delay >= 0 && delay < liverStageMaximumDays);
        return ctx.clock.roundToTSFromDays(delay + latentRelapse);
    }

  private:
    std::unique_ptr<Sampler> sampler;
    double latentRelapse; // days
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
class VivaxBrood {
  public:
    /** Create.
     * 
     * @param host      The host creating this. Not really needed, except to
     * prevent VivaxBrood being default-constructed by its container.
     */
    /// The run this brood belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    VivaxBrood(RunContext& ctx, LocalRng& rng, int origin, WHInterface* host);
    ~VivaxBrood();
    struct UpdResult {
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

    // Whether a clinical / relapse event has been triggered by this brood.
    // Written by update_vivax; no in-model reader — kept for the checkpoint wire.
    bool hadEvent;
    bool hadRelapse;

    /** Equivalent to a blood stage existing. We do not model incidence of
     * gametocytes independently, thus this also tests existance of
     * gametocytes. */
    inline bool isPatent() const { return bloodStageClearDate > ctx.clock.latestTs0(); }

    /** Fully clear blood stage parasites. */
    void treatmentBS();

    /** Fully clear liver stage parasites. */
    void treatmentLS();

  public: // run state; wire in checkpoint/humans.cpp
    // list of times at which the merozoite and hypnozoites release, ordered by
    // time of release, soonest last (i.e. last element is next one to release)
    vector<SimTime> releaseDates;

    // Either never() (no blood stage) or the start of the time step on
    // which the blood stage will clear.
    SimTime bloodStageClearDate = never();

    // Whether the primary blood stage infection has started
    bool primaryHasStarted;

    // Whether the relapse blood stage infection has started
    bool relapseHasStarted;

    //Whether the later category of relapses has started
    bool relapsebHasStarted;

    // Infection origin
    int origin;
};

/**
 * The within-host models. These models encapsulate the infections and related
 * immunity factors of a single human, starting with infection (i.e. assuming
 * successful innoculation), including some drug action code, and outputting
 * parasite densities.
 *
 * True polymorphism (was the Kind-switched single class): one class per
 * model-option variant, selected once per run by ModelOptions
 * (WHInterface::createWithinHostModel below; the kind tag is never
 * serialised):
 *   WHVivax        WHVivax.cpp          the simple vivax model (broods)
 *   WHCommon       CommonWithinHost.cpp  shared P. falciparum model used by
 *                                        the Dummy/Empirical/Molineaux/Penny
 *                                        infection variants
 *   WHDescriptive  DescriptiveWithinHost.cpp  the original descriptive model
 *   WHMock         unittest/WHMock.cpp  unit-test backend owning its state
 * WHCommon and WHDescriptive share the P. falciparum immunity machinery in
 * the intermediate WHFalciparum layer (WHFalciparum.cpp).
 *
 * The base holds the state every variant carries (column binding, infection
 * count, and the immunity/density block that the checkpoint wire writes for
 * every kind); the per-variant approaches are virtual overrides.
 */
class WHInterface {
  public:
    /// The within-host model variant implemented by this instance.
    enum class Kind { Vivax, Common, Descriptive, Mock };

    /** Bind the host human's column row (Human ctor immediately after
     * creation; unit tests bind a mini-column themselves). */
    void bind_columns(HumanColumns& c, size_t r) {
        cols = &c;
        row = r;
    }

    /// The run this model belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    /// @brief Static methods
    //@{
    /// Initialise static parameters (all model variants; selects the variant
    /// via ModelOptions, xml/scenario.cpp).
    static void init(RunContext& ctx, const Parameters& parameters, const scnXml::Scenario& scenario);

    /// Create an instance using the appropriate model (one line per variant).
    static unique_ptr<WHInterface> createWithinHostModel(RunContext& ctx, LocalRng& rng, double comorbidityFactor);
    //@}

    /// Base ctor: no state beyond the tag, no draws — every variant's draws
    /// live in its own ctor at the old ctor* fragment stream positions.
    WHInterface(RunContext& ctx, Kind kind) : ctx(ctx), kind(kind) {}

    virtual ~WHInterface() {}

    WHInterface(const WHInterface&) = delete;
    WHInterface& operator=(const WHInterface&) = delete;
    WHInterface& operator=(WHInterface&&) = delete; // ctx reference member

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
    virtual double probTransmissionToMosquito(vector<double>& probTransGenotype_i,
                                              vector<double>& probTransGenotype_l) const = 0;

    /// @returns true if host has patent parasites
    virtual bool summarize(Human& human) const = 0;

    /// Create a new infection within this human
    virtual void importInfection(LocalRng& rng) = 0;

    /** Carry out the effects of some treatment option, optionally with intervention deployment.
     * 
     * This is equivalent to calling treatSimple, then deploying any included interventions.
     * 
     * This may be used by intervention deployment, thus should use time 
     * ctx.clock.nowOrTs1().
     *
     * Shared by all production variants (treatSimple dispatches virtually;
     * the Vivax variant uses the same code here; the unittest Mock
     * overrides). */
    virtual void treatment(Human& human, TreatmentId treatId);

    /// Conditionally gives Primaquine as a treatment. Reports as appropriate.
    /// (No PQ treatment for falciparum in current models; Vivax overrides.)
    virtual void optionalPqTreatment(Human& human) {}

    /** Treat a patient via the simple treatment model. Return true if any
     * blood-stage treatment is administered. Report any liver-stage treatments. */
    virtual bool treatSimple(Human& human, SimTime timeLiver, SimTime timeBlood) = 0;

    /** Give a patient a course of drugs, via the Pk/Pd model
     * 
     * Note: doses sizes are modified according to age via the dosage
     * table given at the time this function is called.
     *
     * @param schedule Index of a treatment schedule
     * @param dosages Index of a dosage table
     * @param age Age of human in years
     */
    virtual void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d) = 0;

    /** Add new infections and update the parasite densities of existing
     * infections. Also update immune status.
     *
     * @param nNewInfs Number of inoculations this time-step
     * @param genotype_weights For use in selecting infection genotypes. See
     *  documentation of ctx.genotypes.sampleGenotype().
     * @param ageInYears Age of human
     * @param bsvFactor Parasite survival factor for blood-stage vaccines
     */
    virtual void update(Human& human, LocalRng& rng, int& nNewInfs_i, int& nNewInfs_l,
                        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears) = 0;

    /** TODO: this should not be exposed. It is currently used by a
     * severe outcome (pDeath) model inside the EventScheduler "case
     * management" model, and case management diagnostics.
     * (Vivax variant throws.) */
    virtual double getTotalDensity() const { return totalDensity; }

    /** Simulate use of a diagnostic test.
     *
     * Does not report for costing purposes.
     * 
     * Is used both during time step updates and during monitoring.
     * 
     * @returns true when the diagnostic is positive
     */
    virtual bool diagnosticResult(LocalRng& rng, const Diagnostic& diagnostic) const = 0;

    /** Use the pathogenesis model to determine, based on infection status
     * and random draw, this person't morbidity.
     * 
     * This function is called after update() every time step.
     * 
     * @param human A reference to the human, used when reporting
     * @param ageYears Age of human host in years
     * @param isDoomed True if the human is already doomed to die, used for reporting
     */
    virtual StatePair determineMorbidity(Human& human, double ageYears, bool isDoomed) = 0;

    /// Special intervention: clears all immunity.
    virtual void clearImmunity() = 0;

    // TODO(monitoring): these shouldn't have to be exposed (perhaps use summarize to report the data):
    virtual double getCumulative_h() const { return imm_h_ref(); }
    virtual double getCumulative_Y() const { return imm_Y_ref(); }

    /// Vivax/Mock have no per-infection origins (Indigenous); Falciparum
    /// layer derives from its infection list.
    virtual InfectionOrigin getInfectionOrigin() const { return InfectionOrigin::Indigenous; }

    /// Which model this instance is (never serialised; the checkpoint module
    /// switches on it and static_casts to the derived class — the tag is the
    /// proof, no dynamic_cast needed).
    const Kind kind;

  public: // simple data + fragments (test suites and WH models touch these directly)
    /// The host human's SoA columns/row (bind_columns; immunity state lives
    /// in the imm_* / cum_prim_inf columns).
    HumanColumns* cols = nullptr;
    size_t row = 0;

    /// Multiplicity of infection
    int numInfs = 0;

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
    unique_ptr<PathogenesisModel> pathogenesisModel;

    /// End of step on which treatment expires = start of first step after expiry
    SimTime treatExpiryLiver = zero(), treatExpiryBlood = zero();

    /// Length of m_y_lag array. Wouldn't have to be dynamic if Global::interval was known at compile-time.
    /// set by initHumanParameters
    static int y_lag_len;
};

/**
 * Shared P. falciparum immunity machinery, behind the Common and Descriptive
 * variants (bodies: WHFalciparum.cpp; infection-list code shared in
 * CommonWithinHost.cpp / DescriptiveWithinHost.cpp). Ctor draws the innate
 * immunity sample exactly where the old ctorFalciparum() did — always before
 * the derived variant's own ctor body.
 */
class WHFalciparum : public WHInterface {
  public:
    /// Initialise P. falciparum immunity statics (Common/Descriptive variants).
    static void init(RunContext& ctx, const Parameters& parameters, const scnXml::Model& model);
    static void setParams(double cumYStar, double cumHStar, double aM, double dM); // for unit test only

    /// Shared P. falciparum immunity ctor fragment (pathogenesis model +
    /// innate draw + y_lag buffers); draws stream position = old
    /// ctorFalciparum(). Kind comes from the concrete variant.
    WHFalciparum(RunContext& ctx, LocalRng& rng, double comorbidityFactor, Kind kind);
    ~WHFalciparum() override;

    double probTransmissionToMosquito(vector<double>& probTransGenotype_i,
                                      vector<double>& probTransGenotype_l) const override;
    bool diagnosticResult(LocalRng& rng, const Diagnostic& diagnostic) const override;
    bool treatSimple(Human& human, SimTime timeLiver, SimTime timeBlood) override;
    StatePair determineMorbidity(Human& human, double ageYears, bool isDoomed) override;
    void clearImmunity() override;
    InfectionOrigin getInfectionOrigin() const override;

    /** Clear infections of the appropriate stages.
     * 
     * @param stage Which stages of the infection are affected
     */
    void clearInfections(Treatments::Stages stage);

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
    double immunitySurvivalFactor(double ageInYears, double cumulativeExposureJ);

    /** Host-tick-constant part of immunitySurvivalFactor: dH and dA depend on
     * host immunity accumulators and age only, both frozen during an infection
     * update loop (imm_h_ref/imm_Y_ref are mutated outside it). Computing once
     * per host tick saves one exp+div per infection; the per-infection remainder
     * below reproduces the identical IEEE operand sequence ((dY*dH)*dA). */
    struct ImmSurvTick {
        bool active;       ///< imm_h > 1.0 (else dH = dY = 1.0)
        double immY, invY; ///< dY inputs, per infection
        double dH, dA;     ///< hoisted host constants
    };
    /** Host-side guard (ageInYears/imm_h/imm_Y NaN) throws here — same message;
     * a NaN abort ends the run either way. */
    ImmSurvTick immSurvivalTick(double ageInYears);
    /** Per-infection remainder; guards cumulativeExposureJ like the two-double
     * overload above. */
    double immunitySurvivalFactor(const ImmSurvTick& tick, double cumulativeExposureJ);

    /// Reset and refill the m_y_lag_{i,l} infectiousness cache for the current
    /// tick (shared tail of update_common and update_descriptive).
    void cache_y_lag();

  public: // test suites and the checkpoint module touch these directly
    /// Multiplies the mean mass (for this age) as a heterogeneity factor.
    double hetMassMultiplier = 1.0;

    /// Encapsulates drug code for each human
    LSTMModel pkpdModel;

    /** The list of all infections this human has.
     *
     * Since infection models and within host models are very much intertwined,
     * the idea is that each WithinHostModel has its own list of infections. */
    std::list<Infection*> infections;

    bool opt_vaccine_genotype = false;
};

/**
 * Common P. falciparum within-host variant: drives the Dummy/Empirical/
 * Molineaux/Penny infection variants on a daily step (infection selection
 * lives in the Infection hierarchy). Bodies: CommonWithinHost.cpp.
 */
class WHCommon final : public WHFalciparum {
  public:
    /// Initialise Common-variant statics.
    static void init(RunContext& ctx, const scnXml::Scenario& scenario);

    WHCommon(RunContext& ctx, LocalRng& rng, double comorbidityFactor);

    bool summarize(Human& human) const override;
    void importInfection(LocalRng& rng) override;
    void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d) override;
    void update(Human& human, LocalRng& rng, int& nNewInfs_i, int& nNewInfs_l, vector<double>& genotype_weights_i,
                vector<double>& genotype_weights_l, double ageInYears) override;
};

/**
 * Descriptive within-host variant: 5-day step, DescriptiveInfection densities.
 * Bodies: DescriptiveWithinHost.cpp.
 */
class WHDescriptive final : public WHFalciparum {
  public:
    /// Initialise Descriptive-variant statics. Must run after monitoring is set up.
    static void init(RunContext& ctx);

    WHDescriptive(RunContext& ctx, LocalRng& rng, double comorbidityFactor);

    bool summarize(Human& human) const override;
    void importInfection(LocalRng& rng) override;
    void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d) override;
    void update(Human& human, LocalRng& rng, int& nNewInfs_i, int& nNewInfs_l, vector<double>& genotype_weights_i,
                vector<double>& genotype_weights_l, double ageInYears) override;
};

/**
 * Vivax within-host variant: hypnozoite broods, no P. falciparum immunity, no
 * PK/PD. Bodies: WHVivax.cpp.
 */
class WHVivax final : public WHInterface {
  public:
    /// Initialise vivax statics.
    static void init(RunContext& ctx, const Parameters& parameters, const scnXml::Model& model);

    /** Called when health system parameters are loaded.
     * 
     * If no "LiverStageDrug" parameters are present, this is still called but with null pointer.
     */
    static void setHSParameters(const scnXml::LiverStageDrug*);

    WHVivax(RunContext& ctx, LocalRng& rng, double comorbidityFactor);
    ~WHVivax() override;

    double probTransmissionToMosquito(vector<double>& probTransGenotype_i,
                                      vector<double>& probTransGenotype_l) const override;
    bool summarize(Human& human) const override;
    void importInfection(LocalRng& rng) override;
    void optionalPqTreatment(Human& human) override;
    bool treatSimple(Human& human, SimTime timeLiver, SimTime timeBlood) override;
    void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d) override;
    void update(Human& human, LocalRng& rng, int& nNewInfs_i, int& nNewInfs_l, vector<double>& genotype_weights_i,
                vector<double>& genotype_weights_l, double ageInYears) override;
    double getTotalDensity() const override;
    bool diagnosticResult(LocalRng& rng, const Diagnostic& diagnostic) const override;
    StatePair determineMorbidity(Human& human, double ageYears, bool) override;
    void clearImmunity() override;
    double getCumulative_h() const override;
    double getCumulative_Y() const override;

  public: // run state; wire in checkpoint/humans.cpp
    list<VivaxBrood> vivaxInfections;

    /* Is flagged as never getting PQ: this is a heterogeneity factor. Example:
     * Set to zero if everyone can get PQ, 0.5 if females can't get PQ and
     * males aren't tested (i.e. all can get it) or (1+p)/2 where p is the
     * chance a male being tested and found to be G6PD deficient. */
    bool noPQ = false;

    State morbidity = NONE;

    // The number of primary blood stage infections to have started since this
    // human was born: SoA column, see cum_prim_inf_ref().

    // Either never or expiry time of treatment
    SimTime treatExpiryLiverVivax = zero(), treatExpiryBloodVivax = zero();

    // The probability of a clinical event from a primary infection
    double pEvent = numeric_limits<double>::quiet_NaN();
    // The probability of a clinical event from an early relapse
    double pFirstRelapseEvent = numeric_limits<double>::quiet_NaN();
    // The probability of a clinical event from a later relapse
    double pSecondRelapseEvent = numeric_limits<double>::quiet_NaN();

    // Used for reporting; updated by update()
    double pSevere = 0.0;
};

}
