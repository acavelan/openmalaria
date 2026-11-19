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

#ifndef Hmod_Infection
#define Hmod_Infection

#include "Global.h"
#include "core/RunContext.h"
#include "Parameters.h"
#include "Host/WithinHost/Genotypes.h"
#include "util/random.h"

class UnittestUtil;
class MolineauxInfectionSuite;
class PennyInfectionSuite;

namespace OM { namespace WithinHost {

using util::LocalRng;

enum InfectionOrigin
{
    /** Imported infections are directly added to the human population dependent on the importation rate. **/
    Imported,

    /** Introduced infections are locally transmitted infections from mosquitoes who got infected from imported infections. **/
    Introduced,

    /** Indigenous infections are locally transmitted infections from mosquitoes who got infected from introduced or indigenous infections. **/
    Indigenous,
};

// Max duration of sample data for an infection in intervals. Simulated
// infections may last longer; in this case the last sample data is repeated
// until the infection terminates. (Descriptive)
const int32_t maxDurationTS = 83;
const int32_t numDurations = 84;    // Above plus one (include a category for 0)

// The maximum parasite density we allow per DescriptiveInfection. Higher values are set to maxDens.
// NOTE: should this be enforced before or after attenuation (BSV, immunity, IPT)? Was enforced in both places.
const double maxDens = 2000000;

/** A single within-host infection, of one of the modelled types.
 *
 * Was the Infection/CommonInfection hierarchy with Descriptive, Dummy,
 * Empirical, Molineaux and Penny subclasses; now one final class selected by
 * Kind (scenario-derived; Kind is never checkpointed — the checkpointing
 * factory function pointer carries the kind across a resume).
 *
 * Common models (Dummy/Empirical/Molineaux/Penny) use a 1-day time step, with
 * CommonWithinHost handling the conversion when the main simulation uses a
 * 5-day time step; note that for them ctx.clock.ts0(), ts1(), etc. may not be
 * accurate since it is only updated once per main time step. Descriptive uses
 * the 5-day time step directly via determineDensities().
 *
 * Member functions are defined per kind, in the former subclass TUs.
 *
 * Molineaux model: L. Molineaux et al, Parasitology 122 (2001) 379-391.
 * Penny model: Penny et al (2011), within-host dynamics with blood stage
 * vaccines. */
class Infection final {
public:
    enum class Kind { Descriptive, Dummy, Empirical, Molineaux, Penny };

    /// The run this infection belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    /// Pre-erythrocytic latent period (instantiated in WHFalciparum.cpp)
    inline static void init( SimTime latentP ){
        s_latentP = latentP;
    }

    ///@name Static (shared) init — call the one for the model in use
    //@{
    static void initDescriptive(RunContext& ctx, const OM::Parameters& parameters);  // 5-day time step check
    static void initDummy(RunContext& ctx);
    static void initEmpirical(RunContext& ctx);   // 1-day time step check
    static void initMolineaux(RunContext& ctx, const OM::Parameters& parameters);
    static void initPenny(RunContext& ctx);
    //@}

    /// Create a new infection. growthRateMultiplier is only used by Empirical.
    Infection(RunContext& ctx, LocalRng& rng, uint32_t genotype, int origin, Kind kind,
              double growthRateMultiplier = 1.0);
    /// For checkpointing (don't use for anything else)
    Infection(RunContext& ctx, istream& stream, Kind kind);
    ~Infection ();

    /** Update: calculate new density. Call this once per day (Dummy/Empirical/
     * Molineaux/Penny; never for Descriptive).
     *
     * @param survivalFactor Density multiplier to introduce drug & vaccine effects
     * @param now The simulation time. Use this instead of ctx.clock.ts1().
     * @param body_mass Body mass of host in kg
     * @returns True when the infection goes extinct. */
    bool update( LocalRng& rng, double survivalFactor, SimTime now, double body_mass );

    /** Returns true when age reaches the pre-determined duration (i.e. when
     * this infection terminates). Descriptive only; m_duration is never() for
     * the other kinds, which self-terminate in update() instead. */
    bool expired () {
        return ctx.clock.ts0() > m_startDate + m_duration;
    }

    /** Determines parasite density of an individual infection (Descriptive,
     * 5-day time step update).
     *
     * @param cumulativeh Cumulative number of infections
     * @param timeStepMaxDensity (In-out param) Used to return the maximum
     *  parasite density over a 5-day interval.
     * @param innateImmSurvFact Density multiplier for innate immunity.
     * @param bsvFactor Density multiplier for Blood-Stage Vaccine effect.
     */
    void determineDensities(
            LocalRng& rng,
            double cumulativeh,
            double &timeStepMaxDensity,
            double immSurvFact,
            double innateImmSurvFact,
            double bsvFactor);

    /** Decide on an infection duration and return it (Descriptive).
     *
     * Determines infection duration by sampling from the log normal
     * distribution using parameters for 53 patients from Georgia. Mean log
     * duration of an infection values from AJTM p.9 eq.5. Parameters for this
     * model are hard-coded.
     */
    static SimTime infectionDuration(RunContext& ctx, LocalRng& rng);

    /** Get the density of sequestered parasites (Penny). */
    inline double seqDensity(int ageDays){
        size_t todayV = mod_nn(ageDays, delta_V);
        return penny->seqDensities[todayV];
    }

    /** Return true if infection is blood stage.
     * 
     * Note: infections are considered to be liver stage for 5 days. This is
     * hard-coded since it is convenient in a 5-day time step model (and was one
     * of the reasons a 5-day time step was originally used) (TS).
     * 
     * The remainder of the "latentP" (pre-patent) period is blood-stage, where
     * blood-stage drugs do have an effect but parasites are not detectible.
     * 
     * Note 2: this gets called when deciding which infections to clear. If
     * clearing while updating infections (delayed treatment effect), infections
     * are liver-stage on the time step they start and blood-stage on the next
     * update, thus can be cleared the first time step they are considered
     * blood-stage. If clearing immediately (legacy health system and MDA
     * effect), clearance of blood stage infections can only happen after their
     * first update (though due to the latent period densities will still be
     * low). */
    inline bool bloodStage() const{
        return ctx.clock.latestTs0() - m_startDate > sim::fromDays(5);
    }
    
    /// Get the density of the infection as of the last update
    inline double getDensity()const{
        return m_density;
    }
    
    /// Get the cumulative parasite density
    inline double cumulativeExposureJ() {
        return m_cumulativeExposureJ;
    }
    
    /// Get whether the infection is HRP2-deficient
    bool isHrp2Deficient() const {
        return ctx.genotypes.getGenotypes()[m_genotype].hrp2_deficient;
    }
    
    /** Get the infection's genotype. */
    uint32_t genotype()const{ return m_genotype; }

    /** Get the infection's genotype. */
    int origin() const{ return m_origin; }

    /// Resets immunity properties specific to the infection (should only be
    /// called along with clearImmunity() on within-host model).
    inline void clearImmunity(){
        m_cumulativeExposureJ = 0.0;
    }

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }

    /// IC50^slope per drug type, if sampled (PkPd cache). Never serialised;
    /// lazily allocated (phase-5 perf fix — only the PkPd path touches it).
    std::map<size_t, double>& Kn() {
        if (!kn) kn = std::make_unique<std::map<size_t, double>>();
        return *kn;
    }

    // Molineaux constants:
    // v: number of variants per clone (one infection = one new clone); 50 in paper
    //TODO(optimisation): v should have a significant effect on performance, but reducing by,
    // say, half, may not have a big effect on the model. Evaluate.
    static const size_t v = 50;
    // taus: used for the variantTranscending and variantSpecific array, 4 Molineaux time steps = 8 days
    static const size_t taus = 4;

protected:
    void checkpoint (ostream& stream);

    /// Date of inoculation of infection (start of liver stage)
    /// This is the step of inoculation (ts0()).
    SimTime m_startDate = sim::never();
    
    /// Current density of the infection (Pc in the Molineaux/Penny papers)
    double m_density;
    
    /// Cumulative parasite density, since start of this infection
    double m_cumulativeExposureJ;

    /// Pre-erythrocytic latent period (instantiated in WHFalciparum.cpp)
    static SimTime s_latentP;

private:
    /// Genotype of infection (a code; see Genotypes class).
    uint32_t m_genotype;

    int m_origin;

    Kind m_kind;
    friend class ::UnittestUtil;   // reads s_latentP and internals
    friend class ::MolineauxInfectionSuite;
    friend class ::PennyInfectionSuite;

    // ——— ctor bodies and per-kind update/checkpoint, defined in the kind's TU ———
    void ctorBody_Empirical(LocalRng& rng, double growthRateMultiplier);
    void ctorBody_Molineaux(LocalRng& rng);
    void ctorBody_Penny(LocalRng& rng);
    bool updateDensity_Dummy(LocalRng& rng, double survivalFactor, SimTime bsAge, double body_mass);
    bool updateDensity_Empirical(LocalRng& rng, double survivalFactor, SimTime bsAge, double body_mass);
    bool updateDensity_Molineaux(LocalRng& rng, double survivalFactor, SimTime bsAge, double body_mass);
    bool updateDensity_Penny(LocalRng& rng, double survivalFactor, SimTime bsAge, double body_mass);
    void read_Descriptive(istream& stream);
    void read_Empirical(istream& stream);
    void read_Molineaux(istream& stream);
    void read_Penny(istream& stream);
    void checkpoint_Descriptive(ostream& stream);
    void checkpoint_Empirical(ostream& stream);
    void checkpoint_Molineaux(ostream& stream);
    void checkpoint_Penny(ostream& stream);

    // ——— Descriptive state ———
    // Arbitrary predetermined maximum duration of the infection
    SimTime m_duration = sim::never();
    bool notPrintedMDWarning = true;
    /* A triangular matrix: meanLogParasiteCount[i][j] is the
     * Mean Log Parasite Count for age i (in time steps) of an infection which
     * lasts j days. Indices with i>j are unused. */
    static double meanLogParasiteCount[numDurations][numDurations];
    static double sigma0sq;       ///< Sigma0^2 from AJTM p.9 eq. 13
    static double xNuStar;        ///< XNuStar in AJTM p.9 eq. 13

    // ——— Dummy state: none ———

    // ——— Empirical state: in EmpiricalData payload ———

    double getInflatedDensity(LocalRng& rng, double nonInflatedDensity);
    double sigma_noise(int ageDays);
    double sampleSubPatentValue(LocalRng& rng, double mu, double sigma, double upperBound);

    ///@brief Static variables
    ///Set by initEmpirical (and reset by independent parameterisation)
    //@{
    static const int _maximumDurationInDays=418;
    static double _maximumPermittedAmplificationPerCycle;
    static double _subPatentLimit;
    static double _alpha1;
    static double _alpha2;
    static double _alpha3;
    static double _mu1;
    static double _mu2;
    static double _mu3;
    static double _sigma0_res;
    static double _sigmat_res;
    static double _mu_beta1[_maximumDurationInDays];
    static double _sigma_beta1[_maximumDurationInDays];
    static double _mu_beta2[_maximumDurationInDays];
    static double _sigma_beta2[_maximumDurationInDays];
    static double _mu_beta3[_maximumDurationInDays];
    static double _sigma_beta3[_maximumDurationInDays];
    static double _inflationMean;
    static double _inflationVariance;
    static double _extinctionLevel;
    static double _overallMultiplier;
    //@}

    // ——— Molineaux: variant row type (payload below) ———
    struct Variant {
        // Initialise all variables to 0
        Variant ();

        /// Checkpointing
        void operator& (ostream& stream);
        void operator& (istream& stream);

        float Pi1, Pi2;   // Pi(t+1), Pi(t+2): variant's i density (PRBC/μl blood)
        float Si_summation;     // sum in eqn 6
        // index: we use ((bsAge/2) mod 4) for τ = t - δ_v respectively τ = t
        float lagged_Pi[taus];   // Pi(τ) for τ ∈ {t - δ_v, ..., t - 2}
    };

    // ——— Penny statics ———
    // delta_C := delay to clonal antibody response (days) (value 7.2038 round to 7)
    static const int delta_C=7;
    // delta_V := delay to variant specific antibody response in R_V^x (days) (value 6.3572 round to 6)
    static const int delta_V=6;

    static bool immune_threshold_gamma;
    static bool update_density_gamma;

    double getClonalSummation(int ageDays);
    double getVariantSpecificSummation(LocalRng& rng, int ageDays);

    // ——— Per-kind cold state, one heap payload per kind (phase-5 perf fix) ———
    // The profiled hot loop (Descriptive determineDensities/update_descriptive,
    // example@2000) only touches the common fields above; carrying every kind's
    // state inline cost ~4x the cache footprint of the Descriptive hot object
    // (RSS +12% after the 4b merge). Each payload is allocated in the ctor body
    // or the checkpoint reader for its kind only (no RNG draws involved), so
    // the m_rng stream and wire semantics are unchanged.
    struct EmpiricalData {
        double _laggedLogDensities[3];
        double _patentGrowthRateMultiplier;
    };
    struct MolineauxData {
        // m_startDate gives the age; m_density is Pc; m_cumulativeExposureJ is
        // the cumulative density (paper notes in the Molineaux TU)
        float mi[v];        // base multiplication factor per two-day cycle of variant i
        float Sm_summation; // sum in eqn 7
        // index: we use ((bsAge/2) mod 4) for τ = t - δ_v respectively τ = t
        float lagged_Pc[taus];       // Pc(τ) for τ ∈ {t - δ_v, ..., t - 2}
        /* Pc_star, Pm_star: two host-specific critical densities.
         * These two values depend on the first local maximum or the difference
         * between the last positive day and the first positive day. */
        float Pc_star, Pm_star;
        // variant-specific data; variants[i-1] corresponds to variant i in the paper
        vector<Variant> variants;
    };
    struct PennyData {
        // Circulating densities, 1 to delta_C days ago.
        // index (age mod delta_C) corresponds to delta_C days ago.
        double cirDensities[delta_C];
        // as above, but length delta_V
        double seqDensities[delta_V];
        // threshold_*: critical thresholds of the sigmoidal immune functions
        // (N innate, V variant specific, C clonal)
        double threshold_N;
        double threshold_V;
        double threshold_C;
        // tracked summations of densities with decay (sigmoidal functions)
        double variantSpecificSummation;
        double clonalSummation;
    };
    std::unique_ptr<EmpiricalData> empirical;      // set by Kind::Empirical
    std::unique_ptr<MolineauxData> molineaux;      // set by Kind::Molineaux
    std::unique_ptr<PennyData> penny;              // set by Kind::Penny
    std::unique_ptr<std::map<size_t, double>> kn;  // PkPd cache, see Kn()
};

/// Create a dummy infection (used directly by unit test).
Infection* createDummyInfection (RunContext& ctx, LocalRng& rng, uint32_t protID, int origin);

template <typename T>
InfectionOrigin get_infection_origin(const std::list<T*> &infections)
{
    if(infections.empty())
        return InfectionOrigin::Indigenous;

    int nImported = 0, nIntroduced = 0, nIndigenous = 0;
    for( auto inf = infections.begin(); inf != infections.end(); ++inf )
    {
        if((*inf)->origin() == InfectionOrigin::Indigenous) nIndigenous++;
        else if((*inf)->origin() == InfectionOrigin::Introduced) nIntroduced++;
        else nImported++;
    }

    /* The rules are:
    - Imported only if all infections are imported
    - Introduced if at least one Introduced
    - Indigenous otherwise (Imported + Indigenous or just Indigenous infections) */
    if(nIntroduced > 0)
        return InfectionOrigin::Introduced;
    else if(nIndigenous > 0)
        return InfectionOrigin::Indigenous;
    else
        return InfectionOrigin::Imported;
}

} }
#endif
