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

#include "Clock.h"
#include "RunContext.h"
#include "Parameters.h"
#include "Host/WithinHost/Genotypes.h"
#include "util/random.h"

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
 * the 5-day time step directly via determine_densities_descriptive().
 *
 * Per-kind behaviour lives in free functions in the former subclass TUs
 * (phase 7f); shared dispatch (ctor/update/checkpoint) is Infection.cpp.
 *
 * Molineaux model: L. Molineaux et al, Parasitology 122 (2001) 379-391.
 * Penny model: Penny et al (2011), within-host dynamics with blood stage
 * vaccines.
 *
 * Public struct (7.2b precedent: payloads public, friends dead): the per-kind
 * bodies are free functions in the kind's TU (phase 7f). */
struct Infection final {
    enum class Kind { Descriptive, Dummy, Empirical, Molineaux, Penny };

    /// The run this infection belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    /// Which model this infection runs (scenario-derived; never checkpointed —
    /// the checkpointing factory function pointer carries the kind across a
    /// resume).
    Kind kind;

    /// Pre-erythrocytic latent period (instantiated in WHFalciparum.cpp)
    inline static void init( SimTime latentP ){
        s_latentP = latentP;
    }

    /// Shared-field init only: kind-specific construction (which may draw from
    /// rng) is infection_new's / infection_new_checkpointed's job
    /// (Infection.cpp). Pass placeholder genotype/origin when reading a
    /// checkpoint: date/densities/genotype/origin are overwritten from the
    /// stream immediately after.
    Infection(RunContext& ctx, uint32_t genotype, int origin, Kind kind) :
        ctx(ctx),
        kind(kind),
        m_startDate(this->ctx.clock.nowOrTs0()),
        m_density(0.0),
        m_cumulativeExposureJ(0.0),
        m_genotype(genotype),
        m_origin(origin)
    {}

    /** Returns true when age reaches the pre-determined duration (i.e. when
     * this infection terminates). Descriptive only; the other kinds have no
     * Descriptive payload and self-terminate in infection_update() instead. */
    bool expired () {
        return ctx.clock.ts0() > m_startDate + descriptive.duration;
    }

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
        return ctx.genotypes.genotypes[m_genotype].hrp2_deficient;
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
        infection_checkpoint (*this, stream);
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

    /// Date of inoculation of infection (start of liver stage)
    /// This is the step of inoculation (ts0()).
    SimTime m_startDate = sim::never();
    
    /// Current density of the infection (Pc in the Molineaux/Penny papers)
    double m_density;
    
    /// Cumulative parasite density, since start of this infection
    double m_cumulativeExposureJ;

    /// Pre-erythrocytic latent period (instantiated in WHFalciparum.cpp)
    static SimTime s_latentP;

    /// Genotype of infection (a code; see Genotypes class).
    uint32_t m_genotype;

    int m_origin;

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

    // ——— Per-kind cold state, one heap payload per kind (phase-5 perf fix) ———
    // The profiled hot loop (determine_densities_descriptive,
    // example@2000) only touches the common fields above; carrying every kind's
    // state inline cost ~4x the cache footprint of the Descriptive hot object
    // (RSS +12% after the 4b merge). Each payload is allocated in the ctor body
    // or the checkpoint reader for its kind only (no RNG draws involved), so
    // the m_rng stream and wire semantics are unchanged.
    struct DescriptiveData {
        // Arbitrary predetermined maximum duration of the infection
        SimTime duration = sim::never();
        bool notPrintedMDWarning = true;
    };
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
    // Embedded by value (small, hot): the pointer payload cost a malloc per
    // infection plus a dependent load in every tick — bench +2.3%.
    DescriptiveData descriptive;  // set by Kind::Descriptive
    std::unique_ptr<EmpiricalData> empirical;      // set by Kind::Empirical
    std::unique_ptr<MolineauxData> molineaux;      // set by Kind::Molineaux
    std::unique_ptr<PennyData> penny;              // set by Kind::Penny
    std::unique_ptr<std::map<size_t, double>> kn;  // PkPd cache, see Kn()
};

// ——— Shared dispatch (Infection.cpp) ———
/// Create a new infection. growthRateMultiplier is only used by Empirical.
Infection* infection_new(RunContext& ctx, LocalRng& rng, uint32_t genotype, int origin,
                         Infection::Kind kind, double growthRateMultiplier = 1.0);
/// For checkpointing (don't use for anything else)
Infection* infection_new_checkpointed(RunContext& ctx, istream& stream, Infection::Kind kind);

/** Update: calculate new density. Call this once per day (Dummy/Empirical/
 * Molineaux/Penny; never for Descriptive).
 *
 * @param survivalFactor Density multiplier to introduce drug & vaccine effects
 * @param now The simulation time. Use this instead of ctx.clock.ts1().
 * @param body_mass Body mass of host in kg
 * @returns True when the infection goes extinct. */
bool infection_update(Infection& inf, LocalRng& rng, double survivalFactor, SimTime now, double body_mass);

/// Checkpoint write-out (what Infection::operator& forwards to)
void infection_checkpoint(Infection& inf, ostream& stream);

// ——— Descriptive free functions (DescriptiveInfection.cpp) ———
/// 5-day time step check; loads densities.csv matrix + parameters (per run)
void infection_init_descriptive(RunContext& ctx, const OM::Parameters& parameters);

/** Decide on an infection duration and return it.
 *
 * Determines infection duration by sampling from the log normal
 * distribution using parameters for 53 patients from Georgia. Mean log
 * duration of an infection values from AJTM p.9 eq.5. Parameters for this
 * model are hard-coded. */
SimTime infection_duration_descriptive(RunContext& ctx, LocalRng& rng);

/// New-infection hook for Kind::Descriptive (allocates the payload; draws
/// the duration — this is the kind's only RNG draw at construction time).
void infection_ctor_descriptive(Infection& inf, LocalRng& rng);

/** Determines parasite density of an individual infection (Descriptive,
 * 5-day time step update).
 *
 * @param cumulativeh Cumulative number of infections
 * @param timeStepMaxDensity (In-out param) Used to return the maximum
 *  parasite density over a 5-day interval.
 * @param innateImmSurvFact Density multiplier for innate immunity.
 * @param bsvFactor Density multiplier for Blood-Stage Vaccine effect.
 */
void determine_densities_descriptive(
        Infection& inf,
        LocalRng& rng,
        double cumulativeh,
        double &timeStepMaxDensity,
        double immSurvFact,
        double innateImmSurvFact,
        double bsvFactor);

void read_descriptive(Infection& inf, istream& stream);
void checkpoint_descriptive(Infection& inf, ostream& stream);

// ——— Dummy free functions (DummyInfection.cpp) ———
void infection_init_dummy(RunContext& ctx);
/// One update step; @returns true when the infection goes extinct.
bool update_density_dummy(Infection& inf, LocalRng& rng, double survivalFactor, SimTime bsAge, double body_mass);

// ——— Empirical free functions (EmpiricalInfection.cpp) ———
/// 1-day time step check; loads autoRegressionParameters.csv (per run)
void infection_init_empirical(RunContext& ctx);
/// New-infection hook for Kind::Empirical: prepatent-density draws (payload first).
void infection_ctor_empirical(Infection& inf, LocalRng& rng, double growthRateMultiplier);
bool update_density_empirical(Infection& inf, LocalRng& rng, double survivalFactor, SimTime bsAge, double body_mass);
void read_empirical(Infection& inf, istream& stream);
void checkpoint_empirical(Infection& inf, ostream& stream);

// ——— Molineaux free functions (MolineauxInfection.cpp) ———
void infection_init_molineaux(RunContext& ctx, const OM::Parameters& parameters);
/// New-infection hook: mi[] + Pc_star/Pm_star draws (payload first).
void infection_ctor_molineaux(Infection& inf, LocalRng& rng);
bool update_density_molineaux(Infection& inf, LocalRng& rng, double survival_factor, SimTime age_BS, double body_mass);
void read_molineaux(Infection& inf, istream& stream);
void checkpoint_molineaux(Infection& inf, ostream& stream);

// ——— Penny free functions (PennyInfection.cpp) ———
void infection_init_penny(RunContext& ctx);
/// New-infection hook: immune-threshold draws (payload first).
void infection_ctor_penny(Infection& inf, LocalRng& rng);
bool update_density_penny(Infection& inf, LocalRng& rng, double survivalFactor, SimTime bsAge, double body_mass);
void read_penny(Infection& inf, istream& stream);
void checkpoint_penny(Infection& inf, ostream& stream);

/// Create a dummy infection (used directly by unit test).
Infection* createDummyInfection (RunContext& ctx, LocalRng& rng, uint32_t protID, int origin);

template <typename T>
InfectionOrigin get_infection_origin(const std::list<T*> &infections)
{
    if(infections.empty())
        return InfectionOrigin::Indigenous;

    // only >0 tests below matter (all-Imported = neither counter set)
    int nIntroduced = 0, nIndigenous = 0;
    for( auto inf = infections.begin(); inf != infections.end(); ++inf )
    {
        if((*inf)->origin() == InfectionOrigin::Indigenous) nIndigenous++;
        else if((*inf)->origin() == InfectionOrigin::Introduced) nIntroduced++;
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
