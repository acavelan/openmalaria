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
#include "Parameters.h"
#include "Host/WithinHost/Genotypes.h"
#include "util/random.h"
#include <list>

namespace OpenMalaria {

enum InfectionOrigin {
    /** Imported infections are directly added to the human population dependent on the importation rate. **/
    Imported,

    /** Introduced infections are locally transmitted infections from mosquitoes who got infected from imported infections. **/
    Introduced,

    /** Indigenous infections are locally transmitted infections from mosquitoes who got infected from introduced or indigenous infections. **/
    Indigenous,
};

/** A single within-host infection: the base of the per-variant model-option
 * hierarchy below. The base holds ONLY the state every variant shares (start
 * date, density, cumulative exposure, genotype, origin, run context) and the
 * two per-variant operations (update, expired).
 *
 * One class per variant, each owning its own state and implementing its own
 * approach (true polymorphism; was the Kind-switched single class):
 *   DescriptiveInfection  DescriptiveInfection.cpp   5-day step, densities.csv grid
 *   DummyInfection        Infection.cpp              trivial exponential density
 *   EmpiricalInfection    EmpiricalInfection.cpp     1-day step, auto-regression
 *   MolineauxInfection    MolineauxInfection.cpp     1-day step, variant clones
 *   PennyInfection        PennyInfection.cpp         1-day step, circulating/sequestered
 * ModelOptions selects the variant once per run: each variant's init()
 * installs its one-line factory on RunContext::createInfection (Descriptive
 * is created directly by the Descriptive within-host model).
 *
 * Common variants (Dummy/Empirical/Molineaux/Penny) use a 1-day time step,
 * with CommonWithinHost handling the conversion when the main simulation uses
 * a 5-day time step; note that for them ctx.clock.ts0(), ts1(), etc. may not
 * be accurate since it is only updated once per main time step. Descriptive
 * uses the 5-day time step directly via determine_densities().
 *
 * Molineaux model: L. Molineaux et al, Parasitology 122 (2001) 379-391.
 * Penny model: Penny et al (2011), within-host dynamics with blood stage
 * vaccines.
 *
 * kind tag: set by each derived ctor; the checkpoint module (free functions,
 * model/checkpoint/humans.cpp) switches on it and static_casts to the derived
 * class — the tag is the proof, no dynamic_cast. All fields are public
 * (7.2b precedent: payloads public, friends dead). */
struct Infection {
    enum class Kind { Descriptive, Dummy, Empirical, Molineaux, Penny };

    /// The run this infection belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    /// Which variant this infection is (scenario-derived; never checkpointed —
    /// the checkpointing factory RunContext::createInfection carries the kind
    /// across a resume).
    const Kind kind;

    /// Pre-erythrocytic latent period (instantiated in WHFalciparum.cpp)
    inline static void init(SimTime latentP) { s_latentP = latentP; }

    /// Shared-field init only: NO draws here. Variant construction (which may
    /// draw from rng) lives entirely in the derived ctors, at exactly the
    /// stream positions the old infection_new()/infection_ctor_*() held.
    /// Pass placeholder genotype/origin when reading a checkpoint:
    /// date/densities/genotype/origin are overwritten from the stream
    /// immediately after.
    Infection(RunContext& ctx, uint32_t genotype, int origin, Kind kind)
        : ctx(ctx), kind(kind), m_startDate(this->ctx.clock.nowOrTs0()), m_density(0.0), m_cumulativeExposureJ(0.0),
          m_genotype(genotype), m_origin(origin) {}

    virtual ~Infection() {}

    /** Update: calculate new density. Call this once per day (Dummy/Empirical/
     * Molineaux/Penny; never for Descriptive — it self-terminates via
     * expired() and updates via determine_densities()). Each override applies
     * the shared latent-period gate first: nothing happens before
     * now - m_startDate exceeds s_latentP.
     *
     * @param survivalFactor Density multiplier to introduce drug & vaccine effects
     * @param now The simulation time. Use this instead of ctx.clock.ts1().
     * @param body_mass Body mass of host in kg
     * @returns True when the infection goes extinct. */
    virtual bool update(LocalRng& rng, double survivalFactor, SimTime now, double body_mass) {
        assert(false); // no variant updates through the bare base
        return false;
    }

    /** Returns true when age reaches the pre-determined duration (i.e. when
     * this infection terminates). Descriptive variant only (the others
     * self-terminate in update() instead). */
    virtual bool expired() {
        assert(false); // only DescriptiveInfection has a duration
        return false;
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
    inline bool bloodStage() const { return ctx.clock.latestTs0() - m_startDate > fromDays(5); }

    /// Get the density of the infection as of the last update
    inline double getDensity() const { return m_density; }

    /// Get the cumulative parasite density
    inline double cumulativeExposureJ() { return m_cumulativeExposureJ; }

    /// Get whether the infection is HRP2-deficient
    bool isHrp2Deficient() const { return ctx.genotypes.genotypes[m_genotype].hrp2_deficient; }

    /** Get the infection's genotype. */
    uint32_t genotype() const { return m_genotype; }

    /** Get the infection's genotype. */
    int origin() const { return m_origin; }

    /// Resets immunity properties specific to the infection (should only be
    /// called along with clearImmunity() on within-host model).
    inline void clearImmunity() { m_cumulativeExposureJ = 0.0; }

    /// IC50^slope per drug type, if sampled (PkPd cache). Never serialised;
    /// lazily allocated (phase-5 perf fix — only the PkPd path touches it).
    std::map<size_t, double>& Kn() {
        if (!kn)
            kn = std::make_unique<std::map<size_t, double>>();
        return *kn;
    }

    /// Date of inoculation of infection (start of liver stage)
    /// This is the step of inoculation (ts0()).
    SimTime m_startDate = never();

    /// Current density of the infection (Pc in the Molineaux/Penny papers)
    double m_density;

    /// Cumulative parasite density, since start of this infection
    double m_cumulativeExposureJ;

    /// Pre-erythrocytic latent period (instantiated in WHFalciparum.cpp)
    static SimTime s_latentP;

    /// Genotype of infection (a code; see Genotypes class).
    uint32_t m_genotype;

    int m_origin;

    std::unique_ptr<std::map<size_t, double>> kn; // PkPd cache, see Kn()
};

// ———  Descriptive variant (bodies: DescriptiveInfection.cpp)  ———

/** Descriptive model: density read from the meanLogParasiteCount grid
 * (densities.csv) per 5-day time step; duration drawn at construction. The
 * Descriptive within-host model drives it via determine_densities() — this
 * variant has no daily update(). */
struct DescriptiveInfection final : Infection {
    /// 5-day time step check; loads densities.csv matrix + parameters (per run)
    static void init(RunContext& ctx, const Parameters& parameters);

    /// Draws the duration — this variant's only RNG draw at construction time,
    /// at the slot the old infection_ctor_descriptive() held.
    DescriptiveInfection(RunContext& ctx, LocalRng& rng, uint32_t genotype, int origin);

    bool expired() override;

    /** Determines parasite density of this infection (5-day time step update).
     *
     * @param cumulativeh Cumulative number of infections
     * @param timeStepMaxDensity (In-out param) Used to return the maximum
     *  parasite density over a 5-day interval.
     * @param immSurvFact Density multiplier for acquired immunity.
     * @param innateImmSurvFact Density multiplier for innate immunity.
     * @param bsvFactor Density multiplier for Blood-Stage Vaccine effect.
     */
    void determine_densities(LocalRng& rng, double cumulativeh, double& timeStepMaxDensity, double immSurvFact,
                             double innateImmSurvFact, double bsvFactor);

    /** Decide on an infection duration and return it.
     *
     * Determines infection duration by sampling from the log normal
     * distribution using parameters for 53 patients from Georgia. Mean log
     * duration of an infection values from AJTM p.9 eq.5. Parameters for this
     * model are hard-coded. */
    static SimTime sample_duration(RunContext& ctx, LocalRng& rng);

    // Arbitrary predetermined maximum duration of the infection
    SimTime duration = never();
    bool notPrintedMDWarning = true;

    // Max duration of sample data for an infection in intervals. Simulated
    // infections may last longer; in this case the last sample data is repeated
    // until the infection terminates.
    static constexpr int32_t maxDurationTS = 83;
    static constexpr int32_t numDurations = 84; // Above plus one (include a category for 0)

    // The maximum parasite density we allow per DescriptiveInfection. Higher values are set to maxDens.
    // NOTE: should this be enforced before or after attenuation (BSV, immunity, IPT)? Was enforced in both places.
    static constexpr double maxDens = 2000000;
};

// ———  Dummy variant (bodies: Infection.cpp)  ———

/** Dummy model: density grows eightfold per day, modulo 20000, scaled by the
 * survival factor. The simplest within-host infection model. */
struct DummyInfection final : Infection {
    /// Installs the Dummy factory on ctx.createInfection (per run).
    static void init(RunContext& ctx);

    DummyInfection(RunContext& ctx, uint32_t genotype, int origin); // no draws

    bool update(LocalRng& rng, double survivalFactor, SimTime now, double body_mass) override;
};

// ———  Empirical variant (bodies: EmpiricalInfection.cpp)  ———

/** Empirical model: 1-day-step auto-regression on log densities
 * (autoRegressionParameters.csv), fitted to the Malaria Therapy data. */
struct EmpiricalInfection final : Infection {
    /// 1-day time step check; loads autoRegressionParameters.csv (per run)
    static void init(RunContext& ctx);

    /// Prepatent-density draws at the old infection_ctor_empirical() slots;
    /// growthRateMultiplier is only used by this variant.
    EmpiricalInfection(RunContext& ctx, LocalRng& rng, uint32_t genotype, int origin,
                       double growthRateMultiplier = 1.0);

    bool update(LocalRng& rng, double survivalFactor, SimTime now, double body_mass) override;

    double _laggedLogDensities[3];
    double _patentGrowthRateMultiplier;
};

// ———  Molineaux variant (bodies: MolineauxInfection.cpp)  ———

/** Molineaux model (L. Molineaux et al, Parasitology 122 (2001) 379-391):
 * 1-day step with 50 antigenic variants per clone, two-day density cycles
 * under innate/transcending/variant-specific immune functions. */
struct MolineauxInfection final : Infection {
    static void init(RunContext& ctx, const Parameters& parameters);

    /// mi[] + Pc_star/Pm_star draws at the old infection_ctor_molineaux() slots.
    MolineauxInfection(RunContext& ctx, LocalRng& rng, uint32_t genotype, int origin);

    bool update(LocalRng& rng, double survivalFactor, SimTime now, double body_mass) override;

    // Molineaux constants:
    // v: number of variants per clone (one infection = one new clone); 50 in paper
    //TODO(optimisation): v should have a significant effect on performance, but reducing by,
    // say, half, may not have a big effect on the model. Evaluate.
    static constexpr size_t v = 50;
    // taus: used for the variantTranscending and variantSpecific array, 4 Molineaux time steps = 8 days
    static constexpr size_t taus = 4;

    // variant row type
    struct Variant {
        // Initialise all variables to 0
        Variant();

        float Pi1, Pi2;     // Pi(t+1), Pi(t+2): variant's i density (PRBC/μl blood)
        float Si_summation; // sum in eqn 6
        // index: we use ((bsAge/2) mod 4) for τ = t - δ_v respectively τ = t
        float lagged_Pi[taus]; // Pi(τ) for τ ∈ {t - δ_v, ..., t - 2}
    };

    // m_startDate gives the age; m_density is Pc; m_cumulativeExposureJ is
    // the cumulative density (paper notes in the Molineaux TU)
    float mi[v];        // base multiplication factor per two-day cycle of variant i
    float Sm_summation; // sum in eqn 7
    // index: we use ((bsAge/2) mod 4) for τ = t - δ_v respectively τ = t
    float lagged_Pc[taus]; // Pc(τ) for τ ∈ {t - δ_v, ..., t - 2}
    /* Pc_star, Pm_star: two host-specific critical densities.
     * These two values depend on the first local maximum or the difference
     * between the last positive day and the first positive day. */
    float Pc_star, Pm_star;
    // variant-specific data; variants[i-1] corresponds to variant i in the paper
    vector<Variant> variants;
};

// ———  Penny variant (bodies: PennyInfection.cpp)  ———

/** Penny model (Penny et al (2011), within-host dynamics with blood stage
 * vaccines): 1-day step, circulating/sequestered densities with innate,
 * clonal and variant-specific sigmoidal immune responses. */
struct PennyInfection final : Infection {
    static void init(RunContext& ctx);

    /// Immune-threshold draws at the old infection_ctor_penny() slots.
    PennyInfection(RunContext& ctx, LocalRng& rng, uint32_t genotype, int origin);

    bool update(LocalRng& rng, double survivalFactor, SimTime now, double body_mass) override;

    /** Get the density of sequestered parasites. */
    inline double seqDensity(int ageDays) const {
        size_t todayV = mod_nn(ageDays, delta_V);
        return seqDensities[todayV];
    }

    // delta_C := delay to clonal antibody response (days) (value 7.2038 round to 7)
    static constexpr int delta_C = 7;
    // delta_V := delay to variant specific antibody response in R_V^x (days) (value 6.3572 round to 6)
    static constexpr int delta_V = 6;

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

/// Create a dummy infection (used directly by unit tests).
Infection* createDummyInfection(RunContext& ctx, LocalRng& rng, uint32_t protID, int origin);

template <typename T> InfectionOrigin get_infection_origin(const std::list<T*>& infections) {
    if (infections.empty())
        return InfectionOrigin::Indigenous;

    // only >0 tests below matter (all-Imported = neither counter set)
    int nIntroduced = 0, nIndigenous = 0;
    for (auto inf = infections.begin(); inf != infections.end(); ++inf) {
        if ((*inf)->origin() == InfectionOrigin::Indigenous)
            nIndigenous++;
        else if ((*inf)->origin() == InfectionOrigin::Introduced)
            nIntroduced++;
    }

    /* The rules are:
    - Imported only if all infections are imported
    - Introduced if at least one Introduced
    - Indigenous otherwise (Imported + Indigenous or just Indigenous infections) */
    if (nIntroduced > 0)
        return InfectionOrigin::Introduced;
    else if (nIndigenous > 0)
        return InfectionOrigin::Indigenous;
    else
        return InfectionOrigin::Imported;
}

}
