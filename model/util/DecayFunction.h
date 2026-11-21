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

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "Clock.h"
#include "util/sampler.h"
#include "util/UnitParse.h"

namespace scnXml {
class DecayFunction;
}

namespace OM {
class SimClock;
namespace util {

/** A decay function (several also suitible as survival functions): one
 * enum-tagged value type — behaviour is chosen by DATA (Kind switch in
 * compute()/sampleAgeOfDecay(), one case per formula), no class hierarchy.
 *
 * Config vs State: the XML layer builds one configured instance per component
 * parameter (config; hetFactor==0, never eval'd). hetSample() clones it into a
 * per-individual State carrying that individual's hetFactor. Operator/Emax
 * nest child configs; the factor applies to EVERY node of a sampled copy
 * (parity: old makeHetSample recursed hetSample(factor) into children).
 *****************************************************************************/
class DecayFunction {
  public:
    enum class Kind {
        Constant,
        Step,
        Linear,
        Exponential,
        Weibull,
        Hill,
        SmoothCompact,
        Plus,
        Minus,
        Divides,
        Multiplies,
        Emax
    };

    /// Plain-data construction; XML-driven builds go through makeObject.
    DecayFunction(Kind kind, bool increasing, double initialEfficacy, double CV)
        : kind(kind), increasing(increasing), initialEfficacy(initialEfficacy),
          het(Sampler::fromMeanCV(Sampler::Kind::Lognormal, 1.0, CV)) {}

    /** Return a new decay function, constructed from an XML element
     * (definition: model/xml/util.cpp — the only code reading the schema).
     *
     * @param elt XML element specifying which function to use and parameters
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    static unique_ptr<DecayFunction> makeObject(const SimClock& clock, const scnXml::DecayFunction& elt,
                                                const char* eltName);

    DecayFunction(const DecayFunction& other)
        : kind(other.kind), increasing(other.increasing), initialEfficacy(other.initialEfficacy),
          het(other.het ? other.het->clone() : nullptr), hetFactor(other.hetFactor), invL(other.invL), k(other.k),
          Emax(other.Emax), IC50(other.IC50), slope(other.slope), initialConcentration(other.initialConcentration),
          f1(other.f1 ? make_unique<DecayFunction>(*other.f1) : nullptr),
          f2(other.f2 ? make_unique<DecayFunction>(*other.f2) : nullptr), simClock(other.simClock) {}

    /** Sample a DecayFunction value (should be stored per individual).
     *
     * Note that a DecayFunction is needed to call eval() even if heterogeneity
     * is not wanted. If sigma = 0 then the random number stream will not be
     * touched. */
    unique_ptr<DecayFunction> hetSample(LocalRng& rng) const { return hetSample(het->sample(rng)); }

    /** Generate a DecayFunction value from an existing sample. */
    std::unique_ptr<DecayFunction> hetSample(NormalSample sample) const { return hetSample(het->sample(sample)); }

    /** Generate a DecayFunction value from an existing sample: per-individual
     * State = config clone with factor set recursively (Operator/Emax share it). */
    unique_ptr<DecayFunction> hetSample(double factor) const {
        unique_ptr<DecayFunction> sample = make_unique<DecayFunction>(*this);
        sample->setFactor(factor);
        return sample;
    }

    /** Say you have a population of objects which each have two states:
     * decayed and not decayed. If you want to use a DecayFunction to model
     * the proportion of objects which have decayed, you need to work out per
     * object the age of decay. This function does that.
     *
     * This is only valid where mu and sigma parameters are zero. hetFactor is
     * irrelevant here (every kind's draw is of the unscaled age).
     *
     * @returns Age at which an object should decay. */
    SimTime sampleAgeOfDecay(LocalRng& rng) const;

    inline double eval(double ageDays) const {
        if (increasing)
            return 1.0 - compute_with(ageDays, hetFactor) * initialEfficacy;
        else
            return compute_with(ageDays, hetFactor) * initialEfficacy;
    }

    /** eval() with an explicit heterogeneity factor: applies `factor` at
     * EVERY node (like a hetSample clone's recursive setFactor) without
     * cloning. Used by per-host IvRecords (old insecticideDecayHet clones). */
    inline double eval_with(double ageDays, double factor) const {
        if (increasing)
            return 1.0 - compute_with(ageDays, factor) * initialEfficacy;
        else
            return compute_with(ageDays, factor) * initialEfficacy;
    }

    /// Sample just the heterogeneity factor (what hetSample() would clone in).
    double sample_het_factor(LocalRng& rng) const { return het->sample(rng); }
    double sample_het_factor(NormalSample sample) const { return het->sample(sample); }

    /** Checkpoint a sampled copy, recreating it from this configured function.
     * Wire format frozen: bool present + double hetFactor. */

  private:
    /// The formula switch: one case per Kind (was one virtual per subclass).
    /// `factor` is the heterogeneity factor applied at every node.
    double compute_with(double ageDays, double factor) const;

    void setFactor(double factor) {
        hetFactor = factor;
        if (f1)
            f1->setFactor(factor);
        if (f2)
            f2->setFactor(factor);
    }

    Kind kind;
    bool increasing;
    double initialEfficacy;
    std::unique_ptr<util::Sampler> het;

    double hetFactor = 0.0;
    /// Precomputed scale: 1/L days (Step/Linear/Hill/SmoothCompact),
    /// log(2)/L (Exponential), pow(log 2, 1/k)/L (Weibull). Was per-subclass invL.
    double invL = 0.0;
    double k = 0.0; // Weibull/Hill/SmoothCompact shape
    // Emax parameters (0 elsewhere):
    double Emax = 0.0, IC50 = 0.0, slope = 0.0, initialConcentration = 0.0;
    unique_ptr<DecayFunction> f1, f2; // nested children (operator: both; emax: f1)
    /// Set by makeObject; needed by sampleAgeOfDecay (round-to-time-step).
    const SimClock* simClock = nullptr;
};

// ----  formula switch (one case per Kind; origins = old virtual subclasses)  ----

inline double DecayFunction::compute_with(double effectiveAge, double factor) const {
    switch (kind) {
    case Kind::Constant:
        // parity: xml/util.cpp:163 (old ConstantDecayFunction) — note: all decay
        // functions must return 0 when time > 0 and default-constructed.
        if (effectiveAge * factor == std::numeric_limits<double>::infinity())
            return 0.0;
        return 1.0;
    case Kind::Step: // parity: xml/util.cpp:187
        if (effectiveAge * invL * factor < 1.0)
            return 1.0;
        else
            return 0.0;
    case Kind::Linear: // parity: xml/util.cpp:215
        if (effectiveAge * invL * factor < 1.0)
            return 1.0 - effectiveAge * invL * factor;
        else
            return 0.0;
    case Kind::Exponential: // parity: xml/util.cpp:244
        return exp(-effectiveAge * invL * factor);
    case Kind::Weibull: { // parity: xml/util.cpp:270
        double p = -pow(effectiveAge * invL * factor, k);
        if (p < -700.0)
            return 0.0;
        else
            return exp(p);
    }
    case Kind::Hill: // parity: xml/util.cpp:300
        return 1.0 / (1.0 + pow(effectiveAge * invL * factor, k));
    case Kind::SmoothCompact: // parity: xml/util.cpp:326
        if (effectiveAge * invL * factor < 1.0)
            return exp(k - k / (1.0 - pow(effectiveAge * invL * factor, 2.0)));
        else
            return 0.0;
    case Kind::Plus:
    case Kind::Minus:
    case Kind::Divides:
    case Kind::Multiplies: {
        // parity: xml/util.cpp:364 — clamp order std::max(std::min(v, 1.0), 0.0) kept.
        double a = f1->eval_with(effectiveAge, factor), b = f2->eval_with(effectiveAge, factor);
        double v = kind == Kind::Plus ? a + b : kind == Kind::Minus ? a - b : kind == Kind::Divides ? a / b : a * b;
        return std::max(std::min(v, 1.0), 0.0);
    }
    case Kind::Emax: {
        // Emax model: E(C)=Emax*C^slope/(C^slope + IC50^slope), C=C0*f(t).
        // parity: xml/util.cpp:439
        const double concentrationPow = std::pow(initialConcentration * f1->eval_with(effectiveAge, factor), slope);
        const double efficacy = Emax * concentrationPow / (concentrationPow + IC50);
        return std::max(std::min(efficacy, 1.0), 0.0);
    }
    }
    return 0.0; // unreachable: every Kind has a case
}

inline SimTime DecayFunction::sampleAgeOfDecay(LocalRng& rng) const {
    switch (kind) {
    case Kind::Constant:
        return sim::future(); // "in the future" (sim::never() would mean the past)
    case Kind::Step:          // parity: xml/util.cpp:195
        return simClock->roundToTSFromDays(1.0 / invL);
    case Kind::Linear:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:222 (rounds to nearest; may decay instantly or at L)
        return simClock->roundToTSFromDays(rng.uniform_01() / invL);
    case Kind::Exponential:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:248
        return simClock->roundToTSFromDays(-log(rng.uniform_01()) / invL);
    case Kind::Weibull:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:278
        return simClock->roundToTSFromDays(pow(-log(rng.uniform_01()), 1.0 / k) / invL);
    case Kind::Hill:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:304
        return simClock->roundToTSFromDays(pow(1.0 / rng.uniform_01() - 1.0, 1.0 / k) / invL);
    case Kind::SmoothCompact:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:333
        return simClock->roundToTSFromDays(sqrt(1.0 - k / (k - log(rng.uniform_01()))) / invL);
    case Kind::Plus:
    case Kind::Minus:
    case Kind::Divides:
    case Kind::Multiplies:
        // parity: order — the old std::max(f1->…, f2->…) let gcc 16.2 evaluate its
        // arguments RIGHT-TO-LEFT, so f2 draws its uniform first; mirror explicitly.
        {
            const SimTime d2 = f2->sampleAgeOfDecay(rng);
            const SimTime d1 = f1->sampleAgeOfDecay(rng);
            return simClock->roundToTSFromDays(std::max(d1, d2));
        }
    case Kind::Emax: // parity: xml/util.cpp:447
        return simClock->roundToTSFromDays(f1->sampleAgeOfDecay(rng));
    }
    return sim::never(); // unreachable: every Kind has a case
}

}
}
