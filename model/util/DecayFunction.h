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

#ifndef Hmod_DecayFunction
#define Hmod_DecayFunction

#include <limits>
#include <memory>

#include "Global.h"
#include "util/sampler.h"
#include "util/UnitParse.h"

namespace scnXml { class DecayFunction; }

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
class DecayFunction
{
public:
    enum class Kind { Constant, Step, Linear, Exponential, Weibull, Hill,
                      SmoothCompact, Plus, Minus, Divides, Multiplies, Emax };

    /// Plain-data construction; XML-driven builds go through makeObject.
    DecayFunction(Kind kind, bool increasing, double initialEfficacy, double CV) :
        kind(kind), increasing(increasing), initialEfficacy(initialEfficacy),
        het(Sampler::fromMeanCV(Sampler::Kind::Lognormal, 1.0, CV)) {}

    /** Return a new decay function, constructed from an XML element
     * (definition: model/xml/util.cpp — the only code reading the schema).
     *
     * @param elt XML element specifying which function to use and parameters
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    static unique_ptr<DecayFunction> makeObject(
        const SimClock& clock, const scnXml::DecayFunction& elt, const char* eltName
    );

    DecayFunction(const DecayFunction& other) :
        kind(other.kind), increasing(other.increasing),
        initialEfficacy(other.initialEfficacy),
        het(other.het ? other.het->clone() : nullptr),
        hetFactor(other.hetFactor), invL(other.invL), k(other.k),
        Emax(other.Emax), IC50(other.IC50), slope(other.slope),
        initialConcentration(other.initialConcentration),
        f1(other.f1 ? make_unique<DecayFunction>(*other.f1) : nullptr),
        f2(other.f2 ? make_unique<DecayFunction>(*other.f2) : nullptr),
        simClock(other.simClock) {}

    /** Sample a DecayFunction value (should be stored per individual).
     *
     * Note that a DecayFunction is needed to call eval() even if heterogeneity
     * is not wanted. If sigma = 0 then the random number stream will not be
     * touched. */
    unique_ptr<DecayFunction> hetSample(LocalRng& rng) const {
        return hetSample(het->sample(rng));
    }

    /** Generate a DecayFunction value from an existing sample. */
    std::unique_ptr<DecayFunction> hetSample(NormalSample sample) const {
        return hetSample(het->sample(sample));
    }

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
    SimTime sampleAgeOfDecay (LocalRng& rng) const;

    inline double eval(double ageDays) const {
        if(increasing)
            return 1.0 - compute( ageDays ) * initialEfficacy;
        else
            return compute( ageDays ) * initialEfficacy;
    }

    /** Checkpoint a sampled copy, recreating it from this configured function.
     * Wire format frozen: bool present + double hetFactor. */
    void checkpointSample(const unique_ptr<DecayFunction>& sample, ostream& stream) const;
    void checkpointSample(unique_ptr<DecayFunction>& sample, istream& stream) const;

private:
    /// The formula switch: one case per Kind (was one virtual per subclass).
    double compute(double ageDays) const;

    void setFactor(double factor) {
        hetFactor = factor;
        if(f1) f1->setFactor(factor);
        if(f2) f2->setFactor(factor);
    }

    Kind kind;
    bool increasing;
    double initialEfficacy;
    std::unique_ptr<util::Sampler> het;

    double hetFactor = 0.0;
    /// Precomputed scale: 1/L days (Step/Linear/Hill/SmoothCompact),
    /// log(2)/L (Exponential), pow(log 2, 1/k)/L (Weibull). Was per-subclass invL.
    double invL = 0.0;
    double k = 0.0;  // Weibull/Hill/SmoothCompact shape
    // Emax parameters (0 elsewhere):
    double Emax = 0.0, IC50 = 0.0, slope = 0.0, initialConcentration = 0.0;
    unique_ptr<DecayFunction> f1, f2;  // nested children (operator: both; emax: f1)
    /// Set by makeObject; needed by sampleAgeOfDecay (round-to-time-step).
    const SimClock* simClock = nullptr;
};

} }
#endif
