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

#ifndef hzzDecayFunction
#define hzzDecayFunction

#include <limits>
#include <memory>

#include "Global.h"
#include "util/sampler.h"
#include "util/UnitParse.h"

namespace scnXml { class DecayFunction; }

namespace OM {
namespace util {

/** One of a few types of decay function (some of which may also be
 * suitible survival functions).
 *
 * Heterogeneity is implemented by passing a DecayFunction object to the eval
 * function; this should be sampled by the same DecayFunction.
 *****************************************************************************/
class DecayFunction
{
public:
    enum class Kind {
        Constant, Step, Linear, Exponential, Weibull, Hill, SmoothCompact,
        Plus, Minus, Divides, Multiplies, Emax
    };

    DecayFunction(Kind kind, bool increasing = false, double initialEfficacy = 1.0, double CV = 0.0) : 
        kind(kind),
        increasing(increasing),
        initialEfficacy(initialEfficacy),
        het(*LognormalSampler::fromMeanCV(1.0, CV)) {}

    /// Deep copy (including the sub-functions of operator/Emax functions).
    DecayFunction(const DecayFunction& other) : 
        kind(other.kind),
        increasing(other.increasing),
        initialEfficacy(other.initialEfficacy),
        het(other.het),
        invL(other.invL),
        k(other.k),
        constOverLambda(other.constOverLambda),
        Emax(other.Emax),
        IC50(other.IC50),
        slope(other.slope),
        initialConcentration(other.initialConcentration),
        f1(other.f1 ? new DecayFunction(*other.f1) : nullptr),
        f2(other.f2 ? new DecayFunction(*other.f2) : nullptr),
        hetFactor(other.hetFactor) {}

    /** Return a new decay function, constructed from an XML element.
     * 
     * @param elt XML element specifying which function to use and parameters
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    static unique_ptr<DecayFunction> makeObject(
        const scnXml::DecayFunction& elt, const char* eltName
    );

    /** Sample a DecayFunction value (should be stored per individual).
     * 
     * Note that a DecayFunction is needed to call eval() even if heterogeneity
     * is not wanted. If sigma = 0 then the random number stream will not be
     * touched. */
    unique_ptr<DecayFunction> hetSample(LocalRng& rng) const {
        return hetSample(het.sample(rng));
    }
    
    /** Generate a DecayFunction value from an existing sample. */
    std::unique_ptr<DecayFunction> hetSample(NormalSample sample) const {
        return hetSample(het.sample(sample));
    }
    
    /** Generate a DecayFunction value from an existing sample. */
    unique_ptr<DecayFunction> hetSample(double factor) const {
        unique_ptr<DecayFunction> sample = makeHetSample(factor);
        sample->hetFactor = factor;
        return sample;
    }

    /** Say you have a population of objects which each have two states:
     * decayed and not decayed. If you want to use a DecayFunction to model
     * the proportion of objects which have decayed, you need to work out per
     * object the age of decay. This function does that.
     * 
     * This is only valid where mu and sigma parameters are zero.
     * 
     * @returns Age at which an object should decay. */
    SimTime sampleAgeOfDecay (LocalRng& rng) const;
    
    inline double eval(double ageDays) const {
        if(increasing)
            return 1.0 - compute( ageDays ) * initialEfficacy;
        else
            return compute( ageDays ) * initialEfficacy;
    }

    double compute(double ageDays) const;

    /** Checkpoint a sampled copy, recreating it from this configured function. */
    void checkpointSample(const unique_ptr<DecayFunction>& sample, ostream& stream) const;
    void checkpointSample(unique_ptr<DecayFunction>& sample, istream& stream) const;

    /** Deep copy this function; with sub-functions het-sampled at factor. */
    unique_ptr<DecayFunction> makeHetSample(double factor) const;

public:
    Kind kind;
    bool increasing;
    double initialEfficacy;
    util::LognormalSampler het;
    // Parameters, as used by the relevant kind(s):
    double invL = 0.0;          ///< 1 / lifetime (days)
    double k = 0.0;             ///< Weibull / Hill / SmoothCompact shape
    double constOverLambda = 0.0;   ///< Weibull
    double Emax = 0.0;          ///< Emax, IC50^slope, slope, C0:
    double IC50 = 0.0;          ///<   Emax function parameters
    double slope = 0.0;
    double initialConcentration = 0.0;
    unique_ptr<DecayFunction> f1, f2;   ///< operator / Emax sub-functions

public:
    double hetFactor = 0.0;
};

} }
#endif
