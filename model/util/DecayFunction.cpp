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

#include <algorithm>
#include "util/checkpoint.h"
#include <limits>
#include <memory>

#include "Clock.h"
#include "RunContext.h"
#include "util/sampler.h"
#include "util/UnitParse.h"
#include "util/DecayFunction.h"

namespace OM {
namespace util {

// ----  formula switch (one case per Kind; origins = old virtual subclasses)  ----

double DecayFunction::compute_with(double effectiveAge, double factor) const {
    switch(kind) {
    case Kind::Constant:
        // parity: xml/util.cpp:163 (old ConstantDecayFunction) — note: all decay
        // functions must return 0 when time > 0 and default-constructed.
        if( effectiveAge * factor == numeric_limits<double>::infinity() )
            return 0.0;
        return 1.0;
    case Kind::Step:  // parity: xml/util.cpp:187
        if( effectiveAge * invL * factor < 1.0 )
            return 1.0;
        else
            return 0.0;
    case Kind::Linear:  // parity: xml/util.cpp:215
        if( effectiveAge * invL * factor < 1.0 )
            return 1.0 - effectiveAge * invL * factor;
        else
            return 0.0;
    case Kind::Exponential:  // parity: xml/util.cpp:244
        return exp( -effectiveAge * invL * factor );
    case Kind::Weibull: {  // parity: xml/util.cpp:270
        double p = -pow(effectiveAge * invL * factor, k);
        if(p < -700.0)
            return 0.0;
        else
            return exp(p);
    }
    case Kind::Hill:  // parity: xml/util.cpp:300
        return 1.0 / (1.0 + pow(effectiveAge * invL * factor, k));
    case Kind::SmoothCompact:  // parity: xml/util.cpp:326
        if( effectiveAge * invL * factor < 1.0 )
            return exp( k - k / (1.0 - pow(effectiveAge * invL * factor, 2.0)) );
        else
            return 0.0;
    case Kind::Plus: case Kind::Minus: case Kind::Divides: case Kind::Multiplies: {
        // parity: xml/util.cpp:364 — clamp order std::max(std::min(v, 1.0), 0.0) kept.
        double a = f1->eval_with(effectiveAge, factor), b = f2->eval_with(effectiveAge, factor);
        double v = kind == Kind::Plus ? a + b : kind == Kind::Minus ? a - b
               : kind == Kind::Divides ? a / b : a * b;
        return max(min(v, 1.0), 0.0);
    }
    case Kind::Emax: {
        // Emax model: E(C)=Emax*C^slope/(C^slope + IC50^slope), C=C0*f(t).
        // parity: xml/util.cpp:439
        const double concentrationPow =
            std::pow(initialConcentration * f1->eval_with(effectiveAge, factor), slope);
        const double efficacy = Emax * concentrationPow / (concentrationPow + IC50);
        return std::max(std::min(efficacy, 1.0), 0.0);
    }
    }
    return 0.0;  // unreachable: every Kind has a case
}

SimTime DecayFunction::sampleAgeOfDecay(LocalRng& rng) const {
    switch(kind) {
    case Kind::Constant:
        return sim::future();  // "in the future" (sim::never() would mean the past)
    case Kind::Step:  // parity: xml/util.cpp:195
        return simClock->roundToTSFromDays( 1.0 / invL );
    case Kind::Linear:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:222 (rounds to nearest; may decay instantly or at L)
        return simClock->roundToTSFromDays(rng.uniform_01() / invL);
    case Kind::Exponential:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:248
        return simClock->roundToTSFromDays( -log(rng.uniform_01()) / invL );
    case Kind::Weibull:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:278
        return simClock->roundToTSFromDays( pow( -log(rng.uniform_01()), 1.0/k ) / invL );
    case Kind::Hill:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:304
        return simClock->roundToTSFromDays( pow( 1.0 / rng.uniform_01() - 1.0, 1.0/k ) / invL );
    case Kind::SmoothCompact:
        // parity: order — local RNG stream: draw call order fixes stream state
        // parity: xml/util.cpp:333
        return simClock->roundToTSFromDays( sqrt( 1.0 - k / (k - log( rng.uniform_01() )) ) / invL );
    case Kind::Plus: case Kind::Minus: case Kind::Divides: case Kind::Multiplies:
        // parity: order — the old std::max(f1->…, f2->…) let gcc 16.2 evaluate its
        // arguments RIGHT-TO-LEFT, so f2 draws its uniform first; mirror explicitly.
        {   const SimTime d2 = f2->sampleAgeOfDecay(rng);
            const SimTime d1 = f1->sampleAgeOfDecay(rng);
            return simClock->roundToTSFromDays( max(d1, d2) ); }
    case Kind::Emax:  // parity: xml/util.cpp:447
        return simClock->roundToTSFromDays( f1->sampleAgeOfDecay(rng) );
    }
    return sim::never();  // unreachable: every Kind has a case
}

// -----  interface / static functions  -----
void DecayFunction::checkpointSample(
        const unique_ptr<DecayFunction>& sample, ostream& stream
) const {
    bool present = bool(sample);
    present & stream;
    if (present) sample->hetFactor & stream;
}

void DecayFunction::checkpointSample(
        unique_ptr<DecayFunction>& sample, istream& stream
) const {
    bool present;
    present & stream;
    if (present) {
        double factor;
        factor & stream;
        sample = hetSample(factor);
    } else {
        sample.reset();
    }
}


} }
