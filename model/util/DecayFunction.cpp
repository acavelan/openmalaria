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
#include <limits>
#include <memory>

#include "Global.h"
#include "util/sampler.h"
#include "schema/util.h"
#include "util/UnitParse.h"
#include "util/DecayFunction.h"

namespace OM {
namespace util {
    
inline double readLToDays( const scnXml::DecayFunction& elt ){
    if( !elt.getL().present() ){
        throw xml_scenario_error( "decay function: attribute L required" );
    }
    return UnitParse::durationToDays(elt.getL().get(), UnitParse::YEARS);
}

using Kind = DecayFunction::Kind;

double DecayFunction::compute(double effectiveAge) const {
    switch (kind) {
    case Kind::Constant:
        // Note: we now require all decay functions to return 0 when time > 0
        // and the DecayFunction is default-constructed. So const *after deployment*.
        if( effectiveAge * hetFactor == numeric_limits<double>::infinity() )
            return 0.0;
        return 1.0;
    case Kind::Step:
        if( effectiveAge * invL * hetFactor < 1.0 ){
            return 1.0;
        }else{
            return 0.0;
        }
    case Kind::Linear:
        if( effectiveAge * invL * hetFactor < 1.0 )
            return 1.0 - effectiveAge * invL * hetFactor;
        else
            return 0.0;
    case Kind::Exponential:
        return exp( -effectiveAge * invL * hetFactor);
    case Kind::Weibull: {
        double p = -pow(effectiveAge * constOverLambda * hetFactor, k);
        if(p < -700.0)
            return 0.0;
        else
            return exp(p);
    }
    case Kind::Hill:
        return 1.0 / (1.0 + pow(effectiveAge * invL * hetFactor, k));
    case Kind::SmoothCompact:
        if( effectiveAge * invL * hetFactor < 1.0 )
            return exp( k - k / (1.0 - pow(effectiveAge * invL * hetFactor, 2.0)) );
        else
            return 0.0;
    case Kind::Plus: case Kind::Minus: case Kind::Divides: case Kind::Multiplies: {
        const double v1 = f1->eval(effectiveAge), v2 = f2->eval(effectiveAge);
        double v;
        switch (kind) {
        case Kind::Plus: v = v1 + v2; break;
        case Kind::Minus: v = v1 - v2; break;
        case Kind::Divides: v = v1 / v2; break;
        default: v = v1 * v2; break;   // Kind::Multiplies
        }
        return max(min(v, 1.0), 0.0);
    }
    case Kind::Emax: {
        // Emax model: E(C)=Emax*C^slope/(C^slope + IC50^slope), C=C0*f(t).
        const double concentrationPow =
            std::pow(initialConcentration * f1->eval(effectiveAge), slope);
        const double efficacy = Emax * concentrationPow / (concentrationPow + IC50);
        return std::max(std::min(efficacy, 1.0), 0.0);
    }
    }
    return 0.0;   // unreachable: every Kind is handled above
}

SimTime DecayFunction::sampleAgeOfDecay (LocalRng& rng) const {
    switch (kind) {
    case Kind::Constant:
        return sim::future();        // decay occurs "in the future" (don't use sim::never() because that is interpreted as being in the past)
    case Kind::Step:
        return sim::roundToTSFromDays( 1.0 / invL );
    case Kind::Linear:
        // Note: rounds to nearest. Object may decay instantly or at time L.
        return sim::roundToTSFromDays(rng.uniform_01() / invL);
    case Kind::Exponential:
        return sim::roundToTSFromDays( -log(rng.uniform_01()) / invL );
    case Kind::Weibull:
        return sim::roundToTSFromDays( pow( -log(rng.uniform_01()), 1.0/k ) / constOverLambda );
    case Kind::Hill:
        return sim::roundToTSFromDays( pow( 1.0 / rng.uniform_01() - 1.0, 1.0/k ) / invL );
    case Kind::SmoothCompact:
        return sim::roundToTSFromDays( sqrt( 1.0 - k / (k - log( rng.uniform_01() )) ) / invL );
    case Kind::Plus: case Kind::Minus: case Kind::Divides: case Kind::Multiplies:
        return sim::roundToTSFromDays( max(f1->sampleAgeOfDecay(rng), f2->sampleAgeOfDecay(rng)) );
    case Kind::Emax:
        return sim::roundToTSFromDays( f1->sampleAgeOfDecay(rng) );
    }
    return sim::never();   // unreachable: every Kind is handled above
}

unique_ptr<DecayFunction> DecayFunction::makeHetSample(double factor) const {
    unique_ptr<DecayFunction> sample = std::make_unique<DecayFunction>(*this);
    if (sample->f1) sample->f1 = sample->f1->hetSample(factor);
    if (sample->f2) sample->f2 = sample->f2->hetSample(factor);
    return sample;
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

unique_ptr<DecayFunction> DecayFunction::makeObject(
    const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    Kind kind;
    if( func == "constant" )
        kind = Kind::Constant;
    else if( func == "step" )
        kind = Kind::Step;
    else if( func == "linear" )
        kind = Kind::Linear;
    else if( func == "exponential" )
        kind = Kind::Exponential;
    else if( func == "weibull" )
        kind = Kind::Weibull;
    else if( func == "hill" )
        kind = Kind::Hill;
    else if( func == "smooth-compact" )
        kind = Kind::SmoothCompact;
    else if( func == "plus" )
        kind = Kind::Plus;
    else if( func == "minus" )
        kind = Kind::Minus;
    else if( func == "divides" )
        kind = Kind::Divides;
    else if( func == "multiplies" )
        kind = Kind::Multiplies;
    else if( func == "emax" )
        kind = Kind::Emax;
    else
        throw xml_scenario_error("decay function " + string(eltName) + " of type " + string(func) + " unrecognized");

    unique_ptr<DecayFunction> f = std::make_unique<DecayFunction>(
        kind, elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV() );

    switch (kind) {
    case Kind::Step: case Kind::Linear:
        f->invL = 1.0 / readLToDays(elt);
        break;
    case Kind::Hill: case Kind::SmoothCompact:
        f->invL = 1.0 / readLToDays(elt);
        f->k = elt.getK();
        break;
    case Kind::Exponential:
        f->invL = log(2.0) / readLToDays(elt);
        break;
    case Kind::Weibull:
        f->constOverLambda = pow(log(2.0),1.0/elt.getK()) / readLToDays(elt);
        f->k = elt.getK();
        break;
    case Kind::Plus: case Kind::Minus: case Kind::Divides: case Kind::Multiplies: {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 2)
            throw xml_scenario_error("Operator decay function expects two decay functions, " + to_string(decaySequence.size()) +"  were given.");
        f->f1 = DecayFunction::makeObject(decaySequence[0], "Operator::f1");
        f->f2 = DecayFunction::makeObject(decaySequence[1], "Operator::f2");
        break;
    }
    case Kind::Emax: {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 1)
            throw util::xml_scenario_error("The Emax function expects exactly one user function, but " + to_string(decaySequence.size()) +"  were given.");
        if(!elt.getEmax().present())
            throw util::xml_scenario_error("Emax function: the Emax parameter is not declared");
        if(!elt.getIC50().present())
            throw util::xml_scenario_error("Emax function: the IC50 parameter is not declared");
        if(!elt.getSlope().present())
            throw util::xml_scenario_error("Emax function: the Slope parameter is not declared");
        if(!elt.getInitialConcentration().present())
            throw util::xml_scenario_error("Emax function: the InitialConcentration parameter is not declared");
        
        f->f1 = makeObject(decaySequence[0], "Emax::f");

        const double Emax = elt.getEmax().get();
        const double rawIC50 = elt.getIC50().get();
        const double slope = elt.getSlope().get();
        const double initialConcentration = elt.getInitialConcentration().get();

        if (Emax < 0.0 || Emax > 1.0) {
            throw util::xml_scenario_error(
                "Emax function: Emax must be in [0, 1], got " + std::to_string(Emax));
        }

        if (rawIC50 <= 0.0) {
            throw util::xml_scenario_error(
                "Emax function: IC50 must be > 0, got " + std::to_string(rawIC50));
        }

        if (slope <= 0.0) {
            throw util::xml_scenario_error(
                "Emax function: Slope must be > 0, got " + std::to_string(slope));
        }

        if (initialConcentration < 0.0) {
            throw util::xml_scenario_error(
                "Emax function: InitialConcentration must be >= 0, got " + std::to_string(initialConcentration));
        }

        f->Emax = Emax;
        f->slope = slope;
        f->initialConcentration = initialConcentration;
        // Store IC50^slope once since compute() evaluates this repeatedly.
        f->IC50 = std::pow(rawIC50, slope);
        break;
    }
    case Kind::Constant:
        break;
    }
    return f;
}

} }
