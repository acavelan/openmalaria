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

// model/xml/ conversion layer: every scnXml -> plain-value function lives
// under model/xml/; only this directory includes generated schema headers
// (plan §5, §8 F7: XSD codegen + Xerces validation kept). Everything below was
// moved VERBATIM from the noted origin files; declarations stay in the normal
// model headers (scnXml types forward-declared there). Call sites unchanged.

#include "Global.h"
#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/DecayFunction.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "util/sampler.h"
#include <schema/healthSystem.h>
#include <schema/util.h>

namespace OM { namespace util {

// ---- sampler readers (origins: util/sampler.{h,cpp}) ----

void NormalSampler::setParams(const scnXml::SampledValueN& elt){
    mu = elt.getMean();
    if( elt.getDistr() == "const" ){
        if( elt.getSD().present() && elt.getSD().get() != 0.0) {
            throw util::xml_scenario_error( "attribute SD must be zero or omitted when distr=\"const\" or is omitted" );
        }
        sigma = 0.0;
        return;
    }
    
    if( elt.getDistr() != "normal" ){
        throw util::xml_scenario_error( "expected distr to be one of \"const\", \"lognormal\" (note: not all distributions are supported here)" );
    }
    if( !elt.getSD().present() ){
        throw util::xml_scenario_error( "attribute \"SD\" required for sampled normal value when distr is not const" );
    }
    sigma = elt.getSD().get();
}

    /** Create a new <SamplerT> from a xml snippet. */
    template <typename SamplerT>
    inline unique_ptr<SamplerT> createSampler(double mean, const scnXml::SampledValueCV& elt)
    {
        // 1.  Compile‑time identify which concrete sampler we are handling ─────
        constexpr bool isLognormal = std::is_same_v<SamplerT, LognormalSampler>;
        constexpr bool isGamma     = std::is_same_v<SamplerT, GammaSampler>;
        constexpr const char* expectedDistr =
            isLognormal ? "lognormal" :
            isGamma     ? "gamma"
                        : nullptr;

        const std::string& d = elt.getDistr();

        // 2. Special handling for the “const” alias (lognormal with CV = 0)
        if( d == "const" )
        {
            if constexpr (!isLognormal) // only allowed for lognormal
                throw util::xml_scenario_error("\"distr=\"const\": expected \"lognormal\"");

            if( elt.getCV().present() && elt.getCV().get() != 0.0 )
                throw util::xml_scenario_error( "\"distr="+d+"\": attribute \"CV\" must be zero or omitted when distr=\"const\" or is omitted" );
            if( elt.getVariance().present() && elt.getVariance().get() != 0.0 )
                throw util::xml_scenario_error( "\"distr="+d+"\": attribute \"variance\" must be zero or omitted when distr=\"const\" or is omitted" );
            
            return SamplerT::fromMeanCV(mean, 0.0);
        }

        // 3. Check that XML 'distr' matches the compile‑time sampler type
        if (expectedDistr == nullptr)
            throw util::xml_scenario_error("createSampler<SamplerT>: expected SamplerT to be a valid Sampler");

        if (d != expectedDistr)
            throw util::xml_scenario_error("\"distr=" + d + "\": expected \"" + expectedDistr + "\" or \"const\"");

        // 4. Exactly one of CV / variance must be supplied
        const bool hasCV  = elt.getCV().present();
        const bool hasVar = elt.getVariance().present();

        if (hasCV == hasVar)   // both present or both absent
            throw util::xml_scenario_error("\"distr=" + d + "\": exactly one of attributes \"CV\" or \"variance\" must be specified");

        // 5. Check for the optional truncation parameter
        std::optional<double> truncate = elt.getTruncate().present()
            ? std::optional<double>{elt.getTruncate().get()}
            : std::nullopt;

        // 6. Build the sampler
        if(hasCV)
            return SamplerT::fromMeanCV(mean, elt.getCV().get(), truncate);
        else
            return SamplerT::fromMeanVariance(mean, elt.getVariance().get(), truncate);
    }

    /** Create a Sampler from a SampledValueCV xml snippet. */
    std::unique_ptr<util::Sampler> createSampler(double mean, const scnXml::SampledValueCV& elt)
    {
        const std::string& d = elt.getDistr();

        if (d == "const" || d == "lognormal")
            return createSampler<util::LognormalSampler>(mean, elt);
        else if (d == "gamma")
            return createSampler<util::GammaSampler>(mean,  elt);
        else
            throw util::xml_scenario_error("\"distr=" + d + "\": expected \"const\", \"lognormal\" or \"gamma\"");
    }

    /** Create a LognormalSampler from a SampledValueLN xml snippet. */
    std::unique_ptr<LognormalSampler> createSampler(const scnXml::SampledValueLN& elt)
    {
        return createSampler<LognormalSampler>(elt.getMean(), static_cast<const scnXml::SampledValueCV&>(elt));
    }

// ---- BetaSampler / WeibullSampler XML readers (origin: util/sampler.h) ----
void BetaSampler::setParams( const scnXml::BetaMeanSample& elt ){
    setParamsMV( elt.getMean(), elt.getVariance() );
}

void WeibullSampler::setParams( const scnXml::WeibullSample& elt ){
    setScaleShape( elt.getScale(), elt.getShape() );
}

// ---- DecayFunction subclasses + factory (origin: util/DecayFunction.cpp) ----

inline double readLToDays( const scnXml::DecayFunction& elt ){
    if( !elt.getL().present() ){
        throw xml_scenario_error( "decay function: attribute L required" );
    }
    return UnitParse::durationToDays(elt.getL().get(), UnitParse::YEARS);
}

class ConstantDecayFunction : public DecayFunction {
public:
    ConstantDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV())
    {}
    
    double compute(double effectiveAge) const {
        // Note: we now require all decay functions to return 0 when time > 0
        // and the DecayFunction is default-constructed. So const *after deployment*.
        if( effectiveAge * hetFactor == numeric_limits<double>::infinity() )
            return 0.0;
        return 1.0;
    }
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::future();        // decay occurs "in the future" (don't use sim::never() because that is interpreted as being in the past)
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<ConstantDecayFunction> copy = std::make_unique<ConstantDecayFunction>(*this);
        return std::move(copy);
    }
};

class StepDecayFunction : public DecayFunction {
public:
    StepDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) )
    {}
    
    double compute(double effectiveAge) const{
        if( effectiveAge * invL * hetFactor < 1.0 ){
            return 1.0;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( 1.0 / invL );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<StepDecayFunction> copy = make_unique<StepDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double invL;
};

class LinearDecayFunction : public DecayFunction {
public:
    LinearDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) )
    {}
    
    double compute(double effectiveAge) const{
        if( effectiveAge * invL * hetFactor < 1.0 )
            return 1.0 - effectiveAge * invL * hetFactor;
        else
            return 0.0;
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // Note: rounds to nearest. Object may decay instantly or at time L.
        // parity: order — local RNG stream: draw call order fixes stream state
        return sim::roundToTSFromDays(rng.uniform_01() / invL);
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<LinearDecayFunction> copy = make_unique<LinearDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double invL;
};

class ExponentialDecayFunction : public DecayFunction {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( log(2.0) / readLToDays(elt) )
    {}
    
    double compute(double effectiveAge) const{
        return exp( -effectiveAge * invL * hetFactor);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // parity: order — local RNG stream: draw call order fixes stream state
        return sim::roundToTSFromDays( -log(rng.uniform_01()) / invL );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<ExponentialDecayFunction> copy = make_unique<ExponentialDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double invL;
};

class WeibullDecayFunction : public DecayFunction {
public:
    WeibullDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        constOverLambda( pow(log(2.0),1.0/elt.getK()) / readLToDays(elt) ),
        k( elt.getK() )
    {}
    
    double compute(double effectiveAge) const{
        double p = -pow(effectiveAge * constOverLambda * hetFactor, k);
        if(p < -700.0)
            return 0.0;
        else
            return exp(p);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // parity: order — local RNG stream: draw call order fixes stream state
        return sim::roundToTSFromDays( pow( -log(rng.uniform_01()), 1.0/k ) / constOverLambda );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<WeibullDecayFunction> copy = make_unique<WeibullDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double constOverLambda, k;
};

class HillDecayFunction : public DecayFunction {
public:
    HillDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        k( elt.getK() )
    {}
    
    double compute(double effectiveAge) const{
        return 1.0 / (1.0 + pow(effectiveAge * invL * hetFactor, k));
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // parity: order — local RNG stream: draw call order fixes stream state
        return sim::roundToTSFromDays( pow( 1.0 / rng.uniform_01() - 1.0, 1.0/k ) / invL );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<HillDecayFunction> copy = make_unique<HillDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double invL, k;
};

class SmoothCompactDecayFunction : public DecayFunction {
public:
    SmoothCompactDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        k( elt.getK() )
    {}

    double compute(double effectiveAge) const{
        if( effectiveAge * invL * hetFactor < 1.0 )
            return exp( k - k / (1.0 - pow(effectiveAge * invL * hetFactor, 2.0)) );
        else
            return 0.0;
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // parity: order — local RNG stream: draw call order fixes stream state
        return sim::roundToTSFromDays( sqrt( 1.0 - k / (k - log( rng.uniform_01() )) ) / invL );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<SmoothCompactDecayFunction> copy = make_unique<SmoothCompactDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double invL, k;
};

template<class T>
class OperatorDecayFunction : public DecayFunction {
public:
    OperatorDecayFunction( const scnXml::DecayFunction& elt ) : 
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()) {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 2)
            throw xml_scenario_error("Operator decay function expects two decay functions, " + to_string(decaySequence.size()) +"  were given.");

        f1 = DecayFunction::makeObject(decaySequence[0], "Operator::f1");
        f2 = DecayFunction::makeObject(decaySequence[1], "Operator::f2");
    }

    OperatorDecayFunction(const OperatorDecayFunction &copy, unique_ptr<DecayFunction> f1, unique_ptr<DecayFunction> f2) : 
        DecayFunction(copy),
        f1(std::move(f1)), f2(std::move(f2)) {}

    double compute(double effectiveAge) const {
        return max(min(op(f1->eval(effectiveAge), f2->eval(effectiveAge)), 1.0), 0.0);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::roundToTSFromDays( max(f1->sampleAgeOfDecay(rng), f2->sampleAgeOfDecay(rng)) );
    }

    unique_ptr<DecayFunction> makeHetSample(double hetFactor) const {
        unique_ptr<DecayFunction> f1hetSample = f1->hetSample(hetFactor);
        unique_ptr<DecayFunction> f2hetSample = f2->hetSample(hetFactor);
        unique_ptr<OperatorDecayFunction> copy = make_unique<OperatorDecayFunction>(*this, std::move(f1hetSample), std::move(f2hetSample));
        return std::move(copy);
    }

private:
    unique_ptr<DecayFunction> f1, f2;
    T op;
};

class EmaxFunction : public DecayFunction {
public:
    EmaxFunction( const scnXml::DecayFunction& elt ) : 
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()) {
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
        
        f = makeObject(decaySequence[0], "Emax::f");

        Emax = elt.getEmax().get();
        const double rawIC50 = elt.getIC50().get();
        slope = elt.getSlope().get();
        initialConcentration = elt.getInitialConcentration().get();

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

        // Store IC50^slope once since compute() evaluates this repeatedly.
        IC50 = std::pow(rawIC50, slope);
    }

    EmaxFunction(const EmaxFunction &copy, std::unique_ptr<DecayFunction> f) :
        DecayFunction(copy),
        f(std::move(f)),
        Emax(copy.Emax),
        IC50(copy.IC50),
        slope(copy.slope),
        initialConcentration(copy.initialConcentration) {}

    double compute(double effectiveAge) const {
        // Emax model: E(C)=Emax*C^slope/(C^slope + IC50^slope), C=C0*f(t).
        const double concentrationPow =
            std::pow(initialConcentration * f->eval(effectiveAge), slope);
        const double efficacy = Emax * concentrationPow / (concentrationPow + IC50);
        return std::max(std::min(efficacy, 1.0), 0.0);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::roundToTSFromDays( f->sampleAgeOfDecay(rng) );
    }

    std::unique_ptr<DecayFunction> makeHetSample(double hetFactor) const {
        std::unique_ptr<DecayFunction> fhetSample = f->hetSample(hetFactor);
        std::unique_ptr<EmaxFunction> copy =
            std::make_unique<EmaxFunction>(*this, std::move(fhetSample));
        return std::move(copy);
    }

private:
    std::unique_ptr<DecayFunction> f;
    double Emax, IC50, slope, initialConcentration;
};

unique_ptr<DecayFunction> DecayFunction::makeObject(
    const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    if( func == "constant" )
        return make_unique<ConstantDecayFunction>( elt );
    else if( func == "step" )
        return make_unique<StepDecayFunction>( elt );
    else if( func == "linear" )
        return make_unique<LinearDecayFunction>( elt );
    else if( func == "exponential" )
        return make_unique<ExponentialDecayFunction>( elt );
    else if( func == "weibull" )
        return make_unique<WeibullDecayFunction>( elt );
    else if( func == "hill" )
        return make_unique<HillDecayFunction>( elt );
    else if( func == "smooth-compact" )
        return make_unique<SmoothCompactDecayFunction>( elt );
    else if( func == "plus" )
        return make_unique<OperatorDecayFunction<std::plus<double>>>( elt );
    else if( func == "minus" )
        return make_unique<OperatorDecayFunction<std::minus<double>>>( elt );
    else if( func == "divides" )
        return make_unique<OperatorDecayFunction<std::divides<double>>>( elt );
    else if( func == "multiplies" )
        return make_unique<OperatorDecayFunction<std::multiplies<double>>>( elt );
    else if( func == "emax" )
        return make_unique<EmaxFunction>( elt );
    else
        throw xml_scenario_error("decay function " + string(eltName) + " of type " + string(func) + " unrecognized");
}

// ---- AgeGroup interpolation factories (origin: util/AgeGroupInterpolation.cpp) ----

    /** This class gives direct access to input age-group
    * data (discontinuous).
    ********************************************/
    class AgeGroupPiecewiseConstant : public AgeGroupInterpolation
    {
    public:
        AgeGroupPiecewiseConstant ( const scnXml::AgeGroupValues& ageGroups,
                                    const char* eltName )
        {
            if( ageGroups.getGroup().size() == 0 ){
                throw util::xml_scenario_error(string(eltName) + ": at least one age group required" );
            }

            auto pos = dataGroups.begin();
            assert(pos == dataGroups.end());

            double greatestLbound = -1.0;
            for( const scnXml::Group& group : ageGroups.getGroup() ){
                double lbound = group.getLowerbound();
                if( lbound >= greatestLbound ){
                    greatestLbound = lbound;
                    pos = dataGroups.insert( pos, make_pair( lbound, group.getValue() ) );
                } else
                    throw util::xml_scenario_error(string(eltName) + ": lower bound " + to_string(lbound) + " not greater than previous " + to_string(greatestLbound));
            }

            // first lower-bound must be 0
            if( dataGroups.begin()->first != 0.0 )
                throw util::xml_scenario_error(string(eltName) + ": first lower-bound must be 0");

            outputSamples( eltName );
        }

        virtual double eval( double ageYears ) const {
            assert ( ageYears >= 0.0 );
            auto it = dataGroups.upper_bound( ageYears );
            --it;   // previous point when ordered by age (it->first <= ageYears < (it+1)->first)
            return it->second;
        }
        
        virtual void scale( double factor ){
            for( auto it = dataGroups.begin(); it != dataGroups.end(); ++it ){
                it->second *= factor;
            }
        }

        virtual double firstGlobalMaximum() {
            double age = 0.0, max = numeric_limits<double>::min();
            for( auto it = dataGroups.begin(); it != dataGroups.end(); ++it ){
                if( it->second > max ){
                    max = it->second;
                    age = it->first;
                }
            }
            return age;
        }

    protected:
        // All data groups as (lower-age-bound, value) pairs.
        map<double,double> dataGroups;
    };

    
    /** Filter class to convert data groups into points in middle of groups plus
     * stabilization points at ends.
     ********************************************/
    class AgeGroupInterpolationPoints : public AgeGroupInterpolation
    {
    public:
        AgeGroupInterpolationPoints( const scnXml::AgeGroupValues& ageGroups,
                                     const char* eltName )
        {
            // Our read iterator
            auto it = ageGroups.getGroup().begin();

            if( it == ageGroups.getGroup().end() )
                throw util::xml_scenario_error(string(eltName) + ": at least one age group required");

            // first lower-bound must be 0
            if( it->getLowerbound() != 0.0 )
                throw util::xml_scenario_error(string(eltName) + ": first lower-bound must be 0");

            // Our insert iterator (used for performance)
            auto pos = dataPoints.begin();
            assert(pos == dataPoints.end());

            double greatestLbound = 0.0;    // also last lbound: age groups must be ordered lowest to highest
            double lastValue =it->getValue();       // first value is repeated for const start
            for( ; it!=ageGroups.getGroup().end(); ++it ){
                double lbound = it->getLowerbound();
                if( lbound >= greatestLbound ){
                    double insBound = 0.5*(greatestLbound + lbound);
                    pos = dataPoints.insert( pos, make_pair( insBound, lastValue ) );
                    greatestLbound = lbound;
                    lastValue = it->getValue();     // always insert previous value
                } else
                    throw util::xml_scenario_error(string(eltName) + ": lower bound " + to_string(lbound) + " less than previous " + to_string(greatestLbound));
            }

            // add a point in middle of last age group (taking upper bound as max-age-years:
            dataPoints[ 0.5*(greatestLbound + sim::inYears(sim::maxHumanAge())) ] = lastValue;
        }

        virtual void scale( double factor ){
            for( auto it = dataPoints.begin(); it != dataPoints.end(); ++it ){
                it->second *= factor;
            }
        }
        
    protected:
        // Points to interpolate between in the middle of input age groups. Extra
        // points at zero and infinity are added with equal value to first and last
        // points respectively.
        map<double,double> dataPoints;
    };

    
    /** This class gives piecewise linear inpolation on top of input age-group
     * data (continuous but with discontinuous derivative).
     ********************************************/
    class AgeGroupPiecewiseLinear : public AgeGroupInterpolationPoints
    {
    public:
        AgeGroupPiecewiseLinear (
            const scnXml::AgeGroupValues& ageGroups, const char* eltName
        ) : AgeGroupInterpolationPoints(ageGroups, eltName)
        {
            // Add first point at zero
            dataPoints[ 0.0 ] = dataPoints.begin()->second;
            // add a point at infinity to catch remaining values
            dataPoints[ numeric_limits<double>::infinity() ] = dataPoints.rbegin()->second;
            
            outputSamples( eltName );
        }
        
        virtual double eval( double ageYears ) const{
            assert ( ageYears >= 0.0 );
            //TODO: optimise by putting (f1-f0)/(a1-a0) into the array and checking usage
            // Get the first point with age greater than ageYears:
            auto it = dataPoints.upper_bound( ageYears );
            assert( it != dataPoints.end() );
            double a1 = it->first;  // a1 > ageYears
            double f1 = it->second;
            --it;   // previous point when ordered by age
            double a0 = it->first;  // a0 ≤ ageYears
            double f0 = it->second;
            return (ageYears - a0) / (a1 - a0) * (f1 - f0) + f0;
        }

        virtual double firstGlobalMaximum() {
            double age = 0.0, max = numeric_limits<double>::min();
            for( auto it = dataPoints.begin(); it != dataPoints.end(); ++it ){
                if( it->second > max ){
                    max = it->second;
                    age = it->first;
                }
            }
            return age;
        }
    };

    void AgeGroupInterpolator::set(
        const scnXml::AgeGroupValues& ageGroups, const char* eltName
    ){
        reset();
        // Type mostly equivalent to a std::string:
        const scnXml::AgeGroupValues::InterpolationOptional& interp = ageGroups.getInterpolation();
        if( !interp.present() || interp.get() == "linear" ){
            obj = new AgeGroupPiecewiseLinear( ageGroups, eltName );
        }else if( interp.get() == "none" ){
            obj = new AgeGroupPiecewiseConstant( ageGroups, eltName );
        }else
            throw util::xml_scenario_error(string("age group interpolation ") + interp.get() + " not implemented" );
    }

} }
