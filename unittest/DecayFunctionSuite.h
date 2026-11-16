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

#ifndef Hmod_DecayFunctionSuite
#define Hmod_DecayFunctionSuite

#include <cxxtest/TestSuite.h>
#include <bit>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "util/DecayFunction.h"
#include "ScenarioXml.h"

using ::OM::util::DecayFunction;

class DecayFunctionSuite : public CxxTest::TestSuite
{
public:
    DecayFunctionSuite () :
        m_rng(0, 0), dfElt( "" )
    {
        dfElt.setL( "10y" );
        dfElt.setK( 1.6 );
    }
    
    void setUp() {
        m_rng.seed(0, 721347520444481703);
        UnittestUtil::initTime(5);
    }
    
    void testBad() {
        dfElt.setFunction( "unknown" );
        string errMsg = "decay function DecayFunctionSuite of type unknown unrecognized";
        TS_ASSERT_THROWS_EQUALS( df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" ),
                          const std::runtime_error &e, e.what(), errMsg );
    }
    
    void testConstant () {
        dfElt.setFunction( "constant" );
        df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" );
        // Second test: with an appropriately sampled helper value, we should get the results we want
        double dHet = df.heterogeneity.sample(m_rng);
        TS_ASSERT_APPROX( OM::util::eval(df, sim::zero(), dHet) , 1.0 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(10), dHet), 1.0 );
        // Third test: time of decay (age plus now) should always be in the future
        TS_ASSERT_EQUALS( OM::util::sampleAgeOfDecay(df, m_rng, sim::oneTS()), sim::future() );
    }
    
    void testStep () {
        dfElt.setFunction( "step" );
        df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" );
        double dHet = df.heterogeneity.sample(m_rng);
        TS_ASSERT_APPROX( OM::util::eval(df, sim::zero(), dHet), 1.0 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(6), dHet), 1.0 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(20), dHet), 0.0 );
        TS_ASSERT_EQUALS( OM::util::sampleAgeOfDecay(df, m_rng, sim::oneTS()), sim::fromYearsI(10) );
    }
    
    void testLinear () {
        dfElt.setFunction( "linear" );
        df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" );
        double dHet = df.heterogeneity.sample(m_rng);
        TS_ASSERT_APPROX( OM::util::eval(df, sim::zero(), dHet), 1.0 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(6), dHet), 0.4 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(20), dHet), 0.0 );
    }
    
    void testExponential () {
        dfElt.setFunction( "exponential" );
        df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" );
        double dHet = df.heterogeneity.sample(m_rng);
        TS_ASSERT_APPROX( OM::util::eval(df, sim::zero(), dHet), 1.0 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(6), dHet), 0.65975394736842108 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(20), dHet), 0.25 );
    }
    
    void testWeibull () {
        dfElt.setFunction( "weibull" );
        df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" );
        double dHet = df.heterogeneity.sample(m_rng);
        TS_ASSERT_APPROX( OM::util::eval(df, sim::zero(), dHet), 1.0 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(6), dHet), 0.73631084210526321 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(20), dHet), 0.122306 );
    }
    
    void testHill () {
        dfElt.setFunction( "hill" );
        df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" );
        double dHet = df.heterogeneity.sample(m_rng);
        TS_ASSERT_APPROX( OM::util::eval(df, sim::zero(), dHet), 1.0 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(6), dHet), 0.6936673684210527 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(20), dHet), 0.24805074736842106 );
    }
    
    void testSmoothCompact () {
        dfElt.setFunction( "smooth-compact" );
        df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" );
        double dHet = df.heterogeneity.sample(m_rng);
        TS_ASSERT_APPROX( OM::util::eval(df, sim::zero(), dHet), 1.0 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(6), dHet), 0.40656965789473687 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(20), dHet), 0.0 );
    }

    void testEmax () {
        dfElt.setFunction( "emax" );
        dfElt.setEmax( 0.8 );
        dfElt.setIC50( 1.0 );
        dfElt.setSlope( 1.0 );
        dfElt.setInitialConcentration( 2.0 );
        dfElt.getDecay().clear();
        scnXml::DecayFunction inner( "exponential" );
        inner.setL( "10y" );
        dfElt.getDecay().push_back( inner );

        df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" );
        double dHet = df.heterogeneity.sample(m_rng);
        TS_ASSERT_APPROX( OM::util::eval(df, sim::zero(), dHet), 0.5333333333333333 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(10), dHet), 0.4 );
        TS_ASSERT_APPROX( OM::util::eval(df, sim::fromYearsI(20), dHet), 0.26666666666666666 );
    }

    void testEmaxBadIC50 () {
        dfElt.setFunction( "emax" );
        dfElt.setEmax( 0.8 );
        dfElt.setIC50( -1.0 );
        dfElt.setSlope( 0.5 );
        dfElt.setInitialConcentration( 2.0 );
        dfElt.getDecay().clear();
        scnXml::DecayFunction inner( "constant" );
        dfElt.getDecay().push_back( inner );

        string errMsg = "Emax function: IC50 must be > 0, got -1.000000";
        TS_ASSERT_THROWS_EQUALS( df = OM::xml::readDecay( dfElt, "DecayFunctionSuite" ),
                          const std::runtime_error &e, e.what(), errMsg );
    }
    
    void testNativeCompositionAndSharedFactor() {
        using namespace OM::util;
        DecayFunction curve;
        DecayFunction first, second;
        first.initialEfficacy = 0.75;
        second.initialEfficacy = 0.25;
        curve.children = {first, second};
        const DecayKind kinds[] = {DecayKind::Plus, DecayKind::Minus, DecayKind::Divides, DecayKind::Multiplies};
        const double expected[] = {1.0, 0.5, 1.0, 0.1875};
        for (int i = 0; i < 4; ++i) {
            curve.kind = kinds[i];
            TS_ASSERT_EQUALS(eval(curve, 10.0, 1.0), expected[i]);
        }
        curve.kind = DecayKind::Multiplies;
        curve.children[0].kind = DecayKind::Step;
        curve.children[0].rate = 0.1;
        TS_ASSERT_EQUALS(eval(curve, 6.0, 1.0), 0.1875);
        TS_ASSERT_EQUALS(eval(curve, 6.0, 2.0), 0.0);
        curve.increasing = true;
        TS_ASSERT_EQUALS(eval(curve, 6.0, 2.0), 1.0);
    }

    void testCompositePreservesNonFiniteArithmetic() {
        using namespace OM::util;
        DecayFunction first, second, combined;
        first.kind = DecayKind::Exponential;
        first.rate = log(2.0) / 1095.0;
        second.kind = DecayKind::Weibull;
        second.rate = pow(log(2.0), 0.5) / 1825.0;
        second.shape = 2.0;
        combined.children = {first, second};
        for (double age : {0.0, numeric_limits<double>::infinity()}) {
            double factor = age == 0.0 ? numeric_limits<double>::infinity() : 0.0;
            double a = eval(first, age, factor), b = eval(second, age, factor);
            combined.kind = DecayKind::Plus;
            double expected = max(min(std::plus<double>{}(a, b), 1.0), 0.0);
            TS_ASSERT_EQUALS(std::bit_cast<uint64_t>(eval(combined, age, factor)), std::bit_cast<uint64_t>(expected));
            combined.kind = DecayKind::Multiplies;
            expected = max(min(std::multiplies<double>{}(a, b), 1.0), 0.0);
            TS_ASSERT_EQUALS(std::bit_cast<uint64_t>(eval(combined, age, factor)), std::bit_cast<uint64_t>(expected));
        }
    }

    void testLifetimeUsesExplicitIntervalAndPreservesRng() {
        using namespace OM::util;
        DecayFunction curve;
        curve.kind = DecayKind::Step;
        curve.rate = 1.0 / 17.0;
        LocalRng rng(17, 0), reference(17, 0);
        TS_ASSERT_EQUALS(sampleAgeOfDecay(curve, rng, 1), 17);
        TS_ASSERT_EQUALS(sampleAgeOfDecay(curve, rng, 5), 15);
        curve.kind = DecayKind::Linear;
        DecayFunction composed;
        composed.kind = DecayKind::Plus;
        composed.children = {curve, curve};
        sampleAgeOfDecay(composed, rng, 5);
        reference.uniform_01();
        reference.uniform_01();
        std::ostringstream actual, expected;
        rng.checkpoint(actual);
        reference.checkpoint(expected);
        TS_ASSERT_EQUALS(actual.str(), expected.str());
    }

    void testSampleCheckpointKeepsLegacyLayoutAndRestoresAbsence() {
        using namespace OM::util;
        DecaySample sample = 0.75;
        std::ostringstream actual, legacy;
        checkpointDecaySample(sample, actual);
        bool present = true;
        double factor = 0.75;
        present & legacy;
        factor & legacy;
        TS_ASSERT_EQUALS(actual.str(), legacy.str());
        DecaySample restored;
        std::istringstream input(actual.str());
        checkpointDecaySample(restored, input);
        TS_ASSERT(restored.has_value());
        TS_ASSERT_EQUALS(*restored, 0.75);
        std::ostringstream absent;
        checkpointDecaySample(DecaySample{}, absent);
        std::istringstream absentInput(absent.str());
        checkpointDecaySample(restored, absentInput);
        TS_ASSERT(!restored.has_value());
    }

    void testXmlCompositeOperatorsEvaluateTheirNamedOperation() {
        scnXml::DecayFunction first("constant"), second("constant");
        first.setInitialEfficacy(0.25);
        second.setInitialEfficacy(0.5);
        for (const auto& [name, expected] : {
                std::pair{"plus", 0.75}, {"minus", 0.0},
                {"divides", 0.5}, {"multiplies", 0.125}}) {
            scnXml::DecayFunction input(name);
            input.getDecay().push_back(first);
            input.getDecay().push_back(second);
            auto curve = OM::xml::readDecay(input, "test");
            TS_ASSERT_EQUALS(OM::util::eval(curve, 10.0, 1.0), expected);
        }
    }

    void testXmlRejectsMissingScaleAndChildren() {
        scnXml::DecayFunction missingScale("linear");
        TS_ASSERT_THROWS(OM::xml::readDecay(missingScale, "test"), const std::runtime_error&);
        for (const char* kind : {"plus", "minus", "divides", "multiplies", "emax"}) {
            scnXml::DecayFunction missingChildren(kind);
            TS_ASSERT_THROWS(OM::xml::readDecay(missingChildren, "test"), const std::runtime_error&);
        }
    }

private:
    LocalRng m_rng;
    scnXml::DecayFunction dfElt;
    DecayFunction df;
};

#endif
