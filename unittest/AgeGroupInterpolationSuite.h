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

#ifndef Hmod_AgeGroupInterpolationSuite
#define Hmod_AgeGroupInterpolationSuite

#include <cxxtest/TestSuite.h>

#include "ExtraAsserts.h"
#include "util/AgeGroupInterpolation.h"

using OM::util::AgeGroupInterpolator;

class AgeGroupInterpolationSuite : public CxxTest::TestSuite
{
public:
    AgeGroupInterpolationSuite () {
        for (size_t i = 0; i < dataLen; ++i)
            groups.groups.push_back({stdLbounds[i], stdValues[i]});
    }
    ~AgeGroupInterpolationSuite () {
    }
    
    void setUp () {
    }
    void tearDown () {
    }
    
    void testDummy() {
        AgeGroupInterpolator o;
        TS_ASSERT_THROWS( o.eval(5.7), const std::exception& );
    }
    
    void testPiecewiseConst () {
        groups.interpolation = OM::util::AgeInterpolation::Constant;
        AgeGroupInterpolator o;
        o.set( groups, 90.0, "testPiecewiseConst" );
        for( size_t i = 0; i < testLen; ++i ){
            TS_ASSERT_APPROX( o.eval( testAges[ i ] ), piecewiseConstValues[ i ] );
        }
    }
    
    void testLinearInterp () {
        groups.interpolation = OM::util::AgeInterpolation::Linear;
        AgeGroupInterpolator o;
        o.set( groups, 90.0, "testLinearInterp" );
        for( size_t i = 0; i < testLen; ++i ){
            TS_ASSERT_APPROX( o.eval( testAges[ i ] ), linearInterpValues[ i ] );
        }
    }
    
    void testExplicitMaximumAgeAndIndependentCopies() {
        OM::util::AgeGroupValues input{OM::util::AgeInterpolation::Linear, {{0, 2}, {20, 6}}};
        AgeGroupInterpolator a, b;
        a.set(input, 60, "a");
        b.set(input, 100, "b");
        TS_ASSERT_EQUALS(a.eval(40), 6.0);
        TS_ASSERT_EQUALS(b.eval(40), 4.4);
        auto copy = a;
        copy.scale(2.0);
        TS_ASSERT_EQUALS(a.eval(40), 6.0);
        TS_ASSERT_EQUALS(copy.eval(40), 12.0);
        TS_ASSERT_EQUALS(a.firstGlobalMaximum(), 40.0);
    }

    void testValidationAndDuplicateBounds() {
        using namespace OM::util;
        for (auto kind : {AgeInterpolation::Linear, AgeInterpolation::Constant}) {
            AgeGroupInterpolator curve;
            AgeGroupValues input{kind, {}};
            TS_ASSERT_THROWS(curve.set(input, 90, "test"), const xml_scenario_error&);
            input.groups = {{1, 2}};
            TS_ASSERT_THROWS(curve.set(input, 90, "test"), const xml_scenario_error&);
            input.groups = {{0, 2}, {20, 3}, {10, 4}};
            TS_ASSERT_THROWS(curve.set(input, 90, "test"), const xml_scenario_error&);
            TS_ASSERT(!curve.isSet());
            input.groups = {{0, 2}, {0, 9}, {20, 4}};
            curve.set(input, 90, "test");
            TS_ASSERT_EQUALS(curve.eval(0), 2.0);
            curve.reset();
            TS_ASSERT_THROWS(curve.eval(0), const std::exception&);
        }
    }

private:
    static const size_t dataLen = 5;
    static const size_t testLen = 8;
    static const double stdLbounds[dataLen];
    static const double stdValues[dataLen];
    static const double testAges[testLen];
    static const double piecewiseConstValues[testLen];
    static const double linearInterpValues[testLen];
    OM::util::AgeGroupValues groups;
};

const double AgeGroupInterpolationSuite::stdLbounds[] = {
    0,5,10,15,60
};
const double AgeGroupInterpolationSuite::stdValues[] = {
    6.08,3.81,2.62,4.05,5.41
};
const double AgeGroupInterpolationSuite::testAges[] = {
    // various ages, designed to test limits, boundary points and interpolation
    15.2,18.09,7.0,2.5,
    0.0,20.0,900.0,62.0
};
const double AgeGroupInterpolationSuite::piecewiseConstValues[] = {
    4.0499999999999998, 4.0499999999999998, 3.8100000000000001, 6.0800000000000001,
    6.0800000000000001, 4.0499999999999998, 5.4100000000000001, 5.4100000000000001
};
const double AgeGroupInterpolationSuite::linearInterpValues[] = {
    2.7744400000000002, 2.9397479999999998, 4.0369999999999999, 6.0800000000000001,
    6.0800000000000001, 3.0489999999999999, 5.4100000000000001, 4.938533333333333
};

#endif
