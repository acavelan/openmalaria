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

#include <fstream>
#include <limits>
#include <map>

#include "Clock.h"
#include "RunContext.h"
#include "util/CommandLine.h"
#include "util/errors.h"

namespace scnXml {
    class AgeGroupValues;
}
namespace OM { class SimClock; }
namespace OM { namespace util {

/** Deterministic interpolation of data collected by age groups: one
 * enum-tagged value type (Kind switch picks the formula — the old
 * PiecewiseConstant / InterpolationPoints / PiecewiseLinear / Dummy class
 * hierarchy is gone). `data` stays a std::map<double,double>: numeric key
 * order is load-bearing (§3 caveat), insertion order is only an optimisation.
 *
 * Current version does not store an age index, thus an order log(n) lookup
 * must occur each time a value is looked up.
 ********************************************/
struct AgeGroupInterpolator
{
    enum class Kind { Dummy, PiecewiseConstant, Linear };

    /** Create (Dummy: set() must be called before further use). */
    AgeGroupInterpolator() : kind(Kind::Dummy) {}

    /** Set age-group interpolation data from an XML fragment
     * (definition: model/xml/util.cpp — the only code reading the schema).
     *
     * @param ageGroups XML element of per-age-group values
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    void set( const SimClock& clock, const scnXml::AgeGroupValues& ageGroups, const char* eltName );

    /** Set back to initial state. */
    void reset(){ kind = Kind::Dummy; data.clear(); }

    /// Return true if set() was ever called.
    bool isSet(){ return kind != Kind::Dummy; }

    /** Return a value interpolated for age ageYears. */
    inline double eval( double ageYears ) const {
        switch(kind) {
        case Kind::PiecewiseConstant: {  // parity: AgeGroupPiecewiseConstant::eval (xml/util.cpp:531)
            assert ( ageYears >= 0.0 );
            auto it = data.upper_bound( ageYears );
            --it;   // previous point when ordered by age (it->first <= ageYears < (it+1)->first)
            return it->second; }
        case Kind::Linear: {  // parity: AgeGroupPiecewiseLinear::eval (xml/util.cpp:633)
            assert ( ageYears >= 0.0 );
            // Get the first point with age greater than ageYears:
            auto it = data.upper_bound( ageYears );
            assert( it != data.end() );
            double a1 = it->first;  // a1 > ageYears
            double f1 = it->second;
            --it;   // previous point when ordered by age
            double a0 = it->first;  // a0 ≤ ageYears
            double f0 = it->second;
            return (ageYears - a0) / (a1 - a0) * (f1 - f0) + f0; }
        default:
            dummyThrow();
        }
        return 0.0;  // unreachable
    }

    /** Scale function by factor. */
    void scale( double factor ){
        if( kind == Kind::Dummy ) dummyThrow();
        for( auto& kv : data ) kv.second *= factor;
    }

    /** Find the youngest age which is the global maximum (i.e. the age at
     * which individuals are considered adults where all adults have the same
     * same value). */
    double firstGlobalMaximum() const {
        if( kind == Kind::Dummy ) dummyThrow();
        double age = 0.0, max = numeric_limits<double>::min();
        for( const auto& kv : data ){
            if( kv.second > max ){
                max = kv.second;
                age = kv.first;
            }
        }
        return age;
    }

private:
    /** Cold path: old AgeGroupDummy methods threw; keep one throw site. */
    [[noreturn]] void dummyThrow() const {
        throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
    }

    /** Sample interpolator between 0 and max age, outputting to a csv file
     * called name.csv. */
    void outputSamples( const SimClock& clock, const string name ){
        if( !util::CommandLine::option(util::CommandLine::SAMPLE_INTERPOLATIONS) ){
            return;
        }
        ofstream fstream( (name+".csv").c_str() );
        double max = sim::inYears(clock.max_human_age);
        for( double age = 0.0; age < max; age += 0.1 ){
            fstream << age << "," << this->eval( age ) << endl;
        }
    }

    Kind kind;
    /// Breakpoints (lower-age-bound or interpolation points) in numeric key order.
    map<double,double> data;
};

} }
