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

#include <fstream>

#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "schema/healthSystem.h"

namespace OM { namespace util {

using Kind = AgeGroupInterpolator::Kind;

void AgeGroupInterpolator::set(
    const scnXml::AgeGroupValues& ageGroups, const char* eltName
){
    reset();
    // Type mostly equivalent to a std::string:
    const scnXml::AgeGroupValues::InterpolationOptional& interp = ageGroups.getInterpolation();
    if( !interp.present() || interp.get() == "linear" ){
        // Points in the middle of input age groups plus stabilisation points at
        // the ends, interpolated linearly.
        auto it = ageGroups.getGroup().begin();

        if( it == ageGroups.getGroup().end() )
            throw util::xml_scenario_error(string(eltName) + ": at least one age group required");

        // first lower-bound must be 0
        if( it->getLowerbound() != 0.0 )
            throw util::xml_scenario_error(string(eltName) + ": first lower-bound must be 0");

        // Our insert iterator (used for performance)
        auto pos = data.begin();
        assert(pos == data.end());

        double greatestLbound = 0.0;    // also last lbound: age groups must be ordered lowest to highest
        double lastValue = it->getValue();       // first value is repeated for const start
        for( ; it != ageGroups.getGroup().end(); ++it ){
            double lbound = it->getLowerbound();
            if( lbound >= greatestLbound ){
                double insBound = 0.5*(greatestLbound + lbound);
                pos = data.insert( pos, make_pair( insBound, lastValue ) );
                greatestLbound = lbound;
                lastValue = it->getValue();     // always insert previous value
            } else
                throw util::xml_scenario_error(string(eltName) + ": lower bound " + to_string(lbound) + " less than previous " + to_string(greatestLbound));
        }

        // add a point in middle of last age group (taking upper bound as max-age-years:
        data[ 0.5*(greatestLbound + sim::inYears(sim::maxHumanAge())) ] = lastValue;
        // Add first point at zero
        data[ 0.0 ] = data.begin()->second;
        // add a point at infinity to catch remaining values
        data[ numeric_limits<double>::infinity() ] = data.rbegin()->second;

        kind = Kind::PiecewiseLinear;
    }else if( interp.get() == "none" ){
        // Direct access to input age-group data (discontinuous).
        if( ageGroups.getGroup().size() == 0 ){
            throw util::xml_scenario_error(string(eltName) + ": at least one age group required" );
        }

        auto pos = data.begin();
        assert(pos == data.end());

        double greatestLbound = -1.0;
        for( const scnXml::Group& group : ageGroups.getGroup() ){
            double lbound = group.getLowerbound();
            if( lbound >= greatestLbound ){
                greatestLbound = lbound;
                pos = data.insert( pos, make_pair( lbound, group.getValue() ) );
            } else
                throw util::xml_scenario_error(string(eltName) + ": lower bound " + to_string(lbound) + " not greater than previous " + to_string(greatestLbound));
        }

        // first lower-bound must be 0
        if( data.begin()->first != 0.0 )
            throw util::xml_scenario_error(string(eltName) + ": first lower-bound must be 0");

        kind = Kind::PiecewiseConstant;
    }else
        throw util::xml_scenario_error(string("age group interpolation ") + interp.get() + " not implemented" );

    outputSamples( eltName );
}

double AgeGroupInterpolator::eval( double ageYears ) const {
    if( kind == Kind::None )
        throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
    assert ( ageYears >= 0.0 );
    auto it = data.upper_bound( ageYears );
    if( kind == Kind::PiecewiseConstant ){
        --it;   // previous point when ordered by age (it->first <= ageYears < (it+1)->first)
        return it->second;
    }
    // PiecewiseLinear: continuous, but with discontinuous derivative.
    assert( it != data.end() );   // guaranteed by the point at infinity
    //TODO: optimise by putting (f1-f0)/(a1-a0) into the array and checking usage
    // Get the first point with age greater than ageYears:
    double a1 = it->first;  // a1 > ageYears
    double f1 = it->second;
    --it;   // previous point when ordered by age
    double a0 = it->first;  // a0 ≤ ageYears
    double f0 = it->second;
    return (ageYears - a0) / (a1 - a0) * (f1 - f0) + f0;
}

void AgeGroupInterpolator::scale( double factor ){
    if( kind == Kind::None )
        throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
    for( auto it = data.begin(); it != data.end(); ++it ){
        it->second *= factor;
    }
}

double AgeGroupInterpolator::firstGlobalMaximum() const{
    if( kind == Kind::None )
        throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
    double age = 0.0, max = numeric_limits<double>::min();
    for( auto it = data.begin(); it != data.end(); ++it ){
        if( it->second > max ){
            max = it->second;
            age = it->first;
        }
    }
    return age;
}

void AgeGroupInterpolator::outputSamples( const string name ){
    if( !util::CommandLine::option(util::CommandLine::SAMPLE_INTERPOLATIONS) ){
        return;
    }
    ofstream fstream( (name+".csv").c_str() );
    
    double max = sim::inYears(sim::maxHumanAge());
    for( double age = 0.0; age < max; age += 0.1 ){
        fstream << age << "," << this->eval( age ) << endl;
    }
}

void AgeGroupInterpolator::reset(){
    kind = Kind::None;
    data.clear();
}

bool AgeGroupInterpolator::isSet(){
    return kind != Kind::None;
}

} }
