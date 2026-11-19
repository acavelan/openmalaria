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

#include "core/RunContext.h"
#include <fstream>

#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/errors.h"



namespace OM { namespace util {

    void AgeGroupInterpolator::outputSamples( const SimClock& clock, const string name ){
        if( !util::CommandLine::option(util::CommandLine::SAMPLE_INTERPOLATIONS) ){
            return;
        }
        ofstream fstream( (name+".csv").c_str() );
        
        double max = sim::inYears(clock.max_human_age);
        for( double age = 0.0; age < max; age += 0.1 ){
            fstream << age << "," << this->eval( age ) << endl;
        }
    }

    // parity: old AgeGroupDummy methods (this file) — unset interpolator is a bug.
    [[noreturn]] void AgeGroupInterpolator::dummyThrow() const {
        throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
    }

    // parity: identical loop in both old subclasses (xml/util.cpp:538/:601)
    void AgeGroupInterpolator::scale( double factor ){
        if( kind == Kind::Dummy ) dummyThrow();
        for( auto& kv : data ) kv.second *= factor;
    }

    // parity: identical scan in both old subclasses (xml/util.cpp:544/:647);
    // first-wins on ties, min() start value kept (not lowest()).
    double AgeGroupInterpolator::firstGlobalMaximum() const {
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
} }
