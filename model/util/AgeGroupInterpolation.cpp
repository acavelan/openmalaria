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



namespace OM { namespace util {

// ———  AgeGroupInterpolation derivatives  ———
    
    void AgeGroupInterpolation::outputSamples( const string name ){
        if( !util::CommandLine::option(util::CommandLine::SAMPLE_INTERPOLATIONS) ){
            return;
        }
        ofstream fstream( (name+".csv").c_str() );
        
        double max = sim::inYears(sim::maxHumanAge());
        for( double age = 0.0; age < max; age += 0.1 ){
            fstream << age << "," << this->eval( age ) << endl;
        }
    }
    
    /** Dummy object, for initialization. Potentially makes code safer.
     ********************************************/
    class AgeGroupDummy : public AgeGroupInterpolation
    {
    public:
        static AgeGroupDummy singleton;
        
        virtual double eval( double ageYears ) const {
            throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
        }
        virtual void scale( double factor ) {
            throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
        }
        virtual double firstGlobalMaximum() {
            throw TRACED_EXCEPTION_DEFAULT( "AgeGroupDummy" );
        }
    };
    AgeGroupDummy AgeGroupDummy::singleton;
    
    
    // -----  AgeGroupInterpolator  -----
    
    AgeGroupInterpolator::AgeGroupInterpolator() :
        obj(&AgeGroupDummy::singleton) {}
    
    void AgeGroupInterpolator::reset(){
        assert( obj != nullptr );  // should not do that
        if( obj != &AgeGroupDummy::singleton ){
            delete obj;
            obj = &AgeGroupDummy::singleton;
        }
    }
    bool AgeGroupInterpolator::isSet()    {
        return obj != &AgeGroupDummy::singleton;
    }
} }
