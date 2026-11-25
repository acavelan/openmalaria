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

// model/xml/ conversion layer — interventions module (start): imported
// infection rates. Moved verbatim from Host/ImportedInfections.h (the
// inline ImportedInfections::init became an out-of-line definition).

#include "Global.h"
#include "Host/ImportedInfections.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "schema/interventions.h"

#include <algorithm>
#include <string>

using namespace std;

namespace OM {
namespace Host {

void ImportedInfections::init( const scnXml::ImportedInfections& iiElt )
        {
            const scnXml::ImportedInfections::TimedType& tElt = iiElt.getTimed();
            try{
                //NOTE: if changing XSD, this should not have a default unit:
                period = UnitParse::readDuration( tElt.getPeriod(), UnitParse::STEPS );
                if( period < sim::zero() ){
                    throw util::format_error( "cannot be negative" );
                }
            }catch( const util::format_error& e ){
                throw util::xml_scenario_error( string("interventions/importedInfections/timed/period: ").append(e.message()) );
            }
            rate.reserve( tElt.getRate().size() );
            try{
                for( auto it = tElt.getRate().begin(); it != tElt.getRate().end(); ++it ){
                    SimTime date = UnitParse::readDate( it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/ );
                    // convert to per-time-step, per-person
                    double rateVal = it->getValue() * sim::yearsPerStep() * (1.0 / 1000.0);
                    rate.push_back( Rate( date - sim::startDate(), rateVal ) );
                }
                sort( rate.begin(), rate.end() );
                if( rate.size() > 0 ){
                    if( period != sim::zero() && rate[0].time != sim::zero() ){
                        throw util::xml_scenario_error( "must specify rate at time zero when period is not zero" );
                    }
                    // remove useless repeated entries from list
                    double lastRate = rate[0].value;
                    for( size_t i = 1; i < rate.size(); ){
                        if( rate[i].value == lastRate ){
                            rate.erase( rate.begin() + i );
                        }else{
                            lastRate = rate[i].value;
                            ++i;
                        }
                    }
                }
            }catch( const util::format_error& e ){
                throw util::xml_scenario_error( string("interventions/importedInfections/timed/time: ").append(e.message()) );
            }
        }

}
}
