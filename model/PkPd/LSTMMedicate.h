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

/*
 * This file contains code internal to the LSTM PK/PD model to do with
 * dosing/scheduling.
 */
#ifndef Hmod_LSTMMedicate
#define Hmod_LSTMMedicate

#include "Clock.h"
#include "PkPd/LSTMModel.h"
#include "util/errors.h"

#include <vector>
#include <map>

using std::size_t;
using std::vector;
using std::map;

class UnittestUtil;
namespace scnXml { class Treatments; }
namespace OM { namespace PkPd {

struct Schedule {
    vector<MedicateData> medications;
};

struct DosageTable {
    /** Get dosage multilier from age or body mass.
     * 
     * Dosings may be given as a table using age or body mass as the key (first
     * column), or dose may be specified as mg drug / kg body mass. */
    double getMultiplier( double key ){
        if( multMassKg ) return key;
        else{
            auto iter = table.upper_bound( key );
            if( iter == table.end() )
                throw TRACED_EXCEPTION( "bad age/dosage table", util::Error::PkPd );
            return iter->second;
        }
    }
    
    bool useMass;       // false: dosing by age; true: dosing by body mass
    bool multMassKg;    // multiply by mass instead of using table
    map<double,double> table;
};

/** Treatment courses of the current run (was: LSTMTreatments statics). */
inline void invalidModel(){
    throw util::xml_scenario_error( "use of PK/PD treatment requires use of a compatible infection model and specification of pharmacology/treatments" );
}

struct TreatmentRegistry {
    /** Clear previous data. Only needed for testing. */
    void clear(){
        schedules.clear();
        scheduleNames.clear();
        dosages.clear();
        dosagesNames.clear();
    }
    /** Get the index of a named schedule. */
    size_t findSchedule( const std::string& name ){
        if( scheduleNames.size() == 0 ) invalidModel();
        auto iter = scheduleNames.find( name );
        if( iter == scheduleNames.end() ){
            throw util::xml_scenario_error(string("no treatment schedule with this name: ")
                .append(name));
        }
        return iter->second;
    }
    /** Get the index of a named dosage table. */
    size_t findDosages( const std::string& name ){
        if( dosagesNames.size() == 0 ) invalidModel();
        auto iter = dosagesNames.find( name );
        if( iter == dosagesNames.end() ){
            throw util::xml_scenario_error(string("no dosage table with this name: ")
                .append(name));
        }
        return iter->second;
    }
    
    vector<Schedule> schedules;
    map<string,size_t> scheduleNames;
    vector<DosageTable> dosages;
    map<string,size_t> dosagesNames;
};

/** Information about treatment courses: schedules and dosage tables. * 
 * A schedule is a list of times to administrer a drug, along with base drug
 * doses. A dosage table is used to select a multiplier for this base dose
 * according to age (it is assumed for now that age is the only basis on which
 * dosages are selected). */
class LSTMTreatments {
public:
    /** Load treatment data from XML into ctx.treatments.
     * 
     * Load drug data (LSTMDrugType::init()) first. */
    static void init (RunContext& ctx, const scnXml::Treatments& data);
    
private:
    LSTMTreatments();  // not constructible
};

} }
#endif
