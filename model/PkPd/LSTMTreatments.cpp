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

#include "PkPd/LSTMTreatments.h"
#include "PkPd/LSTMMedicate.h"
#include "PkPd/Drug/LSTMDrugType.h"
#include "PkPd/LSTMModel.h"

namespace OM {
namespace PkPd {

void TreatmentRegistry::clear(){
    schedules.clear();
    scheduleNames.clear();
    dosages.clear();
    dosagesNames.clear();
}

inline void invalidModel(){
    throw util::xml_scenario_error( "use of PK/PD treatment requires use of a compatible infection model and specification of pharmacology/treatments" );
}
size_t TreatmentRegistry::findSchedule(const string& name){
    if( scheduleNames.size() == 0 ) invalidModel();
    auto iter = scheduleNames.find( name );
    if( iter == scheduleNames.end() ){
        throw util::xml_scenario_error(string("no treatment schedule with this name: ")
            .append(name));
    }
    return iter->second;
}

size_t TreatmentRegistry::findDosages(const string& name){
    if( dosagesNames.size() == 0 ) invalidModel();
    auto iter = dosagesNames.find( name );
    if( iter == dosagesNames.end() ){
        throw util::xml_scenario_error(string("no dosage table with this name: ")
            .append(name));
    }
    return iter->second;
}

}
}
