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

#ifndef Hmod_LSTMTreatments
#define Hmod_LSTMTreatments

#include "Global.h"
#include "PkPd/Drug/LSTMDrugType.h"
#include "PkPd/LSTMMedicate.h"

using namespace std;

namespace scnXml{
    class Treatments;
    class Pharmacology;
    class Scenario;
}
namespace OM { namespace PkPd {
    
/** Information about treatment courses: schedules and dosage tables.
 * 
 * A schedule is a list of times to administrer a drug, along with base drug
 * doses. A dosage table is used to select a multiplier for this base dose
 * according to age (it is assumed for now that age is the only basis on which
 * dosages are selected). */
/** Treatment schedules and dosage tables read from XML (owned by PkPdState). */
struct TreatmentLibrary {
    vector<Schedule> schedules;
    map<string, size_t> scheduleNames;
    vector<DosageTable> dosages;
    map<string, size_t> dosagesNames;
};

/** Load treatment data from XML (load drug types first). */
void init( TreatmentLibrary& lib, DrugLibrary& drugs, const scnXml::Treatments& data );

/** Get the index of a named schedule. */
size_t findSchedule( TreatmentLibrary& lib, const string& name );

/** Get the index of a named dosage table. */
size_t findDosages( TreatmentLibrary& lib, const string& name );

/** All pharmacology data read from XML: owned by Model, referenced by the
 * per-host LSTMModel. Not checkpointed (rebuilt from scenario data). */
struct PkPdState {
    /// Load drugs and treatments, if the scenario has a pharmacology element.
    void init( const scnXml::Scenario& scenario );
    
    DrugLibrary drugs;
    TreatmentLibrary treatments;
};

} }
#endif
