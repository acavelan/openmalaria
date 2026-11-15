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

#ifndef Hmod_ESCaseManagement
#define Hmod_ESCaseManagement

#include "Global.h"
#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/WHInterface.h"

namespace scnXml {
    class HSEventScheduler;
}

namespace OM { namespace interventions { struct InterventionRegistry; } }
namespace OM { namespace Clinical {
using WithinHost::WHInterface;

struct ClinicalParams;

/** Load health system data from initial data or an intervention's data (both
 * from XML) into params. (Re)loads all data affected by this healthSystem
 * element. */
void loadHSEventScheduler( ClinicalParams& params, const scnXml::HSEventScheduler& esData,
    const interventions::InterventionRegistry& reg );

/** Run the case-management tree for a new event (medicates treatment). */
CMDTOut esExecute( const CMHostData& hostData, const ClinicalParams& params );

} }
#endif
