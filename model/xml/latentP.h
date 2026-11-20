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

// model/xml/: shared latentP read, identical block in scenario.cpp
// (initFalciparum) and vivax.cpp (initVivax); error text/order verbatim.

#ifndef Hmod_xml_latentP
#define Hmod_xml_latentP

#include "Global.h"
#include "core/RunContext.h"
#include "schema/scenario.h"

namespace OM {

/// Default 15d, else model/parameters/latentP; wraps parse errors as
/// "model/parameters/latentP: <msg>" exactly as both callers did.
SimTime read_latentP( const SimClock& clock, const scnXml::Model& model );

} // namespace OM
#endif
