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

#ifndef OM_INTERVENTIONS_ITN
#define OM_INTERVENTIONS_ITN

#include "Transmission/PerHost.h"

#include <map>
#include <memory>
#include <string>

namespace scnXml {
    class ITNDescription;
}

namespace OM {
namespace interventions {

/** Build the ITN parameters from XML. Defined in ITN.cpp (parses the
 * description and the per-species response curves). */
std::unique_ptr<Transmission::HumanVectorInterventionComponent> parseITNComponent(
    ComponentId id, const scnXml::ITNDescription& elt,
    const std::map<std::string, size_t>& species_name_map);

} }

#endif
