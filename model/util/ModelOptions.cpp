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

#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/errors.h"

#include <sstream>
#include <iostream>

namespace OM { namespace util {
    std::bitset<NUM_OPTIONS> ModelOptions::options;
    

    std::bitset<NUM_OPTIONS> ModelOptions::getBaseModelOptions()
    {
        // The "base" model is very simple
        bitset<NUM_OPTIONS> allOptionsOff;

        return allOptionsOff;
    }

    std::bitset<NUM_OPTIONS> ModelOptions::getLegacyDefaultModelOptions()
    {
        bitset<NUM_OPTIONS> defaultOptSet;

        // Kept for legacy reasons.  That is, when an input XML does not specify a model name
        // and does not set these model options to false, they are turned on *for backwards
        // compatibility*.
        defaultOptSet.set (MAX_DENS_CORRECTION);
        defaultOptSet.set (INNATE_MAX_DENS);
        defaultOptSet.set (INDIRECT_MORTALITY_FIX);
        defaultOptSet.set (HEALTH_SYSTEM_MEMORY_FIX);

        return defaultOptSet;
    }

    
} }
