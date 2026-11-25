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

// Leaf vocabulary header: no OM includes beyond Clock.h. Transmission/PerHost.h
// may include this without pulling in the interventions module (phase-5 cycle
// break); everything else should keep using interventions/Interfaces.h.
#pragma once

#include "Clock.h"
#include <limits>

namespace OM { namespace interventions {

/** Enumeration of all available components, in the order that these should be
 * deployed in within a single intervention. */
namespace Component { enum Type {
    RECRUIT_ONLY,     // selection for a sub-population (without other effects)
    SCREEN,     // screening, e.g. as part of MSAT
    TREAT_SIMPLE,       // treatment using the simple model (e.g. for MDA)
    TREAT_PKPD, // treatment using the PK/PD model (e.g. MDA)
    CM_DT,      // execute a "case management" decision tree (e.g. MSAT or more complicated things)
    PEV,        // pre-erythrocytic vaccine
    BSV,        // blood-stage vaccine
    TBV,        // transmission-blocking vaccine
    IPT,        // intermittent preventative treatment
    ITN,        // insecticide treated net
    IRS,        // indoor residual spraying
    GVI,        // generic vector intervention
    CLEAR_IMMUNITY,     // reset accumulated immunity to zero
}; }

/** Essentially just an integer, used as a vector index.
 * 
 * This class wraps the integer to guard against unintended conversions to or
 * from an integer (which could be mis-use). */
struct ComponentId{
    explicit inline ComponentId( size_t id ) : id(id) {}
    inline bool operator== (const ComponentId that) const{ return id == that.id; }
    inline bool operator< (const ComponentId that) const{ return id < that.id; }
    
    /// Special value indicating the whole population
    static ComponentId wholePop() {
        return ComponentId{ numeric_limits<uint32_t>::max() };
    }
    
    size_t id;
};

} }
