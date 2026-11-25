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

// Intervention bases for human-vector interaction, split out of
// Transmission/PerHost.h (phase-5 cycle break): Transmission core now only
// forward-declares these and includes the leaf interventions/ComponentId.h;
// the intervention components (ITN/IRS/GVI) include this header — the
// dependency is one-way: interventions -> Transmission.
#ifndef Hmod_VectorIntervention
#define Hmod_VectorIntervention

#include "Global.h"
#include "interventions/Interfaces.h"
#include "util/random.h"

namespace OM {
namespace Host { class Human; }
namespace Transmission {

/** A base class for human vector intervention parameters (ITN/IRS/GVI).
 *
 * Phase 7.2a: the per-host state half (PerHostInterventionData and the
 * makeHumanPart/asXComponent virtuals) is gone — per-host state is the
 * IvRecord column (humans/HostColumns.h), consumed by interventions/IvState.h
 * free functions. */
class HumanVectorInterventionComponent : public interventions::HumanInterventionComponent {
public:
    virtual ~HumanVectorInterventionComponent() {}
    
    const HumanVectorInterventionComponent* asVectorIntervention() const override{
        return this;
    }
    
protected:
    explicit HumanVectorInterventionComponent(interventions::ComponentId id) :
            HumanInterventionComponent(id) {}
};

}
}
#endif
