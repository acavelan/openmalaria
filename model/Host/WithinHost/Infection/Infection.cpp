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

#include "Context.h"
#include "Host/WithinHost/Infection/Infection.h"

namespace OM {
namespace WithinHost {

// ———  Dummy variant (was DummyInfection.cpp; simplest density model)  ———

Infection* createDummyInfection(Context& ctx, LocalRng& rng, uint32_t protID, int origin) {
    return new DummyInfection(ctx, protID, origin);
}

void DummyInfection::init(Context& ctx) { ctx.createInfection = &createDummyInfection; }

// no draws: the old infection_new Dummy case only set the initial density
DummyInfection::DummyInfection(Context& ctx, uint32_t genotype, int origin)
    : Infection(ctx, genotype, origin, Kind::Dummy) {
    m_density = 16; // increased by DH to avoid zeros in initialKappa
}

bool DummyInfection::update(LocalRng& rng, double survivalFactor, SimTime now, double) {
    SimTime bsAge = now - m_startDate - Infection::s_latentP; // age of post-latent-period blood stage
    if (bsAge < sim::zero())
        return false; // latent period (liver stage) — don't do anything
    const double GROWTH_RATE = 8.0;
    const double PARASITE_THRESHOLD = 1;

    m_density = (mod_nn(int(m_density * GROWTH_RATE), 20000)) * survivalFactor;
    m_cumulativeExposureJ += m_density;

    if (m_density < PARASITE_THRESHOLD) {
        return true;
    } else {
        return false;
    }
}

}
}
