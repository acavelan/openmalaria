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

#include "RunContext.h"
#include "Host/WithinHost/Infection/Infection.h"

namespace OM { namespace WithinHost {

using Kind = Infection::Kind;

// ———  construction  ———

Infection* infection_new(RunContext& ctx, LocalRng& rng, uint32_t genotype, int origin,
        Infection::Kind kind, double growthRateMultiplier) {
    Infection* inf = new Infection(ctx, genotype, origin, kind);
    switch( kind ){
        case Kind::Descriptive:
            // Descriptive samples its duration here; RNG draw point is
            // kind-conditional (was the ctor init list; nothing else draws
            // between construction and this case).
            infection_ctor_descriptive(*inf, rng);
            assert( ctx.clock.one_TS() == 5 );
            break;
        case Kind::Dummy:
            inf->m_density = 16;	// increased by DH to avoid zeros in initialKappa
            break;
        case Kind::Empirical:
            infection_ctor_empirical(*inf, rng, growthRateMultiplier);
            break;
        case Kind::Molineaux:
            infection_ctor_molineaux(*inf, rng);
            break;
        case Kind::Penny:
            infection_ctor_penny(*inf, rng);
            break;
    }
    return inf;
}

// ———  time-step update  ———

bool infection_update( Infection& inf, LocalRng& rng, double survivalFactor, SimTime now, double body_mass ){
    SimTime bsAge = now - inf.m_startDate - Infection::s_latentP;	// age of post-latent-period blood stage
    if( bsAge < sim::zero() )
        return false;	// latent period (liver stage) — don't do anything
    switch( inf.kind ){
        case Kind::Dummy: return update_density_dummy (inf, rng, survivalFactor, bsAge, body_mass);
        case Kind::Empirical: return update_density_empirical (inf, rng, survivalFactor, bsAge, body_mass);
        case Kind::Molineaux: return update_density_molineaux (inf, rng, survivalFactor, bsAge, body_mass);
        case Kind::Penny: return update_density_penny (inf, rng, survivalFactor, bsAge, body_mass);
        default:
            assert( false );  // Descriptive is updated via determine_densities_descriptive()
            return false;
    }
}

// ———  Dummy kind (was DummyInfection.cpp; simplest density model)  ———

Infection* createDummyInfection (RunContext& ctx, LocalRng& rng, uint32_t protID, int origin) {
    return infection_new (ctx, rng, protID, origin, Infection::Kind::Dummy);
}

void infection_init_dummy (RunContext& ctx) {
    ctx.createInfection = &createDummyInfection;
}

bool update_density_dummy( Infection& inf, LocalRng& rng, double survivalFactor, SimTime, double ){
    const double GROWTH_RATE = 8.0;
    const double PARASITE_THRESHOLD = 1;
    
    inf.m_density = (mod_nn(int(inf.m_density*GROWTH_RATE), 20000)) * survivalFactor;
    inf.m_cumulativeExposureJ += inf.m_density;
    
    if (inf.m_density < PARASITE_THRESHOLD) {
        return true;
    } else {
        return false;
    }
}

} }
