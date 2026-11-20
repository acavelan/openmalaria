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

#include "Host/WithinHost/Infection/Infection.h"
#include "util/checkpoint.h"

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

Infection* infection_new_checkpointed(RunContext& ctx, istream& stream, Infection::Kind kind) {
    // the shared fields streamed below are overwritten before any use
    Infection* inf = new Infection(ctx, 0, 0, kind);
    inf->m_startDate & stream;
    inf->m_density & stream;
    inf->m_cumulativeExposureJ & stream;
    inf->m_genotype & stream;
    inf->m_origin & stream;
    switch( kind ){
        case Kind::Descriptive: read_descriptive(*inf, stream); break;
        case Kind::Dummy: break;
        case Kind::Empirical: read_empirical(*inf, stream); break;
        case Kind::Molineaux: read_molineaux(*inf, stream); break;
        case Kind::Penny: read_penny(*inf, stream); break;
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

// ———  checkpointing  ———

void infection_checkpoint (Infection& inf, ostream& stream) {
    inf.m_startDate & stream;
    inf.m_density & stream;
    inf.m_cumulativeExposureJ & stream;
    inf.m_genotype & stream;
    inf.m_origin & stream;
    switch( inf.kind ){
        case Kind::Descriptive: checkpoint_descriptive (inf, stream); break;
        case Kind::Dummy: break;
        case Kind::Empirical: checkpoint_empirical (inf, stream); break;
        case Kind::Molineaux: checkpoint_molineaux (inf, stream); break;
        case Kind::Penny: checkpoint_penny (inf, stream); break;
    }
}

} }
