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

namespace OM { namespace WithinHost {

// ———  construction  ———

Infection::Infection(RunContext& ctx, LocalRng& rng, uint32_t genotype, int origin, Kind kind,
        double growthRateMultiplier) :
    ctx(ctx),
    m_startDate(this->ctx.clock.nowOrTs0()),
    m_density(0.0),
    m_cumulativeExposureJ(0.0),
    m_genotype(genotype),
    m_origin(origin),
    m_kind(kind)
{
    switch( kind ){
        case Kind::Descriptive:
            // Descriptive samples its duration here; RNG draw point is
            // kind-conditional (was the ctor init list; nothing draws between
            // construction and this call).
            infection_ctor_descriptive(*this, rng);
            assert( ctx.clock.one_TS() == 5 );
            break;
        case Kind::Dummy:
            m_density = 16;	// increased by DH to avoid zeros in initialKappa
            break;
        case Kind::Empirical:
            infection_ctor_empirical(*this, rng, growthRateMultiplier);
            break;
        case Kind::Molineaux:
            infection_ctor_molineaux(*this, rng);
            break;
        case Kind::Penny:
            infection_ctor_penny(*this, rng);
            break;
    }
}

Infection::Infection (RunContext& ctx, istream& stream, Kind kind) :
    ctx(ctx),
    m_startDate(sim::never()),
    m_kind(kind)
{
    m_startDate & stream;
    m_density & stream;
    m_cumulativeExposureJ & stream;
    m_genotype & stream;
    m_origin & stream;
    switch( kind ){
        case Kind::Descriptive: read_descriptive(*this, stream); break;
        case Kind::Dummy: break;
        case Kind::Empirical: read_empirical(*this, stream); break;
        case Kind::Molineaux: read_molineaux(*this, stream); break;
        case Kind::Penny: read_penny(*this, stream); break;
    }
}

Infection::~Infection() {}

// ———  time-step update  ———

bool Infection::update( LocalRng& rng, double survivalFactor, SimTime now, double body_mass ){
    SimTime bsAge = now - m_startDate - s_latentP;	// age of post-latent-period blood stage
    if( bsAge < sim::zero() )
        return false;	// latent period (liver stage) — don't do anything
    switch( m_kind ){
        case Kind::Dummy: return update_density_dummy (*this, rng, survivalFactor, bsAge, body_mass);
        case Kind::Empirical: return update_density_empirical (*this, rng, survivalFactor, bsAge, body_mass);
        case Kind::Molineaux: return update_density_molineaux (*this, rng, survivalFactor, bsAge, body_mass);
        case Kind::Penny: return update_density_penny (*this, rng, survivalFactor, bsAge, body_mass);
        default:
            assert( false );  // Descriptive is updated via determine_densities_descriptive()
            return false;
    }
}

// ———  checkpointing  ———

void Infection::checkpoint (ostream& stream) {
    m_startDate & stream;
    m_density & stream;
    m_cumulativeExposureJ & stream;
    m_genotype & stream;
    m_origin & stream;
    switch( m_kind ){
        case Kind::Descriptive: checkpoint_descriptive (*this, stream); break;
        case Kind::Dummy: break;
        case Kind::Empirical: checkpoint_empirical (*this, stream); break;
        case Kind::Molineaux: checkpoint_molineaux (*this, stream); break;
        case Kind::Penny: checkpoint_penny (*this, stream); break;
    }
}

} }
