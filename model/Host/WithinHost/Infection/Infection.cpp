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
    m_kind(kind),
    // Descriptive samples its duration here; RNG draw point is kind-conditional
    m_duration(kind == Kind::Descriptive ? infectionDuration(ctx, rng) : sim::never()),
    notPrintedMDWarning(true)
{
    switch( kind ){
        case Kind::Descriptive:
            assert( ctx.clock.one_TS() == 5 );
            break;
        case Kind::Dummy:
            m_density = 16;	// increased by DH to avoid zeros in initialKappa
            break;
        case Kind::Empirical:
            ctorBody_Empirical(rng, growthRateMultiplier);
            break;
        case Kind::Molineaux:
            ctorBody_Molineaux(rng);
            break;
        case Kind::Penny:
            ctorBody_Penny(rng);
            break;
    }
}

Infection::Infection (RunContext& ctx, istream& stream, Kind kind) :
    ctx(ctx),
    m_startDate(sim::never()),
    m_kind(kind),
    m_duration(sim::never())
{
    m_startDate & stream;
    m_density & stream;
    m_cumulativeExposureJ & stream;
    m_genotype & stream;
    m_origin & stream;
    switch( kind ){
        case Kind::Descriptive: read_Descriptive(stream); break;
        case Kind::Dummy: break;
        case Kind::Empirical: read_Empirical(stream); break;
        case Kind::Molineaux: read_Molineaux(stream); break;
        case Kind::Penny: read_Penny(stream); break;
    }
}

Infection::~Infection() {}

// ———  time-step update  ———

bool Infection::update( LocalRng& rng, double survivalFactor, SimTime now, double body_mass ){
    SimTime bsAge = now - m_startDate - s_latentP;	// age of post-latent-period blood stage
    if( bsAge < sim::zero() )
        return false;	// latent period (liver stage) — don't do anything
    switch( m_kind ){
        case Kind::Dummy: return updateDensity_Dummy (rng, survivalFactor, bsAge, body_mass);
        case Kind::Empirical: return updateDensity_Empirical (rng, survivalFactor, bsAge, body_mass);
        case Kind::Molineaux: return updateDensity_Molineaux (rng, survivalFactor, bsAge, body_mass);
        case Kind::Penny: return updateDensity_Penny (rng, survivalFactor, bsAge, body_mass);
        default:
            assert( false );  // Descriptive is updated via determineDensities()
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
        case Kind::Descriptive: checkpoint_Descriptive (stream); break;
        case Kind::Dummy: break;
        case Kind::Empirical: checkpoint_Empirical (stream); break;
        case Kind::Molineaux: checkpoint_Molineaux (stream); break;
        case Kind::Penny: checkpoint_Penny (stream); break;
    }
}

} }
