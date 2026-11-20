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

#include "Clock.h"
#include "RunContext.h"
#include "interventions/IvState.h"
#include "interventions/Components.h"
#include "interventions/ITN.h"
#include "interventions/IRS.h"
#include "interventions/GVI.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include <cstdlib>

namespace OM { namespace interventions {

void iv_deploy( RunContext& ctx, Host::Human& human, const ComponentConfig& params )
{
    // This adds per-host per-intervention details to the host's data set.
    // Records are never removed: they can contain per-host het samples, and
    // the list order is the effect-product order (parity: c).
    auto& iv = human.cols->iv[human.row];
    for( auto& r : iv ) {
        if( r.id == params.id ) {
            // already have a deployment for that description; just update it
            LocalRng& rng = human.rng();
            switch( r.kind ) {
            case Component::ITN: {
                const ITNParams& p = *ctx.itnComponents[r.id.id];
                r.deploy = ctx.clock.nowOrTs1();
                r.disposal = ctx.clock.nowOrTs1() + p.attritionOfNets->sampleAgeOfDecay(rng);
                r.n_holes = 0;
                r.hole_index = 0.0;
                // this is sampled independently: initial insecticide content
                // doesn't depend on handling
                r.initial_insecticide = p.initialInsecticide.sample(rng);
                if( r.initial_insecticide < 0.0 )
                    r.initial_insecticide = 0.0;      // avoid negative samples
                if( r.initial_insecticide > p.maxInsecticide )
                    r.initial_insecticide = p.maxInsecticide;
                break; }
            case Component::IRS:
                r.deploy = ctx.clock.nowOrTs1();
                r.initial_insecticide = ctx.irsComponents[r.id.id]->sampleInitialInsecticide(rng);
                break;
            default:  // GVI: only the timestamp resets
                r.deploy = ctx.clock.nowOrTs1();
            }
            return;
        }
    }
    // no deployment for that description: append a new record; first-use
    // draws in the old makeHumanPart ctor order (§9 pinned):
    IvRecord r;
    r.id = params.id;
    r.kind = params.kind;
    r.deploy = ctx.clock.now();   // old PerHostInterventionData base ctor
    LocalRng& rng = human.rng();
    switch( r.kind ) {
    case Component::ITN: {
        const ITNParams& p = *ctx.itnComponents[r.id.id];
        // Net rips and insecticide loss are assumed to co-vary dependent on
        // handling of net. They are sampled once per human: human handling
        // is presumed to be the largest cause of variance.
        util::NormalSample x = util::NormalSample::generate(rng);
        r.hole_rate = p.holeRate->sample(x);
        r.rip_rate = p.ripRate->sample(x);
        r.decay = p.insecticideDecay.get();
        r.decay_factor = p.insecticideDecay->sample_het_factor(x);
        // Sample per-deployment variables as in redeploy:
        r.disposal = ctx.clock.now() + p.attritionOfNets->sampleAgeOfDecay(rng);
        // this is sampled independently: initial insecticide content doesn't
        // depend on handling
        r.initial_insecticide = p.initialInsecticide.sample(rng);
        if( r.initial_insecticide < 0.0 )
            r.initial_insecticide = 0.0;      // avoid negative samples
        if( r.initial_insecticide > p.maxInsecticide )
            r.initial_insecticide = p.maxInsecticide;
        break; }
    case Component::IRS: {
        // parity: draw order — old HumanIRS ctor member-init sampled
        // insecticide BEFORE the decay het (both per deployment / per human)
        const IRSParams& p = *ctx.irsComponents[r.id.id];
        r.initial_insecticide = p.sampleInitialInsecticide(rng);
        r.decay = p.insecticideDecay.get();
        r.decay_factor = p.insecticideDecay->sample_het_factor(rng);
        break; }
    case Component::GVI: {
        const GVIParams& p = *ctx.gviComponents[r.id.id];
        // Varience factor of decay is sampled once per human: human is
        // assumed to account for most variance.
        r.decay = p.decay.get();
        r.decay_factor = p.decay->sample_het_factor(rng);
        break; }
    default:
        abort();    // code error: only ITN/IRS/GVI deploy per-host records
    }
    iv.push_back( r );
}

void iv_update( RunContext& ctx, Host::Human& human )
{
    // Only nets mutate per step (old HumanIRS/HumanGVI update() were no-ops).
    for( auto& r : human.cols->iv[human.row] ) {
        if( r.kind != Component::ITN || r.deploy == sim::never() ) continue;
        // First use is at age 0 relative to ts0()
        if( ctx.clock.ts0() >= r.disposal ){
            r.deploy = sim::never();      // DISCARDED state (§9)
            human.removeFromCohort(r.id);
            continue;
        }
        int newHoles = 0;
        if( r.hole_rate > 0)
            // parity: e — per-update hole/rip Poisson draw counts frozen (≤1+1 draws per protected human-day; skip when rates 0) (plan §3e)
            newHoles = human.rng().poisson( r.hole_rate );

        r.n_holes += newHoles;
        r.hole_index += newHoles;

        double lambda = r.n_holes * r.rip_rate;
        if (lambda > 0)
            r.hole_index += ctx.itnComponents[r.id.id]->ripFactor * human.rng().poisson( lambda );
    }
}

// ———  per-record effects (bodies verbatim from the deleted Human* classes)  ———

namespace {
/// ITN/IRS remaining insecticide content (old getInsecticideContent).
inline double insecticide_content( RunContext& ctx, const IvRecord& r ) {
    SimTime age = ctx.clock.nowOrTs1() - r.deploy;  // implies age 1 TS on first use
    return r.initial_insecticide * r.decay->eval_with( age, r.decay_factor );
}
/// GVI effect survival factor (old getEffectSurvival).
inline double gvi_survival( RunContext& ctx, const IvRecord& r ) {
    SimTime age = ctx.clock.nowOrTs1() - r.deploy;  // implies age 1 TS on first use
    return r.decay->eval_with( age, r.decay_factor );
}
}

double iv_rel_attractiveness( RunContext& ctx, const IvRecord& r, size_t speciesIndex ) {
    switch( r.kind ) {
    case Component::ITN: {
        if( r.deploy == sim::never() ) return 1.0;
        const ITNParams& params = *ctx.itnComponents[r.id.id];
        const ITNParams::ITNAnopheles& anoph = params.species[speciesIndex];
        return anoph.relativeAttractiveness( r.hole_index, insecticide_content(ctx, r) );
    }
    case Component::IRS: {
        const IRSParams& params = *ctx.irsComponents[r.id.id];
        const IRSParams::IRSAnopheles& anoph = params.species[speciesIndex];
        double effect = anoph.relativeAttractiveness( insecticide_content(ctx, r) );
        return anoph.byProtection( effect );
    }
    case Component::GVI: {
        const GVIParams& params = *ctx.gviComponents[r.id.id];
        const GVIParams::GVIAnopheles& anoph = params.species[speciesIndex];
        double effect = (1.0 - anoph.deterrency * gvi_survival(ctx, r));
        return anoph.byProtection( effect );
    }
    default: return 1.0;   // unreachable: records are only ever ITN/IRS/GVI
    }
}

double iv_preprandial_survival( RunContext& ctx, const IvRecord& r, size_t speciesIndex ) {
    switch( r.kind ) {
    case Component::ITN: {
        if( r.deploy == sim::never() ) return 1.0;
        const ITNParams& params = *ctx.itnComponents[r.id.id];
        const ITNParams::ITNAnopheles& anoph = params.species[speciesIndex];
        return anoph.preprandialSurvivalFactor( r.hole_index, insecticide_content(ctx, r) );
    }
    case Component::IRS: {
        const IRSParams& params = *ctx.irsComponents[r.id.id];
        const IRSParams::IRSAnopheles& anoph = params.species[speciesIndex];
        double effect = anoph.preprandialSurvivalFactor( insecticide_content(ctx, r) );
        return anoph.byProtection( effect );
    }
    case Component::GVI: {
        const GVIParams& params = *ctx.gviComponents[r.id.id];
        const GVIParams::GVIAnopheles& anoph = params.species[speciesIndex];
        double effect = (1.0 - anoph.preprandialKilling * gvi_survival(ctx, r));
        return anoph.byProtection( effect );
    }
    default: return 1.0;
    }
}

double iv_postprandial_survival( RunContext& ctx, const IvRecord& r, size_t speciesIndex ) {
    switch( r.kind ) {
    case Component::ITN: {
        if( r.deploy == sim::never() ) return 1.0;
        const ITNParams& params = *ctx.itnComponents[r.id.id];
        const ITNParams::ITNAnopheles& anoph = params.species[speciesIndex];
        return anoph.postprandialSurvivalFactor( r.hole_index, insecticide_content(ctx, r) );
    }
    case Component::IRS: {
        const IRSParams& params = *ctx.irsComponents[r.id.id];
        const IRSParams::IRSAnopheles& anoph = params.species[speciesIndex];
        double effect = anoph.postprandialSurvivalFactor( insecticide_content(ctx, r) );
        return anoph.byProtection( effect );
    }
    case Component::GVI: {
        const GVIParams& params = *ctx.gviComponents[r.id.id];
        const GVIParams::GVIAnopheles& anoph = params.species[speciesIndex];
        double effect = (1.0 - anoph.postprandialKilling * gvi_survival(ctx, r));
        return anoph.byProtection( effect );
    }
    default: return 1.0;
    }
}

double iv_rel_fecundity( RunContext& ctx, const IvRecord& r, size_t speciesIndex ) {
    switch( r.kind ) {
    case Component::ITN: {
        if( r.deploy == sim::never() ) return 1.0;
        const ITNParams& params = *ctx.itnComponents[r.id.id];
        const ITNParams::ITNAnopheles& anoph = params.species[speciesIndex];
        return anoph.relFecundity( r.hole_index, insecticide_content(ctx, r) );
    }
    case Component::IRS: {
        const IRSParams& params = *ctx.irsComponents[r.id.id];
        const IRSParams::IRSAnopheles& anoph = params.species[speciesIndex];
        double effect = anoph.fecundityEffect( insecticide_content(ctx, r) );
        return anoph.byProtection( effect );
    }
    case Component::GVI: {
        const GVIParams& params = *ctx.gviComponents[r.id.id];
        const GVIParams::GVIAnopheles& anoph = params.species[speciesIndex];
        double effect = (1.0 - anoph.fecundityReduction * gvi_survival(ctx, r));
        return anoph.byProtection( effect );
    }
    default: return 1.0;
    }
}

} }
