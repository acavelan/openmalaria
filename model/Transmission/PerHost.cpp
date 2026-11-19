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

#include "Transmission/PerHost.h"
#include "Transmission/VectorIntervention.h"
#include "core/RunExtensions.h"
#include "Transmission/TransmissionModel.h"
#include "interventions/InterventionManager.h"
#include "util/errors.h"
#include "util/checkpoint.h"

namespace OM {
namespace Transmission {
using namespace OM::util;

// -----  PerHost static  -----

AgeGroupInterpolator PerHost::relAvailAge;


// -----  PerHost non-static -----

PerHost::PerHost () :
        outsideTransmission(false),
        relativeAvailabilityHet(numeric_limits<double>::signaling_NaN())
{
}

PerHost::~PerHost() = default;
PerHost::PerHost( PerHost&& ) = default;
PerHost& PerHost::operator= ( PerHost&& ) = default;

void PerHost::initialise (RunContext& ctx, LocalRng& rng, double availabilityFactor) {
    relativeAvailabilityHet = availabilityFactor;
    anophEntoAvailabilityRaw.resize(PerHostAnophParams::numSpecies(ctx));
    anophEntoAvailability.resize(PerHostAnophParams::numSpecies(ctx));
    anophProbMosqBiting.resize(PerHostAnophParams::numSpecies(ctx));
    anophProbMosqResting.resize(PerHostAnophParams::numSpecies(ctx));

    for(size_t i = 0; i < PerHostAnophParams::numSpecies(ctx); ++i) {
        const PerHostAnophParams& base = PerHostAnophParams::get(ctx, i);
        anophEntoAvailabilityRaw[i] = base.entoAvailability->sample(rng);
        anophEntoAvailability[i] = base.entoAvailabilityFactor * anophEntoAvailabilityRaw[i] * availabilityFactor;
        anophProbMosqBiting[i] = base.probMosqBiting.sample(rng);
        auto pRest1 = base.probMosqFindRestSite.sample(rng);
        auto pRest2 = base.probMosqSurvivalResting.sample(rng);
        anophProbMosqResting[i] = pRest1 * pRest2;
    }
}

void PerHost::update(Host::Human& human){
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        (*iter)->update(human);
    }
}

void PerHost::deployComponent( RunContext& ctx, LocalRng& rng, const HumanVectorInterventionComponent& params ){
    // This adds per-host per-intervention details to the host's data set.
    // This data is never removed since it can contain per-host heterogeneity samples.
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        if( (*iter)->id() == params.id() ){
            // already have a deployment for that description; just update it
            (*iter)->redeploy( rng, params );
            return;
        }
    }
    // no deployment for that description: must make a new one
    activeComponents.push_back( params.makeHumanPart(ctx, rng) );
}

// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over the time span now - sim::never()
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
// parity: c — history-ordered activeComponents FP product — multiply in vector order, not reassociated
double PerHost::entoAvailabilityHetVecItv (size_t species) const {
    double alpha_i = anophEntoAvailability[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        alpha_i *= (*iter)->relativeAttractiveness( species );
    }
    return alpha_i;
}
double PerHost::probMosqBiting (size_t species) const {
    double P_B_i = anophProbMosqBiting[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        P_B_i *= (*iter)->preprandialSurvivalFactor( species );
    }
    return P_B_i;
}
double PerHost::probMosqResting (size_t species) const {
    double pRest = anophProbMosqResting[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        pRest *= (*iter)->postprandialSurvivalFactor( species );
    }
    return pRest;
}

double PerHost::relMosqFecundity (size_t species) const {
    double relFecundity = 1.0;
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        relFecundity *= (*iter)->relFecundity( species );
    }
    return relFecundity;
}

bool PerHost::hasActiveInterv(interventions::Component::Type type) const{
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        if( (*iter)->isDeployed() ){
            if( interventions::InterventionManager::getComponent( *ctx, (*iter)->id() ).componentType() == type )
                return true;
        }
    }
    return false;
}

void PerHost::checkpointIntervs( ostream& stream ){
    activeComponents.size() & stream;
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        *iter & stream;
    }
}
void PerHost::checkpointIntervs( istream& stream ){
    size_t l;
    l & stream;
    validateListSize(l);
    activeComponents.clear();
    for( size_t i = 0; i < l; ++i ){
        interventions::ComponentId id( stream );
        try{
            const interventions::HumanInterventionComponent& gen_params =
                interventions::InterventionManager::getComponent( *ctx, id );   // may throw
            const HumanVectorInterventionComponent *params =
                gen_params.asVectorIntervention();   // declared seam (was dynamic_cast)
            if( params == 0 )
                throw util::base_exception( "" );       // see catch block below
            activeComponents.push_back( params->makeHumanPart( *ctx, stream, id ) );
        }catch( util::base_exception& e ){
            // two causes, both boil down to index being wrong
            throw util::checkpoint_error( "bad value in checkpoint file" );
        }
    }
}

}
}

// ———  PerHostAnophParams registry (ctx.ext->anophParams)  ———

namespace OM { namespace Transmission {

void PerHostAnophParams::init (RunContext& ctx, const scnXml::Mosq& mosq) {
    ctx.ext->anophParams.emplace_back(mosq);
}

size_t PerHostAnophParams::numSpecies (RunContext& ctx) {
    return ctx.ext->anophParams.size();
}

const PerHostAnophParams& PerHostAnophParams::get (RunContext& ctx, size_t species) {
    return ctx.ext->anophParams[species];
}

void PerHostAnophParams::setEntoAvailabilityFactor(RunContext& ctx, size_t species, double entoAvailability){
    ctx.ext->anophParams[species].entoAvailabilityFactor = entoAvailability;
}

void PerHostAnophParams::calcAvailabilityPercentiles(RunContext& ctx){
    auto& params = ctx.ext->anophParams;
    auto& entoAvailabilityPercentiles = ctx.ext->entoAvailabilityPercentiles;
    int nSamples = 100000;
    vector<double> samples(nSamples, 0.0);

    // Using a different seed to make sure old results will not change
    // Does it really matter?
    // parity: order — MasterRng(0,0) availability-percentile seed is load-bearing for goldens — do not “clean up”
    util::MasterRng seed(0, 0);
    seed.seed(0, 0);

    util::LocalRng rng(seed);

    // Sample
    for(int i=0; i<nSamples; i++)
    {
        double avail = 0.0;
        for(size_t s = 0; s < params.size(); ++s)
            avail += params[s].entoAvailabilityFactor * params[s].entoAvailability->sample(rng);
        samples[i] = avail;
    }

    // Sort samples
    sort(samples.begin(), samples.end());

    entoAvailabilityPercentiles.resize(100);

    // Calc percetiles threshold values
    for(int i=0; i<100; i++)
        entoAvailabilityPercentiles[i] = samples[int(i*nSamples/100)];
    entoAvailabilityPercentiles[0] = 0.0;
}

double PerHostAnophParams::getEntoAvailabilityPercentile(RunContext& ctx, int p){
    auto& entoAvailabilityPercentiles = ctx.ext->entoAvailabilityPercentiles;
    if(p >= 0 && p < 100)
        return entoAvailabilityPercentiles[p];
    else if(p == 100)
        return std::numeric_limits<double>::max();
    else throw util::xml_scenario_error("Availability percentiles in intervention specifications can only be in [0 and 100]");
}

} }
