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

#include "schema/entomology.h"
#include "Transmission/PerHost.h"
#include "Transmission/VectorModel.h"
#include "interventions/InterventionManager.h"
#include "util/errors.h"
#include "util/checkpoint.h"

namespace OM {
namespace Transmission {
using namespace OM::util;

namespace {
// Kind dispatch for PerHostInterventionData (only vector kinds are stored).
using PHID = PerHostInterventionData;

std::unique_ptr<PHID> makeHumanPart(LocalRng& rng, const HumanVectorInterventionComponent& params, SimTime now)
{
    auto d = std::make_unique<PHID>(params.componentType(), params, params.id());
    switch (d->kind) {
    case interventions::Component::ITN: interventions::itnInitNew(*d, rng, params, now); break;
    case interventions::Component::IRS: interventions::irsInitNew(*d, rng, params, now); break;
    case interventions::Component::GVI: interventions::gviInitNew(*d, rng, params, now); break;
    default: break; // not a vector component: never deployed here
    }
    return d;
}
void redeploy(PHID& d, LocalRng& rng, const HumanVectorInterventionComponent& params, SimTime now)
{
    switch (d.kind) {
    case interventions::Component::ITN: interventions::itnRedeploy(d, rng, params, now); break;
    case interventions::Component::IRS: interventions::irsRedeploy(d, rng, params, now); break;
    case interventions::Component::GVI: interventions::gviRedeploy(d, rng, params, now); break;
    default: break;
    }
}
void updateInterv(PHID& d, Host::Human& human, SimTime ts0)
{
    switch (d.kind) {
    case interventions::Component::ITN: interventions::itnUpdate(d, human, ts0); break;
    case interventions::Component::IRS: interventions::irsUpdate(d, human, ts0); break;
    case interventions::Component::GVI: interventions::gviUpdate(d, human, ts0); break;
    default: break;
    }
}
double interv_relativeAttractiveness(const PHID& d, size_t species, SimTime now)
{
    switch (d.kind) {
    case interventions::Component::ITN: return interventions::itnRelativeAttractiveness(d, species, now);
    case interventions::Component::IRS: return interventions::irsRelativeAttractiveness(d, species, now);
    case interventions::Component::GVI: return interventions::gviRelativeAttractiveness(d, species, now);
    }
    return 1.0;
}
double interv_preprandialSurvivalFactor(const PHID& d, size_t species, SimTime now)
{
    switch (d.kind) {
    case interventions::Component::ITN: return interventions::itnPreprandialSurvivalFactor(d, species, now);
    case interventions::Component::IRS: return interventions::irsPreprandialSurvivalFactor(d, species, now);
    case interventions::Component::GVI: return interventions::gviPreprandialSurvivalFactor(d, species, now);
    }
    return 1.0;
}
double interv_postprandialSurvivalFactor(const PHID& d, size_t species, SimTime now)
{
    switch (d.kind) {
    case interventions::Component::ITN: return interventions::itnPostprandialSurvivalFactor(d, species, now);
    case interventions::Component::IRS: return interventions::irsPostprandialSurvivalFactor(d, species, now);
    case interventions::Component::GVI: return interventions::gviPostprandialSurvivalFactor(d, species, now);
    }
    return 1.0;
}
double interv_relFecundity(const PHID& d, size_t species, SimTime now)
{
    switch (d.kind) {
    case interventions::Component::ITN: return interventions::itnRelFecundity(d, species, now);
    case interventions::Component::IRS: return interventions::irsRelFecundity(d, species, now);
    case interventions::Component::GVI: return interventions::gviRelFecundity(d, species, now);
    }
    return 1.0;
}
void checkpointInterv(PHID& d, ostream& stream)
{
    switch (d.kind) {
    case interventions::Component::ITN: interventions::itnCheckpoint(d, stream); break;
    case interventions::Component::IRS: interventions::irsCheckpoint(d, stream); break;
    case interventions::Component::GVI: interventions::gviCheckpoint(d, stream); break;
    default: break;
    }
}
} // anonymous namespace

// -----  PerHost static  -----

PerHostAnophParams::PerHostAnophParams (const scnXml::Mosq& mosq) {
    if(mosq.getAvailability().getDistr() != "const" && (util::ModelOptions::option(NEGATIVE_BINOMIAL_MASS_ACTION) || util::ModelOptions::option(LOGNORMAL_MASS_ACTION) || util::ModelOptions::option(COMORB_TRANS_HET)))
        throw util::xml_scenario_error( "PerHostAnophParams::PerHostAnophParams(): ModelOptions NEGATIVE_BINOMIAL_MASS_ACTION, LOGNORMAL_MASS_ACTION, COMORB_TRANS_HET are not compatible with ento-availability heterogeneity distributions other than \"const\"");

    entoAvailability = util::createSampler(1.0, mosq.getAvailability());
    probMosqBiting.setParams( mosq.getMosqProbBiting() );
    probMosqFindRestSite.setParams( mosq.getMosqProbFindRestSite() );
    probMosqSurvivalResting.setParams( mosq.getMosqProbResting() );
}

// -----  HumanVectorInterventionComponent  -----

void HumanVectorInterventionComponent::deploy(Host::Human& human, mon::Deploy::Method method, interventions::VaccineLimits, SimTime now) const{
    human.perHostTransmission.deployComponent(human.rng, *this, now);
    switch (kind) {
    case interventions::Component::ITN: mon::recordDeploy(mon::measure("nMassITNs"), mon::measure("nEPI_ITNs"), human, method); break;
    case interventions::Component::IRS: mon::recordDeploy(mon::measure("nMassIRS"), mon::measure("nCtsIRS"), human, method); break;
    default: mon::recordDeploy(mon::measure("nMassGVI"), mon::measure("nCtsGVI"), human, method); break; // GVI
    }
}

void HumanVectorInterventionComponent::print_details(std::ostream& out) const{
    out << id().id << "\t" << (kind == interventions::Component::IRS ? "IRS" : (kind == interventions::Component::GVI ? "GVI" : "ITN"));
}

// -----  PerHost non-static -----

PerHost::PerHost () :
        outsideTransmission(false),
        relativeAvailabilityHet(numeric_limits<double>::signaling_NaN())
{
}

void PerHost::initialise (LocalRng& rng, double availabilityFactor, const VectorSpeciesParams& species_params,
                   const AgeGroupInterpolator& age_availability) {
    ageAvailability = &age_availability;
    relativeAvailabilityHet = availabilityFactor;
    anophEntoAvailabilityRaw.resize(species_params.numSpecies());
    anophEntoAvailability.resize(species_params.numSpecies());
    anophProbMosqBiting.resize(species_params.numSpecies());
    anophProbMosqResting.resize(species_params.numSpecies());

    for(size_t i = 0; i < species_params.numSpecies(); ++i) {
        const PerHostAnophParams& base = species_params.get(i);
        anophEntoAvailabilityRaw[i] = base.entoAvailability->sample(rng);
        anophEntoAvailability[i] = base.entoAvailabilityFactor * anophEntoAvailabilityRaw[i] * availabilityFactor;
        anophProbMosqBiting[i] = base.probMosqBiting.sample(rng);
        auto pRest1 = base.probMosqFindRestSite.sample(rng);
        auto pRest2 = base.probMosqSurvivalResting.sample(rng);
        anophProbMosqResting[i] = pRest1 * pRest2;
    }
}

void PerHost::update(Host::Human& human, SimTime ts0){
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        updateInterv(**iter, human, ts0);
    }
}

void PerHost::deployComponent( LocalRng& rng, const HumanVectorInterventionComponent& params, SimTime now ){
    // This adds per-host per-intervention details to the host's data set.
    // This data is never removed since it can contain per-host heterogeneity samples.
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        if( (*iter)->id == params.id() ){
            // already have a deployment for that description; just update it
            redeploy(**iter, rng, params, now );
            return;
        }
    }
    // no deployment for that description: must make a new one
    activeComponents.push_back( makeHumanPart(rng, params, now) );
}

// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over the time span now - sim::never()
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
double PerHost::entoAvailabilityHetVecItv (size_t species, SimTime now) const {
    double alpha_i = anophEntoAvailability[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        alpha_i *= interv_relativeAttractiveness(**iter, species, now );
    }
    return alpha_i;
}
double PerHost::probMosqBiting (size_t species, SimTime now) const {
    double P_B_i = anophProbMosqBiting[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        P_B_i *= interv_preprandialSurvivalFactor(**iter, species, now );
    }
    return P_B_i;
}
double PerHost::probMosqResting (size_t species, SimTime now) const {
    double pRest = anophProbMosqResting[species];
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        pRest *= interv_postprandialSurvivalFactor(**iter, species, now );
    }
    return pRest;
}

double PerHost::relMosqFecundity (size_t species, SimTime now) const {
    double relFecundity = 1.0;
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        relFecundity *= interv_relFecundity(**iter, species, now );
    }
    return relFecundity;
}

bool PerHost::hasActiveInterv(interventions::Component::Type type) const{
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        if( (*iter)->isDeployed() ){
            if( (*iter)->kind == type )
                return true;
        }
    }
    return false;
}

void PerHost::checkpointIntervs( ostream& stream ){
    activeComponents.size() & stream;
    for( auto iter = activeComponents.begin(); iter != activeComponents.end(); ++iter ){
        (*iter)->id & stream; // must be first; read externally so that the
                // correct kind can be found
        checkpointInterv(**iter, stream);
    }
}
void PerHost::load( istream& stream, const interventions::InterventionRegistry& reg ) {
    anophEntoAvailabilityRaw & stream;
    anophEntoAvailability & stream;
    anophProbMosqBiting & stream;
    anophProbMosqResting & stream;
    relativeAvailabilityHet & stream;
    outsideTransmission & stream;
    copulaGaussianSamples & stream;
    checkpointIntervs( stream, reg );
}
void PerHost::checkpointIntervs( istream& stream, const interventions::InterventionRegistry& reg ){
    size_t l;
    l & stream;
    validateListSize(l);
    activeComponents.clear();
    for( size_t i = 0; i < l; ++i ){
        interventions::ComponentId id( stream );
        try{
            const interventions::HumanInterventionComponent& gen_params =
                interventions::getComponent( reg, id );   // may throw
            const HumanVectorInterventionComponent *params =
                dynamic_cast<const HumanVectorInterventionComponent*>( &gen_params );
            if( params == 0 )
                throw util::base_exception( "" );       // see catch block below
            auto d = std::make_unique<PerHostInterventionData>( params->componentType(), *params, id );
            switch (d->kind) {
            case interventions::Component::ITN: interventions::itnRead(stream, *d); break;
            case interventions::Component::IRS: interventions::irsRead(stream, *d); break;
            case interventions::Component::GVI: interventions::gviRead(stream, *d); break;
            default: throw util::base_exception( "" );  // not a vector component
            }
            activeComponents.push_back( std::move(d) );
        }catch( util::base_exception& e ){
            // two causes, both boil down to index being wrong
            throw util::checkpoint_error( "bad value in checkpoint file" );
        }
    }
}

}
}
