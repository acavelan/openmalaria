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

#include "core/RunContext.h"
#include "core/RunExtensions.h"
#include "interventions/ITN.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/SpeciesIndexChecker.h"
#include "Host/Human.h"
#include "R_nmath/qnorm.h"
#include "util/CommandLine.h"
#include <cmath>

namespace OM { namespace interventions {


// —————  utility classes (internal use only)  —————


void factors::SurvivalFactor::init1() {
    useLogitEqns = false;
    a.invBaseSurvival = 1.0;
    a.BF = 0.0;
    a.HF = 0.0;
    a.PF = 0.0;
    a.IF = 0.0;
    a.insecticideScaling = 0.0;
    a.holeScaling = 0.0;
}

double factors::SurvivalFactor::rel_pAtt( double holeIndex, double insecticideContent )const {
    if (!useLogitEqns) {
        const double holeComponent = exp(-holeIndex * a.holeScaling);
        const double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        const double pAtt = a.BF
                + a.HF * holeComponent
                + a.PF * insecticideComponent
                + a.IF * holeComponent * insecticideComponent;
        assert( pAtt >= 0.0 );
        return pAtt / a.BF;
    } else {
        const double holeComponent = min(log(holeIndex + 1.0), b.hMax);
        const double insecticideComponent = log(insecticideContent + 1.0);
        const double x = exp(b.BF
                + b.HF * holeComponent
                + b.PF * insecticideComponent
                + b.IF * holeComponent * insecticideComponent);
        const double pAtt = x / (x + 1.0);
        
        // By construction, this is non-negative:
        return pAtt * b.invBaseEffect;
    }
}

double factors::SurvivalFactor::survivalFactor( double holeIndex, double insecticideContent )const {
    if (!useLogitEqns) {
        double holeComponent = exp(-holeIndex * a.holeScaling);
        double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        double killingEffect = a.BF
                + a.HF * holeComponent
                + a.PF * insecticideComponent
                + a.IF * holeComponent * insecticideComponent;
        
        double survivalFactor = (1.0 - killingEffect) * a.invBaseSurvival;
        assert( killingEffect <= 1.0 );
        // survivalFactor might be out of bounds due to precision error, see #49
        if (survivalFactor < 0.0)
            return 0.0;
        else if (survivalFactor > 1.0)
            return 1.0;
        return survivalFactor;
    } else {
        const double holeComponent = min(log(holeIndex + 1.0), b.hMax);
        const double insecticideComponent = log(insecticideContent + 1.0);
        const double x = exp(b.BF
                + b.HF * holeComponent
                + b.PF * insecticideComponent
                + b.IF * holeComponent * insecticideComponent);
        // K = x/(x+1), so (1-K) = 1/(x+1):
        const double surv = 1.0 / (x + 1.0);
        
        return surv * b.invBaseEffect;
    }
}


double factors::RelativeAttractiveness::relativeAttractiveness (
        double holeIndex, double insecticideContent) const
{
    if (model == SINGLE_STAGE) {
        double holeComponent = exp(-holeIndex * a.holeScaling);
        double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        double relAvail = exp(a.lHF * holeComponent
                + a.lPF * insecticideComponent
                + a.lIF * holeComponent * insecticideComponent );
        assert( relAvail>=0.0 );
        return relAvail;
    }
    
    double factor = numeric_limits<double>::quiet_NaN();
    if (model == TWO_STAGE) {
        // This is essentially a combination of the relative attractiveness as used
        // by IRS and a killing factor.
        
        // Note that an alternative, simpler, model could have been used, but was
        // not for consistency with other models. Alternative (here we don't take
        // the logarithm of PF):
        // pEnt = 1 - PFEntering × insecticideComponent
        
        const double insecticideComponent = 1.0 - exp(-insecticideContent * b.insecticideScalingEntering);
        const double pEnt = exp(b.lPFEntering * insecticideComponent);
        // In this model, effect with 0 insectice, pEnt0 = exp(0) = 1, hence we don't need to
        // divide by a denominator like in the logit model:
        factor = pEnt;
    } else {
        assert (model == TWO_STAGE_LOGIT);
        const double p = log(insecticideContent + 1.0);
        // We directly take the exponential (this is exp(logit.pEnt0):
        const double q = exp(c.entBaseFactor + c.entInsecticideFactor * p);
        const double pEnt = q / (q + 1.0);
        factor = pEnt * c.pEnt0Inv;
    }
    assert( factor >= 0.0 );
    
    // Note: b.pAttacking and c.pAttacking overlap, so we can use either:
    const double rel_pAtt = b.pAttacking.rel_pAtt( holeIndex, insecticideContent );
    return factor * rel_pAtt;
}


// —————  main, public classes  —————

void ITNComponent::deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method, VaccineLimits )const{
    human.perHostTransmission.deployComponent( ctx, human.rng, *this );
    mon::recordDeploy(mon::measure("nMassITNs"), mon::measure("nEPI_ITNs"), human, method);
}

Component::Type ITNComponent::componentType() const{
    return Component::ITN;
}

void ITNComponent::print_details( std::ostream& out )const{
    out << id().id << "\tITN";
}

unique_ptr<PerHostInterventionData> ITNComponent::makeHumanPart(RunContext& ctx, LocalRng& rng) const{
    return unique_ptr<PerHostInterventionData>(new HumanITN( ctx, rng, *this ));
}
unique_ptr<PerHostInterventionData> ITNComponent::makeHumanPart( RunContext& ctx, istream& stream, ComponentId id ) const{
    return unique_ptr<PerHostInterventionData>(new HumanITN( ctx, stream, id ));
}

HumanITN::HumanITN( RunContext& ctx, LocalRng& rng, const ITNComponent& params ) :
        PerHostInterventionData( ctx, params.id() ),
        nHoles( 0 ),
        holeIndex( 0.0 )
{
    // Net rips and insecticide loss are assumed to co-vary dependent on
    // handling of net. They are sampled once per human: human handling is
    // presumed to be the largest cause of variance.
    util::NormalSample x = util::NormalSample::generate(rng);
    holeRate = params.holeRate->sample(x);
    ripRate = params.ripRate->sample(x);
    insecticideDecayHet = params.insecticideDecay->hetSample(x);

    // Sample per-deployment variables as in redeploy:
    disposalTime = ctx.clock.now() + params.attritionOfNets->sampleAgeOfDecay(rng);
    // this is sampled independently: initial insecticide content doesn't depend on handling
    initialInsecticide = params.initialInsecticide.sample(rng);
    if( initialInsecticide < 0.0 )
        initialInsecticide = 0.0;       // avoid negative samples
    if( initialInsecticide > params.maxInsecticide )
        initialInsecticide = params.maxInsecticide;
}

void HumanITN::redeploy(LocalRng& rng, const OM::Transmission::HumanVectorInterventionComponent& params0) {
    const ITNComponent& params = *dynamic_cast<const ITNComponent*>(&params0);
    
    deployTime = ctx.clock.nowOrTs1();
    disposalTime = ctx.clock.nowOrTs1() + params.attritionOfNets->sampleAgeOfDecay(rng);
    nHoles = 0;
    holeIndex = 0.0;
    // this is sampled independently: initial insecticide content doesn't depend on handling
    initialInsecticide = params.initialInsecticide.sample(rng);
    if( initialInsecticide < 0.0 )
        initialInsecticide = 0.0;	// avoid negative samples
    if( initialInsecticide > params.maxInsecticide )
        initialInsecticide = params.maxInsecticide;
}

void HumanITN::update(Host::Human& human){
    const ITNComponent& params = *ctx.ext->itnComponents[m_id.id];
    if( deployTime != sim::never() ){
        // First use is at age 0 relative to ts0()
        if( ctx.clock.ts0() >= disposalTime ){
            deployTime = sim::never();
            human.removeFromCohort(id());
            return;
        }
        
        int newHoles = 0;
        if( holeRate > 0)
            // parity: e — per-update hole/rip Poisson draw counts frozen (≤1+1 draws per protected human-day; skip when rates 0) (plan §3e)
            newHoles = human.rng.poisson( holeRate );

        nHoles += newHoles;
        holeIndex += newHoles;

        double lambda = nHoles * ripRate;
        if (lambda > 0)
            holeIndex += params.ripFactor * human.rng.poisson( lambda );
    }
}

double HumanITN::relativeAttractiveness(size_t speciesIndex) const{
    if( deployTime == sim::never() ) return 1.0;
    const ITNComponent& params = *ctx.ext->itnComponents[m_id.id];
    const ITNComponent::ITNAnopheles& anoph = params.species[speciesIndex];
    return anoph.relativeAttractiveness( holeIndex, getInsecticideContent(params) );
}

double HumanITN::preprandialSurvivalFactor(size_t speciesIndex) const{
    if( deployTime == sim::never() ) return 1.0;
    const ITNComponent& params = *ctx.ext->itnComponents[m_id.id];
    const ITNComponent::ITNAnopheles& anoph = params.species[speciesIndex];
    return anoph.preprandialSurvivalFactor( holeIndex, getInsecticideContent(params) );
}

double HumanITN::postprandialSurvivalFactor(size_t speciesIndex) const{
    if( deployTime == sim::never() ) return 1.0;
    const ITNComponent& params = *ctx.ext->itnComponents[m_id.id];
    const ITNComponent::ITNAnopheles& anoph = params.species[speciesIndex];
    return anoph.postprandialSurvivalFactor( holeIndex, getInsecticideContent(params) );
}
double HumanITN::relFecundity(size_t speciesIndex) const{
    if( deployTime == sim::never() ) return 1.0;
    const ITNComponent& params = *ctx.ext->itnComponents[m_id.id];
    const ITNComponent::ITNAnopheles& anoph = params.species[speciesIndex];
    return anoph.relFecundity( holeIndex, getInsecticideContent(params) );
}

void HumanITN::checkpoint( ostream& stream ){
    deployTime & stream;
    disposalTime & stream;
    nHoles & stream;
    holeIndex & stream;
    initialInsecticide & stream;
    holeRate & stream;
    ripRate & stream;
    ctx.ext->itnComponents[m_id.id]->insecticideDecay->checkpointSample(insecticideDecayHet, stream);
}
HumanITN::HumanITN( RunContext& ctx, istream& stream, ComponentId id ) : PerHostInterventionData( ctx, id )
{
    deployTime & stream;
    disposalTime & stream;
    nHoles & stream;
    holeIndex & stream;
    initialInsecticide & stream;
    holeRate & stream;
    ripRate & stream;
    ctx.ext->itnComponents[m_id.id]->insecticideDecay->checkpointSample(insecticideDecayHet, stream);
}

} }
