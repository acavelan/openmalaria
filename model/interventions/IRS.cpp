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

#include "Global.h"
#include "core/RunExtensions.h"
#include "core/RunContext.h"
#include "interventions/IRS.h"
#include "Host/Human.h"
#include "util/errors.h"
#include "util/SpeciesIndexChecker.h"
#include "util/CommandLine.h"
#include "R_nmath/qnorm.h"
#include <cmath>

namespace OM { namespace interventions {

void IRSComponent::deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method, VaccineLimits )const{
    human.perHostTransmission.deployComponent(ctx, human.rng(), *this);
    mon::recordDeploy(mon::measure("nMassIRS"), mon::measure("nCtsIRS"), human, method);
}

Component::Type IRSComponent::componentType()const{ return Component::IRS; }

void IRSComponent::print_details( std::ostream& out )const{
    out << id().id << "\tIRS";
}

unique_ptr<PerHostInterventionData> IRSComponent::makeHumanPart(RunContext& ctx, LocalRng& rng) const{
    return unique_ptr<PerHostInterventionData>(new HumanIRS( ctx, rng, *this ));
}
unique_ptr<PerHostInterventionData> IRSComponent::makeHumanPart( RunContext& ctx, istream& stream, ComponentId id ) const{
    return unique_ptr<PerHostInterventionData>(new HumanIRS( ctx, stream, id ));
}

IRSComponent::IRSAnopheles::RelativeAttractiveness::RelativeAttractiveness() :
    lPF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{}
IRSComponent::IRSAnopheles::SurvivalFactor::SurvivalFactor() :
    BF( numeric_limits< double >::signaling_NaN() ),
    PF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() ),
    invBaseSurvival( numeric_limits< double >::signaling_NaN() )
{}
void IRSComponent::IRSAnopheles::SurvivalFactor::init1(){
    BF = 0.0;
    PF = 0.0;
    insecticideScaling = 0.0;
    invBaseSurvival = 1.0;
}
double IRSComponent::IRSAnopheles::RelativeAttractiveness::relativeAttractiveness(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lPF*insecticideComponent );
    assert( relAvail>=0.0 );
    return relAvail;
}
double IRSComponent::IRSAnopheles::SurvivalFactor::survivalFactor(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double killingEffect = BF + PF*insecticideComponent;
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    assert( survivalFactor >= 0.0 );
    assert( survivalFactor <= 1.0 );
    return survivalFactor;
}

double IRSComponent::sampleInitialInsecticide(LocalRng& rng) const{
    double sample = initialInsecticide.sample(rng);
    if( sample < 0.0 )
        sample = 0.0;       // avoid negative samples
    if( sample > maxInsecticide )
        sample = maxInsecticide;
    return sample;
}


// ———  per-human data  ———
HumanIRS::HumanIRS( RunContext& ctx, LocalRng& rng, const IRSComponent& params ) :
    PerHostInterventionData( ctx, params.id() ),
    initialInsecticide( params.sampleInitialInsecticide(rng) )
{
    // Varience factor of decay is sampled once per human: human is assumed
    // to account for most variance.
    insecticideDecayHet = params.insecticideDecay->hetSample(rng);
}

void HumanIRS::redeploy( LocalRng& rng, const OM::Transmission::HumanVectorInterventionComponent& params ) {
    deployTime = ctx.clock.nowOrTs1();
    const IRSComponent* irsParams = params.asIRSComponent();  // declared seam (was dynamic_cast)
    assert( irsParams != 0 );   // code error if this fails
    initialInsecticide = irsParams->sampleInitialInsecticide(rng);
}

void HumanIRS::update(Host::Human& human){
}

double HumanIRS::relativeAttractiveness(size_t speciesIndex) const{
    const IRSComponent& params = *ctx.ext->irsComponents[m_id.id];
    const IRSComponent::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.relativeAttractiveness( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

double HumanIRS::preprandialSurvivalFactor(size_t speciesIndex) const{
    const IRSComponent& params = *ctx.ext->irsComponents[m_id.id];
    const IRSComponent::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.preprandialSurvivalFactor( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

double HumanIRS::postprandialSurvivalFactor(size_t speciesIndex) const{
    const IRSComponent& params = *ctx.ext->irsComponents[m_id.id];
    const IRSComponent::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.postprandialSurvivalFactor( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}
double HumanIRS::relFecundity(size_t speciesIndex) const{
    const IRSComponent& params = *ctx.ext->irsComponents[m_id.id];
    const IRSComponent::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.fecundityEffect( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

void HumanIRS::checkpoint( ostream& stream ){
    deployTime & stream;
    initialInsecticide & stream;
    ctx.ext->irsComponents[m_id.id]->insecticideDecay->checkpointSample(insecticideDecayHet, stream);
}
HumanIRS::HumanIRS( RunContext& ctx, istream& stream, ComponentId id ) : PerHostInterventionData( ctx, id )
{
    deployTime & stream;
    initialInsecticide & stream;
    ctx.ext->irsComponents[m_id.id]->insecticideDecay->checkpointSample(insecticideDecayHet, stream);
}

} }
