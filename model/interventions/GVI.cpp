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
#include "core/RunContext.h"
#include "interventions/GVI.h"
#include "Host/Human.h"
#include "util/SpeciesIndexChecker.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include <cmath>

namespace OM { namespace interventions {

vector<GVIComponent*> GVIComponent::componentsByIndex;

void GVIComponent::deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method, VaccineLimits )const{
    human.perHostTransmission.deployComponent(ctx, human.rng, *this);
    mon::recordDeploy(mon::measure("nMassGVI"), mon::measure("nCtsGVI"), human, method);
}

Component::Type GVIComponent::componentType()const{ return Component::GVI; }

void GVIComponent::print_details( std::ostream& out )const{
    out << id().id << "\tGVI";
}

unique_ptr<PerHostInterventionData> GVIComponent::makeHumanPart(RunContext& ctx, LocalRng& rng) const{
    return unique_ptr<PerHostInterventionData>(new HumanGVI( ctx, rng, *this ));
}
unique_ptr<PerHostInterventionData> GVIComponent::makeHumanPart( RunContext& ctx, istream& stream, ComponentId id ) const{
    return unique_ptr<PerHostInterventionData>(new HumanGVI( ctx, stream, id ));
}


// ———  per-human data  ———
HumanGVI::HumanGVI ( RunContext& ctx, LocalRng& rng, const GVIComponent& params ) :
    PerHostInterventionData( ctx, params.id() )
{
    // Varience factor of decay is sampled once per human: human is assumed
    // to account for most variance.
    decayHet = params.decay->hetSample(rng);
}

void HumanGVI::redeploy(LocalRng& rng, const Transmission::HumanVectorInterventionComponent&) {
    deployTime = ctx.clock.nowOrTs1();
}

void HumanGVI::update(Host::Human& human){
}

double HumanGVI::relativeAttractiveness(size_t speciesIndex) const{
    const GVIComponent& params = *GVIComponent::componentsByIndex[m_id.id];
    const GVIComponent::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph.deterrency * getEffectSurvival(params));
    return anoph.byProtection( effect );
}

double HumanGVI::preprandialSurvivalFactor(size_t speciesIndex) const{
    const GVIComponent& params = *GVIComponent::componentsByIndex[m_id.id];
    const GVIComponent::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph.preprandialKilling * getEffectSurvival(params));
    return anoph.byProtection( effect );
}

double HumanGVI::postprandialSurvivalFactor(size_t speciesIndex) const{
    const GVIComponent& params = *GVIComponent::componentsByIndex[m_id.id];
    const GVIComponent::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph.postprandialKilling * getEffectSurvival(params));
    return anoph.byProtection( effect );
}
double HumanGVI::relFecundity(size_t speciesIndex) const{
    const GVIComponent& params = *GVIComponent::componentsByIndex[m_id.id];
    const GVIComponent::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph.fecundityReduction * getEffectSurvival(params));
    return anoph.byProtection( effect );
}

void HumanGVI::checkpoint( ostream& stream ){
    deployTime & stream;
    GVIComponent::componentsByIndex[m_id.id]->decay->checkpointSample(decayHet, stream);
}
HumanGVI::HumanGVI( RunContext& ctx, istream& stream, ComponentId id ) : PerHostInterventionData( ctx, id )
{
    deployTime & stream;
    GVIComponent::componentsByIndex[m_id.id]->decay->checkpointSample(decayHet, stream);
}

} }
