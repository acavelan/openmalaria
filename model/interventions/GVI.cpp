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
#include "interventions/GVI.h"
#include "schema/interventions.h"
#include "Host/Human.h"
#include "util/SpeciesIndexChecker.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include <cmath>

namespace OM { namespace interventions {

GVIComponent::GVIComponent( ComponentId id, const scnXml::GVIDescription& elt,
        const map<string,size_t>& species_name_map ) :
        Transmission::HumanVectorInterventionComponent(id)
{
    decay = DecayFunction::makeObject( elt.getDecay(), "interventions.human.vector.decay" );
    // assume usage modifier is 100% if none is specified
    double propUse;
    if (elt.getUsage().present()) {
        propUse = elt.getUsage().get().getValue();
    }
    else {
        propUse = 1.0;
    }
    if( !( propUse >= 0.0 && propUse <= 1.0 ) ){
        throw util::xml_scenario_error("GVI.description.usage: must be within range [0,1]");
    }
    
    typedef scnXml::GVIDescription::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    species.resize( species_name_map.size() );
    util::SpeciesIndexChecker checker( "GVI intervention", species_name_map );
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ) {
        species[checker.getIndex(it->getMosquito())].init (*it, propUse);
    }
    checker.checkNoneMissed();
}

void GVIComponent::deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits, SimTime now )const{
    human.perHostTransmission.deployComponent(human.rng, *this, now);
    mon::recordDeploy(mon::measure("nMassGVI"), mon::measure("nCtsGVI"), human, method);
}

Component::Type GVIComponent::componentType()const{ return Component::GVI; }

void GVIComponent::print_details( std::ostream& out )const{
    out << id().id << "\tGVI";
}

void GVIComponent::GVIAnopheles::init(const scnXml::GVIDescription::AnophelesParamsType& elt,
                                          double proportionUse)
{
    assert( (std::isnan)(deterrency) ); // double init
    deterrency = elt.getDeterrency().present() ? elt.getDeterrency().get().getValue() : 0.0;
    preprandialKilling = elt.getPreprandialKillingEffect().present() ? elt.getPreprandialKillingEffect().get().getValue() : 0.0;
    postprandialKilling = elt.getPostprandialKillingEffect().present() ? elt.getPostprandialKillingEffect().get().getValue() : 0.0;
    fecundityReduction = elt.getFecundityReduction().present() ? elt.getFecundityReduction().get().getValue() : 0.0;
    // Simpler version of ITN usage/action:
    double propActive = elt.getPropActive();
    if(propActive != 1.0 && propActive != 0.0)
        throw util::xml_scenario_error("GVI.description.propActive: must be either 0 or 1");

    assert( proportionUse >= 0.0 && proportionUse <= 1.0 );
    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = proportionUse * propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}


} } // namespace OM::interventions

namespace OM {
namespace interventions {

// —————  per-host effect state: Transmission::PerHostInterventionData, kind GVI  —————
namespace {
const interventions::GVIComponent& componentOf(const Transmission::PerHostInterventionData& d) {
    return static_cast<const interventions::GVIComponent&>(d.component);
}
/// Survival factor of the effect (implies age 1 TS on first use).
double effectSurvival(const Transmission::PerHostInterventionData& d, SimTime now) {
    return d.decayHet->eval(now - d.deployTime);
}
}

void gviInitNew(Transmission::PerHostInterventionData& d, LocalRng& rng, const Transmission::HumanVectorInterventionComponent&, SimTime /*now*/)
{
    d.deployTime = sim::now();
    // Varience factor of decay is sampled once per human: human is assumed
    // to account for most variance.
    d.decayHet = componentOf(d).decay->hetSample(rng);
}

void gviRedeploy(Transmission::PerHostInterventionData& d, LocalRng&, const Transmission::HumanVectorInterventionComponent&, SimTime now) {
    d.deployTime = now;
}

void gviUpdate(Transmission::PerHostInterventionData&, Host::Human&, SimTime /*ts0*/){
}

double gviRelativeAttractiveness(const Transmission::PerHostInterventionData& d, size_t speciesIndex, SimTime now){
    const interventions::GVIComponent::GVIAnopheles& anoph = componentOf(d).species[speciesIndex];
    double effect = (1.0 - anoph.deterrency * effectSurvival(d, now));
    return anoph.byProtection( effect );
}

double gviPreprandialSurvivalFactor(const Transmission::PerHostInterventionData& d, size_t speciesIndex, SimTime now){
    const interventions::GVIComponent::GVIAnopheles& anoph = componentOf(d).species[speciesIndex];
    double effect = (1.0 - anoph.preprandialKilling * effectSurvival(d, now));
    return anoph.byProtection( effect );
}

double gviPostprandialSurvivalFactor(const Transmission::PerHostInterventionData& d, size_t speciesIndex, SimTime now){
    const interventions::GVIComponent::GVIAnopheles& anoph = componentOf(d).species[speciesIndex];
    double effect = (1.0 - anoph.postprandialKilling * effectSurvival(d, now));
    return anoph.byProtection( effect );
}

double gviRelFecundity(const Transmission::PerHostInterventionData& d, size_t speciesIndex, SimTime now){
    const interventions::GVIComponent::GVIAnopheles& anoph = componentOf(d).species[speciesIndex];
    double effect = (1.0 - anoph.fecundityReduction * effectSurvival(d, now));
    return anoph.byProtection( effect );
}

void gviCheckpoint(Transmission::PerHostInterventionData& d, ostream& stream){
    d.deployTime & stream;
    componentOf(d).decay->checkpointSample(d.decayHet, stream);
}
void gviRead(istream& stream, Transmission::PerHostInterventionData& d){
    d.deployTime & stream;
    componentOf(d).decay->checkpointSample(d.decayHet, stream);
}

} }
