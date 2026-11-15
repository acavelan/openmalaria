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
#include "interventions/IRS.h"
#include "schema/interventions.h"
#include "Host/Human.h"
#include "util/errors.h"
#include "util/SpeciesIndexChecker.h"
#include "util/CommandLine.h"
#include "R_nmath/qnorm.h"
#include <cmath>

namespace OM { namespace interventions {

IRSComponent::IRSComponent( ComponentId id, const scnXml::IRSDescription& elt,
        const map<string,size_t>& species_name_map ) :
        Transmission::HumanVectorInterventionComponent(id)
{
    initialInsecticide.setParams( elt.getInitialInsecticide() );
    const double maxProp = 0.999;       //NOTE: this could be exposed in XML, but probably doesn't need to be
    maxInsecticide = R::qnorm5(maxProp, initialInsecticide.getMu(), initialInsecticide.getSigma(), true, false);
    insecticideDecay = DecayFunction::makeObject( elt.getInsecticideDecay(),
                                                  "interventions.human.IRS.description.insecticideDecay" );
    // assume usage modifier is 100% if none is specified
    double propUse;
    if (elt.getUsage().present()) {
        propUse = elt.getUsage().get().getValue();
    }
    else {
        propUse = 1.0;
    }
    if( !( propUse >= 0.0 && propUse <= 1.0 ) ){
        throw util::xml_scenario_error("ITN.description.usage: must be within range [0,1]");
    }
    
    typedef scnXml::IRSDescription::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    species.resize( species_name_map.size() );
    util::SpeciesIndexChecker checker( "IRS intervention", species_name_map );
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ) {
        species[checker.getIndex(it->getMosquito())].init (*it, propUse, maxInsecticide);
    }
    checker.checkNoneMissed();
}

void IRSComponent::deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits, SimTime now )const{
    human.perHostTransmission.deployComponent(human.rng, *this, now);
    mon::recordDeploy(mon::measure("nMassIRS"), mon::measure("nCtsIRS"), human, method);
}

Component::Type IRSComponent::componentType()const{ return Component::IRS; }

void IRSComponent::print_details( std::ostream& out )const{
    out << id().id << "\tIRS";
}

// The IRS response curves are the shared ITN/IRS factor models with the
// hole terms zero and holeIndex fixed at 0 (see InterventionFactors.h).
// Evaluated at holeIndex=0 the shared formulas reduce bit-exactly to:
//   survival factor   (1−(BF+PF·p))/(1−BF),   p=1−exp(−content·scaling)
//   attractiveness    PF^p
static void initSurvival( factors::SurvivalFactor& sf,
                          const scnXml::IRSKillingEffect& elt,
                          bool postPrandial, double maxInsecticide )
{
    sf.useLogitEqns = false;
    sf.clamp = false;
    sf.a.BF = elt.getBaseFactor();
    sf.a.HF = 0.0;
    sf.a.PF = elt.getInsecticideFactor();
    sf.a.IF = 0.0;
    sf.a.holeScaling = 0.0;
    sf.a.insecticideScaling = elt.getInsecticideScalingFactor();
    sf.a.invBaseSurvival = 1.0 / (1.0 - sf.a.BF);
    if( !( sf.a.BF >= 0.0 && sf.a.BF < 1.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected baseFactor to be in range [0,1]";
        throw util::xml_scenario_error( msg.str() );
    }
    if( !(sf.a.insecticideScaling>=0.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected scaling factor to be non-negative";
        throw util::xml_scenario_error( msg.str() );
    }
    // As for nets: p is in [0,1], so the survival factor stays in [0,1] if
    // PF>=0 and BF+PF·pmax<=1 with pmax the largest possible p.
    double pmax = 1.0-exp(-maxInsecticide*sf.a.insecticideScaling);
    if( !( sf.a.PF >= 0.0 && sf.a.BF+sf.a.PF*pmax <= 1.0 ) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected insecticideFactor≥0, baseFactor+"
            <<pmax<<"×insecticideFactor≤1";
        throw util::xml_scenario_error( msg.str() );
    }
}

static void initDeterrence( factors::RelativeAttractiveness& ra,
                            const scnXml::IRSDeterrency& elt )
{
    ra.model = factors::RelativeAttractiveness::SINGLE_STAGE;
    double PF = elt.getInsecticideFactor();
    ra.a.insecticideScaling = elt.getInsecticideScalingFactor();
    if( !( PF > 0.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams.relativeAvailability: expected insecticideFactor to be positive.";
        throw util::xml_scenario_error( msg.str() );
    }
    ra.a.lHF = 0.0;
    ra.a.lIF = 0.0;
    ra.a.holeScaling = 0.0;
    ra.a.lPF = log( PF );
}

void IRSComponent::IRSAnopheles::init(
    const scnXml::IRSDescription::AnophelesParamsType& elt,
    double proportionUse,
    double maxInsecticide)
{
    initDeterrence( _relativeAttractiveness, elt.getDeterrency() );
    initSurvival( _preprandialKillingEffect, elt.getPreprandialKillingEffect(), false, maxInsecticide );
    initSurvival( _postprandialKillingEffect, elt.getPostprandialKillingEffect(), true, maxInsecticide );
    if (elt.getFecundityReduction().present()) {
        initSurvival( _fecundityEffect, elt.getFecundityReduction().get(), false/*TODO: err msg*/, maxInsecticide );
    } else {
        _fecundityEffect.init1();
    }
    assert( proportionUse >= 0.0 && proportionUse <= 1.0 );
    // Simpler version of ITN usage/action:
    double propActive = elt.getPropActive();
    if(propActive != 1.0 && propActive != 0.0)
        throw util::xml_scenario_error("GVI.description.propActive: must be either 0 or 1");

    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = proportionUse * propActive;
    proportionUnprotected = 1.0 - proportionProtected;
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
} } // namespace OM::interventions

namespace OM {
namespace interventions {

// —————  per-host sprayed-state: Transmission::PerHostInterventionData, kind IRS  —————
namespace {
const interventions::IRSComponent& componentOf(const Transmission::PerHostInterventionData& d) {
    return static_cast<const interventions::IRSComponent&>(d.component);
}
/// Remaining insecticide content (implies age 1 TS on first use).
double insecticideContent(const Transmission::PerHostInterventionData& d, SimTime now) {
    return d.initialInsecticide * d.decayHet->eval(now - d.deployTime);
}
}

void irsInitNew(Transmission::PerHostInterventionData& d, LocalRng& rng, const Transmission::HumanVectorInterventionComponent&, SimTime /*now*/)
{
    const interventions::IRSComponent& params = componentOf(d);
    d.deployTime = sim::now();
    d.initialInsecticide = params.sampleInitialInsecticide(rng);
    // Varience factor of decay is sampled once per human: human is assumed
    // to account for most variance.
    d.decayHet = params.insecticideDecay->hetSample(rng);
}

void irsRedeploy(Transmission::PerHostInterventionData& d, LocalRng& rng, const Transmission::HumanVectorInterventionComponent&, SimTime now) {
    // PerHost::deployComponent only redeploys when the component ids match,
    // so the passed params are always the same component we were built with.
    d.deployTime = now;
    d.initialInsecticide = componentOf(d).sampleInitialInsecticide(rng);
}

void irsUpdate(Transmission::PerHostInterventionData&, Host::Human&, SimTime /*ts0*/){
}

double irsRelativeAttractiveness(const Transmission::PerHostInterventionData& d, size_t speciesIndex, SimTime now){
    const interventions::IRSComponent::IRSAnopheles& anoph = componentOf(d).species[speciesIndex];
    double effect = anoph.relativeAttractiveness( insecticideContent(d, now) );
    return anoph.byProtection( effect );
}

double irsPreprandialSurvivalFactor(const Transmission::PerHostInterventionData& d, size_t speciesIndex, SimTime now){
    const interventions::IRSComponent::IRSAnopheles& anoph = componentOf(d).species[speciesIndex];
    double effect = anoph.preprandialSurvivalFactor( insecticideContent(d, now) );
    return anoph.byProtection( effect );
}

double irsPostprandialSurvivalFactor(const Transmission::PerHostInterventionData& d, size_t speciesIndex, SimTime now){
    const interventions::IRSComponent::IRSAnopheles& anoph = componentOf(d).species[speciesIndex];
    double effect = anoph.postprandialSurvivalFactor( insecticideContent(d, now) );
    return anoph.byProtection( effect );
}

double irsRelFecundity(const Transmission::PerHostInterventionData& d, size_t speciesIndex, SimTime now){
    const interventions::IRSComponent::IRSAnopheles& anoph = componentOf(d).species[speciesIndex];
    double effect = anoph.fecundityEffect( insecticideContent(d, now) );
    return anoph.byProtection( effect );
}

void irsCheckpoint(Transmission::PerHostInterventionData& d, ostream& stream){
    d.deployTime & stream;
    d.initialInsecticide & stream;
    componentOf(d).insecticideDecay->checkpointSample(d.decayHet, stream);
}
void irsRead(istream& stream, Transmission::PerHostInterventionData& d){
    d.deployTime & stream;
    d.initialInsecticide & stream;
    componentOf(d).insecticideDecay->checkpointSample(d.decayHet, stream);
}

} }
