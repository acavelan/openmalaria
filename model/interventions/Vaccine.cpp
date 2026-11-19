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
#include "interventions/Vaccine.h"
#include "interventions/HumanComponents.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include "Host/WithinHost/Genotypes.h"

#include <limits>
#include <cmath>

namespace OM {
namespace interventions {
using namespace OM::util;

vector<VaccineComponent*> VaccineComponent::params;
ComponentId VaccineComponent::reportComponent = ComponentId::wholePop();

void VaccineComponent::deploy(RunContext& ctx, Host::Human& human, mon::Deploy::Method method, VaccineLimits vaccLimits) const
{
    bool administered = human.vaccine.possiblyVaccinate( human, id(), vaccLimits );
    if( administered && VaccineComponent::reportComponent == id() ){
        mon::recordDeploy(mon::measure("nMassVaccinations"), mon::measure("nEPIVaccinations"), human, method);
    }
    assert(type == Vaccine::PEV || type == Vaccine::BSV || type == Vaccine::TBV);
    if (method == mon::Deploy::TREAT)
        mon::recordEvent(mon::measure("nTreatDeployments"), human);
}

Component::Type VaccineComponent::componentType() const
{
    if( type == Vaccine::PEV ) return Component::PEV;
    else if( type == Vaccine::BSV ) return Component::BSV;
    else if( type == Vaccine::TBV ) return Component::TBV;
    else throw SWITCH_DEFAULT_EXCEPTION;
}

void VaccineComponent::print_details(ostream& out) const
{
        out << id().id << "\t" << (type == Vaccine::PEV ? "PEV" : (type == Vaccine::BSV ? "BSV" : "TBV"));
}

double VaccineComponent::getInitialEfficacy (LocalRng& rng, size_t numPrevDoses, uint32_t genotype) const
{
    /* If initialMeanEfficacy.size or more doses have already been given, use
     * the last efficacy. */
    if (numPrevDoses >= perGenotypeInitialMeanEfficacy[genotype].size())
        numPrevDoses = perGenotypeInitialMeanEfficacy[genotype].size() - 1;
    double ime = perGenotypeInitialMeanEfficacy[genotype][numPrevDoses];
    //NOTE(validation): With extra valiadation in random, the first difference is noticed here:
    util::streamValidate(ime);
    util::streamValidate(efficacyB);
    if (ime == 0.0){
        return 0.0;
    } else if (ime < 1.0) {
        // parity: order — local RNG stream: draw call order fixes stream state
        double result = rng.betaWithMean (ime, efficacyB);
        //NOTE(validation):: Without extra validation in random, the first difference is noticed here:
        util::streamValidate(result);
        //TODO(validation):: Why the difference? Bug in/limitation of StreamValidatior?
        // Extra memory allocation due to the extra logging causes some bad
        // memory usage to manifest differently? Am I forgetting to initialise
        // some variable?
        return result;
    } else {
        return 1.0;
    }
}

#if 0
/* R_0: this used to use _the_ vaccine configuration; this is just a check.
 * If reimplementing, use a new vaccine instance; implement some hook to ensure
 * that TBV is given to everyone after enabling the intervention, etc. */
void Vaccine::verifyEnabledForR_0 (){
    if( !types[PEV].active || !types[TBV].active )
        throw util::xml_scenario_error("PEV and TBV vaccines must have a "
                "description to use the insertR_0Case intervention");
}
#endif

// this is only used for checkpointing; loaded values are unimportant
PerEffectPerHumanVaccine::PerEffectPerHumanVaccine() :
    component(ComponentId::wholePop()),
    numDosesAdministered(0)
{
}

PerEffectPerHumanVaccine::PerEffectPerHumanVaccine( LocalRng& rng, ComponentId id, const VaccineComponent& params ) :
    component( id ), numDosesAdministered(0)
{
    hetSample = params.decayFunc->hetSample(rng);
}

void PerEffectPerHumanVaccine::checkpointHet(ostream& stream) {
    VaccineComponent::getParams(component).decayFunc->checkpointSample(hetSample, stream);
}

void PerEffectPerHumanVaccine::checkpointHet(istream& stream) {
    VaccineComponent::getParams(component).decayFunc->checkpointSample(hetSample, stream);
}

double PerHumanVaccine::getFactor( Vaccine::Types type, uint32_t genotype) const {
    double factor = 1.0;
    for( EffectList::const_iterator effect = effects.begin(); effect != effects.end(); ++effect ){
        if( VaccineComponent::getParams(effect->component).type == type ){
            SimTime age = ctx->clock.ts1() - effect->timeLastDeployment;  // implies age 1 TS on first use
            double decayFactor = effect->hetSample->eval( age );
            factor *= 1.0 - effect->perGenotypeInitialEfficacy[genotype] * decayFactor;
        }
    }
    return factor;
}

bool PerHumanVaccine::possiblyVaccinate( Host::Human& human,
        ComponentId componentId, interventions::VaccineLimits vaccLimits )
{
    PerEffectPerHumanVaccine* effect = 0;
    for( EffectList::iterator it = effects.begin(); it != effects.end(); ++it ){
        if( it->component == componentId ){
            effect = &*it;
            break;
        }
    }
    
    uint32_t numDosesAdministered = (effect == 0) ? 0 : effect->numDosesAdministered;
    if( numDosesAdministered < vaccLimits.minPrevDoses ||
        numDosesAdministered >= vaccLimits.maxCumDoses )
        return false;   // no vaccination (this replaces the old schedule for continuous doses)
    
    const VaccineComponent& params = VaccineComponent::getParams( componentId );
    
    if( effect == 0 ){
        effects.push_back( PerEffectPerHumanVaccine( human.rng, componentId, params ) );
        effect = &effects.back();
    }
    
    effect->perGenotypeInitialEfficacy.resize(ctx->genotypes.N());

    for(size_t g=0; g<params.perGenotypeInitialMeanEfficacy.size(); g++)
        effect->perGenotypeInitialEfficacy[g] = params.getInitialEfficacy(human.rng, numDosesAdministered, g);

    util::streamValidate(effect->perGenotypeInitialEfficacy);
    
    effect->numDosesAdministered = numDosesAdministered + 1;
    effect->timeLastDeployment = ctx->clock.nowOrTs1();
    
    return true;
}

}
}
