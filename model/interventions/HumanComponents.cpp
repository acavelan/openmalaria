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

#include "interventions/HumanComponents.h"

namespace OM::interventions {
using namespace util;

void deployVectorIntervention(HumanInterventions& state, LocalRng& rng, const VectorInterventionComponent& params){
    // This adds per-host per-intervention details to the host's data set.
    // This data is never removed since it can contain per-host heterogeneity samples.
    for( auto iter = state.vectorInterventions.begin(); iter != state.vectorInterventions.end(); ++iter ){
        if( (*iter)->id() == params.id() ){
            // already have a deployment for that description; just update it
            (*iter)->redeploy( rng, params );
            return;
        }
    }
    // no deployment for that description: must make a new one
    state.vectorInterventions.push_back( params.makeHumanPart(rng) );
}

bool hasActiveIntervention(const HumanInterventions& state, Component::Type type){
    for( auto iter = state.vectorInterventions.begin(); iter != state.vectorInterventions.end(); ++iter ){
        if( (*iter)->isDeployed() ){
            if( (*iter)->componentType() == type )
                return true;
        }
    }
    return false;
}

void checkpointVectorInterventions(HumanInterventions& state, ostream& stream){
    state.vectorInterventions.size() & stream;
    for( auto iter = state.vectorInterventions.begin(); iter != state.vectorInterventions.end(); ++iter ){
        *iter & stream;
    }
}
void checkpointVectorInterventions(HumanInterventions& state, istream& stream,
        const vector<unique_ptr<interventions::HumanInterventionComponent>>& definitions) {
    size_t l;
    l & stream;
    validateListSize(l);
    state.vectorInterventions.clear();
    for( size_t i = 0; i < l; ++i ){
        interventions::ComponentId id( stream );
        try{
            if (id.id >= definitions.size()) throw util::base_exception("invalid component id");
            const auto& gen_params = *definitions[id.id];
            const VectorInterventionComponent *params =
                dynamic_cast<const VectorInterventionComponent*>( &gen_params );
            if( params == 0 )
                throw util::base_exception( "" );       // see catch block below
            state.vectorInterventions.push_back( params->makeHumanPart( stream, id ) );
        }catch( util::base_exception& e ){
            // two causes, both boil down to index being wrong
            throw util::checkpoint_error( "bad value in checkpoint file" );
        }
    }
}

void updateVectorInterventions(HumanInterventions& state, LocalRng& rng, SimTime stepStart) {
    for (auto& component : state.vectorInterventions) {
        if (component->update(rng, stepStart))
            state.cohortExpiry.erase(component->id());
    }
}

}
