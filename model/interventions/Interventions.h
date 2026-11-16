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

#ifndef OM_INTERVENTIONS_INTERVENTIONS
#define OM_INTERVENTIONS_INTERVENTIONS

#include "Global.h"
#include "interventions/Interfaces.h"
#include "interventions/HumanComponents.h"
#include "Host/ImportedInfections.h"
#include "Transmission/TransmissionModel.h"
#include "schema/interventions.h"

namespace OM {
    struct Population;
    namespace Clinical { struct ClinicalParameters; }
    namespace Host {
        struct Human;
    }

namespace interventions {

class ContinuousHumanDeployment;
class TimedDeployment;

/** Management of interventions deployed on a per-time-step basis. */
struct Interventions {
    Interventions();
    ~Interventions();
    /** Read XML descriptions. */
    void init (const scnXml::Interventions& intervElt, const Population &population, Transmission::TransmissionModel& transmission, Clinical::ClinicalParameters& library);
    
    /// Checkpointing
    template<class S>
    void checkpoint (S& stream) {
        using namespace OM::util::checkpoint;
        // most members are only set from XML,
        // nextTimed varies but is re-set by loadFromCheckpoint
        importedInfections & stream;
    }

    /** Call after loading a checkpoint, passing the intervention-period time.
     * 
     * Serves to replace health-system and EIR where changeHS/changeEIR
     * interventions have been used. */
    void loadFromCheckpoint(
                vector<Host::Human> &population,
                Transmission::TransmissionModel& transmission);
    
    /** @brief Deploy interventions
     *
     * Timed interventions are deployed for this time step.
     * 
     * Continuous interventions are deployed as humans reach the target ages.
     * Unlike with vaccines, missing one schedule doesn't preclude the next. */
    void deploy(vector<Host::Human> &population, Transmission::TransmissionModel& transmission);
    

    /** Get a constant reference to a component class with a certain index.
     * 
     * @throws util::base_exception if the index is out-of-range */
    inline const HumanInterventionComponent& getComponent( ComponentId id ) const{
        if( id.id >= humanComponents.size() )
            throw util::base_exception( "invalid component id" );
        return *humanComponents[id.id];
    }
    
    /** Get a numeric ComponentId from the textual identifier used in the XML.
     * 
     * If textId is unknown, an xml_scenario_error is thrown. */
    ComponentId getComponentId( const std::string textId ) const;
    
    // Map of textual identifiers to numeric identifiers for components
    std::map<std::string,ComponentId> identifierMap;
    // Owns definitions referenced by deployments and personal intervention state.
    vector<unique_ptr<HumanInterventionComponent>> humanComponents;
    // Continuous interventions, sorted by deployment age (weakly increasing)
    vector<ContinuousHumanDeployment> continuous;
    // List of all timed interventions. Should be sorted (time weakly increasing).
    vector<unique_ptr<TimedDeployment>> timed;
    uint32_t nextTimed = 0;  // not checkpointed (see loadFromCheckpoint)
    
    // imported infections are not really interventions, and handled by a separate class
    // (but are grouped here for convenience and due to association in schema)
    OM::Host::ImportedInfections importedInfections;
    CohortRemovalRules removalRules;
};

} }
#endif
