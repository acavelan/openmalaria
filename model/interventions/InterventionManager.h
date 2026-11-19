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
#include "Transmission/TransmissionModel.h"
namespace scnXml { class Interventions; }

namespace OM {
    class Population;
    namespace Host {
        class Human;
    }

namespace interventions {

/** Management of interventions deployed on a per-time-step basis. */
class InterventionManager {
public:
    /** Read XML descriptions. */
    static void init (RunContext& ctx, const scnXml::Interventions& intervElt, const Population &population, Transmission::TransmissionModel& transmission);
    
    /// Checkpointing: most members are only set from XML,
    /// nextTimed varies but is re-set by loadFromCheckpoint; only
    /// imported infections need serialising.
    static void checkpoint (RunContext& ctx, std::istream& stream);
    static void checkpoint (RunContext& ctx, std::ostream& stream);

    /** Call after loading a checkpoint, passing the intervention-period time.
     * 
     * Serves to replace health-system and EIR where changeHS/changeEIR
     * interventions have been used. */
    static void loadFromCheckpoint(
                RunContext& ctx, vector<Host::Human> &population,
                Transmission::TransmissionModel& transmission);
    
    /** @brief Deploy interventions
     *
     * Timed interventions are deployed for this time step.
     * 
     * Continuous interventions are deployed as humans reach the target ages.
     * Unlike with vaccines, missing one schedule doesn't preclude the next. */
    static void deploy(RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel& transmission);
    
    /** Get a constant reference to a component class with a certain index.
     * 
     * @throws util::base_exception if the index is out-of-range */
    static const HumanInterventionComponent& getComponent( RunContext& ctx, ComponentId id );
    
    /** Get a numeric ComponentId from the textual identifier used in the XML.
     * 
     * If textId is unknown, an xml_scenario_error is thrown. */
    static ComponentId getComponentId( RunContext& ctx, const std::string textId );
};

} }
#endif
