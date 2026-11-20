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

#ifndef OM_INTERVENTIONS_INTERFACES
#define OM_INTERVENTIONS_INTERFACES

#include "Clock.h"
#include "RunContext.h"
#include "mon/Monitoring.h"
#include "interventions/ComponentId.h"

namespace scnXml{
    class DeploymentBase;
    class Component;
    class Condition;
    class DTDeploy;
    class Deploy;
    class TriggeredDeployments;
}
namespace xsd{ namespace cxx{ namespace tree{
    // real declaration: template<typename T, bool fund = fundamental_p<T>::r>
    // declarations below spell the bool explicitly so no default is needed
    template<typename T, bool> class sequence;
} } }

namespace OM {
namespace Host {
    class Human;
}
namespace interventions {

/// Component::Type and ComponentId live in the leaf header ComponentId.h
/// (included above) so Transmission/PerHost.h can use them without
/// including this header (phase-5 cycle break).

/** Specifies limits on the number of existing doses when deciding whether to
 * vaccinate a human. */
struct VaccineLimits{
    VaccineLimits() : minPrevDoses( 0 ), maxCumDoses( numeric_limits<uint32_t>::max() ) {}
    void set( const scnXml::DeploymentBase& );
    uint32_t minPrevDoses, maxCumDoses;
};

/* Phase 7.2b: the HumanInterventionComponent class tree is gone. A component
 * is a plain Component {id, kind, duration, variant payload} — see
 * interventions/Components.h; deployment is the deploy_component() switch. */

/** A description of a human intervention (as a list of components). */
class HumanIntervention {
public:
    //NOTE: it would be preferable to use named constructors rather than rely
    // on selection of the correct overloaded constructor here, but doing so
    // would result in an extra copy without C++11's move semantics.
    
    /** Create from a list of XML elements: list of intervention components
     * with unconditional deployment (for triggered deployments). **/
    explicit HumanIntervention(
        RunContext& ctx,
        const xsd::cxx::tree::sequence<scnXml::Component, false>& componentList );
    /** Create from a list of XML elements: list of intervention components
     * and list of conditions of deployment. **/
    explicit HumanIntervention(
        RunContext& ctx,
        const xsd::cxx::tree::sequence<scnXml::Component, false>& componentList,
        const xsd::cxx::tree::sequence<scnXml::Condition, false>& conditionList );
    /** Create from a list of XML elements: list of intervention components
     * with unconditional deployment (for triggered deployments). **/
    explicit HumanIntervention(
        RunContext& ctx,
        const xsd::cxx::tree::sequence<scnXml::DTDeploy, false>& componentList );
    
    /** Deploy all components to a pre-selected human. */
    void deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method,
        VaccineLimits vaccLimits ) const;
    
    void print_details( std::ostream& out )const;
    
protected:
    // Components to deploy, in type-rank order (contract b). Resolved by id
    // against the run registry at deploy time (InterventionManager owns them).
    vector<ComponentId> components;
    // List of conditions of deployment. All must be satisfied to deploy.
    vector<size_t> conditions;
};

//TODO: this class is likely more complicated than it needs to be (e.g. SubList
// could just store ComponentId instead of inheriting HumanIntervention)
class TriggeredDeployments {
public:
    TriggeredDeployments( RunContext& ctx, const scnXml::TriggeredDeployments& elt );
    
    void deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method,
            VaccineLimits vaccLimits )const;
    
private:
    struct SubList : protected HumanIntervention{
        SubList( RunContext& ctx, const scnXml::Deploy& elt );
        
        void deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method,
            VaccineLimits vaccLimits )const;
        
        // Deployment restrictions:
        SimTime minAge, maxAge;
        double coverage;
    };
    
    vector<SubList> lists;
};

} }
#endif
