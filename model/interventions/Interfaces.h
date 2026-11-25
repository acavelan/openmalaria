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

#include "Global.h"
#include "core/RunContext.h"
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
namespace Transmission {
    class HumanVectorInterventionComponent;
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

/** A description of one component of a human intervention.
 * 
 * Note that one component can have several "actions", but that deployment and
 * decay of these "actions" is usually related.
 * 
 * This is a base class. */
class HumanInterventionComponent {
public:
    virtual ~HumanInterventionComponent() {}
    
    /** Deploy the component to a pre-selected human.
     * 
     * This may be called between updates (usual intervention deployment time)
     * or from case management. It should use time ctx.clock.nowOrTs1().
     * 
     * @param human Individual receiving the intervention
     * @param method Channel of deployment (mass, continuous)
     */
    virtual void deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method,
        VaccineLimits vaccLimits ) const =0;
    
    /** Get the component identifier. */
    inline ComponentId id()const{ return m_id; }
    
    /** Type-safe downcast seam to the vector-intervention base (phase-5:
     * replaces the dynamic_cast in PerHost's checkpoint restore; overridden
     * by Transmission::HumanVectorInterventionComponent, null otherwise). */
    virtual const Transmission::HumanVectorInterventionComponent* asVectorIntervention() const{
        return nullptr;
    }
    
    /** Get the duration. */
    inline SimTime duration()const{ return m_duration; }
    
    /** Returns the appropriate descriptor from the Component::Type enum.
     * 
     * This is only used a small number of times during setup, so doesn't need
     * to be fast. */
    virtual Component::Type componentType() const=0;
    
    virtual void print_details( std::ostream& out )const =0;
    
    // Only for use by InterventionManager:
    inline void setExpireAfter( SimTime duration ){
        m_duration = duration;
    }
    
protected:
    /** Construct (from a derived class).
     * 
     * @param id Component identifier; used as an identifier for cumulative
     *  deployment as well as to match human-specific components to general
     *  parameters (i.e. objects of the class extending this one).
     */
    explicit HumanInterventionComponent(ComponentId id) : m_id(id) {}
    
private:
    /** Don't copy (this may be possible but shouldn't be needed). */
    HumanInterventionComponent( const HumanInterventionComponent& ) = delete;
    HumanInterventionComponent& operator= ( const HumanInterventionComponent& ) = delete;
    
    ComponentId m_id;
    SimTime m_duration = sim::never();
};

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
    // List of pointers to components. Does not manage memory (InterventionManager::humanComponents does that).
    vector<const HumanInterventionComponent*> components;
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
