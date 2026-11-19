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

// Intervention bases for human-vector interaction, split out of
// Transmission/PerHost.h (phase-5 cycle break): Transmission core now only
// forward-declares these and includes the leaf interventions/ComponentId.h;
// the intervention components (ITN/IRS/GVI) include this header — the
// dependency is one-way: interventions -> Transmission.
#ifndef Hmod_VectorIntervention
#define Hmod_VectorIntervention

#include "Global.h"
#include "interventions/Interfaces.h"
#include "util/random.h"

namespace OM {
namespace Host { class Human; }
namespace interventions {
    class ITNComponent;
    class IRSComponent;
}
namespace Transmission {

/**
 * A base class for interventions affecting human-vector interaction.
 * 
 * The constructor should initialise the data to represent an intervention
 * deployed at this time (ctx.clock.now()).
 * 
 * redeploy() should reset the intervention to a freshly deployed state. If
 * necessary, PerHost::deployComponent can be updated to make it create a new
 * instance instead of calling redeploy.
 */
class PerHostInterventionData {
public:
    /// The run this data belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    virtual ~PerHostInterventionData() {}
    
    /** Deploy an intervention. */
    virtual void redeploy( LocalRng& rng, const HumanVectorInterventionComponent& params ) =0;
    
    /** Per time step update. Used by ITNs to update hole decay. */
    virtual void update(Host::Human& human) =0;
    
    /** Get effect of deterrencies of interventions, as an attractiveness multiplier.
     * 
     * @return a value describing effect on attractiveness. Must not be
     * negative. 0 means mosquitoes are fully deterred, 1 that the intervention
     * has no effect, 2 that the intervention attracts twice as many mosquitoes
     * as would otherwise come. */
    virtual double relativeAttractiveness(size_t species) const =0;
    /** Get the killing effect on mosquitoes before they've eaten as a survival
     * multiplier. */
    virtual double preprandialSurvivalFactor(size_t species) const =0;
    /** Get the killing effect on mosquitoes after they've eaten as a survival
     * multiplier. */
    virtual double postprandialSurvivalFactor(size_t species) const =0;
    /// Get the mosquito fecundity multiplier (1 for no effect).
    virtual double relFecundity(size_t species) const =0;
    
    /// Index of effect describing the intervention
    inline interventions::ComponentId id() const { return m_id; }
    
    /// Return true if this component is deployed (i.e. currently active)
    inline bool isDeployed() const{ return deployTime != sim::never(); }
    
    /// Checkpointing (write only)
    void operator& (ostream& stream) {
        m_id & stream; // must be first; read externally so that the correct
                // makeHumanPart function can be found
        checkpoint( stream );
    }
    
protected:
    /// Set the component id
    explicit PerHostInterventionData( RunContext& ctx, interventions::ComponentId id ) :
            ctx(ctx),
            deployTime( this->ctx.clock.now() ),
            m_id(id) {}
    
    /// Checkpointing: write
    virtual void checkpoint( ostream& stream ) =0;
    
    SimTime deployTime = sim::never();        // time of deployment or sim::never()
    interventions::ComponentId m_id;       // component id; don't change
};

/** A base class for human vector intervention parameters. */
class HumanVectorInterventionComponent : public interventions::HumanInterventionComponent {
public:
    virtual ~HumanVectorInterventionComponent() {}
    
    /** Create a new object to store human-specific details of deployment. */
    virtual unique_ptr<PerHostInterventionData> makeHumanPart(RunContext& ctx, LocalRng&) const =0;
    virtual unique_ptr<PerHostInterventionData> makeHumanPart( RunContext& ctx, istream& stream,
            interventions::ComponentId id ) const =0;
    
    // Type-safe downcast seams (phase-5: replace the dynamic_casts in ITN.cpp
    // and IRS.cpp redeploy; owner-class API, null unless the component is one):
    virtual const interventions::ITNComponent* asITNComponent() const{ return nullptr; }
    virtual const interventions::IRSComponent* asIRSComponent() const{ return nullptr; }
    
    const HumanVectorInterventionComponent* asVectorIntervention() const override{
        return this;
    }
    
protected:
    explicit HumanVectorInterventionComponent(interventions::ComponentId id) :
            HumanInterventionComponent(id) {}
};

}
}
#endif
