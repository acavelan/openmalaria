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

#ifndef Hmod_interventions_HumanComponents
#define Hmod_interventions_HumanComponents

#include "Global.h"
#include "interventions/Interfaces.h"
#include "util/sampler.h"
#include "util/DecayFunction.h"
#include "util/checkpoint_containers.h"
#include <map>

namespace OM {
namespace Host {
    struct Human;
}
namespace interventions {
    class VaccineComponent;
    using util::LocalRng;

class VectorInterventionComponent;

class VectorInterventionState {
public:
    virtual ~VectorInterventionState() {}
    
    /** Deploy an intervention. */
    virtual void redeploy( LocalRng& rng, const VectorInterventionComponent& params ) =0;
    
    /** Per time step update. Used by ITNs to update hole decay. */
    // Returns true when this update removes the component's cohort membership.
    virtual bool update(LocalRng&, SimTime) { return false; }
    virtual interventions::Component::Type componentType() const =0;
    
    /** Get effect of deterrencies of interventions, as an attractiveness multiplier.
     * 
     * @return a value describing effect on attractiveness. Must not be
     * negative. 0 means mosquitoes are fully deterred, 1 that the intervention
     * has no effect, 2 that the intervention attracts twice as many mosquitoes
     * as would otherwise come. */
    virtual double relativeAttractiveness(size_t species, SimTime effectTime) const =0;
    /** Get the killing effect on mosquitoes before they've eaten as a survival
     * multiplier. */
    virtual double preprandialSurvivalFactor(size_t species, SimTime effectTime) const =0;
    /** Get the killing effect on mosquitoes after they've eaten as a survival
     * multiplier. */
    virtual double postprandialSurvivalFactor(size_t species, SimTime effectTime) const =0;
    /// Get the mosquito fecundity multiplier (1 for no effect).
    virtual double relFecundity(size_t species, SimTime effectTime) const =0;
    
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
    explicit VectorInterventionState( interventions::ComponentId id ) :
            deployTime( sim::now() ),
            m_id(id) {}
    
    /// Checkpointing: write
    virtual void checkpoint( ostream& stream ) =0;
    
    SimTime deployTime = sim::never();        // time of deployment or sim::never()
    interventions::ComponentId m_id;       // component id; don't change
};

/** A base class for human vector intervention parameters. */
class VectorInterventionComponent : public interventions::HumanInterventionComponent {
public:
    virtual ~VectorInterventionComponent() {}
    
    /** Create a new object to store human-specific details of deployment. */
    virtual unique_ptr<VectorInterventionState> makeHumanPart(LocalRng&) const =0;
    virtual unique_ptr<VectorInterventionState> makeHumanPart( istream& stream,
            interventions::ComponentId id ) const =0;
protected:
    explicit VectorInterventionComponent(interventions::ComponentId id) :
            HumanInterventionComponent(id) {}
};


// ———  sub-population removal  ———

namespace SubPopRemove{
    /** Identifiers for types of sub-population removal option, followed by NUM
     * in last place (as a counter).
     * 
     * ON_FIRST_BOUT: remove humans from some sub-populations at the start of
     * their first clinical events following inclusion in the sub-population.
     * 
     * ON_FIRST_INFECTION: remove patent humans from some sub-population(s)
     * during each survey.
     * 
     * ON_FIRST_TREATMENT: remove humans from some sub-populations when
     * assigning each course of treatment.
     */
    enum RemoveAtCode {
        ON_FIRST_BOUT,
        ON_FIRST_INFECTION,
        ON_FIRST_TREATMENT,
        NUM };
}

/** For each RemoveAtCode (excluding NUM), this is a list of all sub-population
 * identifiers for which the option is enabled. */
extern vector<ComponentId> removeAtIds[SubPopRemove::NUM];

// ———  vaccines  ———

namespace Vaccine{
    enum Types { PEV, BSV, TBV, NumVaccineTypes };
}

/** Per vaccine effect (type), per human details. */
class PerEffectPerHumanVaccine {
public:
    //Note: this constructor is only for checkpointing
    PerEffectPerHumanVaccine();
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        component & stream;
        numDosesAdministered & stream;
        timeLastDeployment & stream;
        perGenotypeInitialEfficacy & stream;
        checkpointHet(stream);
    }
    
private:
    PerEffectPerHumanVaccine( LocalRng& rng, ComponentId id, const VaccineComponent& params );
    
    ComponentId component;      // id of component (for finding parameters)
    // Type of vaccine
    /** Number of vaccine doses this individual has received.
     *
     * If an individual misses one EPI (continuous) vaccine dose, it's
     * intentional that they also miss following EPI doses (unless a timed mass
     * vaccination reintroduces them to the EPI schedule). */
    uint32_t numDosesAdministered;
    /// Time of last vaccination with this vaccine type
    SimTime timeLastDeployment = sim::never();
    /// Efficacy at last deployment (undecayed)
    vector<double> perGenotypeInitialEfficacy;
    unique_ptr<util::DecayFunction> hetSample;

    void checkpointHet(ostream& stream);
    void checkpointHet(istream& stream);
    
    friend class PerHumanVaccine;
};

/** Per-human vaccine code. */
class PerHumanVaccine {
public:
    PerHumanVaccine() {}
    
    /** Get one minus the efficacy of the vaccine (1 for no effect, 0 for full effect). */
    double getFactor( Vaccine::Types type, uint32_t genotype = 0) const;
    
    /** Vaccinate unless the passed VaccineLimits specify not to.
     * 
     * @returns true when the vaccine is administered */
    bool possiblyVaccinate( Host::Human& human, ComponentId componentId,
                           interventions::VaccineLimits vaccLimits );

#if 0
    /// Hack for R_0 experiment: make current human the infection source
    inline void specialR_0(){
        // At this point special TBV should have been given to everyone but the PEV to no-one
        effects[R_0_PEV].initialEfficacy = 1.0;
        effects[R_0_TBV].initialEfficacy = 0.0;
    }
#endif
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        effects & stream;
    }

private:
    /// Details for each deployed vaccine for this human
    typedef std::vector<PerEffectPerHumanVaccine> EffectList;
    EffectList effects;
};

struct HumanInterventions {
    PerHumanVaccine vaccine;
    std::map<ComponentId, SimTime> cohortExpiry;
    uint32_t cohortSet = 0;
    uint32_t nextContinuousDeployment = 0;
    std::map<ComponentId, double> copulaSamples;
    // Retain expired instances: redeployment reuses personal heterogeneity samples.
    vector<unique_ptr<VectorInterventionState>> vectorInterventions;
};

void deployVectorIntervention(HumanInterventions&, LocalRng&, const VectorInterventionComponent&);
void updateVectorInterventions(HumanInterventions&, LocalRng&, SimTime stepStart);
bool hasActiveIntervention(const HumanInterventions&, Component::Type);
void checkpointVectorInterventions(HumanInterventions&, ostream&);
void checkpointVectorInterventions(HumanInterventions&, istream&,
    const vector<unique_ptr<HumanInterventionComponent>>& definitions);

}
}
#endif
