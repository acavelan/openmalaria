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

// Note: this file is included by exactly one source (Interventions.cpp)
// and contains definitions as well as declarations.

// The includes here are more for documentation than required.
#include "Global.h"
#include "util/errors.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/ClinicalXml.h"
#include "Population.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/VectorModel.h"
#include "Host/HumanTransmission.h"
#include "util/random.h"
#include "util/errors.h"
#include <schema/interventions.h>

#include "Host/Human.h"

namespace OM { namespace interventions {

struct ByDeployTime;    // forward decl for friend

// ———  TimedDeployment and derivatives  ———

/** Interface for timed deployment of an intervention. */
class TimedDeployment {
public:
    /// Create, passing time of deployment
    explicit TimedDeployment(SimTime deployDate) :
            date( deployDate )
    {
        if( deployDate < sim::startDate() ){
            throw util::xml_scenario_error("timed intervention deployment: may not be negative");
        }else if( deployDate >= sim::endDate() ){
            cerr << "Warning: timed intervention deployment at date "<< deployDate;
            cerr << " happens after last survey" << endl;
        }
    }
    virtual ~TimedDeployment() {}
    
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission, const Interventions& definitions) =0;
    
    virtual void print_details( std::ostream& out )const =0;
    
    // Read access required in this file; don't really need protection:
    SimTime date = sim::never();
};

class TimedChangeHSDeployment : public TimedDeployment {
public:
    TimedChangeHSDeployment(SimTime date, unique_ptr<scnXml::HealthSystem> hs, Clinical::ClinicalParameters& library) :
        TimedDeployment( date ),
        newHS(std::move(hs)), library(library)
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission, const Interventions& definitions) {
        xml::setHealthSystem( *newHS, definitions, library );
        newHS.reset();
    }
    virtual void print_details( std::ostream& out )const{
        out << "ChangeHS";
    }
    
private:
    unique_ptr<scnXml::HealthSystem> newHS;
    Clinical::ClinicalParameters& library;
};

class TimedChangeEIRDeployment : public TimedDeployment {
public:
    TimedChangeEIRDeployment(SimTime date, vector<double> daily, int eipDuration) :
        TimedDeployment( date ),
        newEIR(std::move(daily)),
        eipDuration(eipDuration)
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission, const Interventions& definitions) {
        transmission.changeEIRIntervention(newEIR, eipDuration);
        vector<double>().swap(newEIR);
    }
    virtual void print_details( std::ostream& out )const{
        out << "ChangeEIR";
    }
    
private:
    vector<double> newEIR;
    int eipDuration;
};

class TimedUninfectVectorsDeployment : public TimedDeployment {
public:
    TimedUninfectVectorsDeployment( SimTime date ) :
        TimedDeployment( date )
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission, const Interventions& definitions) {
        transmission.uninfectVectors();
    }
    virtual void print_details( std::ostream& out )const{
        out << "UninfectVectors";
    }
};

void VaccineLimits::set( const scnXml::DeploymentBase& deploy ){
    if( deploy.getVaccMinPrevDoses().present() ){
        int v = deploy.getVaccMinPrevDoses().get();
        if( v < 0 ) throw util::xml_scenario_error( "vaccMinPrevDoses: min value is 0" );
        minPrevDoses = v;
    }
    if( deploy.getVaccMaxCumDoses().present() ){
        int v = deploy.getVaccMaxCumDoses().get();
        if( v < 0 ) throw util::xml_scenario_error( "vaccMaxCumDoses: min value is 0" );
        maxCumDoses = v;
    }
}

/// Base class for TimedHumanDeployment and ContinuousHumanDeployment
class HumanDeploymentBase {
protected:
    /**
     * @param deploy XML element describing deployment
     * @param intervention The intervention to deploy (list of components)
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     * @param complement Whether to take the complement of the sub-population
     *  to which deployment will be restricted
     */
    HumanDeploymentBase( const scnXml::DeploymentBase& deploy,
                         shared_ptr<const HumanIntervention> intervention,
                         ComponentId subPop, bool complement ) :
            coverage( deploy.getCoverage() ),
            subPop( subPop ),
            complement( complement ),
            intervention( intervention )
    {
        if( !(coverage >= 0.0 && coverage <= 1.0) ){
            throw util::xml_scenario_error("intervention deployment coverage must be in range [0,1]");
        }
        // If coverageCorr or coverageVar is present, then the user is asking for the copula
        if(deploy.getCoverageCorr().present() || deploy.getCoverageVar().present())
        {
            if(!(deploy.getCoverageCorr().present() && deploy.getCoverageVar().present()))
                throw std::runtime_error("[Gaussian copula]: both coverageCorr and coverageVar must be present");

            if(deploy.getCoverageVar().get() < 0)
                throw std::runtime_error("[Gaussian copula]: negative variance is not allowed");
            
            // if variance is 0, then no copula and fall back to probability = coverage
            if(deploy.getCoverageVar().get() > 0)
            {
                coverageCorr = deploy.getCoverageCorr().get();
                coverageVar = deploy.getCoverageVar().get();
                copula = true;
            }
        }

        vaccLimits.set( deploy );
    }
    
    inline void deployToHuman( Host::Human& human, mon::Deploy::Method method ) const{
        intervention->deploy( human, method, vaccLimits );
    }
    
    bool copula = false;
    double coverage;    // proportion coverage within group meeting above restrictions
    double coverageCorr = 0.0, coverageVar = 0.0;
    VaccineLimits vaccLimits;
    ComponentId subPop;      // ComponentId::wholePop() if deployment is not restricted to a sub-population
    bool complement;
    shared_ptr<const HumanIntervention> intervention;
};


/// Timed deployment of human-specific interventions
class TimedHumanDeployment : public TimedDeployment, protected HumanDeploymentBase {
public:
    /** 
     * @param date Time of deployment
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     */
    TimedHumanDeployment( SimTime date,
                           const scnXml::MassDeployment& mass,
                           shared_ptr<const HumanIntervention> intervention,
                           ComponentId subPop, bool complement ) :
        TimedDeployment( date ),
        HumanDeploymentBase( mass, intervention, subPop, complement ),
        minAge( sim::fromYearsN( mass.getMinAge() ) ),
        maxAge( sim::future() ),
        minAvailability( 0 ),
        maxAvailability( 100 )
    {
        if( mass.getMaxAge().present() )
            maxAge = sim::fromYearsN( mass.getMaxAge().get() );
            
        if( minAge < sim::zero() || maxAge < minAge )
            throw util::xml_scenario_error("timed intervention must have 0 <= minAvailability <= maxAvailability <= 100");

        if( mass.getMaxAvailability().present() )
            maxAvailability = mass.getMaxAvailability().get();

        if( mass.getMinAvailability().present() )
            minAvailability = mass.getMinAvailability().get();
            
        if( minAvailability < sim::zero() || maxAvailability < minAvailability )
            throw util::xml_scenario_error("timed intervention must have 0 <= minAvailability <= maxAvailability <= 100");
    }
    
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission, const Interventions& definitions) {
        for(Human& human : population) {
            const SimTime age = human.age(sim::now());
            if( age >= minAge && age < maxAge ){
                double availability = 0.0;
                for(size_t i = 0; i < transmission.hostParameters.species.size(); ++i)
                    availability += human.transmission.species[i].availability;

                double probability = coverage;

                if(copula && coverage > 0)
                {
                    try 
                    {
                        // Beta distribution for intervention (Beta distributed)
                        const double beta_mean = coverage;  // Assuming this value is given
                        const double beta_var = coverageVar;       // Given as well
                        const double alpha = ((1.0 - beta_mean) / beta_var - 1.0 / beta_mean) * (beta_mean * beta_mean);
                        const double beta = alpha * (1.0 / beta_mean - 1.0);

                        if(transmission.hostParameters.species.size() > 1)
                            throw std::runtime_error("only supports one mosquito species");

                        const double ux = transmission.hostParameters.species[0].entoAvailability->cdf(human.transmission.species[0].sampledAvailability);
                        if(ux < 1)
                        {
                            double xx = gsl_cdf_ugaussian_Pinv(ux);                      // Unit interval to Normal

                            // Apply the Gaussian copula transformation for correlation
                            const ComponentId cid = subPop;

                            /** human.rng.gauss(0.0, 1.0) is calculated once per component and per human. Re-deployment of
                             * the same component on the same human must use the same gaussian sample. Therefore we store
                             * this value in the human the first time it is calculated. */
                            double g = 0.0;
                            const auto it = human.interventions.copulaSamples.find(cid);
                            if (it != human.interventions.copulaSamples.end())
                                g = it->second;
                            else
                            {
                                g = human.rng.gauss(0.0, 1.0);
                                human.interventions.copulaSamples[cid] = g;
                            }
                            const double yy = coverageCorr * xx + g * sqrt(1 - coverageCorr * coverageCorr);

                            // Transform back to unit interval
                            const double uy = gsl_cdf_ugaussian_P(yy);

                            if (alpha < 0 || beta < 0)
                                throw std::runtime_error("resulting alpha and beta parameters must be positive");

                            // Get the final probability from Beta distribution
                            probability = gsl_cdf_beta_Pinv(uy, alpha, beta);
                        }
                        else // Gamma with CV=0
                            probability = human.rng.beta(alpha, beta);
                    }
                    catch (const std::exception& e) {
                        std::ostringstream oss;
                        oss << "[Gaussian copula]: computing correlated intervention probability using Gaussian copula and Beta distribution: "
                            << e.what()
                            << ". Possible reason: coverage might be too high and/or population size might be too small.";
                        throw std::runtime_error(oss.str());
                    }
                }

                if( availability >= transmission.hostParameters.availabilityPercentile(minAvailability) && availability <= transmission.hostParameters.availabilityPercentile(maxAvailability) ) {
                    if( subPop == ComponentId::wholePop() || (human.isInSubPop( subPop ) != complement) ){
                        if( human.rng.bernoulli( probability ) ){
                            deployToHuman( human, mon::Deploy::TIMED );
                        }
                    }
                }
            }
        }
    }
    
    virtual void print_details( std::ostream& out )const{
        out << "Human:" << endl;
        out << "\tage: " << sim::inYears(minAge) << "y\tmax age: " << sim::inYears(maxAge) << "y\tsubpop: ";
        if( subPop == ComponentId::wholePop() ) out << "whole";
        else out << subPop.id;
        out << "\tcomplement: " << complement << "\tcoverage: " << coverage << "\t";
        intervention->print_details( out );
    }
    
protected:
    // restrictions on deployment
    SimTime minAge, maxAge;
    SimTime minAvailability, maxAvailability;
};

/// Timed deployment of human-specific interventions in cumulative mode
class TimedCumulativeHumanDeployment : public TimedHumanDeployment {
public:
    /** 
     * @param date Time of deployment
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     * @param cumCuvId Id of component to test coverage for
     */
    TimedCumulativeHumanDeployment( SimTime date,
                           const scnXml::MassDeployment& mass,
                           shared_ptr<const HumanIntervention> intervention,
                           ComponentId subPop, bool complement,
                           ComponentId cumCuvId ) :
        TimedHumanDeployment( date, mass, intervention, subPop, complement ),
        cumCovInd( cumCuvId )
    {
    }
    
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission, const Interventions& definitions) {
        // Cumulative case: bring target group's coverage up to target coverage
        vector<Host::Human*> unprotected;
        size_t total = 0;       // number of humans within age bound and optionally subPop
        for(Host::Human &human : population) {
            SimTime age = human.age(sim::now());
            if( age >= minAge && age < maxAge ){
                double availability = 0.0;
                for(size_t i = 0; i < transmission.hostParameters.species.size(); ++i)
                    availability += human.transmission.species[i].availability;
                if( availability >= transmission.hostParameters.availabilityPercentile(minAvailability) && availability <= transmission.hostParameters.availabilityPercentile(maxAvailability) ) {
                    if( subPop == ComponentId::wholePop() || (human.isInSubPop( subPop ) != complement) ){
                        total+=1;
                        if( !human.isInSubPop(cumCovInd) )
                            unprotected.push_back( &human );
                    }
                }
            }
        }
        
        if( total == 0 ) return;        // no humans to deploy to; avoid divide by zero
        double propProtected = static_cast<double>( total - unprotected.size() ) / static_cast<double>( total );
        if( propProtected < coverage ){
            // Proportion propProtected are already covered, so need to
            // additionally cover the proportion (coverage - propProtected),
            // selected from the list unprotected.
            double additionalCoverage = (coverage - propProtected) / (1.0 - propProtected);
            cerr << "cum deployment: prop protected " << propProtected << "; additionalCoverage " << additionalCoverage << "; total " << total << endl;
            for(Human* human : unprotected) {
                if( human->rng.uniform_01() < additionalCoverage ){
                    deployToHuman( *human, mon::Deploy::TIMED );
                }
            }
        }
    }
    
protected:
    ComponentId cumCovInd;
};

class TimedVectorDeployment : public TimedDeployment {
public:
    enum class Action { Population, Trap, AddAnimalHosts, ModifyAnimalHosts };
    TimedVectorDeployment(SimTime date, Action action, size_t instance = 0) :
        TimedDeployment(date), action(action), instance(instance) {}

    Action action;
    size_t instance;
    string name;
    double ratio = 0.0;
    SimTime lifespan = sim::never();
    vector<Transmission::Anopheles::NhhDeploymentParameters> nhhParams;

    void deploy(vector<Host::Human>& population, Transmission::TransmissionModel& transmission,
            const Interventions&) override {
        auto* model = dynamic_cast<Transmission::VectorModel*>(&transmission);
        if (!model) return;
        switch (action) {
        case Action::Population:
            model->deployVectorPopInterv(instance);
            break;
        case Action::Trap:
            model->deployVectorTrap(instance, population.size() * ratio, lifespan);
            break;
        case Action::ModifyAnimalHosts:
            for (size_t i = 0; i < model->speciesIndex.size(); ++i) {
                auto& species = *model->species[i];
                if (species.nhhInstances.count(name) == 0)
                    throw util::xml_scenario_error("non human hosts type " + name + " not deployed during non human hosts intervention deployment");
                Transmission::Anopheles::deploy(species.nhhInterventions[name][instance], model->m_rng, sim::now());
            }
            break;
        case Action::AddAnimalHosts: {
            double popSize = population.size();
            if (model->interventionMode != Transmission::SimulationMode::dynamicEIR)
                throw util::xml_scenario_error("vector interventions can only be used in dynamic transmission mode (mode=\"dynamic\")");
            for (size_t i = 0; i < model->speciesIndex.size(); ++i) {
                auto& species = *model->species[i];
                if (species.nhhInstances.count(name) != 0)
                    throw util::xml_scenario_error("non human hosts type " + name + " already deployed during non human hosts deployment");
                const auto& p = nhhParams[i];
                const double adultAvail = transmission.hostParameters.species[i].entoAvailability->mean();
                const double availability = transmission.hostParameters.species[i].entoAvailabilityFactor * popSize * adultAvail * p.relativeAvailabilityHuman;
                species.nhhInstances[name] = {availability, p.biting, p.findingRestSite,
                    p.resting, p.fecundity, sim::now() + lifespan};
            }
            break;
        }
        }
    }

    void print_details(std::ostream& out) const override {
        switch (action) {
        case Action::Population: out << "Vector"; break;
        case Action::Trap:
            out << "VectorTrap:" << endl << "\tratio: " << ratio << "\tlifespan: " << lifespan;
            break;
        case Action::ModifyAnimalHosts: out << "NonHumanHosts(" << name << ")"; break;
        case Action::AddAnimalHosts:
            out << "AddNonHumanHosts(" << name << "):" << endl << "\tlifespan: " << lifespan;
            break;
        }
    }
};

// ———  ContinuousHumanDeployment  ———

/** Interface for continuous deployment of an intervention. */
class ContinuousHumanDeployment : protected HumanDeploymentBase {
public:
    /// Create, passing deployment age
    ContinuousHumanDeployment( SimTime begin, SimTime end,
                                 const ::scnXml::ContinuousDeployment& elt,
                                 shared_ptr<const HumanIntervention> intervention,
                                 ComponentId subPop, bool complement ) :
            HumanDeploymentBase( elt, intervention, subPop, complement ),
            begin( begin ), end( end ),
            deployAge( sim::fromYearsN( elt.getTargetAgeYrs() ) ),
            minAvailability( 0 ),
            maxAvailability( 100 )
    {
        if( begin < sim::startDate() || end < begin ){
            throw util::xml_scenario_error("continuous intervention must have startDate <= begin <= end");
        }
        if( deployAge <= sim::zero() ){
            ostringstream msg;
            msg << "continuous intervention with target age "<<elt.getTargetAgeYrs();
            msg << " years corresponds to time step " << sim::inSteps(deployAge);
            msg << "; must be at least 1.";
            throw util::xml_scenario_error( msg.str() );
        }
        if( deployAge > sim::maxHumanAge() ){
            ostringstream msg;
            msg << "continuous intervention must have target age no greater than ";
            msg << sim::inYears(sim::maxHumanAge());
            throw util::xml_scenario_error( msg.str() );
        }
        if( elt.getMaxAvailability().present() )
            maxAvailability = elt.getMaxAvailability().get();

        if( elt.getMinAvailability().present() )
            minAvailability = elt.getMinAvailability().get();
            
        if( minAvailability < sim::zero() || maxAvailability < minAvailability )
            throw util::xml_scenario_error("timed intervention must have 0 <= minAvailability <= maxAvailability <= 100");
    }
    
    /** Apply filters and potentially deploy.
     * 
     * @returns false iff this deployment (and thus all later ones in the
     *  ordered list) happens in the future. */
    bool filterAndDeploy( Host::Human& human, const Transmission::HostTransmissionParameters& hostParameters ) const{
        SimTime age = human.age(sim::now());
        if( deployAge > age ){
            // stop processing continuous deployments for this
            // human for now because remaining ones happen in the future
            return false;
        }else if( deployAge == age ){
            double availability = 0.0;
            for(size_t i = 0; i < hostParameters.species.size(); ++i)
                availability += human.transmission.species[i].sampledAvailability;

            if( availability >= hostParameters.availabilityPercentile(minAvailability) && availability <= hostParameters.availabilityPercentile(maxAvailability) )
            {
                auto now = sim::intervDate();
                if( begin <= now && now < end &&
                    ( subPop == ComponentId::wholePop() ||
                        (human.isInSubPop( subPop ) != complement)
                    ) &&
                    human.rng.uniform_01() < coverage )     // RNG call should be last test
                {
                    deployToHuman( human, mon::Deploy::CTS );
                }
            }
        }//else: for some reason, a deployment age was missed; ignore it
        return true;
    }
    
    inline void print_details( std::ostream& out )const{
        out << "Human:" << endl;
        out << "\tbegin: " << begin << "\tend: ";
        if( end == sim::future() ) out << "(none)";
        else out << end;
        out << "\tdeployAge: " << sim::inYears(deployAge) << "y\tsubpop: ";
        if( subPop == ComponentId::wholePop() ) out << "whole";
        else out << subPop.id;
        out << "\tcomplement: " << complement << "\tcoverage: " << coverage << "\t";
        intervention->print_details( out );
    }
    
protected:
    SimTime begin, end;    // first time step active and first time step no-longer active
    SimTime deployAge = sim::never();
    SimTime minAvailability, maxAvailability;

    friend ByDeployTime;
};

struct ByDeployTime {
    bool operator() (const unique_ptr<TimedDeployment>& a, const unique_ptr<TimedDeployment>& b) const{
        return a->date < b->date;
    }
    
    bool operator()(const ContinuousHumanDeployment& a, const ContinuousHumanDeployment& b) const{
        return a.deployAge < b.deployAge;
    }
} byDeployTime;

} }
