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

// Note: this file is included by several sources (InterventionManager.cpp,
// xml/interventions.cpp and core/RunExtensions.h). All definitions here must
// stay inline/template or header-only.

// The includes here are more for documentation than required.
#ifndef OM_INTERVENTIONS_DEPLOYMENTS_H
#define OM_INTERVENTIONS_DEPLOYMENTS_H

#include "Global.h"
#include "core/RunContext.h"
#include "util/errors.h"
#include "Population.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/PerHost.h"
#include "util/random.h"
#include "util/errors.h"
namespace scnXml
{
    class DeploymentBase;
    class MassDeployment;
    class ContinuousDeployment;
}

#include "util/SpeciesIndexChecker.h"
#include "Host/Human.h"

namespace OM { namespace interventions {

using Host::Human;

struct ByDeployTime;    // forward decl for friend

// ———  TimedDeployment and derivatives  ———

/** Interface for timed deployment of an intervention. */
class TimedDeployment {
public:
    /// Create, passing time of deployment
    explicit TimedDeployment(RunContext& ctx, SimTime deployDate) :
            date( deployDate )
    {
        if( deployDate < ctx.clock.start_date ){
            throw util::xml_scenario_error("timed intervention deployment: may not be negative");
        }else if( deployDate >= ctx.clock.end_date ){
            cerr << "Warning: timed intervention deployment at date "<< deployDate;
            cerr << " happens after last survey" << endl;
        }
    }
    virtual ~TimedDeployment() {}
    
    virtual void deploy (RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel& transmission) =0;
    
    /** Re-applied at checkpoint restore (true only for changeHS/changeEIR —
     * they mutate model state, unlike per-host deployments). Phase-5: the
     * declared-API replacement for dynamic_cast in loadFromCheckpoint. */
    virtual bool replayedOnRestore() const { return false; }
    
    virtual void print_details( std::ostream& out )const =0;
    
    // Read access required in this file; don't really need protection:
    SimTime date = sim::never();
};

class DummyTimedDeployment : public TimedDeployment {
public:
    DummyTimedDeployment(RunContext& ctx) :
        TimedDeployment( ctx, ctx.clock.start_date /* hack */ )
    {
        // TimedDeployment's ctor checks that the deployment time-step is
        // within the intervention period. We want this time to be after the
        // last time-step, so set the time here after TimedDeployment's ctor
        // check has been done (hacky).
        date = sim::future();
    }
    virtual void deploy (RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel& transmission) {}
    virtual void print_details( std::ostream& out )const{
        out << "Dummy";
    }
};

class TimedUninfectVectorsDeployment : public TimedDeployment {
public:
    TimedUninfectVectorsDeployment( RunContext& ctx, SimTime date ) :
        TimedDeployment( ctx, date )
    {}
    virtual void deploy (RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel& transmission) {
        transmission.uninfectVectors();
    }
    virtual void print_details( std::ostream& out )const{
        out << "UninfectVectors";
    }
};

class HumanDeploymentBase {
protected:
    RunContext* ctx;  /**< owning deployment ctor supplied; pointer so the class stays copyable */

    /**
     * @param deploy XML element describing deployment
     * @param intervention The intervention to deploy (list of components)
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     * @param complement Whether to take the complement of the sub-population
     *  to which deployment will be restricted
     */
    HumanDeploymentBase( RunContext& ctx, const scnXml::DeploymentBase& deploy,
                         shared_ptr<const HumanIntervention> intervention,
                         ComponentId subPop, bool complement );
    
    inline void deployToHuman( Host::Human& human, mon::Deploy::Method method ) const{
        intervention->deploy( *ctx, human, method, vaccLimits );
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
    TimedHumanDeployment( RunContext& ctx, SimTime date,
                           const scnXml::MassDeployment& mass,
                           shared_ptr<const HumanIntervention> intervention,
                           ComponentId subPop, bool complement );
    
    virtual void deploy (RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel& transmission) {
        for(Human& human : population) {
            const SimTime age = human.age(ctx.clock.now());
            if( age >= minAge && age < maxAge ){
                double availability = 0.0;
                for(size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(ctx); ++i)
                    availability += human.perHostTransmission.anophEntoAvailability[i];

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

                        if(Transmission::PerHostAnophParams::numSpecies(ctx) > 1)
                            throw std::runtime_error("only supports one mosquito species");

                        const double ux = Transmission::PerHostAnophParams::get(ctx, 0).entoAvailability->cdf(human.perHostTransmission.anophEntoAvailabilityRaw[0]);
                        if(ux < 1)
                        {
                            // parity: a — per-human heterogeneity lognormal/beta draws (gsl_cdf pure quantiles, no stream cost) then bernoulli-last; per-candidate draw order/count frozen (plan §3a)
                            double xx = gsl_cdf_ugaussian_Pinv(ux);                      // Unit interval to Normal

                            // Apply the Gaussian copula transformation for correlation
                            const ComponentId cid = subPop;

                            /** human.rng().gauss(0.0, 1.0) is calculated once per component and per human. Re-deployment of
                             * the same component on the same human must use the same gaussian sample. Therefore we store
                             * this value in the human the first time it is calculated. */
                            double g = 0.0;
                            const auto it = human.perHostTransmission.copulaGaussianSamples.find(cid);
                            if (it != human.perHostTransmission.copulaGaussianSamples.end())
                                g = it->second;
                            else
                            {
                                g = human.rng().gauss(0.0, 1.0);
                                human.perHostTransmission.copulaGaussianSamples[cid] = g;
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
                            probability = human.rng().beta(alpha, beta);
                    }
                    catch (const std::exception& e) {
                        std::ostringstream oss;
                        oss << "[Gaussian copula]: computing correlated intervention probability using Gaussian copula and Beta distribution: "
                            << e.what()
                            << ". Possible reason: coverage might be too high and/or population size might be too small.";
                        throw std::runtime_error(oss.str());
                    }
                }

                if( availability >= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, minAvailability) && availability <= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, maxAvailability) ) {
                    if( subPop == ComponentId::wholePop() || (human.isInSubPop( subPop ) != complement) ){
                        // parity: a — deploy short-circuit: RNG draw is the LAST test of the deploy condition; filter/age/cohort checks must precede it (plan §3a)
                        if( human.rng().bernoulli( probability ) ){
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
    TimedCumulativeHumanDeployment( RunContext& ctx, SimTime date,
                           const scnXml::MassDeployment& mass,
                           shared_ptr<const HumanIntervention> intervention,
                           ComponentId subPop, bool complement,
                           ComponentId cumCuvId ) :
        TimedHumanDeployment( ctx, date, mass, intervention, subPop, complement ),
        cumCovInd( cumCuvId )
    {
    }
    
    virtual void deploy (RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel& transmission) {
        // Cumulative case: bring target group's coverage up to target coverage
        vector<Host::Human*> unprotected;
        size_t total = 0;       // number of humans within age bound and optionally subPop
        for(Host::Human &human : population) {
            SimTime age = human.age(ctx.clock.now());
            if( age >= minAge && age < maxAge ){
                double availability = 0.0;
                for(size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(ctx); ++i)
                    availability += human.perHostTransmission.anophEntoAvailability[i];
                if( availability >= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, minAvailability) && availability <= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, maxAvailability) ) {
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
                // parity: a — deploy short-circuit: RNG draw is the LAST test of the deploy condition; filter/age/cohort checks must precede it (plan §3a)
                if( human->rng().uniform_01() < additionalCoverage ){
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
    TimedVectorDeployment( RunContext& ctx, SimTime date, size_t instance ) :
        TimedDeployment( ctx, date ),
        inst(instance)
    {}
    virtual void deploy (RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel& transmission) {
        Transmission::TransmissionModel *vectorModel = transmission.isVector() ? &transmission : nullptr;
        if(vectorModel)
            vectorModel->deployVectorPopInterv(inst);
    }
    virtual void print_details( std::ostream& out )const{
        out << "Vector";
    }
    
private:
    size_t inst;
};

class TimedTrapDeployment : public TimedDeployment {
public:
    TimedTrapDeployment( RunContext& ctx, SimTime date, size_t instance, double ratio, SimTime lifespan ) :
        TimedDeployment( ctx, date ), inst(instance), ratio(ratio), lifespan(lifespan)
    {}
    virtual void deploy (RunContext& ctx, vector<Host::Human>& population, Transmission::TransmissionModel& transmission) {
        Transmission::TransmissionModel *vectorModel = transmission.isVector() ? &transmission : nullptr;
        if(vectorModel)
        {
            double number = population.size() * ratio;
            vectorModel->deployVectorTrap(inst, number, lifespan);
        }
    }
    virtual void print_details( std::ostream& out )const{
        out << "VectorTrap:" << endl;
        out << "\tratio: " << ratio << "\tlifespan: " << lifespan; 
    }
    
private:
    size_t inst;
    double ratio;
    SimTime lifespan = sim::never();
};

// ———  ContinuousHumanDeployment  ———

/** Interface for continuous deployment of an intervention. */
class ContinuousHumanDeployment : protected HumanDeploymentBase {
public:
    /// Create, passing deployment age
    ContinuousHumanDeployment( RunContext& ctx, SimTime begin, SimTime end,
                                 const ::scnXml::ContinuousDeployment& elt,
                                 shared_ptr<const HumanIntervention> intervention,
                                 ComponentId subPop, bool complement );
    
    /** Apply filters and potentially deploy.
     * 
     * @returns false iff this deployment (and thus all later ones in the
     *  ordered list) happens in the future. */
    bool filterAndDeploy( RunContext& ctx, Host::Human& human ) const{
        SimTime age = human.age(ctx.clock.now());
        if( deployAge > age ){
            // stop processing continuous deployments for this
            // human for now because remaining ones happen in the future
            return false;
        }else if( deployAge == age ){
            double availability = 0.0;
            for(size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(ctx); ++i)
                availability += human.perHostTransmission.anophEntoAvailabilityRaw[i];

            if( availability >= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, minAvailability) && availability <= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, maxAvailability) )
            {
                auto now = ctx.clock.interv_date();
                if( begin <= now && now < end &&
                    ( subPop == ComponentId::wholePop() ||
                        (human.isInSubPop( subPop ) != complement)
                    ) &&
                    // parity: a — deploy short-circuit: RNG draw is the LAST test of the deploy condition; filter/age/cohort checks must precede it (plan §3a)
                    human.rng().uniform_01() < coverage )     // RNG call should be last test
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
};
// inline variable: Deployments.h is now included by more than one source
inline const ByDeployTime byDeployTime{};

} }

#endif
