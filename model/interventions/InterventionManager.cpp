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

#include "interventions/InterventionManager.h"
#include "schema/interventions.h"
#include "Population.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "interventions/GVI.h"
#include "interventions/IRS.h"
#include "interventions/ITN.h"
#include "interventions/Vaccine.h"
#include "util/ModelOptions.h"
#include "Clinical/ESCaseManagement.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Host/Human.h"
#include "Transmission/TransmissionModel.h"
#include "Host/WithinHost/Diagnostic.h"
#include "PkPd/LSTMTreatments.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include "interventions/Deployments.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Transmission/VectorModel.h"
#include "Clinical/ClinicalModel.h"
#include "Model.h"
#include "Transmission/PerHost.h"
#include "util/SpeciesIndexChecker.h"

namespace OM
{
namespace interventions
{
// ———  component classes and per-deployment helpers  ———
using Host::Human;
using WithinHost::Diagnostic;
using WithinHost::diagnostics;

// ———  HumanIntervention  ———

bool componentCmp(const HumanInterventionComponent *a, const HumanInterventionComponent *b){
    return a->componentType() < b->componentType();
}
HumanIntervention::HumanIntervention(
    const xsd::cxx::tree::sequence<scnXml::Component>& componentList,
    const InterventionRegistry& reg )
{
    components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it ) {
        ComponentId id = getComponentId( reg, it->getId() );
        const HumanInterventionComponent* component = &getComponent( reg, id );
        components.push_back( component );
    }
    
    /** Sort components according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort( components.begin(), components.end(), componentCmp );
}
HumanIntervention::HumanIntervention(
    const xsd::cxx::tree::sequence<scnXml::Component>& componentList,
    const xsd::cxx::tree::sequence<scnXml::Condition>& conditionList,
    const InterventionRegistry& reg )
{
    components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it ) {
        ComponentId id = getComponentId( reg, it->getId() );
        const HumanInterventionComponent* component = &getComponent( reg, id );
        components.push_back( component );
    }
    
    /** Sort components according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort( components.begin(), components.end(), componentCmp );
    
    for( auto it = conditionList.begin(), end = conditionList.end(); it != end; ++it ) {
        double min = it->getMinValue().present() ? it->getMinValue().get() : -numeric_limits<double>::infinity();
        double max = it->getMaxValue().present() ? it->getMaxValue().get() : numeric_limits<double>::infinity();
        conditions.push_back(mon::setupCondition( it->getMeasure(), min, max, it->getInitialState() ));
    }
}
HumanIntervention::HumanIntervention( const xsd::cxx::tree::sequence<scnXml::DTDeploy>& componentList,
    const InterventionRegistry& reg ){
    components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it ) {
        ComponentId id = getComponentId( reg, it->getComponent() );
        const HumanInterventionComponent* component = &getComponent( reg, id );
        components.push_back( component );
    }
    
    /** Sort components according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort( components.begin(), components.end(), componentCmp );
}

void HumanIntervention::deploy( Human& human, mon::Deploy::Method method,
    VaccineLimits vaccLimits ) const
{
    for( size_t cond : conditions ){
        // Abort deployment if any condition is false.
        if( !mon::checkCondition(cond) ) return;
    }
    
    for( auto it = components.begin(); it != components.end(); ++it ) {
        const interventions::HumanInterventionComponent& component = **it;
        // we must report first, since it can change cohort and sub-population
        // which may affect what deployment does (at least in the case of reporting deployments)
        human.addToCohort(component.id(), component.duration());
        component.deploy( human, method, vaccLimits );
    }
}

void HumanIntervention::print_details( std::ostream& out )const{
    out << "human:";
    for( auto it = components.begin(); it != components.end(); ++it ){
        out << '\t' << (*it)->id().id;
    }
}

// ———  Utilities  ———

TriggeredDeployments::TriggeredDeployments(const scnXml::TriggeredDeployments& elt, const InterventionRegistry& reg){
    lists.reserve( elt.getDeploy().size() );
    for( auto it = elt.getDeploy().begin(), end = elt.getDeploy().end(); it != end; ++it ) {
        lists.push_back( SubList( *it, reg ) );
    }
}
void TriggeredDeployments::deploy(Human& human,
        mon::Deploy::Method method, VaccineLimits vaccLimits) const
{
    for( auto it = lists.begin(); it != lists.end(); ++it ){
        it->deploy( human, method, vaccLimits );
    }
}

TriggeredDeployments::SubList::SubList( const scnXml::TriggeredDeployments::DeployType& elt, const InterventionRegistry& reg ) :
        HumanIntervention( elt.getComponent(), reg ),
        minAge( sim::fromYearsN( elt.getMinAge() ) ),
        maxAge( sim::future() ),
        coverage( elt.getP() )
{
    if( elt.getMaxAge().present() )
        maxAge = sim::fromYearsN( elt.getMaxAge().get() );
    
    if( minAge < sim::zero() || maxAge < minAge ){
        throw util::xml_scenario_error("triggered intervention must have 0 <= minAge <= maxAge");
    }
    
    if( coverage < 0.0 || coverage > 1.0 ){
        throw util::xml_scenario_error( "triggered intervention must have 0 <= coverage <= 1" );
    }
    
    // Zero coverage: optimise
    if( coverage <= 0.0 || minAge >= maxAge ) components.clear();
}
void TriggeredDeployments::SubList::deploy( Host::Human& human,
        mon::Deploy::Method method, VaccineLimits vaccLimits )const
{
    // This may be used within a time step; in that case we use human age at
    // the beginning of the step since this gives the expected result with age
    // limits of 0 to 1 time step.
    //TODO: why does the above contradict what we do? Is this due to a date being shifted later?
    SimTime age = human.age(sim::nowOrTs1());
    if( age >= minAge && age < maxAge ){
        if( coverage >= 1.0 || human.rng.bernoulli( coverage ) ){
            HumanIntervention::deploy( human, method, vaccLimits );
        }
    }
}


// ———  Derivatives of HumanInterventionComponent  ———

class RecruitmentOnlyComponent : public HumanInterventionComponent {
public:
    RecruitmentOnlyComponent( ComponentId id ) :
        HumanInterventionComponent(id)
    {}
    virtual ~RecruitmentOnlyComponent() {}
    
    /// Reports to monitoring, nothing else
    virtual void deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits ) const{
        mon::recordDeploy(mon::measure("nMassRecruitOnly"), mon::measure("nCtsRecruitOnly"), human, method);
    }
    
    virtual Component::Type componentType() const{
        return Component::RECRUIT_ONLY;
    }
    
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tRecruit only";
    }
};

class ScreenComponent : public HumanInterventionComponent {
public:
    ScreenComponent( ComponentId id, const scnXml::Screen& elt, const InterventionRegistry& reg ) :
        HumanInterventionComponent(id),
        diagnostic(diagnostics::get(elt.getDiagnostic())),
        positive( elt.getPositive(), reg ),
        negative( elt.getNegative(), reg )
    {}
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits vaccLimits ) const{
        mon::recordDeploy(mon::measure("nMassScreenings"), mon::measure("nCtsScreenings"), human, method);
        if( human.withinHostModel->diagnosticResult(human.rng, diagnostic) ){
            positive.deploy( human, method, vaccLimits );
        }else{
            negative.deploy( human, method, vaccLimits );
        }
    }
    
    virtual Component::Type componentType() const{ return Component::SCREEN; }
    
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tScreen";
    }
    
public:
    const Diagnostic& diagnostic;
    HumanIntervention positive, negative;
};

/// Simple treatment: no PK/PD, just remove parasites
class TreatSimpleComponent : public HumanInterventionComponent {
public:
    TreatSimpleComponent( ComponentId id, const scnXml::DTTreatSimple& elt ) :
        HumanInterventionComponent(id)
    {
        //NOTE: this code is currently identical to that in CMDTTreatSimple
        try{
            SimTime durL = UnitParse::readShortDuration( elt.getDurationLiver(), UnitParse::NONE ),
                durB = UnitParse::readShortDuration( elt.getDurationBlood(), UnitParse::NONE );
            SimTime neg1 = -sim::oneTS();
            if( durL < neg1 || durB < neg1 ){
                throw util::xml_scenario_error( "treatSimple: cannot have durationBlood or durationLiver less than -1" );
            }
            if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
                if( durL != sim::zero() || durB != neg1 )
                    throw util::unimplemented_exception( "vivax model only supports timestepsLiver=0, timestepsBlood=-1" );
                // Actually, the model ignores these parameters; we just don't want somebody thinking it doesn't.
            }
            timeLiver = durL;
            timeBlood = durB;
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("treatSimple: ").append(e.message()) );
        }
    }
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits ) const{
        mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        human.withinHostModel->treatSimple( human, timeLiver, timeBlood );
    }
    
    virtual Component::Type componentType() const{ return Component::TREAT_SIMPLE; }
    
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\ttreatSimple";
    }
    
public:
    SimTime timeLiver, timeBlood;
};

//NOTE: this is similar to CMDTTreatPKPD, but supporting only a single "treatment"
class TreatPKPDComponent : public HumanInterventionComponent {
public:
    TreatPKPDComponent( ComponentId id, const scnXml::DTTreatPKPD& elt ) :
        HumanInterventionComponent(id),
            schedule(PkPd::LSTMTreatments::findSchedule(elt.getSchedule())),
            dosage(PkPd::LSTMTreatments::findDosages(elt.getDosage())),
            delay_d(elt.getDelay_h() / 24.0 /* convert to days */)
    {
    }
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits ) const{
        mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        double age = sim::inYears(human.age(sim::nowOrTs1()));
        human.withinHostModel->treatPkPd( schedule, dosage, age, delay_d );
    }
    
    virtual Component::Type componentType() const{ return Component::TREAT_PKPD; }
    
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\ttreatPKPD";
    }
    
public:
    size_t schedule;        // index of the schedule
    size_t dosage;          // name of the dosage table
    double delay_d;         // delay in days
};

class DecisionTreeComponent : public HumanInterventionComponent {
public:
    DecisionTreeComponent( ComponentId id, const scnXml::DecisionTree& elt, const InterventionRegistry& reg ) :
        HumanInterventionComponent(id),
        tree(Clinical::CMDecisionTree::create(elt, false, reg))
    {}
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits )const{
        Clinical::CMDTOut out = tree.exec( Clinical::CMHostData(human,
                sim::inYears(human.age(sim::nowOrTs1())),
                human.clinicalModel->getLatestState()) );
        if( out.treated ){
            mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        }
        if( out.screened ){
            mon::recordDeploy(mon::measure("nMassScreenings"), mon::measure("nCtsScreenings"), human, method);
        }
    }
    
    virtual Component::Type componentType() const{ return Component::CM_DT; }
    
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tDecision tree";
    }
    
public:
    const Clinical::CMDecisionTree& tree;
};

class ClearImmunityComponent : public HumanInterventionComponent {
public:
    ClearImmunityComponent( ComponentId id ) :
        HumanInterventionComponent(id) {}
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits )const{
        human.withinHostModel->clearImmunity();
        //NOTE: not reported
    }
    
    virtual Component::Type componentType() const{ return Component::CLEAR_IMMUNITY; }
    
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tclear immunity";
    }
};



// ———  InterventionManager  ———

// ———  Deployments: free functions over the data structs  ———

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

SimTime checkDeployDate(SimTime date) {
    if( date < sim::startDate() ){
        throw util::xml_scenario_error("timed intervention deployment: may not be negative");
    }else if( date >= sim::endDate() ){
        cerr << "Warning: timed intervention deployment at date "<< date;
        cerr << " happens after last survey" << endl;
    }
    return date;
}

HumanDeployment::HumanDeployment(const scnXml::DeploymentBase& deploy,
                                 std::shared_ptr<const HumanIntervention> intervention,
                                 ComponentId subPop, bool complement) :
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

TimedDeployment::Human::Human(const scnXml::MassDeployment& mass,
                              std::shared_ptr<const HumanIntervention> intervention,
                              ComponentId subPop, bool complement) :
    base( mass, intervention, subPop, complement ),
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

TimedDeployment::CumulativeHuman::CumulativeHuman(const scnXml::MassDeployment& mass,
                              std::shared_ptr<const HumanIntervention> intervention,
                              ComponentId subPop, bool complement, ComponentId cumCovId) :
    human( mass, intervention, subPop, complement ), cumCovId( cumCovId )
{
}

ContinuousHumanDeployment::ContinuousHumanDeployment(SimTime begin, SimTime end,
                                 const ::scnXml::ContinuousDeployment& elt,
                                 std::shared_ptr<const HumanIntervention> intervention,
                                 ComponentId subPop, bool complement) :
        base( elt, intervention, subPop, complement ),
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

const string vec_mode_err = "vector interventions can only be used in dynamic transmission mode (mode=\"dynamic\")";

namespace {
/** Is the human's total mosquito availability within [min,max] percentiles?
 * Timed deployments test the current (intervention-adjusted) availability,
 * continuous deployments the raw baseline value (behaviour kept as-is). */
bool inAvailabilityRange(const Host::Human& human, bool raw, SimTime minAvail, SimTime maxAvail,
                         const Transmission::VectorSpeciesParams& species_params)
{
    const vector<double>& values = raw ? human.perHostTransmission.anophEntoAvailabilityRaw
                                       : human.perHostTransmission.anophEntoAvailability;
    double availability = 0.0;
    for (size_t i = 0; i < species_params.numSpecies(); ++i)
        availability += values[i];
    return availability >= species_params.getEntoAvailabilityPercentile(minAvail)
        && availability <= species_params.getEntoAvailabilityPercentile(maxAvail);
}
}

void apply(TimedDeployment& dep, vector<Host::Human>& population, Transmission::TransmissionModel& transmission,
    const InterventionRegistry& reg)
{
    struct Action {
        vector<Host::Human>& population;
        Transmission::TransmissionModel& transmission;
        const InterventionRegistry& reg;

        void operator() (std::monostate) const {}  // list-end marker
        void operator() (TimedDeployment::ChangeHS& a) const {
            Clinical::ClinicalModel::setHS( *a.newHS, reg );
            a.newHS.reset();
        }
        void operator() (TimedDeployment::ChangeEIR& a) const {
            transmission.changeEIRIntervention( *a.newEIR );
            a.newEIR.reset();
        }
        void operator() (const TimedDeployment::UninfectVectors&) const {
            transmission.uninfectVectors();
        }
        void operator() (const TimedDeployment::Human& a) const {
            for(Host::Human& human : population) {
                const SimTime age = human.age(sim::now());
                if( age >= a.minAge && age < a.maxAge ){
                    double probability = a.base.coverage;

                    if(a.base.copula && a.base.coverage > 0)
                    {
                        try
                        {
                            // Beta distribution for intervention (Beta distributed)
                            const double beta_mean = a.base.coverage;  // Assuming this value is given
                            const double beta_var = a.base.coverageVar;       // Given as well
                            const double alpha = ((1.0 - beta_mean) / beta_var - 1.0 / beta_mean) * (beta_mean * beta_mean);
                            const double beta = alpha * (1.0 / beta_mean - 1.0);

                            if(transmission.speciesParams.numSpecies() > 1)
                                throw std::runtime_error("only supports one mosquito species");

                            const double ux = transmission.speciesParams.get(0).entoAvailability->cdf(human.perHostTransmission.anophEntoAvailabilityRaw[0]);
                            if(ux < 1)
                            {
                                double xx = gsl_cdf_ugaussian_Pinv(ux);                      // Unit interval to Normal

                                // Apply the Gaussian copula transformation for correlation
                                const ComponentId cid = a.base.subPop;

                                /** human.rng.gauss(0.0, 1.0) is calculated once per component and per human. Re-deployment of
                                 * the same component on the same human must use the same gaussian sample. Therefore we store
                                 * this value in the human the first time it is calculated. */
                                double g = 0.0;
                                const auto it = human.perHostTransmission.copulaGaussianSamples.find(cid);
                                if (it != human.perHostTransmission.copulaGaussianSamples.end())
                                    g = it->second;
                                else
                                {
                                    g = human.rng.gauss(0.0, 1.0);
                                    human.perHostTransmission.copulaGaussianSamples[cid] = g;
                                }
                                const double yy = a.base.coverageCorr * xx + g * sqrt(1 - a.base.coverageCorr * a.base.coverageCorr);

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

                    if( inAvailabilityRange(human, false, a.minAvailability, a.maxAvailability, transmission.speciesParams) ) {
                        if( a.base.subPop == ComponentId::wholePop() || (human.isInSubPop( a.base.subPop ) != a.base.complement) ){
                            if( human.rng.bernoulli( probability ) ){
                                a.base.intervention->deploy( human, mon::Deploy::TIMED, a.base.vaccLimits );
                            }
                        }
                    }
                }
            }
        }
        void operator() (const TimedDeployment::CumulativeHuman& a) const {
            const TimedDeployment::Human& h = a.human;
            // Cumulative case: bring target group's coverage up to target coverage
            vector<Host::Human*> unprotected;
            size_t total = 0;       // number of humans within age bound and optionally subPop
            for(Host::Human &human : population) {
                SimTime age = human.age(sim::now());
                if( age >= h.minAge && age < h.maxAge ){
                    if( inAvailabilityRange(human, false, h.minAvailability, h.maxAvailability, transmission.speciesParams) ) {
                        if( h.base.subPop == ComponentId::wholePop() || (human.isInSubPop( h.base.subPop ) != h.base.complement) ){
                            total+=1;
                            if( !human.isInSubPop(a.cumCovId) )
                                unprotected.push_back( &human );
                        }
                    }
                }
            }

            if( total == 0 ) return;        // no humans to deploy to; avoid divide by zero
            double propProtected = static_cast<double>( total - unprotected.size() ) / static_cast<double>( total );
            if( propProtected < h.base.coverage ){
                // Proportion propProtected are already covered, so need to
                // additionally cover the proportion (coverage - propProtected),
                // selected from the list unprotected.
                double additionalCoverage = (h.base.coverage - propProtected) / (1.0 - propProtected);
                cerr << "cum deployment: prop protected " << propProtected << "; additionalCoverage " << additionalCoverage << "; total " << total << endl;
                for(Host::Human* human : unprotected) {
                    if( human->rng.uniform_01() < additionalCoverage ){
                        h.base.intervention->deploy( *human, mon::Deploy::TIMED, h.base.vaccLimits );
                    }
                }
            }
        }
        void operator() (const TimedDeployment::VectorPop& a) const {
            Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
            if(vectorModel)
                vectorModel->deployVectorPopInterv(a.instance);
        }
        void operator() (const TimedDeployment::Trap& a) const {
            Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
            if(vectorModel)
            {
                double number = population.size() * a.ratio;
                vectorModel->deployVectorTrap(a.instance, number, a.lifespan);
            }
        }
        void operator() (const TimedDeployment::ReduceNhh& a) const {
            Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
            if(vectorModel)
            {
                for (size_t i = 0; i < vectorModel->speciesIndex.size(); ++i)
                {
                    Transmission::Anopheles::AnophelesModel *anophModel = vectorModel->species[i].get();

                    if (anophModel->nhhInstances.count(a.intervName) == 0)
                        throw util::xml_scenario_error("non human hosts type " + a.intervName + " not deployed during non human hosts intervention deployment");

                    anophModel->reduceNhhAvailability[a.intervName][a.instance].deploy(vectorModel->m_rng, sim::now());
                    anophModel->reduceP_B_I[a.intervName][a.instance].deploy(vectorModel->m_rng, sim::now());
                    anophModel->reduceP_C_I[a.intervName][a.instance].deploy(vectorModel->m_rng, sim::now());
                    anophModel->reduceP_D_I[a.intervName][a.instance].deploy(vectorModel->m_rng, sim::now());
                    anophModel->reduceFecundity[a.intervName][a.instance].deploy(vectorModel->m_rng, sim::now());
                }
            }
        }
        void operator() (TimedDeployment::AddNhh& a) const {
            Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
            if(vectorModel)
            {
                double popSize = population.size();
                if (vectorModel->interventionMode != Transmission::SimulationMode::dynamicEIR) { throw util::xml_scenario_error(vec_mode_err); }
                for (size_t i = 0; i < vectorModel->speciesIndex.size(); ++i)
                {
                    Transmission::Anopheles::AnophelesModel *anophModel = vectorModel->species[i].get();

                    if (anophModel->nhhInstances.count(a.intervName) != 0)
                        throw util::xml_scenario_error("non human hosts type " + a.intervName + " already deployed during non human hosts deployment");

                    const TimedDeployment::AddNhh::SpeciesParams &p = a.params[anophModel->mosq.name];
                    const double adultAvail = transmission.speciesParams.get(i).entoAvailability->mean();
                    const double avail_i = transmission.speciesParams.get(i).entoAvailabilityFactor * popSize * adultAvail * p.mosqRelativeAvailabilityHuman;
                    Transmission::Anopheles::Nhh nhh;

                    nhh.avail_i = avail_i;
                    nhh.P_B_I = p.mosqProbBiting;
                    nhh.P_C_I = p.mosqProbFindingRestSite;
                    nhh.P_D_I = p.mosqProbResting;
                    nhh.rel_fecundity = p.hostFecundityFactor;
                    nhh.expiry = sim::now() + a.lifespan;

                    // add the nhh to the active nhh instances
                    anophModel->nhhInstances[a.intervName] = nhh;
                }
            }
        }
    };
    std::visit(Action{population, transmission, reg}, dep.action);
}

static void initReduceNhh(size_t instance, const std::string &intervName,
                   const scnXml::Description2::AnophelesSequence list,
                   const scnXml::DecayFunction &decay, Transmission::TransmissionModel& transmission)
{
    Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
    if(vectorModel)
    {
        util::SpeciesIndexChecker checker(intervName, vectorModel->speciesIndex);
        for (const scnXml::NonHumanHostsSpeciesIntervention &anoph : list)
        {
            const string &mosqName = anoph.getMosquito();
            size_t i = checker.getIndex(mosqName);

            Transmission::Anopheles::AnophelesModel *anophModel = vectorModel->species[i].get();

            if (anophModel->reduceNhhAvailability[intervName].size() <= instance) anophModel->reduceNhhAvailability[intervName].resize(instance + 1);
            if (anophModel->reduceP_B_I[intervName].size() <= instance) anophModel->reduceP_B_I[intervName].resize(instance + 1);
            if (anophModel->reduceP_C_I[intervName].size() <= instance) anophModel->reduceP_C_I[intervName].resize(instance + 1);
            if (anophModel->reduceP_D_I[intervName].size() <= instance) anophModel->reduceP_D_I[intervName].resize(instance + 1);
            if (anophModel->reduceFecundity[intervName].size() <= instance) anophModel->reduceFecundity[intervName].resize(instance + 1);

            if (anoph.getAvailabilityReduction().present())
            {
                const scnXml::AvailabilityReduction &decayFunc = anoph.getAvailabilityReduction().get();
                if (decayFunc.getInitial() > 1.0) throw util::xml_scenario_error("availabilityReduction intervention: initial effect must be <= 1");
                anophModel->reduceNhhAvailability[intervName][instance].set(decayFunc.getInitial(), decay, "availabilityReduction");
            }
            if (anoph.getPreprandialKillingEffect().present())
            {
                const scnXml::PreprandialKillingEffect &decayFunc = anoph.getPreprandialKillingEffect().get();
                if (decayFunc.getInitial() < 0 || decayFunc.getInitial() > 1)
                    throw util::xml_scenario_error("PreprandialKillingEffect intervention: initial effect must be between 0 and 1");
                anophModel->reduceP_B_I[intervName][instance].set(decayFunc.getInitial(), decay, "reduceP_B_I");
            }
            if (anoph.getPostprandialKillingEffect().present())
            {
                const scnXml::PostprandialKillingEffect &decayFunc = anoph.getPostprandialKillingEffect().get();
                if (decayFunc.getInitial() < 0 || decayFunc.getInitial() > 1)
                    throw util::xml_scenario_error("PostprandialKillingEffect intervention: initial effect must be between 0 and 1");
                anophModel->reduceP_C_I[intervName][instance].set(decayFunc.getInitial(), decay, "reduceP_C_I");
            }
            if (anoph.getRestingKillingEffect().present())
            {
                const scnXml::RestingKillingEffect &decayFunc = anoph.getRestingKillingEffect().get();
                if (decayFunc.getInitial() < 0 || decayFunc.getInitial() > 1)
                    throw util::xml_scenario_error("RestingKillingEffect intervention: initial effect must be be between 0 and 1");
                anophModel->reduceP_D_I[intervName][instance].set(decayFunc.getInitial(), decay, "reduceP_D_I");
            }
            if (anoph.getFecundityReduction().present())
            {
                const scnXml::FecundityReduction &decayFunc = anoph.getFecundityReduction().get();
                if (decayFunc.getInitial() < 0 || decayFunc.getInitial() > 1)
                    throw util::xml_scenario_error("FecundityReduction intervention: initial effect must be be between 0 and 1");
                anophModel->reduceFecundity[intervName][instance].set(decayFunc.getInitial(), decay, "reduceFecundity");
            }
        }
        checker.checkNoneMissed();
    }
}

static TimedDeployment::AddNhh readAddNhh(const std::string &intervName, SimTime lifespan,
                                   const scnXml::Description3::AnophelesSequence list,
                                   Transmission::TransmissionModel& transmission)
{
    TimedDeployment::AddNhh result;
    result.intervName = intervName;
    result.lifespan = lifespan;
    Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
    if(vectorModel)
    {
        util::SpeciesIndexChecker checker(intervName, vectorModel->speciesIndex);
        for (const scnXml::NonHumanHostsVectorSpecies &anoph : list)
        {
            const string &mosqName = anoph.getMosquito();
            checker.getIndex(mosqName);

            TimedDeployment::AddNhh::SpeciesParams &p = result.params[mosqName];
            p.mosqRelativeAvailabilityHuman = anoph.getMosqRelativeAvailabilityHuman().getValue();
            p.mosqProbBiting = anoph.getMosqProbBiting().getValue();
            p.mosqProbFindingRestSite = anoph.getMosqProbFindRestSite().getValue();
            p.mosqProbResting = anoph.getMosqProbResting().getValue();
            p.hostFecundityFactor = anoph.getHostFecundityFactor().getValue();
        }
        checker.checkNoneMissed();
    }
    return result;
}

bool filterAndDeploy(const ContinuousHumanDeployment& dep, Host::Human& human,
                     const Transmission::VectorSpeciesParams& species_params)
{
    SimTime age = human.age(sim::now());
    if( dep.deployAge > age ){
        // stop processing continuous deployments for this
        // human for now because remaining ones happen in the future
        return false;
    }else if( dep.deployAge == age ){
        if( inAvailabilityRange(human, true, dep.minAvailability, dep.maxAvailability, species_params) )
        {
            auto now = sim::intervDate();
            if( dep.begin <= now && now < dep.end &&
                ( dep.base.subPop == ComponentId::wholePop() ||
                    (human.isInSubPop( dep.base.subPop ) != dep.base.complement)
                ) &&
                human.rng.uniform_01() < dep.base.coverage )     // RNG call should be last test
            {
                dep.base.intervention->deploy( human, mon::Deploy::CTS, dep.base.vaccLimits );
            }
        }
    }//else: for some reason, a deployment age was missed; ignore it
    return true;
}

void print_details(std::ostream& out, const TimedDeployment& dep)
{
    struct Print {
        std::ostream& out;
        void operator() (std::monostate) const { out << "Dummy"; }
        void operator() (const TimedDeployment::ChangeHS&) const { out << "ChangeHS"; }
        void operator() (const TimedDeployment::ChangeEIR&) const { out << "ChangeEIR"; }
        void operator() (const TimedDeployment::UninfectVectors&) const { out << "UninfectVectors"; }
        void operator() (const TimedDeployment::Human& a) const {
            out << "Human:" << endl;
            out << "	age: " << sim::inYears(a.minAge) << "y	max age: " << sim::inYears(a.maxAge) << "y	subpop: ";
            if( a.base.subPop == ComponentId::wholePop() ) out << "whole";
            else out << a.base.subPop.id;
            out << "	complement: " << a.base.complement << "	coverage: " << a.base.coverage << "	";
            a.base.intervention->print_details( out );
        }
        void operator() (const TimedDeployment::CumulativeHuman& a) const { (*this)(a.human); }
        void operator() (const TimedDeployment::VectorPop&) const { out << "Vector"; }
        void operator() (const TimedDeployment::Trap& a) const {
            out << "VectorTrap:" << endl;
            out << "	ratio: " << a.ratio << "	lifespan: " << a.lifespan;
        }
        void operator() (const TimedDeployment::ReduceNhh& a) const {
            out << "NonHumanHosts(" << a.intervName << ")";
        }
        void operator() (const TimedDeployment::AddNhh& a) const {
            out << "AddNonHumanHosts(" << a.intervName << "):" << endl;
            out << "	lifespan: " << a.lifespan;
        }
    };
    std::visit(Print{out}, dep.action);
}

void print_details(std::ostream& out, const ContinuousHumanDeployment& dep)
{
    out << "Human:" << endl;
    out << "	begin: " << dep.begin << "	end: ";
    if( dep.end == sim::future() ) out << "(none)";
    else out << dep.end;
    out << "	deployAge: " << sim::inYears(dep.deployAge) << "y	subpop: ";
    if( dep.base.subPop == ComponentId::wholePop() ) out << "whole";
    else out << dep.base.subPop.id;
    out << "	complement: " << dep.base.complement << "	coverage: " << dep.base.coverage << "	";
    dep.base.intervention->print_details( out );
}


// declared in HumanComponents.h:
vector<ComponentId> removeAtIds[SubPopRemove::NUM];

// static functions:

void init(InterventionManager& mgr, const scnXml::Interventions &intervElt, const Population &population,
          Transmission::TransmissionModel &transmission)
{
    auto& continuous = mgr.continuous;
    auto& timed = mgr.timed;
    uint32_t& nextTimed = mgr.nextTimed;
    auto& importedInfections = mgr.importedInfections;
    nextTimed = 0;

    if (intervElt.getChangeHS().present())
    {
        const scnXml::ChangeHS &chs = intervElt.getChangeHS().get();
        if (chs.getTimedDeployment().size() > 0)
        {
            // timed deployments:
            for (auto it = chs.getTimedDeployment().begin(); it != chs.getTimedDeployment().end(); ++it)
            {
                try
                {
                    SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                    timed.push_back(TimedDeployment{checkDeployDate(date), TimedDeployment::ChangeHS{
                        unique_ptr<scnXml::HealthSystem>(it->_clone())}});
                }
                catch (const util::format_error &e)
                {
                    throw util::xml_scenario_error(string("interventions/changeHS/timedDeployment/time: ").append(e.message()));
                }
            }
        }
    }
    if (intervElt.getChangeEIR().present())
    {
        const scnXml::ChangeEIR &eir = intervElt.getChangeEIR().get();
        if (eir.getTimedDeployment().size() > 0)
        {
            // timed deployments:
            for (auto it = eir.getTimedDeployment().begin(); it != eir.getTimedDeployment().end(); ++it)
            {
                try
                {
                    SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                    timed.push_back(TimedDeployment{checkDeployDate(date), TimedDeployment::ChangeEIR{
                        unique_ptr<scnXml::NonVector>(it->_clone())}});
                }
                catch (const util::format_error &e)
                {
                    throw util::xml_scenario_error(string("interventions/changeEIR/timedDeployment/time: ").append(e.message()));
                }
            }
        }
    }
    // species_index_map is not available with the non-vector model or
    // non-dynamic mode, so setting it (lazily) also checks sim mode:
    map<string, size_t> species_index_map;

    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel != nullptr) species_index_map = vectorModel->getSpeciesIndexMap();
    }

    if (intervElt.getHuman().present())
    {
        const scnXml::HumanInterventions &human = intervElt.getHuman().get();

        // 1. Read components
        for (auto it = human.getComponent().begin(), end = human.getComponent().end(); it != end; ++it)
        {
            const scnXml::HumanInterventionComponent &component = *it;
            if (mgr.registry.identifiers.count(component.getId()) > 0)
            {
                ostringstream msg;
                msg << "The id attribute of interventions.human.component elements must be unique; found \"" << component.getId()
                    << "\" twice.";
                throw util::xml_scenario_error(msg.str());
            }
            ComponentId id(mgr.registry.components.size()); // i.e. index of next item
            mgr.registry.identifiers.insert(make_pair(component.getId(), id));

            SimTime expireAfter = sim::future();
            if (component.getSubPopRemoval().present())
            {
                const scnXml::SubPopRemoval &removeOpts = component.getSubPopRemoval().get();
                if (removeOpts.getOnFirstBout()) { removeAtIds[SubPopRemove::ON_FIRST_BOUT].push_back(id); }
                if (removeOpts.getOnFirstInfection()) { removeAtIds[SubPopRemove::ON_FIRST_INFECTION].push_back(id); }
                if (removeOpts.getOnFirstTreatment()) { removeAtIds[SubPopRemove::ON_FIRST_TREATMENT].push_back(id); }
                if (removeOpts.getAfterYears().present()) { expireAfter = sim::fromYearsN(removeOpts.getAfterYears().get()); }
            }

            HumanInterventionComponent *hiComponent;
            if (component.getScreen().present()) { hiComponent = new ScreenComponent(id, component.getScreen().get(), mgr.registry); }
            else if (component.getTreatSimple().present())
            {
                hiComponent = new TreatSimpleComponent(id, component.getTreatSimple().get());
            }
            else if (component.getTreatPKPD().present())
            {
                hiComponent = new TreatPKPDComponent(id, component.getTreatPKPD().get());
            }
            else if (component.getDecisionTree().present())
            {
                hiComponent = new DecisionTreeComponent(id, component.getDecisionTree().get(), mgr.registry);
            }
            else if (component.getPEV().present())
            {
                hiComponent = new VaccineComponent(id, component.getPEV().get(), Vaccine::PEV);
            }
            else if (component.getBSV().present())
            {
                hiComponent = new VaccineComponent(id, component.getBSV().get(), Vaccine::BSV);
            }
            else if (component.getTBV().present())
            {
                hiComponent = new VaccineComponent(id, component.getTBV().get(), Vaccine::TBV);
            }
            else if (component.getITN().present())
            {
                hiComponent = new ITNComponent(id, component.getITN().get(), species_index_map);
            }
            else if (component.getIRS().present())
            {
                hiComponent = new IRSComponent(id, component.getIRS().get(), species_index_map);
            }
            else if (component.getGVI().present())
            {
                hiComponent = new GVIComponent(id, component.getGVI().get(), species_index_map);
            }
            else if (component.getRecruitmentOnly().present())
            {
                hiComponent = new RecruitmentOnlyComponent(id);
            }
            else if (component.getClearImmunity().present())
            {
                hiComponent = new ClearImmunityComponent(id);
            }
            else
            {
                throw util::xml_scenario_error("expected intervention.human.component element to have a "
                                               "child, didn't find it (perhaps I need updating)");
            }
            hiComponent->setExpireAfter(expireAfter);
            mgr.registry.components.push_back(unique_ptr<HumanInterventionComponent>(hiComponent));
        }

        // 2. Read the list of deployments
        for (auto it = human.getDeployment().begin(), end = human.getDeployment().end(); it != end; ++it)
        {
            const scnXml::Deployment &elt = *it;
            // 2.a intervention components
            auto intervention = shared_ptr<HumanIntervention>(new HumanIntervention(elt.getComponent(), elt.getCondition(), mgr.registry));

            // 2.b intervention deployments
            for (auto ctsIt = elt.getContinuous().begin(); ctsIt != elt.getContinuous().end(); ++ctsIt)
            {
                ComponentId subPop = ComponentId::wholePop();
                bool complement = false;
                if (ctsIt->getRestrictToSubPop().present())
                {
                    const string &subPopStr = ctsIt->getRestrictToSubPop().get().getId();
                    subPop = getComponentId(mgr.registry, subPopStr);
                    complement = ctsIt->getRestrictToSubPop().get().getComplement();
                }
                const scnXml::ContinuousList::DeploySequence &ctsSeq = ctsIt->getDeploy();
                for (auto it2 = ctsSeq.begin(), end2 = ctsSeq.end(); it2 != end2; ++it2)
                {
                    try
                    {
                        SimTime begin = sim::startDate();
                        if (it2->getBegin().present())
                        {
                            begin =
                                UnitParse::readDate(it2->getBegin().get(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        }
                        SimTime end = sim::future();
                        if (it2->getEnd().present())
                        {
                            end = UnitParse::readDate(it2->getEnd().get(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        }
                        continuous.push_back(ContinuousHumanDeployment(begin, end, *it2, intervention, subPop, complement));
                    }
                    catch (const util::format_error &e)
                    {
                        throw util::xml_scenario_error(string("interventions/human/deployment/continuous/deploy: ").append(e.message()));
                    }
                }
            }
            for (auto timedIt = elt.getTimed().begin(); timedIt != elt.getTimed().end(); ++timedIt)
            {
                ComponentId subPop = ComponentId::wholePop();
                bool complement = false;
                if (timedIt->getRestrictToSubPop().present())
                {
                    const string &subPopStr = timedIt->getRestrictToSubPop().get().getId();
                    subPop = getComponentId(mgr.registry, subPopStr);
                    complement = timedIt->getRestrictToSubPop().get().getComplement();
                }
                try
                {
                    multimap<SimTime, const scnXml::MassDeployment *> deployTimes;
                    for (auto it2 = timedIt->getDeploy().begin(), end2 = timedIt->getDeploy().end(); it2 != end2; ++it2)
                    {
                        SimTime date = UnitParse::readDate(it2->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);

                        if (it2->getRepeatStep().present() != it2->getRepeatEnd().present())
                        {
                            throw util::xml_scenario_error("deploy: use of repeatStep or repeatEnd without other");
                        }
                        if (it2->getRepeatStep().present())
                        {
                            SimTime step = UnitParse::readDuration(it2->getRepeatStep().get(), UnitParse::NONE);
                            if (step < sim::oneTS()) { throw util::xml_scenario_error("deploy: repeatStep must be >= 1"); }
                            SimTime end = UnitParse::readDate(it2->getRepeatEnd().get(), UnitParse::NONE);
                            while (date < end)
                            {
                                deployTimes.insert(make_pair(date, &*it2));
                                date = date + step;
                            }
                        }
                        else
                        {
                            deployTimes.insert(make_pair(date, &*it2));
                        }
                    }
                    SimTime lastDate = sim::never();
                    for (auto deploy = deployTimes.begin(), dpEnd = deployTimes.end(); deploy != dpEnd; ++deploy)
                    {
                        if (deploy->first == lastDate)
                        {
                            ostringstream msg;
                            msg << "Timed deployment of components ";
                            bool first = true;
                            for (auto cp = elt.getComponent().begin(), cpEnd = elt.getComponent().end(); cp != cpEnd; ++cp)
                            {
                                if (!first) { msg << ", "; }
                                else
                                {
                                    first = false;
                                }
                                msg << cp->getId();
                            }
                            msg << " has multiple deplyoments at date " << lastDate << " (step " << sim::inSteps(lastDate - sim::startDate())
                                << ')';
                            throw util::xml_scenario_error(msg.str());
                        }
                        lastDate = deploy->first;
                    }

                    if (timedIt->getCumulativeCoverage().present())
                    {
                        const scnXml::CumulativeCoverage &cumCov = timedIt->getCumulativeCoverage().get();
                        ComponentId cumCovComponent = getComponentId(mgr.registry, cumCov.getComponent());
                        for (auto deploy = deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy)
                        {
                            timed.push_back(TimedDeployment{checkDeployDate(deploy->first), TimedDeployment::CumulativeHuman(
                                *deploy->second, intervention, subPop, complement, cumCovComponent)});
                        }
                    }
                    else
                    {
                        for (auto deploy = deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy)
                        {
                            timed.push_back(TimedDeployment{checkDeployDate(deploy->first), TimedDeployment::Human(
                                *deploy->second, intervention, subPop, complement)});
                        }
                    }
                }
                catch (const util::format_error &e)
                {
                    throw util::xml_scenario_error(string("interventions/human/deployment/timed/deploy/time: ").append(e.message()));
                }
            }
        }
    }
    if (intervElt.getImportedInfections().present())
    {
        const scnXml::ImportedInfections &ii = intervElt.getImportedInfections().get();
        importedInfections.init(ii);
    }
    // Must come after vaccines are initialised:
    if (intervElt.getInsertR_0Case().present())
    {
        throw util::xml_scenario_error("R_0 code is disabled to reduce "
                                       "maintenance. If you need it, either use an old OpenMalaria "
                                       "version (pre 32) or request it be reenabled.");
    }
    if (intervElt.getUninfectVectors().present())
    {
        const scnXml::UninfectVectors &elt = intervElt.getUninfectVectors().get();
        if (elt.getTimedDeployment().size() > 0)
        {
            // timed deployments:
            for (auto it = elt.getTimedDeployment().begin(); it != elt.getTimedDeployment().end(); ++it)
            {
                SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                timed.push_back(TimedDeployment{checkDeployDate(date), TimedDeployment::UninfectVectors{}});
            }
        }
    }
    if (intervElt.getVectorPop().present())
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel)
        {
            typedef scnXml::VectorPop::InterventionSequence SeqT;
            const SeqT &seq = intervElt.getVectorPop().get().getIntervention();
            size_t instance = 0;
            for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
            {
                const scnXml::VectorIntervention &elt = *it;
                if (elt.getTimed().present())
                {
                    vectorModel->initVectorInterv(elt.getDescription().getAnopheles(), instance, elt.getName());

                    const scnXml::TimedBaseList::DeploySequence &seq = elt.getTimed().get().getDeploy();
                    for (auto it = seq.begin(); it != seq.end(); ++it)
                    {
                        SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        timed.push_back(TimedDeployment{checkDeployDate(date), TimedDeployment::VectorPop{instance}});
                    }
                    instance++;
                }
            }
        }
        else
            throw util::xml_scenario_error("vector model interventions can not be used with the non-vector model");
    }
    if (intervElt.getAddNonHumanHosts().present())
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel)
        {
            typedef scnXml::AddNonHumanHosts::NonHumanHostsSequence SeqT;
            const SeqT &seq = intervElt.getAddNonHumanHosts().get().getNonHumanHosts();
            for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
            {
                const scnXml::NonHumanHosts2 &elt = *it;
                if (elt.getTimed().present())
                {
                    for (const scnXml::Deploy2 &deploy : elt.getTimed().get().getDeploy())
                    {
                        SimTime date =
                            UnitParse::readDate(deploy.getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        SimTime lifespan = UnitParse::readDuration(deploy.getLifespan(), UnitParse::NONE);
                        timed.push_back(TimedDeployment{checkDeployDate(date), readAddNhh(elt.getName(), lifespan,
                            elt.getDescription().getAnopheles(), transmission)});
                    }
                }
            }
        }
        else
            throw util::xml_scenario_error("vector model interventions can not be used with the non-vector model");
    }
    if (intervElt.getNonHumanHostsModifications().present())
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel)
        {
            typedef scnXml::NonHumanHostsModifications::InterventionSequence SeqT;
            const SeqT &seq = intervElt.getNonHumanHostsModifications().get().getIntervention();
            size_t instance = 0;
            for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
            {
                const scnXml::NonHumanHostsIntervention &elt = *it;
                if (elt.getTimed().present())
                {
                    const scnXml::DecayFunction &decay = elt.getDecay();
                    const scnXml::TimedBaseList::DeploySequence &seq = elt.getTimed().get().getDeploy();
                    for (auto it = seq.begin(); it != seq.end(); ++it)
                    {
                        SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        const string& nhhName = elt.getNonHumanHostsName();
                        date = checkDeployDate(date);
                        initReduceNhh(instance, nhhName, elt.getDescription().getAnopheles(), decay, transmission);
                        timed.push_back(TimedDeployment{date, TimedDeployment::ReduceNhh{instance, nhhName}});
                    }
                    instance++;
                }
            }
        }
        else
            throw util::xml_scenario_error("vector model interventions can not be used with the non-vector model");
    }
    if (intervElt.getVectorTrap().present())
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if (vectorModel)
        {
            size_t instance = 0;
            for (const scnXml::VectorTrap &trap : intervElt.getVectorTrap().get().getIntervention())
            {
                vectorModel->initVectorTrap(trap.getDescription(), instance, trap.getName());
                if (trap.getTimed().present())
                {
                    for (const scnXml::Deploy1 &deploy : trap.getTimed().get().getDeploy())
                    {
                        SimTime date = UnitParse::readDate(deploy.getTime(), UnitParse::STEPS);
                        double ratio = deploy.getRatioToHumans();
                        SimTime lifespan = UnitParse::readDuration(deploy.getLifespan(), UnitParse::NONE);
                        timed.push_back(TimedDeployment{checkDeployDate(date), TimedDeployment::Trap{instance, ratio, lifespan}});
                    }
                }
                instance += 1;
            }
        }
        else
            throw util::xml_scenario_error("vector model interventions can not be used with the non-vector model");
    }

    // lists must be sorted, increasing
    // For reproducability, we need to use stable_sort, not sort.
    stable_sort(continuous.begin(), continuous.end(),
        [](const ContinuousHumanDeployment& a, const ContinuousHumanDeployment& b) { return a.deployAge < b.deployAge; });
    stable_sort(timed.begin(), timed.end(),
        [](const TimedDeployment& a, const TimedDeployment& b) { return a.date < b.date; });

    // make sure the list ends with something always in the future, so we don't
    // have to check nextTimed is within range:
    timed.push_back(TimedDeployment{sim::future(), std::monostate{}});

    if (util::CommandLine::option(util::CommandLine::PRINT_INTERVENTIONS))
    {
        cout << "Continuous deployments:" << endl << "begin\tend\tage\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for (auto it = continuous.begin(); it != continuous.end(); ++it)
        {
            print_details(std::cout, *it);
            cout << endl;
        }
        cout << "Timed deployments:" << endl << "time\tmin age\tmax age\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for (auto &deploy : timed)
        {
            print_details(std::cout, deploy);
            cout << endl;
        }
        cout << "Human components:" << endl;
        for (auto &component : mgr.registry.components)
        {
            component->print_details(cout);
            cout << endl;
        }
    }
}

const HumanInterventionComponent& getComponent( const InterventionRegistry& reg, ComponentId id ){
    if( id.id >= reg.components.size() )
        throw util::base_exception( "invalid component id" );
    return *reg.components[id.id];
}

ComponentId getComponentId(const InterventionRegistry& reg, const string textId)
{
    auto it = reg.identifiers.find(textId);
    if (it == reg.identifiers.end())
    {
        ostringstream msg;
        msg << "unable to find an intervention component with id \"" << textId
            << "\" (wrong name, no definition or used before definition?)";
        throw util::xml_scenario_error(msg.str());
    }
    return it->second;
}

void loadFromCheckpoint(Model& model)
{
    InterventionManager& mgr = model.interventions;
    vector<Host::Human>& population = model.humans.humans;
    Transmission::TransmissionModel& transmission = *model.transmission;
    auto& timed = mgr.timed;
    uint32_t& nextTimed = mgr.nextTimed;

    SimTime date = sim::intervDate();
    // We need to re-deploy changeHS and changeEIR interventions, but nothing
    // else. nextTimed should be zero so we can go through all past interventions.
    // Only redeploy those which happened before this time step.
    assert(nextTimed == 0);
    while (timed[nextTimed].date < date)
    {
        const auto& action = timed[nextTimed].action;
        if (std::holds_alternative<TimedDeployment::ChangeHS>(action) || std::holds_alternative<TimedDeployment::ChangeEIR>(action))
        {
            // Note: neither changeHS nor changeEIR interventions care what the
            // current time step is when they are deployed, so we don't need to
            // tell them the deployment date.
            apply(timed[nextTimed], population, transmission, mgr.registry);
        }
        nextTimed += 1;
    }
}

void deploy(Model& model)
{
    InterventionManager& mgr = model.interventions;
    vector<Host::Human>& population = model.humans.humans;
    Transmission::TransmissionModel& transmission = *model.transmission;
    auto& continuous = mgr.continuous;
    auto& timed = mgr.timed;
    uint32_t& nextTimed = mgr.nextTimed;
    auto& importedInfections = mgr.importedInfections;

    if (sim::intervTime() < sim::zero()) return;

    // deploy imported infections (not strictly speaking an intervention)
    importedInfections.import(population);

    // deploy timed interventions
    SimTime now = sim::intervDate();
    while (timed[nextTimed].date <= now)
    {
        if (util::CommandLine::option(util::CommandLine::VERBOSE))
        {
            cout << "deploy::" << now << "::Timed";
            print_details(cout, timed[nextTimed]);
            cout << endl;
        }

        apply(timed[nextTimed], population, transmission, mgr.registry);
        nextTimed += 1;
    }

    // deploy continuous interventions
    for (Host::Human &human : population)
    {
        uint32_t nextCtsDist = human.nextCtsDist;
        // deploy continuous interventions
        while (nextCtsDist < continuous.size())
        {
            if (util::CommandLine::option(util::CommandLine::VERBOSE))
            {
                cout << "deploy::" << now << "::Continuous";
                print_details(cout, continuous[nextCtsDist]);
                cout << endl;
            }
            if (!filterAndDeploy(continuous[nextCtsDist], human, transmission.speciesParams)) break; // deployment (and all remaining) happens in the future
            nextCtsDist = ++human.nextCtsDist;
        }
    }
}


} // namespace interventions
} // namespace OM
