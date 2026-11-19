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

// model/xml/ conversion layer — interventions module (start): imported
// infection rates. Moved verbatim from Host/ImportedInfections.h (the
// inline ImportedInfections::init became an out-of-line definition).

#include "Global.h"
#include "Host/ImportedInfections.h"
#include "Population.h"
#include "interventions/InterventionManager.h"
#include "interventions/HumanComponents.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "interventions/GVI.h"
#include "interventions/IRS.h"
#include "interventions/ITN.h"
#include "interventions/Vaccine.h"
#include "xml/HumanInterventionComponents.h"
#include "xml/entomology.h"
#include "interventions/Deployments.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Transmission/VectorModel.h"
#include "Transmission/PerHost.h"
#include "Clinical/ClinicalModel.h"
#include "mon/Monitoring.h"
#include "util/errors.h"
#include "schema/interventions.h"
#include "schema/healthSystem.h"

#include <algorithm>
#include <string>

using namespace std;

namespace OM {
namespace Host {

void ImportedInfections::init( const scnXml::ImportedInfections& iiElt )
        {
            const scnXml::ImportedInfections::TimedType& tElt = iiElt.getTimed();
            try{
                //NOTE: if changing XSD, this should not have a default unit:
                period = UnitParse::readDuration( tElt.getPeriod(), UnitParse::STEPS );
                if( period < sim::zero() ){
                    throw util::format_error( "cannot be negative" );
                }
            }catch( const util::format_error& e ){
                throw util::xml_scenario_error( string("interventions/importedInfections/timed/period: ").append(e.message()) );
            }
            rate.reserve( tElt.getRate().size() );
            try{
                for( auto it = tElt.getRate().begin(); it != tElt.getRate().end(); ++it ){
                    SimTime date = UnitParse::readDate( it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/ );
                    // convert to per-time-step, per-person
                    double rateVal = it->getValue() * sim::yearsPerStep() * (1.0 / 1000.0);
                    rate.push_back( Rate( date - sim::startDate(), rateVal ) );
                }
                sort( rate.begin(), rate.end() );
                if( rate.size() > 0 ){
                    if( period != sim::zero() && rate[0].time != sim::zero() ){
                        throw util::xml_scenario_error( "must specify rate at time zero when period is not zero" );
                    }
                    // remove useless repeated entries from list
                    double lastRate = rate[0].value;
                    for( size_t i = 1; i < rate.size(); ){
                        if( rate[i].value == lastRate ){
                            rate.erase( rate.begin() + i );
                        }else{
                            lastRate = rate[i].value;
                            ++i;
                        }
                    }
                }
            }catch( const util::format_error& e ){
                throw util::xml_scenario_error( string("interventions/importedInfections/timed/time: ").append(e.message()) );
            }
        }

}
}

namespace OM {
namespace interventions {

class TimedChangeHSDeployment : public TimedDeployment {
public:
    TimedChangeHSDeployment( SimTime date, const scnXml::ChangeHS::TimedDeploymentType& hs ) :
        TimedDeployment( date ),
        newHS( hs._clone() )
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        Clinical::ClinicalModel::setHS( *newHS );
        delete newHS;
        newHS = 0;
    }
    virtual void print_details( std::ostream& out )const{
        out << "ChangeHS";
    }
    
private:
    scnXml::HealthSystem *newHS;
};

class TimedChangeEIRDeployment : public TimedDeployment {
public:
    TimedChangeEIRDeployment( SimTime date, const scnXml::ChangeEIR::TimedDeploymentType& nv ) :
        TimedDeployment( date ),
        newEIR( nv._clone() )
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        transmission.changeEIRIntervention( *newEIR );
        delete newEIR;
        newEIR = 0;
    }
    virtual void print_details( std::ostream& out )const{
        out << "ChangeEIR";
    }
    
private:
    scnXml::NonVector *newEIR;
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

    HumanDeploymentBase::HumanDeploymentBase( const scnXml::DeploymentBase& deploy,
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

    TimedHumanDeployment::TimedHumanDeployment( SimTime date,
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

    ContinuousHumanDeployment::ContinuousHumanDeployment( SimTime begin, SimTime end,
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

class TimedNonHumanHostsDeployment : public TimedDeployment {
public:
    TimedNonHumanHostsDeployment( SimTime date, size_t instance, string intervName, const scnXml::Description2::AnophelesSequence list, const scnXml::DecayFunction &decay, Transmission::TransmissionModel& transmission) :
        TimedDeployment( date ),
        instance(instance),
        intervName(intervName)
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
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission)
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
        {
            for (size_t i = 0; i < vectorModel->speciesIndex.size(); ++i)
            {
                Transmission::Anopheles::AnophelesModel *anophModel = vectorModel->species[i].get();

                if (anophModel->nhhInstances.count(intervName) == 0)
                    throw util::xml_scenario_error("non human hosts type " + intervName + " not deployed during non human hosts intervention deployment");

                anophModel->reduceNhhAvailability[intervName][instance].deploy(vectorModel->m_rng, sim::now());
                anophModel->reduceP_B_I[intervName][instance].deploy(vectorModel->m_rng, sim::now());
                anophModel->reduceP_C_I[intervName][instance].deploy(vectorModel->m_rng, sim::now());
                anophModel->reduceP_D_I[intervName][instance].deploy(vectorModel->m_rng, sim::now());
                anophModel->reduceFecundity[intervName][instance].deploy(vectorModel->m_rng, sim::now());
            }
        }
    }
    virtual void print_details( std::ostream& out )const{
        out << "NonHumanHosts(" << intervName << ")";
    }
    
private:
    size_t instance;
    string intervName;
};
const string vec_mode_err = "vector interventions can only be used in dynamic transmission mode (mode=\"dynamic\")";

class TimedAddNonHumanHostsDeployment : public TimedDeployment {
public:
    TimedAddNonHumanHostsDeployment( SimTime date, const string &intervName, SimTime lifespan, const scnXml::Description3::AnophelesSequence list, Transmission::TransmissionModel& transmission) :
        TimedDeployment( date ), intervName(intervName), lifespan(lifespan)
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
        {
            util::SpeciesIndexChecker checker(intervName, vectorModel->speciesIndex);
            for (const scnXml::NonHumanHostsVectorSpecies &anoph : list)
            {
                const string &mosqName = anoph.getMosquito();
                checker.getIndex(mosqName);

                NhhParamsInterv nhh;
                nhh.mosqRelativeAvailabilityHuman = anoph.getMosqRelativeAvailabilityHuman().getValue();
                nhh.mosqProbBiting = anoph.getMosqProbBiting().getValue();
                nhh.mosqProbFindingRestSite = anoph.getMosqProbFindRestSite().getValue();
                nhh.mosqProbResting = anoph.getMosqProbResting().getValue();
                nhh.hostFecundityFactor = anoph.getHostFecundityFactor().getValue();
                nhhParams[mosqName] = nhh;
            }
            checker.checkNoneMissed();
        }
    }

    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission)
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
        {
            double popSize = population.size();
            if (vectorModel->interventionMode != Transmission::SimulationMode::dynamicEIR) { throw util::xml_scenario_error(vec_mode_err); }
            for (size_t i = 0; i < vectorModel->speciesIndex.size(); ++i)
            {
                Transmission::Anopheles::AnophelesModel *anophModel = vectorModel->species[i].get();
        
                if (anophModel->nhhInstances.count(intervName) != 0)
                    throw util::xml_scenario_error("non human hosts type " + intervName + " already deployed during non human hosts deployment");

                const NhhParamsInterv &p = nhhParams[anophModel->mosq.name];
                const double adultAvail = Transmission::PerHostAnophParams::get(i).entoAvailability->mean();
                const double avail_i = Transmission::PerHostAnophParams::get(i).entoAvailabilityFactor * popSize * adultAvail * p.mosqRelativeAvailabilityHuman;
                Transmission::Anopheles::Nhh nhh;

                nhh.avail_i = avail_i;
                nhh.P_B_I = p.mosqProbBiting;
                nhh.P_C_I = p.mosqProbFindingRestSite;
                nhh.P_D_I = p.mosqProbResting;
                nhh.rel_fecundity = p.hostFecundityFactor;
                nhh.expiry = sim::now() + lifespan;

                // add the nhh to the active nhh instances
                anophModel->nhhInstances[intervName] = nhh;
            }
        }
    }
    virtual void print_details( std::ostream& out )const{
        out << "AddNonHumanHosts(" << intervName << "):" << endl;
        out << "\tlifespan: " << lifespan;
    }
    
private:
    struct NhhParamsInterv {
        double mosqRelativeAvailabilityHuman;
        double mosqProbBiting;
        double mosqProbFindingRestSite;
        double mosqProbResting;
        double hostFecundityFactor;
    };
    string intervName;
    SimTime lifespan = sim::never();
    std::map<string, NhhParamsInterv> nhhParams; 
};

// ———  InterventionManager  ———

void InterventionManager::init(const scnXml::Interventions &intervElt, const Population &population,
                               Transmission::TransmissionModel &transmission)
{
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
                    timed.push_back(unique_ptr<TimedDeployment>(new TimedChangeHSDeployment(date, *it)));
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
                    timed.push_back(unique_ptr<TimedDeployment>(new TimedChangeEIRDeployment(date, *it)));
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
            if (identifierMap.count(component.getId()) > 0)
            {
                ostringstream msg;
                msg << "The id attribute of interventions.human.component elements must be unique; found \"" << component.getId()
                    << "\" twice.";
                throw util::xml_scenario_error(msg.str());
            }
            ComponentId id(humanComponents.size()); // i.e. index of next item
            identifierMap.insert(make_pair(component.getId(), id));

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
            if (component.getScreen().present()) { hiComponent = new ScreenComponent(id, component.getScreen().get()); }
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
                hiComponent = new DecisionTreeComponent(id, component.getDecisionTree().get());
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
            humanComponents.push_back(unique_ptr<HumanInterventionComponent>(hiComponent));
        }

        // 2. Read the list of deployments
        for (auto it = human.getDeployment().begin(), end = human.getDeployment().end(); it != end; ++it)
        {
            const scnXml::Deployment &elt = *it;
            // 2.a intervention components
            auto intervention = shared_ptr<HumanIntervention>(new HumanIntervention(elt.getComponent(), elt.getCondition()));

            // 2.b intervention deployments
            for (auto ctsIt = elt.getContinuous().begin(); ctsIt != elt.getContinuous().end(); ++ctsIt)
            {
                ComponentId subPop = ComponentId::wholePop();
                bool complement = false;
                if (ctsIt->getRestrictToSubPop().present())
                {
                    const string &subPopStr = ctsIt->getRestrictToSubPop().get().getId();
                    subPop = getComponentId(subPopStr);
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
                    subPop = getComponentId(subPopStr);
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
                        ComponentId cumCovComponent = getComponentId(cumCov.getComponent());
                        for (auto deploy = deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy)
                        {
                            timed.push_back(unique_ptr<TimedDeployment>(new TimedCumulativeHumanDeployment(
                                deploy->first, *deploy->second, intervention, subPop, complement, cumCovComponent)));
                        }
                    }
                    else
                    {
                        for (auto deploy = deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy)
                        {
                            timed.push_back(unique_ptr<TimedDeployment>(
                                new TimedHumanDeployment(deploy->first, *deploy->second, intervention, subPop, complement)));
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
        // code disabled; search for R_0 to find more blocks like this:
#if 0
        const scnXml::InsertR_0Case& elt = intervElt.getInsertR_0Case().get();
        if( elt.getTimedDeployment().size() > 0 ){
            Vaccine::verifyEnabledForR_0();
            // timed deployments:
            for( auto it = elt.getTimedDeployment().begin(); it != elt.getTimedDeployment().end(); ++it ){
                SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                timed.push_back( new TimedR_0Deployment( date ) );
            }
        }
#endif
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
                timed.push_back(unique_ptr<TimedDeployment>(new TimedUninfectVectorsDeployment(date)));
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
                    initVectorInterv(*vectorModel, elt.getDescription().getAnopheles(), instance, elt.getName());

                    const scnXml::TimedBaseList::DeploySequence &seq = elt.getTimed().get().getDeploy();
                    for (auto it = seq.begin(); it != seq.end(); ++it)
                    {
                        SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        timed.push_back(unique_ptr<TimedDeployment>(new TimedVectorDeployment(date, instance)));
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
            size_t instance = 0;
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
                        timed.push_back(unique_ptr<TimedDeployment>(new TimedAddNonHumanHostsDeployment(date, elt.getName(), lifespan, elt.getDescription().getAnopheles(), transmission)));
                    }
                    instance++;
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
                        timed.push_back(
                            unique_ptr<TimedDeployment>(new TimedNonHumanHostsDeployment(date, instance, elt.getNonHumanHostsName(), elt.getDescription().getAnopheles(), decay, transmission)));
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
                initVectorTrap(*vectorModel, trap.getDescription(), instance, trap.getName());
                if (trap.getTimed().present())
                {
                    for (const scnXml::Deploy1 &deploy : trap.getTimed().get().getDeploy())
                    {
                        SimTime date = UnitParse::readDate(deploy.getTime(), UnitParse::STEPS);
                        double ratio = deploy.getRatioToHumans();
                        SimTime lifespan = UnitParse::readDuration(deploy.getLifespan(), UnitParse::NONE);
                        timed.push_back(unique_ptr<TimedDeployment>(new TimedTrapDeployment(date, instance, ratio, lifespan)));
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
    // parity: b — stable_sort order = reproducibility contract; comparator ties resolve by insertion order
    stable_sort(continuous.begin(), continuous.end(), byDeployTime);
    stable_sort(timed.begin(), timed.end(), byDeployTime);

    // make sure the list ends with something always in the future, so we don't
    // have to check nextTimed is within range:
    timed.push_back(unique_ptr<TimedDeployment>(new DummyTimedDeployment()));

    if (util::CommandLine::option(util::CommandLine::PRINT_INTERVENTIONS))
    {
        cout << "Continuous deployments:" << endl << "begin\tend\tage\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for (auto it = continuous.begin(); it != continuous.end(); ++it)
        {
            it->print_details(std::cout);
            cout << endl;
        }
        cout << "Timed deployments:" << endl << "time\tmin age\tmax age\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for (auto &deploy : timed)
        {
            deploy->print_details(std::cout);
            cout << endl;
        }
        cout << "Human components:" << endl;
        for (auto &component : humanComponents)
        {
            component->print_details(cout);
            cout << endl;
        }
    }
}

void InterventionManager::loadFromCheckpoint(vector<Host::Human> &population, Transmission::TransmissionModel &transmission)
{
    SimTime date = sim::intervDate();
    // We need to re-deploy changeHS and changeEIR interventions, but nothing
    // else. nextTimed should be zero so we can go through all past interventions.
    // Only redeploy those which happened before this time step.
    assert(nextTimed == 0);
    while (timed[nextTimed]->date < date)
    {
        TimedDeployment *deployment = &*timed[nextTimed];
        if (dynamic_cast<TimedChangeHSDeployment *>(deployment) != 0 || dynamic_cast<TimedChangeEIRDeployment *>(deployment) != 0)
        {
            // Note: neither changeHS nor changeEIR interventions care what the
            // current time step is when they are deployed, so we don't need to
            // tell them the deployment date.
            deployment->deploy(population, transmission);
        }
        nextTimed += 1;
    }
}

} }
