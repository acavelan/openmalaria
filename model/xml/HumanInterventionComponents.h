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
#include "core/RunContext.h"
#include "interventions/Components.h"
#include "interventions/InterventionManager.h"
#include "interventions/IvState.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "Clinical/ESCaseManagement.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Host/Human.h"
#include "Transmission/TransmissionModel.h"
#include "Host/WithinHost/Diagnostic.h"
#include "PkPd/LSTMTreatments.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include <schema/healthSystem.h>
#include <schema/interventions.h>

namespace OM { namespace interventions {
using Host::Human;
using WithinHost::Diagnostic;
using WithinHost::diagnostics;

// ———  HumanIntervention  ———

// Components sort into standard type-rank order (contract b).
static void sort_components( RunContext& ctx, std::vector<ComponentId>& components ){
    auto kind = [&ctx]( ComponentId id ){ return ctx.ext->humanComponents[id.id].kind; };
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
    std::stable_sort( components.begin(), components.end(),
        [&kind]( ComponentId a, ComponentId b ){ return kind(a) < kind(b); } );
}
static ComponentId component_id( RunContext& ctx, const std::string& textId ){
    ComponentId id = InterventionManager::getComponentId( ctx, textId );
    if( id.id >= ctx.ext->humanComponents.size() )
        throw util::base_exception( "invalid component id" );  // used before defined
    return id;
}
HumanIntervention::HumanIntervention(
    RunContext& ctx,
    const xsd::cxx::tree::sequence<scnXml::Component>& componentList )
{
    components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it )
        components.push_back( component_id( ctx, it->getId() ) );
    sort_components( ctx, components );
}
HumanIntervention::HumanIntervention(
    RunContext& ctx,
    const xsd::cxx::tree::sequence<scnXml::Component>& componentList,
    const xsd::cxx::tree::sequence<scnXml::Condition>& conditionList )
{
    components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it )
        components.push_back( component_id( ctx, it->getId() ) );
    sort_components( ctx, components );
    
    for( auto it = conditionList.begin(), end = conditionList.end(); it != end; ++it ) {
        double min = it->getMinValue().present() ? it->getMinValue().get() : -numeric_limits<double>::infinity();
        double max = it->getMaxValue().present() ? it->getMaxValue().get() : numeric_limits<double>::infinity();
        conditions.push_back(mon::setupCondition( ctx, it->getMeasure(), min, max, it->getInitialState() ));
    }
}
HumanIntervention::HumanIntervention(
    RunContext& ctx, const xsd::cxx::tree::sequence<scnXml::DTDeploy>& componentList ){
    components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it )
        components.push_back( component_id( ctx, it->getComponent() ) );
    sort_components( ctx, components );
}

void HumanIntervention::deploy( RunContext& ctx, Human& human, mon::Deploy::Method method,
    VaccineLimits vaccLimits ) const
{
    for( size_t cond : conditions ){
        // Abort deployment if any condition is false.
        if( !mon::checkCondition(ctx, cond) ) return;
    }
    
    for( ComponentId id : components ) {
        const ComponentConfig& component = ctx.ext->humanComponents[id.id];
        // we must report first, since it can change cohort and sub-population
        // which may affect what deployment does (at least in the case of reporting deployments)
        human.addToCohort(id, component.duration);
        deploy_component( ctx, component, human, method, vaccLimits );
    }
}

void HumanIntervention::print_details( std::ostream& out )const{
    out << "human:";
    for( ComponentId id : components ){
        out << '\t' << id.id;
    }
}

// ———  Utilities  ———

TriggeredDeployments::TriggeredDeployments(RunContext& ctx, const scnXml::TriggeredDeployments& elt){
    lists.reserve( elt.getDeploy().size() );
    for( auto it = elt.getDeploy().begin(), end = elt.getDeploy().end(); it != end; ++it ) {
        lists.push_back( SubList( ctx, *it ) );
    }
}
void TriggeredDeployments::deploy( RunContext& ctx, Human& human,
        mon::Deploy::Method method, VaccineLimits vaccLimits) const
{
    for( auto it = lists.begin(); it != lists.end(); ++it ){
        it->deploy( ctx, human, method, vaccLimits );
    }
}

TriggeredDeployments::SubList::SubList( RunContext& ctx, const scnXml::Deploy& elt ) :
        HumanIntervention( ctx, elt.getComponent() ),
        minAge( ctx.clock.from_years_n( elt.getMinAge() ) ),
        maxAge( sim::future() ),
        coverage( elt.getP() )
{
    if( elt.getMaxAge().present() )
        maxAge = ctx.clock.from_years_n( elt.getMaxAge().get() );
    
    if( minAge < sim::zero() || maxAge < minAge ){
        throw util::xml_scenario_error("triggered intervention must have 0 <= minAge <= maxAge");
    }
    
    if( coverage < 0.0 || coverage > 1.0 ){
        throw util::xml_scenario_error( "triggered intervention must have 0 <= coverage <= 1" );
    }
    
    // Zero coverage: optimise
    if( coverage <= 0.0 || minAge >= maxAge ) components.clear();
}
void TriggeredDeployments::SubList::deploy( RunContext& ctx, Host::Human& human,
        mon::Deploy::Method method, VaccineLimits vaccLimits )const
{
    // This may be used within a time step; in that case we use human age at
    // the beginning of the step since this gives the expected result with age
    // limits of 0 to 1 time step.
    //TODO: why does the above contradict what we do? Is this due to a date being shifted later?
    SimTime age = human.age(ctx.clock.nowOrTs1());
    if( age >= minAge && age < maxAge ){
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if( coverage >= 1.0 || human.rng().bernoulli( coverage ) ){
            HumanIntervention::deploy( ctx, human, method, vaccLimits );
        }
    }
}


// ———  Component payload constructors  ———

ScreenParams::ScreenParams( RunContext& ctx, const scnXml::Screen& elt ) :
    diagnostic(diagnostics::get(ctx, elt.getDiagnostic())),
    positive( ctx, elt.getPositive() ),
    negative( ctx, elt.getNegative() )
{}

/// Simple treatment: no PK/PD, just remove parasites
TreatSimpleParams::TreatSimpleParams( RunContext& ctx, const scnXml::DTTreatSimple& elt )
{
    //NOTE: this code is currently identical to that in CMDTTreatSimple
    try{
        SimTime durL = UnitParse::readShortDuration(ctx.clock,  elt.getDurationLiver(), UnitParse::NONE ),
            durB = UnitParse::readShortDuration(ctx.clock,  elt.getDurationBlood(), UnitParse::NONE );
        SimTime neg1 = -ctx.clock.one_TS();
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

//NOTE: this is similar to CMDTTreatPKPD, but supporting only a single "treatment"
TreatPKPDParams::TreatPKPDParams( RunContext& ctx, const scnXml::DTTreatPKPD& elt ) :
    schedule(ctx.treatments.findSchedule(elt.getSchedule())),
    dosage(ctx.treatments.findDosages(elt.getDosage())),
    delay_d(elt.getDelay_h() / 24.0 /* convert to days */)
{
}

DecisionTreeParams::DecisionTreeParams( RunContext& ctx, const scnXml::DecisionTree& elt ) :
    tree(Clinical::compile_decision_tree(ctx, elt, false))
{}

// ———  Component deployment / printing (replaces the class tree's virtuals;
//      reports and draws happen in exactly the order of the deleted bodies)  ———

void deploy_component( RunContext& ctx, const ComponentConfig& c, Human& human,
    mon::Deploy::Method method, VaccineLimits vaccLimits )
{
    switch( c.kind ){
    case Component::RECRUIT_ONLY:  // reports to monitoring, nothing else
        mon::recordDeploy(mon::measure("nMassRecruitOnly"), mon::measure("nCtsRecruitOnly"), human, method);
        break;
    case Component::SCREEN: { const ScreenParams& p = std::get<ScreenParams>(c.data);
        mon::recordDeploy(mon::measure("nMassScreenings"), mon::measure("nCtsScreenings"), human, method);
        if( human.withinHostModel->diagnosticResult(human.rng(), p.diagnostic) )
            p.positive.deploy( ctx, human, method, vaccLimits );
        else
            p.negative.deploy( ctx, human, method, vaccLimits );
        break; }
    case Component::TREAT_SIMPLE: { const TreatSimpleParams& p = std::get<TreatSimpleParams>(c.data);
        mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        human.withinHostModel->treatSimple( human, p.timeLiver, p.timeBlood );
        break; }
    case Component::TREAT_PKPD: { const TreatPKPDParams& p = std::get<TreatPKPDParams>(c.data);
        mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        double age = sim::inYears(human.age(ctx.clock.nowOrTs1()));
        human.withinHostModel->treatPkPd( p.schedule, p.dosage, age, p.delay_d );
        break; }
    case Component::CM_DT: { const DecisionTreeParams& p = std::get<DecisionTreeParams>(c.data);
        Clinical::DTExecCtx ec( ctx, human,
                sim::inYears(human.age(ctx.clock.nowOrTs1())),
                human.clinicalModel->getLatestState() );
        Clinical::CMDTOut out = Clinical::exec_decision_tree(p.tree, ec);
        if( out.treated ){
            mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        }
        if( out.screened ){
            mon::recordDeploy(mon::measure("nMassScreenings"), mon::measure("nCtsScreenings"), human, method);
        }
        break; }
    case Component::PEV: case Component::BSV: case Component::TBV:
        if( vaccinate( ctx, human, c.id, vaccLimits ) && ctx.ext->vaccineReportComponent == c.id ){
            mon::recordDeploy(mon::measure("nMassVaccinations"), mon::measure("nEPIVaccinations"), human, method);
        }
        if (method == mon::Deploy::TREAT)
            mon::recordEvent(mon::measure("nTreatDeployments"), human);
        break;
    case Component::ITN:
        iv_deploy( ctx, human, c );
        mon::recordDeploy(mon::measure("nMassITNs"), mon::measure("nEPI_ITNs"), human, method);
        break;
    case Component::IRS:
        iv_deploy( ctx, human, c );
        mon::recordDeploy(mon::measure("nMassIRS"), mon::measure("nCtsIRS"), human, method);
        break;
    case Component::GVI:
        iv_deploy( ctx, human, c );
        mon::recordDeploy(mon::measure("nMassGVI"), mon::measure("nCtsGVI"), human, method);
        break;
    case Component::CLEAR_IMMUNITY:
        human.withinHostModel->clearImmunity();
        //NOTE: not reported
        break;
    default:  // Component::IPT: no implementation exists
        throw SWITCH_DEFAULT_EXCEPTION;
    }
}

void print_component( std::ostream& out, const ComponentConfig& c ){
    out << c.id.id << '\t' << (
        c.kind == Component::RECRUIT_ONLY ? "Recruit only" :
        c.kind == Component::SCREEN ? "Screen" :
        c.kind == Component::TREAT_SIMPLE ? "treatSimple" :
        c.kind == Component::TREAT_PKPD ? "treatPKPD" :
        c.kind == Component::CM_DT ? "Decision tree" :
        c.kind == Component::PEV ? "PEV" :
        c.kind == Component::BSV ? "BSV" :
        c.kind == Component::TBV ? "TBV" :
        c.kind == Component::ITN ? "ITN" :
        c.kind == Component::IRS ? "IRS" :
        c.kind == Component::GVI ? "GVI" :
        /* Component::CLEAR_IMMUNITY */ "clear immunity");
}

} }
