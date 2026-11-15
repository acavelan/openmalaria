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

#include "schema/interventions.h"
#include "Clinical/CMDecisionTree.h"
#include "Clinical/ClinicalModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Diagnostic.h"
#include "PkPd/LSTMTreatments.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/UnitParse.h"
#include "interventions/Interfaces.h"
#include "interventions/InterventionManager.h"

#include <limits>
#include <list>
#include <vector>
#include <map>

namespace OM { namespace Clinical {

using WithinHost::Diagnostic;
using WithinHost::diagnostics;
using namespace OM::util;
using std::map;
using std::vector;
using Kind = CMDecisionTree::Kind;

// ———  execution  ———

CMDTOut CMDecisionTree::exec( CMHostData hostData ) const {
    switch (kind) {
    case Kind::Multiple: {
        // Branch out to multiple descendants.
        CMDTOut result(false);
        for( auto it = children.begin(), end = children.end(); it != end; ++it ) {
            CMDTOut r2 = (*it)->exec( hostData );
            result.treated = result.treated || r2.treated;
        }
        return result;
    }
    case Kind::CaseType: {
        // Should the patient recieve first or second line treatment?
        // Uses of this in complicated cases should trigger an exception during initialisation.
        assert( (hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED) );
        
        if( hostData.pgState & Episode::SECOND_CASE ) return secondLine->exec( hostData );
        else return firstLine->exec( hostData );
    }
    case Kind::InfectionOrigin: {
        // Should the patient recieve a different treatment based on its infection type?
        WithinHost::InfectionOrigin InfectionOrigin = hostData.withinHost().getInfectionOrigin();
        if(InfectionOrigin == WithinHost::InfectionOrigin::Imported)
            return imported->exec( hostData );
        else if (InfectionOrigin == WithinHost::InfectionOrigin::Introduced)
            return introduced->exec( hostData );
        else
            return indigenous->exec( hostData );
    }
    case Kind::Diagnostic: {
        // Use a diagnostic, and call another decision tree based on the outcome.
        CMDTOut result;
        
        if( hostData.withinHost().diagnosticResult( hostData.human.rng, *diagnostic ) )
            result = positive->exec( hostData );
        else
            result = negative->exec( hostData );
        result.screened = true;

        return result;
    }
    case Kind::Uncomplicated: {
        const Clinical::Episode &latest = hostData.human.clinicalModel->getLatestReport();

        if ( ((hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED)) || (hostData.pgState & Episode::MALARIA) )
        {
            if (latest.time + memory >= sim::nowOrTs0())
                return positive->exec( hostData );
            else
                return negative->exec( hostData );
        }
        else
            return negative->exec( hostData );
    }
    case Kind::Severe: {
        CMDTOut result;

        if (hostData.pgState & Episode::COMPLICATED)
            result = positive->exec( hostData );
        else
            result = negative->exec( hostData );

        // result.screened = true;
        return result;
    }
    case Kind::Random: {
        // Choose a branch randomly.
        auto it = branches.upper_bound( hostData.human.rng.uniform_01() );
        assert( it != branches.end() );
        return it->second->exec( hostData );
    }
    case Kind::Age: {
        // Choose what to do based on the patient's age.
        // age is that of human at start of time step (i.e. may be as low as 0)
        auto it = branches.upper_bound( hostData.ageYears );
        if( it == branches.end() )
            throw TRACED_EXCEPTION( "bad age-based decision tree switch", util::Error::PkPd );
        return it->second->exec( hostData );
    }
    case Kind::Cohort: {
        CMDTOut result;

        // Rely on the health system memory to not count the same episode twice
        if(hostData.human.isInSubPop(*component))
            result = positive->exec( hostData );
        else
            result = negative->exec( hostData );

        // result.screened = true;
        return result;
    }
    case Kind::NoTreatment:
        return CMDTOut(false);
    case Kind::Report: {
        for( const int outId : outIds ){
            mon::recordEvent(mon::measure("nCMDTReport"), hostData.human, 1, outId);
        }
        return CMDTOut(false);
    }
    case Kind::TreatFailure:
        // Report treatment without affecting parasites.
        return CMDTOut(true /*report treatment*/);
    case Kind::TreatPKPD: {
        // Deliver one or more treatments via the PK/PD model.
        for( const TreatInfo& treatment : treatments ){
            hostData.withinHost().treatPkPd( treatment.schedule, treatment.dosage, hostData.ageYears, 0.0 );
        }
        return CMDTOut(true);
    }
    case Kind::TreatSimple: {
        // Deliver one or more simple treatments.
        bool bsTreatment = false;
        for(size_t i=0; i<timeLiver.size(); i++)
        {
            bsTreatment = hostData.withinHost().treatSimple( hostData.human, timeLiver[i], timeBlood[i] );
        }
        return CMDTOut(bsTreatment);
    }
    case Kind::Deploy: {
        // Deploy one or more interventions.
        deployment->deploy( hostData.human,
                  mon::Deploy::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/),
                  sim::nowOrTs1() );
        
        //NOTE: it's not intuitively obvious what value should be returned here
        // in the case of intervention deployment. This at least means that
        // repeat seekers get second-line treatment.
        return CMDTOut(false);
    }
    }
    return CMDTOut(false);   // unreachable: every Kind is handled above
}

// ———  construction  ———

// Memory management: lists all decisions and frees memory at program exit.
// We store pointers into this list, so elements must not move.
vector<unique_ptr<CMDecisionTree>> decision_library;

// Takes ownership of a new decision and returns a const reference to it.
// Trees are also created at runtime by changeHS interventions; old ones are
// kept forever (a few per event, never freed — negligible).
const CMDecisionTree& save_decision( CMDecisionTree* decision ){
    decision_library.push_back( unique_ptr<CMDecisionTree>(decision) );
    return *decision_library.back();
}

CMDecisionTree::TreatInfo::TreatInfo( const string& s, const string& d, double h ) :
    schedule(PkPd::LSTMTreatments::findSchedule(s)),
    dosage(PkPd::LSTMTreatments::findDosages(d)),
    delay_h(h) {}

namespace {

using Reg = const interventions::InterventionRegistry&;

unique_ptr<CMDecisionTree> make( Kind kind ){
    return std::make_unique<CMDecisionTree>(kind);
}

const CMDecisionTree& createCaseType( const scnXml::DTCaseType& node, bool isUC, Reg reg);
const CMDecisionTree& createInfectionOrigin( const scnXml::DTInfectionOrigin& node, bool isUC, Reg reg);
const CMDecisionTree& createDiagnostic( const scnXml::DTDiagnostic& node, bool isUC, Reg reg);
const CMDecisionTree& createUncomplicated( const scnXml::DTUncomplicated& node, bool isUC, Reg reg);
const CMDecisionTree& createSevere( const scnXml::DTSevere& node, bool isUC, Reg reg);
const CMDecisionTree& createRandom( const scnXml::DTRandom& node, bool isUC, Reg reg);
const CMDecisionTree& createAge( const scnXml::DTAge& node, bool isUC, Reg reg);
const CMDecisionTree& createCohort( const scnXml::DTCohort& node, bool isUC, Reg reg);
const CMDecisionTree& createReport( const scnXml::DecisionTree::ReportSequence& seq, bool isUC, Reg reg);
const CMDecisionTree& createTreatPKPD( const scnXml::DecisionTree::TreatPKPDSequence& seq, bool isUC, Reg reg);
const CMDecisionTree& createTreatSimple( const scnXml::DecisionTree::TreatSimpleSequence& seq, bool isUC, Reg reg);
const CMDecisionTree& createDeploy( const scnXml::DecisionTree::DeploySequence& seq, bool isUC, Reg reg);

const CMDecisionTree& createMultiple( const scnXml::DTMultiple& node, bool isUC, Reg reg){
    unique_ptr<CMDecisionTree> self = make(Kind::Multiple);
    for( const scnXml::DTCaseType& sn : node.getCaseType() ){
        self->children.push_back( &createCaseType(sn, isUC, reg) );
    }
    for( const scnXml::DTInfectionOrigin& sn : node.getInfectionOrigin() ){
        self->children.push_back( &createInfectionOrigin(sn, isUC, reg) );
    }
    for( const scnXml::DTDiagnostic& sn : node.getDiagnostic() ){
        self->children.push_back( &createDiagnostic(sn, isUC, reg) );
    }
    for( const scnXml::DTUncomplicated& sn : node.getUncomplicated() ){
        self->children.push_back( &createUncomplicated(sn, isUC, reg) );
    }
    for( const scnXml::DTSevere& sn : node.getSevere() ){
        self->children.push_back( &createSevere(sn, isUC, reg) );
    }
    for( const scnXml::DTRandom& sn : node.getRandom() ){
        self->children.push_back( &createRandom(sn, isUC, reg) );
    }
    for( const scnXml::DTAge& sn : node.getAge() ){
        self->children.push_back( &createAge(sn, isUC, reg) );
    }
    for( const scnXml::DTCohort& sn : node.getCohort() ){
        self->children.push_back( &createCohort(sn, isUC, reg) );
    }
    if( node.getTreatPKPD().size() ){
        self->children.push_back( &createTreatPKPD(node.getTreatPKPD(), isUC, reg) );
    }
    if( node.getTreatSimple().size() ){
        self->children.push_back( &createTreatSimple(node.getTreatSimple(), isUC, reg) );
    }
    if( node.getReport().size() ){
        self->children.push_back( &createReport(node.getReport(), isUC, reg) );
    }
    if( node.getDeploy().size() ){
        self->children.push_back( &createDeploy(node.getDeploy(), isUC, reg) );
    }
    return save_decision( self.release() );
}

const CMDecisionTree& createCaseType( const scnXml::DTCaseType& node, bool isUC, Reg reg){
    if( !isUC ){
        throw util::xml_scenario_error( "decision tree: caseType can only be used for uncomplicated cases" );
    }
    unique_ptr<CMDecisionTree> self = make(Kind::CaseType);
    self->firstLine = &CMDecisionTree::create( node.getFirstLine(), isUC, reg);
    self->secondLine = &CMDecisionTree::create( node.getSecondLine(), isUC, reg);
    return save_decision( self.release() );
}

const CMDecisionTree& createInfectionOrigin( const scnXml::DTInfectionOrigin& node, bool isUC, Reg reg){
    unique_ptr<CMDecisionTree> self = make(Kind::InfectionOrigin);
    self->imported = &CMDecisionTree::create( node.getImported(), isUC, reg);
    self->introduced = &CMDecisionTree::create( node.getIntroduced(), isUC, reg);
    self->indigenous = &CMDecisionTree::create( node.getIndigenous(), isUC, reg);
    return save_decision( self.release() );
}

const CMDecisionTree& createDiagnostic( const scnXml::DTDiagnostic& node, bool isUC, Reg reg){
    unique_ptr<CMDecisionTree> self = make(Kind::Diagnostic);
    self->diagnostic = &diagnostics::get( node.getDiagnostic() );
    self->positive = &CMDecisionTree::create( node.getPositive(), isUC, reg);
    self->negative = &CMDecisionTree::create( node.getNegative(), isUC, reg);
    return save_decision( self.release() );
}

const CMDecisionTree& createUncomplicated( const scnXml::DTUncomplicated& node, bool isUC, Reg reg){
    unique_ptr<CMDecisionTree> self = make(Kind::Uncomplicated);
    self->positive = &CMDecisionTree::create( node.getPositive(), isUC, reg);
    self->negative = &CMDecisionTree::create( node.getNegative(), isUC, reg);
    self->memory = UnitParse::readShortDuration(node.getMemory(), UnitParse::STEPS);
    if(self->memory > ClinicalModel::hsMemory())
        throw util::xml_scenario_error( "<uncomplicated> memory parameter must be less than or equal to the healthsystem memory (hsmemory parameter)" );
    return save_decision( self.release() );
}

const CMDecisionTree& createSevere( const scnXml::DTSevere& node, bool isUC, Reg reg){
    unique_ptr<CMDecisionTree> self = make(Kind::Severe);
    self->positive = &CMDecisionTree::create( node.getPositive(), isUC, reg);
    self->negative = &CMDecisionTree::create( node.getNegative(), isUC, reg);
    return save_decision( self.release() );
}

const CMDecisionTree& createRandom( const scnXml::DTRandom& node, bool isUC, Reg reg){
    unique_ptr<CMDecisionTree> result = make(Kind::Random);
    
    double cum_p = 0.0;
    for( const scnXml::Outcome& outcome : node.getOutcome() ){
        cum_p += outcome.getP();
        result->branches.insert( make_pair(cum_p, &CMDecisionTree::create( outcome, isUC, reg)) );
    }
    
    // Test cum_p is approx. 1.0 in case the input tree is wrong. We require no
    // less than one to make sure generated random numbers are not greater than
    // the last option.
    if (cum_p < 1.0 || cum_p > 1.001){
        throw util::xml_scenario_error (
            "decision tree (random node): expected probability sum to be 1.0 but found " + to_string(cum_p));
    }
    
    return save_decision( result.release() );
}

const CMDecisionTree& createAge( const scnXml::DTAge& node, bool isUC, Reg reg){
    unique_ptr<CMDecisionTree> result = make(Kind::Age);
    
    double lastAge = numeric_limits<double>::quiet_NaN();
    const CMDecisionTree* lastNode = nullptr;
    for( const scnXml::Age& age : node.getAge() ){
        if( lastAge != lastAge ){
            if( age.getLb() != 0.0 )
                throw util::xml_scenario_error( "decision tree age switch must have first lower bound equal 0" );
        }else{
            if( age.getLb() <= lastAge ){
                throw util::xml_scenario_error( "decision tree age switch must list age groups in increasing order" );
            }
            double lb = age.getLb();
            result->branches.insert( make_pair(lb, lastNode) );
        }
        lastNode = &CMDecisionTree::create(age, isUC, reg);
        lastAge = age.getLb();
    }
    double noLb = numeric_limits<double>::infinity();
    result->branches.insert( make_pair(noLb, lastNode) );
    
    return save_decision( result.release() );
}

const CMDecisionTree& createCohort( const scnXml::DTCohort& node, bool isUC, Reg reg){
    unique_ptr<CMDecisionTree> self = make(Kind::Cohort);
    self->positive = &CMDecisionTree::create( node.getPositive(), isUC, reg);
    self->negative = &CMDecisionTree::create( node.getNegative(), isUC, reg);
    self->component = interventions::getComponentId(reg, node.getComponent());
    return save_decision( self.release() );
}

const CMDecisionTree& createReport( const scnXml::DecisionTree::ReportSequence& seq, bool /*isUC*/, Reg /*reg*/){
    unique_ptr<CMDecisionTree> self = make(Kind::Report);
    self->outIds.reserve( seq.size() );
    for( const scnXml::DTReport& reportElt : seq ){
        self->outIds.push_back( reportElt.getOutputNumber() );
    }
    assert( self->outIds.size() > 0 );        // a Report node should not be empty
    return save_decision( self.release() );
}

const CMDecisionTree& createTreatPKPD( const scnXml::DecisionTree::TreatPKPDSequence& seq, bool /*isUC*/, Reg /*reg*/){
    unique_ptr<CMDecisionTree> self = make(Kind::TreatPKPD);
    self->treatments.reserve( seq.size() );
    for( const scnXml::DTTreatPKPD& treatElt : seq ){
        self->treatments.push_back( CMDecisionTree::TreatInfo(
            treatElt.getSchedule(),
            treatElt.getDosage(),
            treatElt.getDelay_h()
        ) );
    }
    assert( self->treatments.size() > 0 );        // a TreatPKPD node should not be empty
    return save_decision( self.release() );
}

const CMDecisionTree& createTreatSimple( const scnXml::DecisionTree::TreatSimpleSequence& seq, bool /*isUC*/, Reg /*reg*/){
    unique_ptr<CMDecisionTree> self = make(Kind::TreatSimple);
    for (const auto &treatElt : seq) {
        //NOTE: this code is currently identical to that in SimpleTreatComponent
        try{
            SimTime durL = UnitParse::readShortDuration( treatElt.getDurationLiver(), UnitParse::NONE ),
                durB = UnitParse::readShortDuration( treatElt.getDurationBlood(), UnitParse::NONE );
            SimTime neg1 = -sim::oneTS();
            if( durL < neg1 || durB < neg1 ){
                throw util::xml_scenario_error( "treatSimple: cannot have durationBlood or durationLiver less than -1" );
            }
            self->timeLiver.push_back(durL);
            self->timeBlood.push_back(durB);
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("treatSimple: ").append(e.message()) );
        }
    }
    return save_decision( self.release() );
}

/*NOTE: we use interventions::HumanIntervention to read and sort the list of
 * interventions to be deployed, though we don't use pointers it finds. */
const CMDecisionTree& createDeploy( const scnXml::DecisionTree::DeploySequence& seq, bool /*isUC*/, Reg reg){
    unique_ptr<CMDecisionTree> self = make(Kind::Deploy);
    self->deployment.emplace( seq, reg );
    return save_decision( self.release() );
}

} // anonymous namespace

const CMDecisionTree& CMDecisionTree::create( const scnXml::DecisionTree& node, bool isUC, const interventions::InterventionRegistry& reg){
    if( node.getMultiple().present() ) return createMultiple( node.getMultiple().get(), isUC, reg);
    // branching nodes
    if( node.getCaseType().present() ) return createCaseType( node.getCaseType().get(), isUC, reg);
    if( node.getInfectionOrigin().present() ) return createInfectionOrigin( node.getInfectionOrigin().get(), isUC, reg);
    if( node.getDiagnostic().present() ) return createDiagnostic( node.getDiagnostic().get(), isUC, reg);
    if( node.getUncomplicated().present() ) return createUncomplicated( node.getUncomplicated().get(), isUC, reg);
    if( node.getSevere().present() ) return createSevere( node.getSevere().get(), true, reg);
    if( node.getRandom().present() ) return createRandom( node.getRandom().get(), isUC, reg);
    if( node.getAge().present() ) return createAge( node.getAge().get(), isUC, reg);
    if( node.getCohort().present() ) return createCohort( node.getCohort().get(), isUC, reg);

    // action nodes
    if( node.getNoTreatment().present() ) return save_decision( make(Kind::NoTreatment).release() );
    if( node.getReport().size() ) return createReport( node.getReport(), isUC, reg );
    if( node.getTreatFailure().present() ) return save_decision( make(Kind::TreatFailure).release() );
    if( node.getTreatPKPD().size() ) return createTreatPKPD( node.getTreatPKPD(), isUC, reg );
    if( node.getTreatSimple().size() ) return createTreatSimple( node.getTreatSimple(), isUC, reg );
    if( node.getDeploy().size() ) return createDeploy( node.getDeploy(), isUC, reg );
    throw xml_scenario_error( "unterminated decision tree" );
}

} }
