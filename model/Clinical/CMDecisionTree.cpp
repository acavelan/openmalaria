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

#include "Clinical/CMDecisionTree.h"
#include "schema/healthSystem.h"
#include "Clinical/ClinicalModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Diagnostic.h"
#include "PkPd/LSTMTreatments.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/UnitParse.h"
#include "interventions/Interfaces.h"
#include "interventions/Interventions.h"

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


// ———  special 'multiple' node  ———

/**
 * Branch out to multiple descendants.
 */
class CMDTMultiple : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTMultiple& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library);
    
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTMultiple*>(&that);
        return p && children == p->children;
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        CMDTOut result(false);
        for( auto it = children.begin(), end = children.end(); it != end; ++it ) {
            CMDTOut r2 = (*it)->exec( hostData );
            result.treated = result.treated || r2.treated;
        }
        return result;
    }
    
private:
    CMDTMultiple( /*size_t capacity*/ ){
//         children.reserve( capacity );
    }
    
    typedef vector<const CMDecisionTree*> Children_t;
    Children_t children;
};

// ———  branching nodes  ———

// Branch construction keeps the original argument evaluation order for node interning.
static const CMDecisionTree& createBranch(const scnXml::DTCaseType&, bool isUC, const interventions::Interventions&, ClinicalParameters&);
static const CMDecisionTree& createBranch(const scnXml::DTInfectionOrigin&, bool isUC, const interventions::Interventions&, ClinicalParameters&);
static const CMDecisionTree& createBranch(const scnXml::DTDiagnostic&, bool isUC, const interventions::Interventions&, ClinicalParameters&);
static const CMDecisionTree& createBranch(const scnXml::DTUncomplicated&, bool isUC, const interventions::Interventions&, ClinicalParameters&);
static const CMDecisionTree& createBranch(const scnXml::DTSevere&, bool isUC, const interventions::Interventions&, ClinicalParameters&);
static const CMDecisionTree& createBranch(const scnXml::DTCohort&, bool isUC, const interventions::Interventions&, ClinicalParameters&);

class CMDTBranch : public CMDecisionTree {
public:
    enum class Kind { CaseType, InfectionOrigin, Diagnostic, Uncomplicated, Severe, Cohort };
    CMDTBranch(Kind kind, const CMDecisionTree& first, const CMDecisionTree& second) :
        kind(kind), first(first), second(second) {}
    CMDTBranch(Kind kind, const CMDecisionTree& first, const CMDecisionTree& second,
            const CMDecisionTree& third) : CMDTBranch(kind, first, second) { this->third = &third; }
    CMDTBranch(const Diagnostic& diagnostic, const CMDecisionTree& first, const CMDecisionTree& second) :
        CMDTBranch(Kind::Diagnostic, first, second) { this->diagnostic = &diagnostic; }
    CMDTBranch(SimTime memory, const CMDecisionTree& first, const CMDecisionTree& second, SimTime hsMemory) :
        CMDTBranch(Kind::Uncomplicated, first, second) {
        this->memory = memory;
        if (memory > hsMemory)
            throw util::xml_scenario_error("<uncomplicated> memory parameter must be less than or equal to the healthsystem memory (hsmemory parameter)");
    }
    CMDTBranch(interventions::ComponentId component, const CMDecisionTree& first, const CMDecisionTree& second) :
        CMDTBranch(Kind::Cohort, first, second) { this->component = component; }

    Kind kind;
    const CMDecisionTree& first;
    const CMDecisionTree& second;
    const CMDecisionTree* third = nullptr;
    const Diagnostic* diagnostic = nullptr;
    SimTime memory = 0;
    interventions::ComponentId component = interventions::ComponentId::wholePop();

    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTBranch*>(&that);
        if (!p || kind != p->kind) return false;
        if (kind == Kind::Diagnostic && *diagnostic != *p->diagnostic) return false;
        // Legacy interning ignores uncomplicated memory and cohort IDs.
        return first == p->first && second == p->second &&
            (kind != Kind::InfectionOrigin || *third == *p->third);
    }

    CMDTOut exec(CMHostData hostData) const override {
        switch (kind) {
        case Kind::CaseType:
            assert((hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED));
            return (hostData.pgState & Episode::SECOND_CASE ? second : first).exec(hostData);
        case Kind::InfectionOrigin: {
            const auto origin = hostData.withinHost().getInfectionOrigin();
            if (origin == WithinHost::InfectionOrigin::Imported) return first.exec(hostData);
            if (origin == WithinHost::InfectionOrigin::Introduced) return second.exec(hostData);
            return third->exec(hostData);
        }
        case Kind::Diagnostic: {
            auto result = (hostData.withinHost().diagnosticResult(hostData.human.rng, *diagnostic)
                ? first : second).exec(hostData);
            result.screened = true;
            return result;
        }
        case Kind::Uncomplicated: {
            const auto& latest = hostData.human.clinicalModel->getLatestReport();
            if (((hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED)) ||
                    (hostData.pgState & Episode::MALARIA)) {
                if (latest.time + memory >= sim::nowOrTs0()) return first.exec(hostData);
            }
            return second.exec(hostData);
        }
        case Kind::Severe:
            return (hostData.pgState & Episode::COMPLICATED ? first : second).exec(hostData);
        case Kind::Cohort:
            return (hostData.human.isInSubPop(component) ? first : second).exec(hostData);
        }
        throw SWITCH_DEFAULT_EXCEPTION;
    }
};

/**
 * Choose a branch randomly.
 * 
 * TODO(optimisation): can we re-order the trees such that this can be abstracted out or
 * used only once?
 */
class CMDTRandom : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTRandom& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library);
    
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTRandom*>(&that);
        return p && branches == p->branches;
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        auto it = branches.upper_bound( hostData.human.rng.uniform_01() );
        assert( it != branches.end() );
        return it->second->exec( hostData );
    }
    
private:
    CMDTRandom(){}
    
    // keys are cumulative probabilities; last entry should equal 1
    typedef map<double,const CMDecisionTree*> Branches_t;
    Branches_t branches;
};

/**
 * Choose what to do based on the patient's age.
 */
class CMDTAge : public CMDecisionTree {
public:
    static const CMDecisionTree& create( const ::scnXml::DTAge& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library);
    
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTAge*>(&that);
        return p && branches == p->branches;
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        // age is that of human at start of time step (i.e. may be as low as 0)
        auto it = branches.upper_bound( hostData.ageYears );
        if( it == branches.end() )
            throw TRACED_EXCEPTION( "bad age-based decision tree switch", util::Error::PkPd );
        return it->second->exec( hostData );
    }
    
private:
    CMDTAge() {}
    
    // keys are upper bounds of age categories
    typedef map<double,const CMDecisionTree*> Branches_t;
    Branches_t branches;
};


// ———  action nodes  ———

/** Do nothing. **/
class CMDTNoTreatment : public CMDecisionTree {
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTNoTreatment*>(&that);
        return p != nullptr;
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        return CMDTOut(false);
    }
};

class CMDTReport : public CMDecisionTree {
public:
    CMDTReport( const scnXml::DecisionTree::ReportSequence& seq ){
        outIds.reserve( seq.size() );
        for( const scnXml::DTReport& reportElt : seq ){
            outIds.push_back( reportElt.getOutputNumber() );
        }
        assert( outIds.size() > 0 );        // CMDTTreatPKPD should not be used in this case
    }
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTReport*>(&that);
        return p && outIds == p->outIds;
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        for( const int outId : outIds ){
            mon::recordEvent(mon::measure("nCMDTReport"), hostData.human, 1, outId);
        }
        return CMDTOut(false);
    }
private:
    vector<int> outIds;
};

/** Report treament without affecting parasites. **/
class CMDTTreatFailure : public CMDecisionTree {
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTTreatFailure*>(&that);
        return p != nullptr;
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        return CMDTOut(true /*report treatment*/);
    }
};

/**
 * Deliver one or more treatments via the PK/PD model.
 */
class CMDTTreatPKPD : public CMDecisionTree {
public:
    CMDTTreatPKPD( const scnXml::DecisionTree::TreatPKPDSequence& seq ){
        treatments.reserve( seq.size() );
        for( const scnXml::DTTreatPKPD& treatElt : seq ){
            treatments.push_back( TreatInfo(
                treatElt.getSchedule(),
                treatElt.getDosage(),
                treatElt.getDelay_h()
            ) );
        }
        assert( treatments.size() > 0 );        // CMDTTreatPKPD should not be used in this case
    }
    
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTTreatPKPD*>(&that);
        return p && treatments == p->treatments;
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        for( const TreatInfo& treatment : treatments ){
            hostData.withinHost().treatPkPd( treatment.schedule, treatment.dosage, hostData.ageYears, 0.0 );
        }
        return CMDTOut(true);
    }
    
private:
    struct TreatInfo{
        TreatInfo( const string& s, const string& d, double h ) :
            schedule(PkPd::LSTMTreatments::findSchedule(s)),
            dosage(PkPd::LSTMTreatments::findDosages(d)),
            delay_h(h) {}
        bool operator==(const TreatInfo&) const = default;
        size_t schedule;        // index of the schedule
        size_t dosage;          // name of the dosage table
        double delay_h;         // delay in hours
    };
    vector<TreatInfo> treatments;
};

/**
 * Deliver one or more simple treatments
 */
class CMDTTreatSimple : public CMDecisionTree {
public:
    CMDTTreatSimple( const scnXml::DecisionTree::TreatSimpleSequence& elt ) :
        timeLiver(sim::zero()), timeBlood(sim::zero())
    {
        for (const auto &treatElt : elt) {
            //NOTE: this code is currently identical to that in SimpleTreatComponent
            try{
                SimTime durL = UnitParse::readShortDuration( treatElt.getDurationLiver(), UnitParse::NONE ),
                    durB = UnitParse::readShortDuration( treatElt.getDurationBlood(), UnitParse::NONE );
                SimTime neg1 = -sim::oneTS();
                if( durL < neg1 || durB < neg1 ){
                    throw util::xml_scenario_error( "treatSimple: cannot have durationBlood or durationLiver less than -1" );
                }
                timeLiver.push_back(durL);
                timeBlood.push_back(durB);
            }catch( const util::format_error& e ){
                throw util::xml_scenario_error( string("treatSimple: ").append(e.message()) );
            }
        }
    }
    
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTTreatSimple*>(&that);
        // Legacy interning compares blood-duration count, not values.
        return p && timeBlood.size() == p->timeBlood.size() && timeLiver == p->timeLiver;
    }

    
    virtual CMDTOut exec( CMHostData hostData ) const{
        bool bsTreatment = false;
        for(size_t i=0; i<timeLiver.size(); i++)
        {
            bsTreatment = hostData.withinHost().treatSimple( hostData.human, timeLiver[i], timeBlood[i] );
        }
        return CMDTOut(bsTreatment);
    }
    
private:
    vector<SimTime> timeLiver, timeBlood;
};

/**
 * Deploy one or more interventions.
 */
/*NOTE: we use interventions::HumanIntervention to read and sort the list of
 * interventions to be deployed, though we don't use pointers it finds. */
class CMDTDeploy : public CMDecisionTree, interventions::HumanIntervention {
public:
    CMDTDeploy( const scnXml::DecisionTree::DeploySequence& seq, const interventions::Interventions& definitions ) :
        HumanIntervention(seq, definitions) {}
    
protected:
    bool operator==(const CMDecisionTree& that) const override {
        if (this == &that) return true;
        const auto* p = dynamic_cast<const CMDTDeploy*>(&that);
        return p && components == p->components;
    }
    
    virtual CMDTOut exec( CMHostData hostData ) const{
        deploy( hostData.human,
                  mon::Deploy::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/) );
        
        //NOTE: it's not intuitively obvious what value should be returned here
        // in the case of intervention deployment. This at least means that
        // repeat seekers get second-line treatment.
        return CMDTOut(false);
    }
};

// ———  static functions  ———

// Memory management: lists all decisions and frees memory at program exit.
// We store pointers into this list, so elements must not move.
// Saves a decision to library.nodes, making it const.
// Also optimises away duplicates
const CMDecisionTree& save_decision(ClinicalParameters& library, CMDecisionTree* rawDecision) {
    unique_ptr<CMDecisionTree> decision(rawDecision);
    // We search the library for a duplicate, and delete this one if there is
    // a duplicate. Note that this is not implemented efficiently, but a little
    // wasted time at start up is hardly a concern.
    for( auto& d : library.nodes ){
        if( *d == *decision ){
            return *d;
        }
    }
    
    // No match: add to the library.
    library.nodes.push_back(std::move(decision));
    return *library.nodes.back();
}

const CMDecisionTree& CMDecisionTree::create( const scnXml::DecisionTree& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    if( node.getMultiple().present() ) return CMDTMultiple::create( node.getMultiple().get(), isUC, definitions, library);
    // branching nodes
    if( node.getCaseType().present() ) return createBranch( node.getCaseType().get(), isUC, definitions, library);
    if( node.getInfectionOrigin().present() ) return createBranch( node.getInfectionOrigin().get(), isUC, definitions, library);
    if( node.getDiagnostic().present() ) return createBranch( node.getDiagnostic().get(), isUC, definitions, library);
    if( node.getUncomplicated().present() ) return createBranch( node.getUncomplicated().get(), isUC, definitions, library);
    if( node.getSevere().present() ) return createBranch( node.getSevere().get(), true, definitions, library);
    if( node.getRandom().present() ) return CMDTRandom::create( node.getRandom().get(), isUC, definitions, library);
    if( node.getAge().present() ) return CMDTAge::create( node.getAge().get(), isUC, definitions, library);
    if( node.getCohort().present() ) return createBranch( node.getCohort().get(), isUC, definitions, library);

    // action nodes
    if( node.getNoTreatment().present() ) return save_decision(library, new CMDTNoTreatment() );
    if( node.getReport().size() ) return save_decision(library, new CMDTReport(node.getReport()) );
    if( node.getTreatFailure().present() ) return save_decision(library, new CMDTTreatFailure() );
    if( node.getTreatPKPD().size() ) return save_decision(library,
        new CMDTTreatPKPD( node.getTreatPKPD() ) );
    if( node.getTreatSimple().size() ) return save_decision(library,
        new CMDTTreatSimple( node.getTreatSimple() ) );
    if( node.getDeploy().size() ) return save_decision(library,
        new CMDTDeploy( node.getDeploy(), definitions ) );
    throw xml_scenario_error( "unterminated decision tree" );
}

const CMDecisionTree& CMDTMultiple::create( const scnXml::DTMultiple& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    unique_ptr<CMDTMultiple> self(new CMDTMultiple());
    for( const scnXml::DTCaseType& sn : node.getCaseType() ){
        self->children.push_back( &createBranch(sn, isUC, definitions, library) );
    }
    for( const scnXml::DTInfectionOrigin& sn : node.getInfectionOrigin() ){
        self->children.push_back( &createBranch(sn, isUC, definitions, library) );
    }
    for( const scnXml::DTDiagnostic& sn : node.getDiagnostic() ){
        self->children.push_back( &createBranch(sn, isUC, definitions, library) );
    }
    for( const scnXml::DTUncomplicated& sn : node.getUncomplicated() ){
        self->children.push_back( &createBranch(sn, isUC, definitions, library) );
    }
    for( const scnXml::DTSevere& sn : node.getSevere() ){
        self->children.push_back( &createBranch(sn, isUC, definitions, library) );
    }
    for( const scnXml::DTRandom& sn : node.getRandom() ){
        self->children.push_back( &CMDTRandom::create(sn, isUC, definitions, library) );
    }
    for( const scnXml::DTAge& sn : node.getAge() ){
        self->children.push_back( &CMDTAge::create(sn, isUC, definitions, library) );
    }
    for( const scnXml::DTCohort& sn : node.getCohort() ){
        self->children.push_back( &createBranch(sn, isUC, definitions, library) );
    }
    if( node.getTreatPKPD().size() ){
        self->children.push_back( &save_decision(library, new CMDTTreatPKPD(node.getTreatPKPD())) );
    }
    if( node.getTreatSimple().size() ){
        self->children.push_back( &save_decision(library, new CMDTTreatSimple(node.getTreatSimple())) );
    }
    if( node.getReport().size() ){
        self->children.push_back( &save_decision(library, new CMDTReport(node.getReport())) );
    }
    if( node.getDeploy().size() ){
        self->children.push_back( &save_decision(library, new CMDTDeploy(node.getDeploy(), definitions)) );
    }
    return save_decision(library, self.release());
}

const CMDecisionTree& createBranch( const scnXml::DTCaseType& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    if( !isUC ){
        throw util::xml_scenario_error( "decision tree: caseType can only be used for uncomplicated cases" );
    }
    return save_decision(library, new CMDTBranch(CMDTBranch::Kind::CaseType,
        CMDecisionTree::create( node.getFirstLine(), isUC, definitions, library),
        CMDecisionTree::create( node.getSecondLine(), isUC, definitions, library)
    ) );
}

const CMDecisionTree& createBranch( const scnXml::DTInfectionOrigin& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    return save_decision(library, new CMDTBranch(CMDTBranch::Kind::InfectionOrigin,
        CMDecisionTree::create( node.getImported(), isUC, definitions, library),
        CMDecisionTree::create( node.getIntroduced(), isUC, definitions, library),
        CMDecisionTree::create( node.getIndigenous(), isUC, definitions, library)
    ) );
}

const CMDecisionTree& createBranch( const scnXml::DTDiagnostic& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    return save_decision(library, new CMDTBranch(
        diagnostics::get( node.getDiagnostic() ),
        CMDecisionTree::create( node.getPositive(), isUC, definitions, library),
        CMDecisionTree::create( node.getNegative(), isUC, definitions, library)
    ) );
}

const CMDecisionTree& createBranch( const scnXml::DTUncomplicated& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    return save_decision(library, new CMDTBranch(
        UnitParse::readShortDuration(node.getMemory(), UnitParse::STEPS),
        CMDecisionTree::create( node.getPositive(), isUC, definitions, library),
        CMDecisionTree::create( node.getNegative(), isUC, definitions, library), library.healthSystemMemory
    ) );
}

const CMDecisionTree& createBranch( const scnXml::DTSevere& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    return save_decision(library, new CMDTBranch(CMDTBranch::Kind::Severe,
        CMDecisionTree::create( node.getPositive(), isUC, definitions, library),
        CMDecisionTree::create( node.getNegative(), isUC, definitions, library)
    ) );
}

const CMDecisionTree& createBranch( const scnXml::DTCohort& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    return save_decision(library, new CMDTBranch(
        definitions.getComponentId(node.getComponent()),
        CMDecisionTree::create( node.getPositive(), isUC, definitions, library),
        CMDecisionTree::create( node.getNegative(), isUC, definitions, library)
    ) );
}


const CMDecisionTree& CMDTRandom::create(
    const scnXml::DTRandom& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library)
{
    unique_ptr<CMDTRandom> result(new CMDTRandom());
    
    double cum_p = 0.0;
    for( const scnXml::Outcome& outcome : node.getOutcome() ){
        cum_p += outcome.getP();
        result->branches.insert( make_pair(cum_p, &CMDecisionTree::create( outcome, isUC, definitions, library)) );
    }
    
    // Test cum_p is approx. 1.0 in case the input tree is wrong. We require no
    // less than one to make sure generated random numbers are not greater than
    // the last option.
    if (cum_p < 1.0 || cum_p > 1.001){
        throw util::xml_scenario_error (
            "decision tree (random node): expected probability sum to be 1.0 but found " + to_string(cum_p));
    }
    
    return save_decision(library, result.release());
}

const CMDecisionTree& CMDTAge::create(const scnXml::DTAge& node, bool isUC, const interventions::Interventions& definitions, ClinicalParameters& library){
    unique_ptr<CMDTAge> result(new CMDTAge());
    
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
        lastNode = &CMDecisionTree::create(age, isUC, definitions, library);
        lastAge = age.getLb();
    }
    double noLb = numeric_limits<double>::infinity();
    result->branches.insert( make_pair(noLb, lastNode) );
    
    return save_decision(library, result.release());
}

CMDTOut executeCaseManagement(const CMHostData& hostData, const ClinicalParameters& config) {
    assert(hostData.pgState & Episode::SICK);
    const auto* tree = (hostData.pgState & Episode::COMPLICATED)
        ? config.escm_complicated : config.escm_uncomplicated;
    return tree->exec(hostData);
}

} }
