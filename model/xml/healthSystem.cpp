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

// model/xml/ conversion layer — health system / clinical side of the host
// modules: non-malaria fevers (PathogenesisModel), neonatal mortality,
// diagnostic selection, within-host treatment options and the whole
// decision-tree / case-management construction from scnXml.
// Moved verbatim from Host/WithinHost/Pathogenesis/PathogenesisModel.cpp,
// Host/NeonatalMortality.cpp and Host/WithinHost/Treatments.cpp. State stays
// defined in the original files (extern re-declared below).

#include "core/RunExtensions.h"
#include "Global.h"
#include "core/RunContext.h"
#include "Parameters.h"
#include "Host/NeonatalMortality.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "Host/WithinHost/Pathogenesis/Submodels.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "Clinical/CM5DayCommon.h"
#include "Clinical/CMDecisionTree.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/DecisionTree5Day.h"
#include "Clinical/ESCaseManagement.h"
#include "Clinical/Episode.h"
#include "Clinical/EventScheduler.h"
#include "Clinical/ImmediateOutcomes.h"
#include "PkPd/LSTMTreatments.h"
#include "interventions/Interfaces.h"
#include "interventions/InterventionManager.h"
#include "mon/Monitoring.h"
#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "schema/healthSystem.h"
#include "schema/scenario.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace std;

namespace OM {
namespace WithinHost {
using namespace OM::util;
namespace Pathogenesis {

// State defined in PathogenesisModel.cpp (external linkage):
extern double pg_indirRiskCoFactor, pg_severeMalThreshold, pg_comorbIntercept,
    pg_inv_critAgeComorb;
extern AgeGroupInterpolator pg_NMF_incidence;
extern bool opt_predetermined_episodes, opt_mueller_pres_model;

void PathogenesisModel::init( RunContext& ctx, const Parameters& parameters, const scnXml::Clinical& clinical, bool nmfOnly ){
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc = clinical.getNonMalariaFevers().get();
        pg_NMF_incidence.set( ctx.clock, nmfDesc.getIncidence(), "incidence" );
    }
    if( nmfOnly ) return;
    
    pg_indirRiskCoFactor = 1 - exp(-parameters[Parameter::INDIRECT_RISK_COFACTOR]);
    pg_severeMalThreshold = parameters[Parameter::SEVERE_MALARIA_THRESHHOLD] + 1;
    pg_comorbIntercept = 1 - exp(-parameters[Parameter::COMORBIDITY_INTERCEPT]);
    pg_inv_critAgeComorb = 1 / parameters[Parameter::CRITICAL_AGE_FOR_COMORBIDITY];

    if (util::ModelOptions::option (util::PREDETERMINED_EPISODES)) {
        opt_predetermined_episodes = true;
        //no separate init:
        PyrogenPathogenesis::init( ctx, parameters );
    } else {
        if (util::ModelOptions::option (util::MUELLER_PRESENTATION_MODEL)){
            opt_mueller_pres_model = true;
            MuellerPathogenesis::init( ctx, parameters );
        }else{
            PyrogenPathogenesis::init( ctx, parameters );
        }
    }
}

}

TreatmentId Treatments::addTreatment( RunContext& ctx, const scnXml::TreatmentOption& desc ){
    Treatments treatment( ctx, desc );
    TreatmentId id( ctx.ext->whTreatments.size() );
    ctx.ext->whTreatments.push_back( std::move(treatment) );
    return id;
}

Treatments::Stages stageFromString( const std::string& str ){
    if( str == "liver" ) return Treatments::LIVER;
    if( str == "blood" ) return Treatments::BLOOD;
    if( str == "both" ) return Treatments::BOTH;
    throw util::format_error( std::string("stage must be liver, blood or both, not ").append(str) );
}
Treatments::Treatments( RunContext& ctx, const scnXml::TreatmentOption& elt ) :
    TriggeredDeployments(ctx, elt), timeLiver(sim::zero()), timeBlood(sim::zero())
{
    for( auto it = elt.getClearInfections().begin(), end = elt.getClearInfections().end();
        it != end; ++it )
    {
        try{
            //NOTE: if changing XSD, this should not be called "timesteps" or have a default unit
            SimTime len = UnitParse::readShortDuration(ctx.clock,  it->getTimesteps(), UnitParse::STEPS );
            if( len < -ctx.clock.one_TS() || len == sim::zero() ){ 
                throw util::format_error( "timesteps must be ≥ 1 or have special value -1" );
            }
            Stages stage = stageFromString( it->getStage() );
            if( stage & LIVER ){
                if( timeLiver != sim::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of liver stage effect" );
                timeLiver = len;
            }
            if( stage & BLOOD ){
                if( timeBlood != sim::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of blood stage effect" );
                timeBlood = len;
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string(".../clearInfections/timesteps: ").append(e.message()) );
        }
    }
}

}

namespace Host {
using namespace OM::util;

// State defined in NeonatalMortality.cpp (external linkage):
extern std::vector<double> prevByGestationalAge;
extern const WithinHost::Diagnostic* neonatalDiagnostic;

void NeonatalMortality::init( RunContext& ctx, const scnXml::Clinical& clinical ){
    SimTime fiveMonths = sim::fromDays( 5 * 30 );
    prevByGestationalAge.assign( ctx.clock.inSteps(fiveMonths), 0.0 );
    
    if( clinical.getNeonatalMortality().present() ){
        neonatalDiagnostic = &WithinHost::diagnostics::get( ctx,
            clinical.getNeonatalMortality().get().getDiagnostic() );
    }else{
        //NOTE: this is a compatibility option for older scenarios
        if( CommandLine::option( CommandLine::DEPRECATION_WARNINGS ) ){
            cerr << "Deprecation warning: specification of the diagnostic "
                "used by the Neonatal Mortality model is recommended "
                "(model/clinical/neonatalMortality)" << endl;
        }
        neonatalDiagnostic = &WithinHost::diagnostics::monitoringDiagnostic( ctx );
    }
}
}
}

// ———  clinical / health-system case management  ———

namespace OM {
namespace Clinical {

// NOTE: the following block (decision-tree node classes, the decision
// library and the create() factories, plus the setHealthSystem()/init()
// readers of the clinical models) is moved verbatim from
// Clinical/CMDecisionTree.cpp, EventScheduler.cpp, ESCaseManagement.cpp,
// ImmediateOutcomes.cpp, DecisionTree5Day.cpp and ClinicalModel.cpp.
using WithinHost::Diagnostic;
using WithinHost::diagnostics;
using namespace OM::util;
using std::map;
using std::vector;

// State defined in ClinicalModel.cpp (external linkage):
extern bool opt_event_scheduler, opt_imm_outcomes;
extern double oddsRatioThreshold;
// State defined in EventScheduler.cpp (external linkage):
extern bool opt_non_malaria_fevers;
// State defined in ESCaseManagement.cpp (external linkage):
extern const CMDecisionTree *escm_uncomplicated, *escm_complicated;
// State defined in DecisionTree5Day.cpp (external linkage):
extern const CMDecisionTree *treeUCOfficial, *treeUCSelfTreat;

// ———  special 'multiple' node  ———

/**
 * Branch out to multiple descendants.
 */
class CMDTMultiple : public CMDecisionTree {
public:
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTMultiple& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTMultiple* p = dynamic_cast<const CMDTMultiple*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( children.size() != p->children.size() ) return false;
        for( size_t i = 0; i < children.size(); ++i ){
            if( children[i] != p->children[i] ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        CMDTOut result(false);
        for( auto it = children.begin(), end = children.end(); it != end; ++it ) {
            CMDTOut r2 = (*it)->exec( ctx, hostData );
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

/**
 * Should the patient recieve first or second line treatment?
 */
class CMDTCaseType : public CMDecisionTree {
public:
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTCaseType& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTCaseType* p = dynamic_cast<const CMDTCaseType*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( firstLine != p->firstLine ) return false;
        if( secondLine != p->secondLine ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        // Uses of this in complicated cases should trigger an exception during initialisation.
        assert( (hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED) );
        
        if( hostData.pgState & Episode::SECOND_CASE ) return secondLine.exec( ctx, hostData );
        else return firstLine.exec( ctx, hostData );
    }
    
private:
    CMDTCaseType( const CMDecisionTree& firstLine,
                  const CMDecisionTree& secondLine ) :
                  firstLine(firstLine), secondLine(secondLine)
                  {}
    
    const CMDecisionTree& firstLine;
    const CMDecisionTree& secondLine;
};

/**
 * Should the patient recieve a different treatment based on its infection type?
 */
class CMDTInfectionOrigin : public CMDecisionTree {
public:
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTInfectionOrigin& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTInfectionOrigin* p = dynamic_cast<const CMDTInfectionOrigin*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( imported != p->imported ) return false;
        if( introduced != p->introduced ) return false;
        if( indigenous != p->indigenous ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        WithinHost::InfectionOrigin InfectionOrigin = hostData.withinHost().getInfectionOrigin();
        if(InfectionOrigin == WithinHost::InfectionOrigin::Imported)
            return imported.exec( ctx, hostData );
        else if (InfectionOrigin == WithinHost::InfectionOrigin::Introduced)
            return introduced.exec( ctx, hostData );
        else
            return indigenous.exec( ctx, hostData );
    }
    
private:
    CMDTInfectionOrigin( const CMDecisionTree& imported,
                  const CMDecisionTree& introduced,
                  const CMDecisionTree& indigenous ) :
                  imported(imported), introduced(introduced), indigenous(indigenous)
                  {}
    
    const CMDecisionTree& imported, &introduced, &indigenous;
};

/**
 * Use a diagnostic, and call another decision tree based on the outcome.
 */
class CMDTDiagnostic : public CMDecisionTree {
public:
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTDiagnostic& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTDiagnostic* p = dynamic_cast<const CMDTDiagnostic*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( diagnostic != p->diagnostic ) return false;
        if( positive != p->positive ) return false;
        if( negative != p->negative ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        CMDTOut result;
        
        if( hostData.withinHost().diagnosticResult( hostData.human.rng, diagnostic ) )
            result = positive.exec( ctx, hostData );
        else
            result = negative.exec( ctx, hostData );
        result.screened = true;

        return result;
    }
    
private:
    CMDTDiagnostic( const Diagnostic& diagnostic,
        const CMDecisionTree& positive,
        const CMDecisionTree& negative ) :
        diagnostic(diagnostic),
        positive(positive), negative(negative)
    {}
    
    const Diagnostic& diagnostic;
    const CMDecisionTree& positive;
    const CMDecisionTree& negative;
};

class CMDTUncomplicated : public CMDecisionTree {
public:
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTUncomplicated& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTUncomplicated* p = dynamic_cast<const CMDTUncomplicated*>( &that );
        if( p == 0 ) return false;      // different type of node
        // if( diagnostic != p->diagnostic ) return false;
        if( positive != p->positive ) return false;
        if( negative != p->negative ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        const Clinical::Episode &latest = hostData.human.clinicalModel->getLatestReport();

        if ( ((hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED)) || (hostData.pgState & Episode::MALARIA) )
        {
            if (latest.time + memory >= ctx.clock.nowOrTs0())
                return positive.exec( ctx, hostData );
            else
                return negative.exec( ctx, hostData );
        }
        else
            return negative.exec( ctx, hostData );
    }
    
private:
    CMDTUncomplicated(
        const SimTime memory,
        const CMDecisionTree& positive,
        const CMDecisionTree& negative ) :
        memory(memory), positive(positive), negative(negative)
    {
        if(memory > ClinicalModel::hsMemory())
            throw util::xml_scenario_error( "<uncomplicated> memory parameter must be less than or equal to the healthsystem memory (hsmemory parameter)" );
    }
    
    const SimTime memory;
    const CMDecisionTree& positive;
    const CMDecisionTree& negative;
};

class CMDTSevere : public CMDecisionTree {
public:
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTSevere& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTSevere* p = dynamic_cast<const CMDTSevere*>( &that );
        if( p == 0 ) return false;      // different type of node
        // if( diagnostic != p->diagnostic ) return false;
        if( positive != p->positive ) return false;
        if( negative != p->negative ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        CMDTOut result;

        if (hostData.pgState & Episode::COMPLICATED)
            result = positive.exec( ctx, hostData );
        else
            result = negative.exec( ctx, hostData );

        // result.screened = true;
        return result;
    }
    
private:
    CMDTSevere(
        const CMDecisionTree& positive,
        const CMDecisionTree& negative ) :
        positive(positive), negative(negative)
    {}
    
    const CMDecisionTree& positive;
    const CMDecisionTree& negative;
};

/**
 * Choose a branch randomly.
 * 
 * TODO(optimisation): can we re-order the trees such that this can be abstracted out or
 * used only once?
 */
class CMDTRandom : public CMDecisionTree {
public:
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTRandom& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTRandom* p = dynamic_cast<const CMDTRandom*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( branches.size() != p->branches.size() ) return false;
        for( auto it1 = branches.begin(), it2 = p->branches.begin();
            it1 != branches.end(); ++it1, ++it2 )
        {
            if( it1->first != it2->first ) return false;
            if( it1->second != it2->second ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        // parity: d — decision-tree contract: child order = XSD sequence, cumulative keys, one uniform_01 per Random node (plan §3d)
        auto it = branches.upper_bound( hostData.human.rng.uniform_01() );
        assert( it != branches.end() );
        return it->second->exec( ctx, hostData );
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
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTAge& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTAge* p = dynamic_cast<const CMDTAge*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( branches.size() != p->branches.size() ) return false;
        for( auto it1 = branches.begin(), it2 = p->branches.begin();
            it1 != branches.end(); ++it1, ++it2 )
        {
            if( it1->first != it2->first ) return false;
            if( it1->second != it2->second ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        // age is that of human at start of time step (i.e. may be as low as 0)
        auto it = branches.upper_bound( hostData.ageYears );
        if( it == branches.end() )
            throw TRACED_EXCEPTION( "bad age-based decision tree switch", util::Error::PkPd );
        return it->second->exec( ctx, hostData );
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
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTNoTreatment* p = dynamic_cast<const CMDTNoTreatment*>( &that );
        return p != 0;  // same type: is equivalent
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
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
    virtual bool operator==( const CMDecisionTree& that ) const{
     if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTReport* p = dynamic_cast<const CMDTReport*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( outIds.size() != p->outIds.size() ) return false;
        for( auto it1 = outIds.begin(), it2 = p->outIds.begin();
            it1 != outIds.end(); ++it1, ++it2 )
        {
            if( *it1 != *it2 ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
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
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTTreatFailure* p = dynamic_cast<const CMDTTreatFailure*>( &that );
        return p != 0;  // same type: is equivalent
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        return CMDTOut(true /*report treatment*/);
    }
};

/**
 * Deliver one or more treatments via the PK/PD model.
 */
class CMDTTreatPKPD : public CMDecisionTree {
public:
    CMDTTreatPKPD( RunContext& ctx, const scnXml::DecisionTree::TreatPKPDSequence& seq ){
        treatments.reserve( seq.size() );
        for( const scnXml::DTTreatPKPD& treatElt : seq ){
            treatments.push_back( TreatInfo(
                ctx.treatments,
                treatElt.getSchedule(),
                treatElt.getDosage(),
                treatElt.getDelay_h()
            ) );
        }
        assert( treatments.size() > 0 );        // CMDTTreatPKPD should not be used in this case
    }
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTTreatPKPD* p = dynamic_cast<const CMDTTreatPKPD*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( treatments.size() != p->treatments.size() ) return false;
        for( auto it1 = treatments.begin(), it2 = p->treatments.begin();
            it1 != treatments.end(); ++it1, ++it2 )
        {
            if( *it1 != *it2 ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        for( const TreatInfo& treatment : treatments ){
            hostData.withinHost().treatPkPd( treatment.schedule, treatment.dosage, hostData.ageYears, 0.0 );
        }
        return CMDTOut(true);
    }
    
private:
    struct TreatInfo{
        TreatInfo( PkPd::TreatmentRegistry& tr, const string& s, const string& d, double h ) :
            schedule(tr.findSchedule(s)),
            dosage(tr.findDosages(d)),
            delay_h(h) {}
        inline bool operator!=( const TreatInfo& that )const{
            return schedule != that.schedule ||
                dosage != that.dosage ||
                delay_h != that.delay_h;
        }
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
    CMDTTreatSimple( RunContext& ctx, const scnXml::DecisionTree::TreatSimpleSequence& elt ) :
        timeLiver(sim::zero()), timeBlood(sim::zero())
    {
        for (const auto &treatElt : elt) {
            //NOTE: this code is currently identical to that in SimpleTreatComponent
            try{
                SimTime durL = UnitParse::readShortDuration(ctx.clock,  treatElt.getDurationLiver(), UnitParse::NONE ),
                    durB = UnitParse::readShortDuration(ctx.clock,  treatElt.getDurationBlood(), UnitParse::NONE );
                SimTime neg1 = -ctx.clock.one_TS();
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
    virtual bool operator==( const CMDecisionTree& that ) const{
     if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTTreatSimple* p = dynamic_cast<const CMDTTreatSimple*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( timeLiver.size() != p->timeLiver.size() ) return false;
        if( timeBlood.size() != p->timeBlood.size() ) return false;
        for( auto it1 = timeLiver.begin(), it2 = p->timeLiver.begin();
            it1 != timeLiver.end(); ++it1, ++it2 )
        {
            if( *it1 != *it2 ) return false;
        }
        return true;    // no tests failed; must be the same
    }

    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
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
    CMDTDeploy( RunContext& ctx, const scnXml::DecisionTree::DeploySequence& seq ) :
        HumanIntervention(seq) {}
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTDeploy* p = dynamic_cast<const CMDTDeploy*>( &that );
        if( p == 0 ) return false;      // different type of node
        if( components.size() != p->components.size() ) return false;
        for( auto it1 = components.begin(), it2 = p->components.begin();
            it1 != components.end(); ++it1, ++it2 )
        {
            // here we can compare pointers since components are always de-duplicated
            if( *it1 != *it2 ) return false;
        }
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        deploy( ctx, hostData.human,
                  mon::Deploy::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/) );
        
        //NOTE: it's not intuitively obvious what value should be returned here
        // in the case of intervention deployment. This at least means that
        // repeat seekers get second-line treatment.
        return CMDTOut(false);
    }
};

class CMDTCohort : public CMDecisionTree {
public:
    static const CMDecisionTree& create( RunContext& ctx, const ::scnXml::DTCohort& node, bool isUC );
    
protected:
    virtual bool operator==( const CMDecisionTree& that ) const{
        if( this == &that ) return true; // short cut: same object thus equivalent
        const CMDTCohort* p = dynamic_cast<const CMDTCohort*>( &that );
        if( p == 0 ) return false;      // different type of node
        // if( diagnostic != p->diagnostic ) return false;
        if( positive != p->positive ) return false;
        if( negative != p->negative ) return false;
        return true;    // no tests failed; must be the same
    }
    
    virtual CMDTOut exec( RunContext& ctx, CMHostData hostData ) const{
        CMDTOut result;

        // Rely on the health system memory to not count the same episode twice
        if(hostData.human.isInSubPop(component))
            result = positive.exec( ctx, hostData );
        else
            result = negative.exec( ctx, hostData );

        // result.screened = true;
        return result;
    }
    
private:
    CMDTCohort(
        const std::string component,
        const CMDecisionTree& positive,
        const CMDecisionTree& negative ) :
        component(interventions::InterventionManager::getComponentId(component)), positive(positive), negative(negative)
    {
    }
    
    interventions::ComponentId component;
    const CMDecisionTree& positive;
    const CMDecisionTree& negative;
};

// ———  static functions  ———

// Saves a decision to the run's decision library (ctx.ext->decision_library;
// frees memory at end of run). We store pointers into this list, so elements
// must not move. Also optimises away duplicates.
const CMDecisionTree& save_decision( RunContext& ctx, CMDecisionTree* decision ){
    auto& decision_library = ctx.ext->decision_library;
    // We search the library for a duplicate, and delete this one if there is
    // a duplicate. Note that this is not implemented efficiently, but a little
    // wasted time at start up is hardly a concern.
    for( auto& d : decision_library ){
        if( *d == *decision ){
            return *d;
        }
    }
    
    // No match: add to the library.
    decision_library.push_back( unique_ptr<CMDecisionTree>(decision) );
    return *decision_library.back();
}

const CMDecisionTree& CMDecisionTree::create( RunContext& ctx, const scnXml::DecisionTree& node, bool isUC){
    if( node.getMultiple().present() ) return CMDTMultiple::create( ctx, node.getMultiple().get(), isUC);
    // branching nodes
    if( node.getCaseType().present() ) return CMDTCaseType::create( ctx, node.getCaseType().get(), isUC);
    if( node.getInfectionOrigin().present() ) return CMDTInfectionOrigin::create( ctx, node.getInfectionOrigin().get(), isUC);
    if( node.getDiagnostic().present() ) return CMDTDiagnostic::create( ctx, node.getDiagnostic().get(), isUC);
    if( node.getUncomplicated().present() ) return CMDTUncomplicated::create( ctx, node.getUncomplicated().get(), isUC);
    if( node.getSevere().present() ) return CMDTSevere::create( ctx, node.getSevere().get(), true);
    if( node.getRandom().present() ) return CMDTRandom::create( ctx, node.getRandom().get(), isUC);
    if( node.getAge().present() ) return CMDTAge::create( ctx, node.getAge().get(), isUC);
    if( node.getCohort().present() ) return CMDTCohort::create( ctx, node.getCohort().get(), isUC);

    // action nodes
    if( node.getNoTreatment().present() ) return save_decision( ctx, new CMDTNoTreatment() );
    if( node.getReport().size() ) return save_decision( ctx, new CMDTReport(node.getReport()) );
    if( node.getTreatFailure().present() ) return save_decision( ctx, new CMDTTreatFailure() );
    if( node.getTreatPKPD().size() ) return save_decision( ctx, new CMDTTreatPKPD( ctx, node.getTreatPKPD() ) );
    if( node.getTreatSimple().size() ) return save_decision( ctx, new CMDTTreatSimple( ctx, node.getTreatSimple() ) );
    if( node.getDeploy().size() ) return save_decision( ctx, new CMDTDeploy( ctx, node.getDeploy() ) );
    throw xml_scenario_error( "unterminated decision tree" );
}

const CMDecisionTree& CMDTMultiple::create( RunContext& ctx, const scnXml::DTMultiple& node, bool isUC ){
    CMDTMultiple* self = new CMDTMultiple();
    for( const scnXml::DTCaseType& sn : node.getCaseType() ){
        // parity: d — child chain built in XSD element-sequence order (plan §3d)
        self->children.push_back( &CMDTCaseType::create( ctx, sn, isUC) );
    }
    for( const scnXml::DTInfectionOrigin& sn : node.getInfectionOrigin() ){
        self->children.push_back( &CMDTInfectionOrigin::create( ctx, sn, isUC) );
    }
    for( const scnXml::DTDiagnostic& sn : node.getDiagnostic() ){
        self->children.push_back( &CMDTDiagnostic::create( ctx, sn, isUC) );
    }
    for( const scnXml::DTUncomplicated& sn : node.getUncomplicated() ){
        self->children.push_back( &CMDTUncomplicated::create( ctx, sn, isUC) );
    }
    for( const scnXml::DTSevere& sn : node.getSevere() ){
        self->children.push_back( &CMDTSevere::create( ctx, sn, isUC) );
    }
    for( const scnXml::DTRandom& sn : node.getRandom() ){
        self->children.push_back( &CMDTRandom::create( ctx, sn, isUC) );
    }
    for( const scnXml::DTAge& sn : node.getAge() ){
        self->children.push_back( &CMDTAge::create( ctx, sn, isUC) );
    }
    for( const scnXml::DTCohort& sn : node.getCohort() ){
        self->children.push_back( &CMDTCohort::create( ctx, sn, isUC) );
    }
    if( node.getTreatPKPD().size() ){
        self->children.push_back( &save_decision( ctx, new CMDTTreatPKPD(ctx, node.getTreatPKPD())) );
    }
    if( node.getTreatSimple().size() ){
        self->children.push_back( &save_decision( ctx, new CMDTTreatSimple(ctx, node.getTreatSimple())) );
    }
    if( node.getReport().size() ){
        self->children.push_back( &save_decision( ctx, new CMDTReport(node.getReport())) );
    }
    if( node.getDeploy().size() ){
        self->children.push_back( &save_decision( ctx, new CMDTDeploy(ctx, node.getDeploy())) );
    }
    return save_decision( ctx, self );
}

const CMDecisionTree& CMDTCaseType::create( RunContext& ctx, const scnXml::DTCaseType& node, bool isUC ){
    if( !isUC ){
        throw util::xml_scenario_error( "decision tree: caseType can only be used for uncomplicated cases" );
    }
    return save_decision( ctx, new CMDTCaseType(
        CMDecisionTree::create( ctx, node.getFirstLine(), isUC ),
        CMDecisionTree::create( ctx, node.getSecondLine(), isUC )
    ) );
}

const CMDecisionTree& CMDTInfectionOrigin::create( RunContext& ctx, const scnXml::DTInfectionOrigin& node, bool isUC ){
    return save_decision( ctx, new CMDTInfectionOrigin(
        CMDecisionTree::create( ctx, node.getImported(), isUC ),
        CMDecisionTree::create( ctx, node.getIntroduced(), isUC ),
        CMDecisionTree::create( ctx, node.getIndigenous(), isUC )
    ) );
}

const CMDecisionTree& CMDTDiagnostic::create( RunContext& ctx, const scnXml::DTDiagnostic& node, bool isUC ){
    return save_decision( ctx, new CMDTDiagnostic(
        diagnostics::get( ctx, node.getDiagnostic() ),
        CMDecisionTree::create( ctx, node.getPositive(), isUC ),
        CMDecisionTree::create( ctx, node.getNegative(), isUC )
    ) );
}

const CMDecisionTree& CMDTUncomplicated::create( RunContext& ctx, const scnXml::DTUncomplicated& node, bool isUC ){
    return save_decision( ctx, new CMDTUncomplicated(
        UnitParse::readShortDuration(ctx.clock, node.getMemory(), UnitParse::STEPS),
        CMDecisionTree::create( ctx, node.getPositive(), isUC ),
        CMDecisionTree::create( ctx, node.getNegative(), isUC )
    ) );
}

const CMDecisionTree& CMDTSevere::create( RunContext& ctx, const scnXml::DTSevere& node, bool isUC ){
    return save_decision( ctx, new CMDTSevere(
        CMDecisionTree::create( ctx, node.getPositive(), isUC ),
        CMDecisionTree::create( ctx, node.getNegative(), isUC )
    ) );
}

const CMDecisionTree& CMDTCohort::create( RunContext& ctx, const scnXml::DTCohort& node, bool isUC ){
    return save_decision( ctx, new CMDTCohort(
        node.getComponent(),
        CMDecisionTree::create( ctx, node.getPositive(), isUC ),
        CMDecisionTree::create( ctx, node.getNegative(), isUC )
    ) );
}


const CMDecisionTree& CMDTRandom::create(RunContext& ctx,
    const scnXml::DTRandom& node, bool isUC )
{
    CMDTRandom* result = new CMDTRandom();
    
    double cum_p = 0.0;
    for( const scnXml::Outcome& outcome : node.getOutcome() ){
        // parity: d — cumulative keys summed in XSD outcome order; map<double> numeric key order; last key must be 1 (plan §3d)
        cum_p += outcome.getP();
        result->branches.insert( make_pair(cum_p, &CMDecisionTree::create( ctx, outcome, isUC )) );
    }
    
    // Test cum_p is approx. 1.0 in case the input tree is wrong. We require no
    // less than one to make sure generated random numbers are not greater than
    // the last option.
    if (cum_p < 1.0 || cum_p > 1.001){
        throw util::xml_scenario_error (
            "decision tree (random node): expected probability sum to be 1.0 but found " + to_string(cum_p));
    }
    
    return save_decision( ctx, result );
}

const CMDecisionTree& CMDTAge::create(RunContext& ctx, const scnXml::DTAge& node, bool isUC){
    CMDTAge* result = new CMDTAge();
    
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
        lastNode = &CMDecisionTree::create( ctx, age, isUC);
        lastAge = age.getLb();
    }
    double noLb = numeric_limits<double>::infinity();
    result->branches.insert( make_pair(noLb, lastNode) );
    
    return save_decision( ctx, result );
}
void ClinicalEventScheduler::init( RunContext& ctx, const Parameters& parameters, const scnXml::Clinical& clinical )
{
    if( ctx.clock.one_TS() != sim::oneDay() )
        throw util::xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day time step.");
    
    opt_non_malaria_fevers = util::ModelOptions::option( util::NON_MALARIA_FEVERS );
    
    alpha = exp( -parameters[Parameter::CFR_NEG_LOG_ALPHA] );
    if( !(0.0<=alpha && alpha<=1.0) ){
        throw util::xml_scenario_error(
            "Clinical outcomes: propDeathsFirstDay should be within range [0,1]");
    }
    neg_v = -parameters[Parameter::CFR_SCALE_FACTOR];
    
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc2 =
            clinical.getNonMalariaFevers().get();
        if( !nmfDesc2.getPrNeedTreatmentMF().present() ||
            !nmfDesc2.getPrNeedTreatmentNMF().present() ){
            throw util::xml_scenario_error( "prNeedTreatmentMF and "
                "prNeedTreatmentNMF elements required in model->clinical->"
                "NonMalariaFevers" );
        }
        NMF_need_antibiotic.set( ctx.clock, nmfDesc2.getPrNeedTreatmentNMF().get(), "prNeedTreatmentNMF" );
        MF_need_antibiotic.set( ctx.clock, nmfDesc2.getPrNeedTreatmentMF().get(), "prNeedTreatmentMF" );
    }
}

void ClinicalEventScheduler::setParameters(const SimClock& clock, const scnXml::HSEventScheduler& esData) {
    const scnXml::ClinicalOutcomes& coData = esData.getClinicalOutcomes();
    
    maxUCSeekingMemory = sim::fromDays(coData.getMaxUCSeekingMemory());
    uncomplicatedCaseDuration = sim::fromDays(coData.getUncomplicatedCaseDuration());
    complicatedCaseDuration = sim::fromDays(coData.getComplicatedCaseDuration());
    extraDaysAtRisk = sim::fromDays(coData.getComplicatedRiskDuration()) - complicatedCaseDuration;
    if( uncomplicatedCaseDuration < sim::fromDays(1)
	|| complicatedCaseDuration < sim::fromDays(1)
	|| maxUCSeekingMemory < sim::zero()
	|| extraDaysAtRisk + complicatedCaseDuration < sim::fromDays(1) // at risk at least 1 day
	|| extraDaysAtRisk > sim::zero()        // at risk longer than case duration
    ){
	throw util::xml_scenario_error(
            "Clinical outcomes: constraints on case/risk/memory duration not met (see documentation)");
    }
    
    cumDailyPrImmUCTS.reserve( coData.getDailyPrImmUCTS().size() );
    double cumP = 0.0;
    for( auto it = coData.getDailyPrImmUCTS().begin(); it != coData.getDailyPrImmUCTS().end(); ++it ){
        cumP += *it;
        cumDailyPrImmUCTS.push_back( cumP );
    }
    if( cumP < 0.999 || cumP > 1.001 ){
        throw util::xml_scenario_error( "Event scheduler: dailyPrImmUCTS seq must add up to 1" );
    }
    cumDailyPrImmUCTS.back() = 1.0;
    
    caseFatalityRate.scale( alpha );
    
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !esData.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of healthSystem->EventScheduler required");
        }
        const scnXml::HSESNMF& nmfDesc = esData.getNonMalariaFevers().get();
        
        double pT = nmfDesc.getPrTreatment();
        logOddsAbBase = log( pT / (1.0 - pT) );
        logOddsAbNegTest = log(nmfDesc.getEffectNegativeTest());
        logOddsAbPosTest = log(nmfDesc.getEffectPositiveTest());
        logOddsAbNeed = log(nmfDesc.getEffectNeed());
        logOddsAbInformal = log(nmfDesc.getEffectInformal());
        oneMinusEfficacyAb = 1.0 - nmfDesc.getTreatmentEfficacy();
        severeNmfMortality.set( clock, nmfDesc.getCFR(), "CFR" );
    }
}
void ESCaseManagement::setHealthSystem(RunContext& ctx, const scnXml::HSEventScheduler& esData){
    escm_uncomplicated = &CMDecisionTree::create( ctx, esData.getUncomplicated(), true );
    escm_complicated = &CMDecisionTree::create( ctx, esData.getComplicated(), false );
    
    // Calling our parent class like this is messy. Changing this would require
    // moving change-of-health-system handling into ClinicalModel.
    ClinicalEventScheduler::setParameters( ctx.clock, esData );
}
// ———  static, utility functions  ———

double getHealthSystemACRByName (const scnXml::TreatmentDetails& td, const string& drug)
{
    double val;
    if (drug == "CQ")
        val = td.getCQ().present() ? td.getCQ().get().getValue() : 0;
    else if (drug == "SP")
        val = td.getSP().present() ? td.getSP().get().getValue() : 0;
    else if (drug == "AQ")
        val = td.getAQ().present() ? td.getAQ().get().getValue() : 0;
    else if (drug == "SPAQ")
        val = td.getSPAQ().present() ? td.getSPAQ().get().getValue() : 0;
    else if (drug == "ACT")
        val = td.getACT().present() ? td.getACT().get().getValue() : 0;
    else if (drug == "QN")
        val = td.getQN().present() ? td.getQN().get().getValue() : 0;
    else if (drug == "selfTreatment")
        val = td.getSelfTreatment().getValue();
    else {
        throw util::xml_scenario_error ("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
    }
    
    if( !(val >= 0.0 && val <= 1.0) ){
        throw util::xml_scenario_error (string("healthSystem initialACR/compliance/nonCompliersEffective: ").append(drug).append(" value must be in range [0,1]"));
    }
    return val;
}

// note: "CQOptional" is the same type as for other drugs
const scnXml::TreatmentActions::CQOptional& getHealthSystemTreatmentOptional(
    const scnXml::TreatmentActions& treatActs, const string& drug )
{
    if( drug == "CQ" ) return treatActs.getCQ();
    if( drug == "SP" ) return treatActs.getSP();
    if( drug == "AQ" ) return treatActs.getAQ();
    if( drug == "SPAQ" ) return treatActs.getSPAQ();
    if( drug == "ACT" ) return treatActs.getACT();
    if( drug == "QN" ) return treatActs.getQN();
    throw util::xml_scenario_error ("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
}
WithinHost::TreatmentId getHealthSystemTreatmentByName( RunContext& ctx, const scnXml::TreatmentActions& treatActs, const string& drug )
{
    const scnXml::TreatmentActions::CQOptional& elt = getHealthSystemTreatmentOptional( treatActs, drug );
    if( !elt.present() )
        throw util::xml_scenario_error( string("healthSystem.treatmentActions: description required for ").append(drug) );
    return WithinHost::WHInterface::addTreatment( ctx, elt.get() );
}


// ———  static, init  ———

void ImmediateOutcomes::setHealthSystem (RunContext& ctx, const scnXml::HSImmediateOutcomes& hsDescription) {
    const string &firstLine = hsDescription.getDrugRegimen().getFirstLine(),
        &secondLine = hsDescription.getDrugRegimen().getSecondLine(),
        &inpatient = hsDescription.getDrugRegimen().getInpatient();
    
    // ———  calculate probability of getting treatment  ———

    double accessUCOfficial1 = hsDescription.getPSeekOfficialCareUncomplicated1().getValue();
    double accessUCOfficial2 = hsDescription.getPSeekOfficialCareUncomplicated2().getValue();
    // Note: this asymmetry is historical, and probably matters little:
    accessUCSelfTreat[FirstLine] = hsDescription.getPSelfTreatUncomplicated().getValue();
    accessUCSelfTreat[SecondLine] = 0.0;
    accessUCAny[FirstLine] = accessUCOfficial1 + accessUCSelfTreat[FirstLine];
    accessUCAny[SecondLine] = accessUCOfficial1 + accessUCSelfTreat[SecondLine];
    accessSevere = hsDescription.getPSeekOfficialCareSevere().getValue();
    
    if( accessUCOfficial1 < 0.0 || accessUCSelfTreat[FirstLine] < 0.0 || accessUCAny[FirstLine] > 1.0 ||
        accessUCOfficial2 < 0.0 || accessUCSelfTreat[SecondLine] < 0.0 || accessUCAny[SecondLine] > 1.0 ||
        accessSevere < 0.0 || accessSevere > 1.0 )
    {
        throw util::xml_scenario_error ("healthSystem: "
            "pSeekOfficialCareUncomplicated1 and pSelfTreatUncomplicated must be"
            " at least 0 and their sum must be at most 1, and"
            "pSeekOfficialCareUncomplicated2 and pSeekOfficialCareSevere must "
            "be in range [0,1]");
    }
    
    // ———  calculate probability of clearing parasites  ———
    
    const double complianceFirstLine = getHealthSystemACRByName (hsDescription.getCompliance(), firstLine);
    const double complianceSecondLine = getHealthSystemACRByName (hsDescription.getCompliance(), secondLine);

    const double cureRateFirstLine = getHealthSystemACRByName (hsDescription.getInitialACR(), firstLine);
    const double cureRateSecondLine = getHealthSystemACRByName (hsDescription.getInitialACR(), secondLine);

    const double nonCompliersEffectiveFirstLine = getHealthSystemACRByName (hsDescription.getNonCompliersEffective(), firstLine);
    const double nonCompliersEffectiveSecondLine = getHealthSystemACRByName (hsDescription.getNonCompliersEffective(), secondLine);

    const double complianceSelfTreatment = hsDescription.getCompliance().getSelfTreatment().getValue();
    const double cureRateSelfTreatment = hsDescription.getInitialACR().getSelfTreatment().getValue();
    if( !(
        complianceSelfTreatment >= 0.0 && complianceSelfTreatment <= 1.0
        && cureRateSelfTreatment >= 0.0 && cureRateSelfTreatment <= 1.0
    ) ){
        throw util::xml_scenario_error ("healthSystem initialACR/compliance/nonCompliersEffective: pSelfTreatment must be in range [0,1]");
    }
    
    
    cureRateUCOfficial[FirstLine] = complianceFirstLine * cureRateFirstLine
            + (1 - complianceFirstLine) * nonCompliersEffectiveFirstLine;
    cureRateUCOfficial[SecondLine] = complianceSecondLine * cureRateSecondLine
            + (1 - complianceSecondLine) * nonCompliersEffectiveSecondLine;
    cureRateUCSelfTreat[FirstLine] = complianceSelfTreatment * cureRateSelfTreatment
            + (1 - complianceSelfTreatment) * nonCompliersEffectiveFirstLine;
    cureRateUCSelfTreat[SecondLine] = 0.0;      // as with access, this is a historical asymmetry
    
    cureRateSevere = getHealthSystemACRByName (hsDescription.getInitialACR(), inpatient);
    
    treatmentUC[FirstLine] = getHealthSystemTreatmentByName(ctx, hsDescription.getTreatmentActions(), firstLine);
    
    if( secondLine == firstLine ) treatmentUC[SecondLine] = treatmentUC[FirstLine];
    else treatmentUC[SecondLine] = getHealthSystemTreatmentByName(ctx, hsDescription.getTreatmentActions(), secondLine);
    
    if( inpatient == firstLine ) treatmentSevere = treatmentUC[FirstLine];
    else if( inpatient == secondLine ) treatmentSevere = treatmentUC[SecondLine];
    else treatmentSevere = getHealthSystemTreatmentByName(ctx, hsDescription.getTreatmentActions(), inpatient);
    
    useDiagnosticUC = hsDescription.getUseDiagnosticUC();
    
    if( ModelOptions::option(util::VIVAX_SIMPLE_MODEL) ){
        WithinHost::WHVivax::setHSParameters(
            hsDescription.getLiverStageDrug().present() ?
            &hsDescription.getLiverStageDrug().get() : 0 );
    }else if( hsDescription.getLiverStageDrug().present() ){
        throw util::xml_scenario_error( "health-system's liverStageDrug element only supported by vivax" );
    }
}
void DecisionTree5Day::setHealthSystem(RunContext& ctx, const scnXml::HSDT5Day& hsDescription){
    double accessUCOfficial1 = hsDescription.getPSeekOfficialCareUncomplicated1().getValue();
    double accessUCOfficial2 = hsDescription.getPSeekOfficialCareUncomplicated2().getValue();
    // Note: this asymmetry is historical, and probably matters little:
    accessUCSelfTreat[FirstLine] = hsDescription.getPSelfTreatUncomplicated().getValue();
    accessUCSelfTreat[SecondLine] = 0.0;
    accessUCAny[FirstLine] = accessUCOfficial1 + accessUCSelfTreat[FirstLine];
    accessUCAny[SecondLine] = accessUCOfficial1 + accessUCSelfTreat[SecondLine];
    accessSevere = hsDescription.getPSeekOfficialCareSevere().getValue();
    
    if( accessUCOfficial1 < 0.0 || accessUCSelfTreat[FirstLine] < 0.0 || accessUCAny[FirstLine] > 1.0 ||
        accessUCOfficial2 < 0.0 || accessUCSelfTreat[SecondLine] < 0.0 || accessUCAny[SecondLine] > 1.0 ||
        accessSevere < 0.0 || accessSevere > 1.0 )
    {
        throw util::xml_scenario_error ("healthSystem: "
            "pSeekOfficialCareUncomplicated1 and pSelfTreatUncomplicated must be"
            " at least 0 and their sum must be at most 1, and"
            "pSeekOfficialCareUncomplicated2 and pSeekOfficialCareSevere must "
            "be in range [0,1]");
    }
    
    treeUCOfficial = &CMDecisionTree::create( ctx, hsDescription.getTreeUCOfficial(), true );
    treeUCSelfTreat = &CMDecisionTree::create( ctx, hsDescription.getTreeUCSelfTreat(), true );
    
    cureRateSevere = hsDescription.getCureRateSevere().getValue();
    treatmentSevere = WHInterface::addTreatment( ctx, hsDescription.getTreatmentSevere() );
    
    if( ModelOptions::option(util::VIVAX_SIMPLE_MODEL) ){
        WithinHost::WHVivax::setHSParameters(
            hsDescription.getLiverStageDrug().present() ?
            &hsDescription.getLiverStageDrug().get() : 0 );
    }else if( hsDescription.getLiverStageDrug().present() ){
        throw util::xml_scenario_error( "health-system's liverStageDrug element only supported by vivax" );
    }
}
void ClinicalModel::init( RunContext& ctx, const Parameters& parameters, const scnXml::Scenario& scenario ) {
    const scnXml::Clinical& clinical = scenario.getModel().getClinical();
    try{
        indirectMortBugfix = util::ModelOptions::option (util::INDIRECT_MORTALITY_FIX);
        //NOTE: if changing XSD, this should not have a default unit:
        healthSystemMemory = UnitParse::readShortDuration(ctx.clock, clinical.getHealthSystemMemory(), UnitParse::STEPS);
        oddsRatioThreshold = exp( parameters[Parameter::LOG_ODDS_RATIO_CF_COMMUNITY] );
        InfantMortality::init( ctx, parameters );
    }catch( const util::format_error& e ){
        throw util::xml_scenario_error( string("model/clinical/healthSystemMemory: ").append(e.message()) );
    }
    
    if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER)){
        opt_event_scheduler = true;
        ClinicalEventScheduler::init( ctx, parameters, clinical );
    }else{
        if( scenario.getHealthSystem().getImmediateOutcomes().present() ){
            opt_imm_outcomes = true;
            
            if( util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS) ){
                cerr << "Deprecation warning: healthSystem: use of ImmediateOutcomes can be replaced by the more flexible DecisionTree5Day (optional)" << endl;
            }
        }
        // else: decision tree 5 day
        
        CM5DayCommon::init();
    }
}

void ClinicalModel::setHS( RunContext& ctx, const scnXml::HealthSystem& healthSystem ){
    caseFatalityRate.set( ctx.clock, healthSystem.getCFR(), "CFR" );
    pSequelaeInpatient.set( ctx.clock, healthSystem.getPSequelaeInpatient(), "pSequelaeInpatient" );
    if( opt_event_scheduler ){
        if( !healthSystem.getEventScheduler().present() ){
            throw util::xml_scenario_error ("Expected EventScheduler "
                "section in healthSystem data (initial or intervention)");
        }
        ESCaseManagement::setHealthSystem(ctx, healthSystem.getEventScheduler().get());
    }else if( opt_imm_outcomes ){
        if ( !healthSystem.getImmediateOutcomes().present() ){
            throw util::xml_scenario_error ("Expected ImmediateOutcomes "
                "section in healthSystem data (initial or intervention)");
        }
        ImmediateOutcomes::setHealthSystem(ctx, healthSystem.getImmediateOutcomes().get());
    }else{
        if( !healthSystem.getDecisionTree5Day().present() ){
            throw util::xml_scenario_error ("Expected DecisionTree5Day "
                "section in healthSystem data (initial or intervention)");
        }
        DecisionTree5Day::setHealthSystem(ctx, healthSystem.getDecisionTree5Day().get());
    }
}
}
}
