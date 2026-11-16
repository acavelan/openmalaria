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

// Unittest for the EventScheduler case management

#ifndef Hmod_CMDecisionTreeSuite
#define Hmod_CMDecisionTreeSuite

#include "interventions/Interventions.h"
#include <cxxtest/TestSuite.h>
#include "Clinical/CMDecisionTree.h"
#include "Clinical/ClinicalXml.h"
#include "Clinical/DecisionTree5Day.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/EventScheduler.h"
#include "util/random.h"
#include "UnittestUtil.h"
#include "WHMock.h"
#include <limits>

using namespace OM::Clinical;
using namespace OM::WithinHost;
using UnitTest::WHMock;

// Internal interning entry point: verify ownership without exposing it in the model API.
namespace OM::Clinical { const CMDecisionTree& save_decision(ClinicalParameters&, CMDecisionTree*); }

class CountedDecision : public CMDecisionTree {
public:
    inline static int live = 0;
    CountedDecision() { ++live; }
    ~CountedDecision() override { --live; }
    bool operator==(const CMDecisionTree& other) const override {
        return dynamic_cast<const CountedDecision*>(&other) != nullptr;
    }
    CMDTOut exec(CMHostData) const override { return CMDTOut(false); }
};

class CMDecisionTreeSuite : public CxxTest::TestSuite
{
    OM::ModelParameters parameters;
    OM::interventions::Interventions definitions;
    util::MasterRng masterRng{0, 0};
public:
    CMDecisionTreeSuite () {}
    ~CMDecisionTreeSuite () {}
    
    void setUp () {
        // Note: cannot create whm in constructor, since it uses random number
        // generator which is initialized after constructor runs.
        UnittestUtil::initTime(5);
        UnittestUtil::setDiagnostics();

        UnittestUtil::EmpiricalWHM_setup(parameters.withinHost);

        human.reset( UnittestUtil::createHuman(sim::nowOrTs0(), parameters, masterRng).release() );

        ETS_ASSERT( human.get() != 0 );
        whm = dynamic_cast<WHMock*>(UnittestUtil::setHumanWH( *human,
                unique_ptr<WithinHost::WHInterface>(new WHMock()) ));
        ETS_ASSERT( whm != 0 );

        hd.reset( new CMHostData( *human, 21 /* any real age will do for most tests */,
                                      Episode::NONE ) );

        // could seed random-number-generator, but shouldn't affect outcomes
        UnittestUtil::PkPdSuiteSetup ();
    }
    void tearDown () {
        human.reset();
        whm = nullptr;
        PkPd::LSTMDrugType::clear();
        PkPd::LSTMTreatments::clear();
    }
    
    /* Runs the decision tree N times
     * returns the proportion of these runs where the output was any treatment
     */
    double propTreatmentsNReps (int N, const scnXml::DecisionTree& dt) {
        const CMDecisionTree& cmdt = CMDecisionTree::create( dt, true, definitions, parameters.clinical );
        
	whm->nTreatments = 0;
        int secondCounter = 0;
	for (int i = 0; i < N; ++i) {
	    secondCounter += (cmdt.exec( *hd ).treated ? 1 : 0);
	}
	TS_ASSERT_EQUALS( whm->nTreatments, secondCounter );
	return double(whm->nTreatments) / double(N);
    }
    
    void testTreatmentValuesAndValidation() {
        scnXml::TreatmentOption input;
        input.getClearInfections().push_back(scnXml::ClearInfections("-1t", "liver"));
        input.getClearInfections().push_back(scnXml::ClearInfections("2t", "blood"));
        auto treatment = OM::xml::readTreatment(input, definitions);
        TS_ASSERT_EQUALS(treatment.timeLiver, -sim::oneTS());
        TS_ASSERT_EQUALS(treatment.timeBlood, 2 * sim::oneTS());
        ClinicalParameters other;
        other.treatmentSevere = treatment;
        treatment.timeBlood = 3 * sim::oneTS();
        TS_ASSERT_EQUALS(other.treatmentSevere.timeBlood, 2 * sim::oneTS());
        whm->nTreatments = 0;
        whm->treatment(*human, other.treatmentSevere);
        TS_ASSERT_EQUALS(whm->nTreatments, 1);
        TS_ASSERT_EQUALS(whm->lastTimeLiver, -sim::oneTS());
        TS_ASSERT_EQUALS(whm->lastTimeBlood, 2 * sim::oneTS());

        input.getClearInfections().push_back(scnXml::ClearInfections("1t", "both"));
        TS_ASSERT_THROWS(OM::xml::readTreatment(input, definitions), const OM::util::xml_scenario_error&);
        for (const auto& bad : {scnXml::ClearInfections("0t", "blood"),
                               scnXml::ClearInfections("-2t", "liver"),
                               scnXml::ClearInfections("1t", "unknown")}) {
            input.getClearInfections().clear();
            input.getClearInfections().push_back(bad);
            TS_ASSERT_THROWS(OM::xml::readTreatment(input, definitions), const OM::util::xml_scenario_error&);
        }
    }

    void testClinicalConfigurationAndActiveTreesAreIndependent() {
        ClinicalParameters first, second;
        first.healthSystemMemory = 10;
        second.healthSystemMemory = 25;
        second.kind = ClinicalParameters::Kind::ImmediateOutcomes;
        auto a = ClinicalModel::createClinicalModel(1.0, first);
        auto b = ClinicalModel::createClinicalModel(1.0, second);
        TS_ASSERT(dynamic_cast<DecisionTree5Day*>(a.get()));
        TS_ASSERT(dynamic_cast<ImmediateOutcomes*>(b.get()));
        TS_ASSERT_EQUALS(a->hsMemory(), 10);
        TS_ASSERT_EQUALS(b->hsMemory(), 25);
        first.healthSystemMemory = 15;
        TS_ASSERT_EQUALS(a->hsMemory(), 15);
        TS_ASSERT_EQUALS(b->hsMemory(), 25);

        scnXml::DecisionTree treated, untreated;
        treated.setTreatFailure(scnXml::DTTreatFailure());
        untreated.setNoTreatment(scnXml::DTNoTreatment());
        first.escm_uncomplicated = &CMDecisionTree::create(treated, true, definitions, first);
        first.escm_complicated = &CMDecisionTree::create(untreated, false, definitions, first);
        second.escm_uncomplicated = &CMDecisionTree::create(untreated, true, definitions, second);
        second.escm_complicated = &CMDecisionTree::create(treated, false, definitions, second);
        hd->pgState = Episode::SICK;
        TS_ASSERT(executeCaseManagement(*hd, first).treated);
        TS_ASSERT(!executeCaseManagement(*hd, second).treated);
        hd->pgState = static_cast<Episode::State>(Episode::SICK | Episode::COMPLICATED);
        TS_ASSERT(!executeCaseManagement(*hd, first).treated);
        TS_ASSERT(executeCaseManagement(*hd, second).treated);

        ClinicalParameters daily;
        daily.kind = ClinicalParameters::Kind::EventScheduler;
        auto c = ClinicalModel::createClinicalModel(1.0, daily);
        TS_ASSERT(dynamic_cast<ClinicalEventScheduler*>(c.get()));
    }

    void testDecisionLibrariesHaveIndependentLifetimes() {
        TS_ASSERT_EQUALS(CountedDecision::live, 0);
        {
            OM::ModelParameters second;
            const auto& survivor = save_decision(second.clinical, new CountedDecision());
            scnXml::DecisionTree input;
            input.setNoTreatment(scnXml::DTNoTreatment());
            const auto& secondTree = CMDecisionTree::create(input, true, definitions, second.clinical);
            {
                OM::ModelParameters first;
                const auto& independent = save_decision(first.clinical, new CountedDecision());
                TS_ASSERT_DIFFERS(&independent, &survivor);
                TS_ASSERT_EQUALS(CountedDecision::live, 2);
                const auto& firstTree = CMDecisionTree::create(input, true, definitions, first.clinical);
                TS_ASSERT_DIFFERS(&firstTree, &secondTree);
                TS_ASSERT_EQUALS(&firstTree,
                    &CMDecisionTree::create(input, true, definitions, first.clinical));
            }
            TS_ASSERT_EQUALS(CountedDecision::live, 1);
            TS_ASSERT(!survivor.exec(*hd).treated);
            TS_ASSERT(!secondTree.exec(*hd).treated);
        }
        TS_ASSERT_EQUALS(CountedDecision::live, 0);
    }

    void testDuplicateNodesAreDestroyed() {
        ClinicalParameters library;
        const auto& first = save_decision(library, new CountedDecision());
        for (int i = 0; i < 10; ++i)
            TS_ASSERT_EQUALS(&save_decision(library, new CountedDecision()), &first);
        TS_ASSERT_EQUALS(CountedDecision::live, 1);
    }

    void testLegacySimpleTreatmentInterningIgnoresBloodValues() {
        scnXml::DecisionTree first, second, different;
        first.getTreatSimple().push_back(scnXml::DTTreatSimple("35d", "5d"));
        second.getTreatSimple().push_back(scnXml::DTTreatSimple("35d", "10d"));
        different.getTreatSimple().push_back(scnXml::DTTreatSimple("40d", "5d"));
        const auto& saved = CMDecisionTree::create(first, true, definitions, parameters.clinical);
        const auto& reused = CMDecisionTree::create(second, true, definitions, parameters.clinical);
        TS_ASSERT_EQUALS(&saved, &reused);
        TS_ASSERT_DIFFERS(&saved, &CMDecisionTree::create(different, true, definitions, parameters.clinical));
        reused.exec(*hd);
        TS_ASSERT_EQUALS(whm->lastTimeBlood, 5);
    }

    void testLegacyCohortInterningIgnoresComponent() {
        using OM::interventions::ComponentId;
        definitions.identifierMap.insert_or_assign("comparison-first", ComponentId(101));
        definitions.identifierMap.insert_or_assign("comparison-second", ComponentId(102));
        scnXml::DecisionTree positive, negative, first, second;
        scnXml::DTReport report;
        report.setOutputNumber(10051);
        positive.getReport().push_back(report);
        report.setOutputNumber(10052);
        negative.getReport().push_back(report);
        first.setCohort(scnXml::DTCohort(positive, negative, "comparison-first"));
        second.setCohort(scnXml::DTCohort(positive, negative, "comparison-second"));
        TS_ASSERT_EQUALS(&CMDecisionTree::create(first, true, definitions, parameters.clinical),
                         &CMDecisionTree::create(second, true, definitions, parameters.clinical));
    }

    void testBranchKindsRemainDistinctAndPreserveOutcomes() {
        scnXml::DecisionTree treated, untreated, firstLine, severe, diagnostic, origin;
        treated.setTreatFailure(scnXml::DTTreatFailure());
        untreated.setNoTreatment(scnXml::DTNoTreatment());
        firstLine.setCaseType(scnXml::DTCaseType(treated, untreated));
        severe.setSevere(scnXml::DTSevere(treated, untreated));
        diagnostic.setDiagnostic(scnXml::DTDiagnostic(treated, untreated, "microscopy"));
        origin.setInfectionOrigin(scnXml::DTInfectionOrigin(untreated, untreated, treated));
        const auto& first = CMDecisionTree::create(firstLine, true, definitions, parameters.clinical);
        const auto& severity = CMDecisionTree::create(severe, true, definitions, parameters.clinical);
        TS_ASSERT_DIFFERS(&first, &severity);
        hd->pgState = Episode::SICK;
        TS_ASSERT(first.exec(*hd).treated);
        TS_ASSERT(!severity.exec(*hd).treated);
        hd->pgState = static_cast<Episode::State>(Episode::SICK | Episode::SECOND_CASE);
        TS_ASSERT(!first.exec(*hd).treated);
        hd->pgState = static_cast<Episode::State>(Episode::SICK | Episode::COMPLICATED);
        TS_ASSERT(severity.exec(*hd).treated);
        whm->totalDensity = 0.0;
        TS_ASSERT(CMDecisionTree::create(diagnostic, true, definitions, parameters.clinical).exec(*hd).screened);
        TS_ASSERT(CMDecisionTree::create(origin, true, definitions, parameters.clinical).exec(*hd).treated);
        TS_ASSERT_THROWS(CMDecisionTree::create(firstLine, false, definitions, parameters.clinical), const OM::util::xml_scenario_error&);
    }

    void testLegacyUncomplicatedInterningIgnoresMemory() {
        struct Restore {
            SimTime& value;
            SimTime old;
            ~Restore() { value = old; }
        } restore{parameters.clinical.healthSystemMemory, parameters.clinical.healthSystemMemory};
        parameters.clinical.healthSystemMemory = 30;
        scnXml::DecisionTree positive, negative, first, second;
        positive.setTreatFailure(scnXml::DTTreatFailure());
        negative.setNoTreatment(scnXml::DTNoTreatment());
        scnXml::DTUncomplicated branch(positive, negative);
        branch.setMemory("10d"); first.setUncomplicated(branch);
        branch.setMemory("20d"); second.setUncomplicated(branch);
        TS_ASSERT_EQUALS(&CMDecisionTree::create(first, true, definitions, parameters.clinical),
                         &CMDecisionTree::create(second, true, definitions, parameters.clinical));
        branch.setMemory("35d"); second.setUncomplicated(branch);
        TS_ASSERT_THROWS(CMDecisionTree::create(second, true, definitions, parameters.clinical), const OM::util::xml_scenario_error&);
    }

    void testReportOrderDistinguishesNodes() {
        scnXml::DecisionTree first, second;
        scnXml::DTReport a, b;
        a.setOutputNumber(10061); b.setOutputNumber(10062);
        first.getReport().push_back(a); first.getReport().push_back(b);
        second.getReport().push_back(b); second.getReport().push_back(a);
        const auto& saved = CMDecisionTree::create(first, true, definitions, parameters.clinical);
        TS_ASSERT_EQUALS(&saved, &CMDecisionTree::create(first, true, definitions, parameters.clinical));
        TS_ASSERT_DIFFERS(&saved, &CMDecisionTree::create(second, true, definitions, parameters.clinical));
    }

    void testRandomP () {
	// random decision
        // option (a) is to treat, option (b) is to do nothing
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        
        scnXml::Outcome o1r2( 0.9 ), o2r2( 0.1 );
        o1r2.getTreatPKPD().push_back( treat1 );
        o2r2.setNoTreatment( scnXml::DTNoTreatment() );
        scnXml::DTRandom r2;
        r2.getOutcome().push_back( o1r2 );
        r2.getOutcome().push_back( o2r2 );
        
        scnXml::Outcome o1r3( 0.7 ), o2r3( 0.3 );
        o1r3.getTreatPKPD().push_back( treat1 );
        o2r3.setNoTreatment( scnXml::DTNoTreatment() );
        scnXml::DTRandom r3;
        r3.getOutcome().push_back( o1r3 );
        r3.getOutcome().push_back( o2r3 );
        
        scnXml::Outcome o1r1( 0.5 ), o2r1( 0.5 );
        o1r1.setRandom( r2 );
        o2r1.setRandom( r3 );
        
        scnXml::DTRandom r1;
        r1.getOutcome().push_back( o1r1 );
        r1.getOutcome().push_back( o2r1 );
        
        scnXml::DecisionTree dt;
        dt.setRandom( r1 );
	
	const int N = 10000;
	const double LIM = .02;
	
        // test that dt.exec chooses to treat 80% and no action 20% of the time:
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt ), 0.8, LIM );
    }

    void testUC2Test () {
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back( treat1 );
        scnXml::DecisionTree noAction;
        noAction.setNoTreatment( scnXml::DTNoTreatment() );

        scnXml::DTCaseType ct( simpleTreat,  // first line: simple treatment
                            noAction );     // second line: no action
        scnXml::DecisionTree dt;
        dt.setCaseType( ct );

        using Underlying = std::underlying_type_t<Pathogenesis::State>;
        hd->pgState = static_cast<Episode::State>( static_cast<Underlying>(Pathogenesis::STATE_MALARIA) );
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );

        hd->pgState = static_cast<Episode::State>( static_cast<Underlying>(Pathogenesis::STATE_MALARIA) | Episode::SECOND_CASE );

        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 0 );
    }

    void testParasiteTest () {
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back( treat1 );
        scnXml::DecisionTree noAction;
        noAction.setNoTreatment( scnXml::DTNoTreatment() );
        
        scnXml::DTDiagnostic microscopy( simpleTreat,  // positive: simple treatment
                            noAction,     // negative: no action
                            "microscopy" );     // type of diagnostic
        scnXml::DecisionTree dt_mic;
        dt_mic.setDiagnostic( microscopy );
        
        scnXml::DTDiagnostic rdt( simpleTreat,  // positive: simple treatment
                            noAction,     // negative: no action
                            "RDT" );     // type of diagnostic
        scnXml::DecisionTree dt_rdt;
        dt_rdt.setDiagnostic( rdt );
        
        hd->pgState = static_cast<Episode::State>( Pathogenesis::STATE_MALARIA );
	const int N = 20000;
	const double LIM = .02;
	
	whm->totalDensity = 0.0;	// no parasites (so we test specificity)
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic ), 1 - 0.75, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt ), 1 - 0.942, LIM );
	
        whm->totalDensity = 80.0;	// a few parasites
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic ), 0.96815, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt ), 0.90599, LIM );
	
        whm->totalDensity = 2000.0;	// lots of parasites
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic ), 0.99867, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt ), 0.99584, LIM );
    }
    
    void testAgeSwitch(){
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DTAge ageSwitch;
        
        scnXml::Age treatYoung( 0.0 );  // ages 0 to 2.5
        treatYoung.getTreatPKPD().push_back( treat1 );
        ageSwitch.getAge().push_back( treatYoung );
        
        scnXml::Age noTreat( 2.5 );     // ages 2.5 to 50
        noTreat.setNoTreatment( scnXml::DTNoTreatment() );
        ageSwitch.getAge().push_back( noTreat );
        
        scnXml::Age treatOlder( 50.0 );  // ages 50+
        treatOlder.getTreatPKPD().push_back( treat1 );
        ageSwitch.getAge().push_back( treatOlder );
        
        scnXml::DecisionTree dt;
        dt.setAge( ageSwitch );
        
        hd->ageYears = 1;
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
        hd->ageYears = 2.5;
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 0 );
        hd->ageYears = 50;
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
        hd->ageYears = 1e6;     // there is no upper bound, so this should work!
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
    }
    
    void testSimpleTreat(){
        TS_ASSERT_EQUALS( whm->lastTimeLiver, sim::never() );
        TS_ASSERT_EQUALS( whm->lastTimeBlood, sim::never() );
        
        scnXml::DTTreatSimple treat1( "0t", "1t" );   // 0 time steps liver, 1 blood (using 1 day TS)
        scnXml::DecisionTree dt1;
        dt1.getTreatSimple().push_back( treat1 );
        
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt1 ), 1 );
        TS_ASSERT_EQUALS( whm->lastTimeLiver, 0 );
        TS_ASSERT_EQUALS( whm->lastTimeBlood, 1*5 );
        
        scnXml::DTTreatSimple treat2( "15d", "-1t" );
        scnXml::DecisionTree dt2;
        dt2.getTreatSimple().push_back( treat2 );
        
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt2 ), 1 );
        TS_ASSERT_EQUALS( whm->lastTimeLiver, 3*5 );
        TS_ASSERT_EQUALS( whm->lastTimeBlood, -1*5 );
    }
    
    double runAndGetMgPrescribed( scnXml::DecisionTree& dt, double age ){
        hd->ageYears = age;
        UnittestUtil::clearMedicateQueue( whm->pkpd );
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
        return UnittestUtil::getPrescribedMg( whm->pkpd );
    }
    
    void testDosing(){
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree dt1;
        dt1.getTreatPKPD().push_back( treat1 );
        
        // Test our dosing table. Set with a multiplier of 1 below 5 and 5 from 5.
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt1, 0 ), 6, 1e-8 );
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt1, 4.9 ), 6, 1e-8 );
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt1, 5 ), 30, 1e-8 );
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt1, 99 ), 30, 1e-8 );
        
        scnXml::DTTreatPKPD treat2( "sched2", "dosage1" );
        scnXml::DecisionTree dt2;
        dt2.getTreatPKPD().push_back( treat2 );
        
        // Test our dosing table. Set with a multiplier of 1 below 5 and 5 from 5.
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt2, 0 ), 7, 1e-8 );
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt2, 99 ), 35, 1e-8 );
    }
    
private:
    unique_ptr<Host::Human> human;
    WHMock* whm;
    unique_ptr<CMHostData> hd;
};

#endif
