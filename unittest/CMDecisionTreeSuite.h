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

#pragma once

#include <cxxtest/TestSuite.h>
#include "Clinical/CMDecisionTree.h"
#include "util/random.h"
#include "UnittestUtil.h"
#include "WHMock.h"
#include <limits>

using UnitTest::WHMockState;

class CMDecisionTreeSuite : public CxxTest::TestSuite {
  public:
    CMDecisionTreeSuite() {}
    ~CMDecisionTreeSuite() {}

    void setUp() {
        // Note: cannot create whm in constructor, since it uses random number
        // generator which is initialized after constructor runs.
        UnittestUtil::initTime(5);
        UnittestUtil::setDiagnostics();

        UnittestUtil::EmpiricalWHM_setup();

        human.reset(UnittestUtil::createHuman(UnittestUtil::clock().nowOrTs0()).release());

        ETS_ASSERT(human.get() != 0);
        wh = UnittestUtil::setHumanWH(*human, UnitTest::makeWHMock(UnittestUtil::ctx()));
        whm = &UnitTest::whMockState(*wh);
        ETS_ASSERT(whm != 0);

        ageYears = 21; // any real age will do for most tests
        pgState = NONE;

        // could seed random-number-generator, but shouldn't affect outcomes
        UnittestUtil::PkPdSuiteSetup();
    }
    void tearDown() {
        human.reset();
        wh = nullptr;
        whm = nullptr;
        UnittestUtil::ctx().drugs.clear();
        UnittestUtil::ctx().treatments.clear();
    }

    /* Runs the decision tree N times
     * returns the proportion of these runs where the output was any treatment
     */
    double propTreatmentsNReps(int N, const scnXml::DecisionTree& dt) {
        const FlatTree& cmdt = compile_decision_tree(UnittestUtil::ctx(), dt, true);

        whm->nTreatments = 0;
        int secondCounter = 0;
        for (int i = 0; i < N; ++i) {
            DTExecCtx ec(UnittestUtil::ctx(), *human, ageYears, pgState);
            secondCounter += (exec_decision_tree(cmdt, ec).treated ? 1 : 0);
        }
        TS_ASSERT_EQUALS(whm->nTreatments, secondCounter);
        return double(whm->nTreatments) / double(N);
    }

    void testRandomP() {
        // random decision
        // option (a) is to treat, option (b) is to do nothing
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");

        scnXml::Outcome o1r2(0.9), o2r2(0.1);
        o1r2.getTreatPKPD().push_back(treat1);
        o2r2.setNoTreatment(scnXml::DTNoTreatment());
        scnXml::DTRandom r2;
        r2.getOutcome().push_back(o1r2);
        r2.getOutcome().push_back(o2r2);

        scnXml::Outcome o1r3(0.7), o2r3(0.3);
        o1r3.getTreatPKPD().push_back(treat1);
        o2r3.setNoTreatment(scnXml::DTNoTreatment());
        scnXml::DTRandom r3;
        r3.getOutcome().push_back(o1r3);
        r3.getOutcome().push_back(o2r3);

        scnXml::Outcome o1r1(0.5), o2r1(0.5);
        o1r1.setRandom(r2);
        o2r1.setRandom(r3);

        scnXml::DTRandom r1;
        r1.getOutcome().push_back(o1r1);
        r1.getOutcome().push_back(o2r1);

        scnXml::DecisionTree dt;
        dt.setRandom(r1);

        const int N = 10000;
        const double LIM = .02;

        // test that dt.exec chooses to treat 80% and no action 20% of the time:
        TS_ASSERT_DELTA(propTreatmentsNReps(N, dt), 0.8, LIM);
    }

    void testUC2Test() {
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back(treat1);
        scnXml::DecisionTree noAction;
        noAction.setNoTreatment(scnXml::DTNoTreatment());

        scnXml::DTCaseType ct(simpleTreat, // first line: simple treatment
                              noAction);   // second line: no action
        scnXml::DecisionTree dt;
        dt.setCaseType(ct);

        using Underlying = std::underlying_type_t<State>;
        pgState = static_cast<Episode::State>(static_cast<Underlying>(STATE_MALARIA));
        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt), 1);

        pgState = static_cast<Episode::State>(static_cast<Underlying>(STATE_MALARIA) | SECOND_CASE);

        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt), 0);
    }

    void testParasiteTest() {
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back(treat1);
        scnXml::DecisionTree noAction;
        noAction.setNoTreatment(scnXml::DTNoTreatment());

        scnXml::DTDiagnostic microscopy(simpleTreat,   // positive: simple treatment
                                        noAction,      // negative: no action
                                        "microscopy"); // type of diagnostic
        scnXml::DecisionTree dt_mic;
        dt_mic.setDiagnostic(microscopy);

        scnXml::DTDiagnostic rdt(simpleTreat, // positive: simple treatment
                                 noAction,    // negative: no action
                                 "RDT");      // type of diagnostic
        scnXml::DecisionTree dt_rdt;
        dt_rdt.setDiagnostic(rdt);

        pgState = static_cast<Episode::State>(STATE_MALARIA);
        const int N = 20000;
        const double LIM = .02;

        whm->totalDensity = 0.0; // no parasites (so we test specificity)
        TS_ASSERT_DELTA(propTreatmentsNReps(N, dt_mic), 1 - 0.75, LIM);
        TS_ASSERT_DELTA(propTreatmentsNReps(N, dt_rdt), 1 - 0.942, LIM);

        whm->totalDensity = 80.0; // a few parasites
        TS_ASSERT_DELTA(propTreatmentsNReps(N, dt_mic), 0.96815, LIM);
        TS_ASSERT_DELTA(propTreatmentsNReps(N, dt_rdt), 0.90599, LIM);

        whm->totalDensity = 2000.0; // lots of parasites
        TS_ASSERT_DELTA(propTreatmentsNReps(N, dt_mic), 0.99867, LIM);
        TS_ASSERT_DELTA(propTreatmentsNReps(N, dt_rdt), 0.99584, LIM);
    }

    void testAgeSwitch() {
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");
        scnXml::DTAge ageSwitch;

        scnXml::Age treatYoung(0.0); // ages 0 to 2.5
        treatYoung.getTreatPKPD().push_back(treat1);
        ageSwitch.getAge().push_back(treatYoung);

        scnXml::Age noTreat(2.5); // ages 2.5 to 50
        noTreat.setNoTreatment(scnXml::DTNoTreatment());
        ageSwitch.getAge().push_back(noTreat);

        scnXml::Age treatOlder(50.0); // ages 50+
        treatOlder.getTreatPKPD().push_back(treat1);
        ageSwitch.getAge().push_back(treatOlder);

        scnXml::DecisionTree dt;
        dt.setAge(ageSwitch);

        ageYears = 1;
        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt), 1);
        ageYears = 2.5;
        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt), 0);
        ageYears = 50;
        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt), 1);
        ageYears = 1e6; // there is no upper bound, so this should work!
        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt), 1);
    }

    void testSimpleTreat() {
        TS_ASSERT_EQUALS(whm->lastTimeLiver, never());
        TS_ASSERT_EQUALS(whm->lastTimeBlood, never());

        scnXml::DTTreatSimple treat1("0t", "1t"); // 0 time steps liver, 1 blood (using 1 day TS)
        scnXml::DecisionTree dt1;
        dt1.getTreatSimple().push_back(treat1);

        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt1), 1);
        TS_ASSERT_EQUALS(whm->lastTimeLiver, 0);
        TS_ASSERT_EQUALS(whm->lastTimeBlood, 1 * 5);

        scnXml::DTTreatSimple treat2("15d", "-1t");
        scnXml::DecisionTree dt2;
        dt2.getTreatSimple().push_back(treat2);

        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt2), 1);
        TS_ASSERT_EQUALS(whm->lastTimeLiver, 3 * 5);
        TS_ASSERT_EQUALS(whm->lastTimeBlood, -1 * 5);
    }

    double runAndGetMgPrescribed(scnXml::DecisionTree& dt, double age) {
        ageYears = age;
        UnittestUtil::clearMedicateQueue(whm->pkpd);
        TS_ASSERT_EQUALS(propTreatmentsNReps(1, dt), 1);
        return UnittestUtil::getPrescribedMg(whm->pkpd);
    }

    // ———  boundary/compile-canonicalization goldens (phase 4a)  ———
    // Golden values below were captured from the pre-rewrite polymorphic
    // decision-tree code (p3-green). They pin: age-switch exact-boundary
    // routing, cumulative-key ties (p=0 duplicate keys are dropped by
    // map-insert: the earlier child wins), leaf-run fold seams (adjacent
    // outcomes with identical subtrees must select identically), and the
    // CMDTMultiple screened-swallowing quirk. Golden == today's behavior.

    /** Run the tree once; return (treated, screened) with treatment-marker
     * observation via WHMock simple-treat timestamps. */
    CMDTOut execOnce(const scnXml::DecisionTree& dt) {
        const FlatTree& cmdt = compile_decision_tree(UnittestUtil::ctx(), dt, true);
        whm->lastTimeLiver = whm->lastTimeBlood = never();
        DTExecCtx ec(UnittestUtil::ctx(), *human, ageYears, pgState);
        return exec_decision_tree(cmdt, ec);
    }

    // treatSimple markers distinguish branch identity without RNG prediction:
    static constexpr SimTime MARK_A = 0;     // "0t"/"1t" → liver 0
    static constexpr SimTime MARK_B = 1 * 5; // "1t"/"1t" → liver 5 (5-day TS)

    void testAgeFoldSeams() {
        // age groups [0,2.5)→treat, [2.5,50)→treat (identical subtree),
        // [50,∞)→none: adjacent identical leaves are a fold seam.
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");
        scnXml::DTAge ageSwitch;
        scnXml::Age g0(0.0);
        g0.getTreatPKPD().push_back(treat1);
        scnXml::Age g1(2.5);
        g1.getTreatPKPD().push_back(treat1);
        scnXml::Age g2(50.0);
        g2.setNoTreatment(scnXml::DTNoTreatment());
        ageSwitch.getAge().push_back(g0);
        ageSwitch.getAge().push_back(g1);
        ageSwitch.getAge().push_back(g2);
        scnXml::DecisionTree dt;
        dt.setAge(ageSwitch);

        // golden (pre-rewrite): key order = {2.5→g0, 50→g1, +inf→g2},
        // upper_bound strictly-greater selection at each boundary:
        const double ages[] = {2.4, 2.5, 3.0, 50.0, 1e6};
        const int golden[] = {1, 1, 1, 0, 0};
        for (size_t i = 0; i < 5; ++i) {
            ageYears = ages[i];
            TS_ASSERT_EQUALS(execOnce(dt).treated, golden[i] != 0);
        }
        ageYears = 21;
    }

    void testRandomTieDrop() {
        // p=0 middle outcome creates a cumulative-key tie with the first
        // outcome (0.5). map-insert drops the duplicate key: the SECOND
        // outcome's subtree (marker B) must never be selected.
        scnXml::DTTreatSimple a("0t", "1t"), b("1t", "1t");
        scnXml::Outcome oA(0.5), oB(0.0), oC(0.5);
        oA.getTreatSimple().push_back(a);
        oB.getTreatSimple().push_back(b);
        oC.setNoTreatment(scnXml::DTNoTreatment());
        scnXml::DTRandom r;
        r.getOutcome().push_back(oA);
        r.getOutcome().push_back(oB);
        r.getOutcome().push_back(oC);
        scnXml::DecisionTree dt;
        dt.setRandom(r);

        const int N = 20000;
        int nTreated = 0, nMarkB = 0, nMarkA = 0;
        for (int i = 0; i < N; ++i) {
            CMDTOut out = execOnce(dt);
            nTreated += out.treated ? 1 : 0;
            nMarkA += (whm->lastTimeLiver == MARK_A) ? 1 : 0;
            nMarkB +=
                (whm->lastTimeLiver == MARK_B && whm->lastTimeBlood == 1 * 5 && whm->lastTimeLiver != MARK_A) ? 1 : 0;
        }
        // golden (pre-rewrite, DTGOLDEN capture): dropped key never fires
        TS_ASSERT_EQUALS(nMarkB, 0);
        TS_ASSERT_EQUALS(nTreated, nMarkA);
        TS_ASSERT_DELTA(double(nTreated) / N, 0.5, 0.02);
    }

    void testRandomAllFold() {
        // All outcomes share one identical subtree: keys {0.5,1.0} fold to
        // {1.0}; every draw must treat (exact, no statistics).
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");
        scnXml::Outcome o1(0.5), o2(0.5);
        o1.getTreatPKPD().push_back(treat1);
        o2.getTreatPKPD().push_back(treat1);
        scnXml::DTRandom r;
        r.getOutcome().push_back(o1);
        r.getOutcome().push_back(o2);
        scnXml::DecisionTree dt;
        dt.setRandom(r);

        const int N = 1000;
        int nTreated = 0;
        for (int i = 0; i < N; ++i)
            nTreated += execOnce(dt).treated ? 1 : 0;
        TS_ASSERT_EQUALS(nTreated, 1000); // golden: uniform_01 < 0.5 or not, target is identical
    }

    void testRandomZeroPFirst() {
        // p=0 first outcome: key 0.0 is never upper-bounded above by
        // uniform_01() in [0,1) — actually key 0.0 is skipped even for
        // x==0.0 (strictly-greater semantics) — never selected.
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");
        scnXml::Outcome oT(0.0), oN(1.0);
        oT.getTreatPKPD().push_back(treat1);
        oN.setNoTreatment(scnXml::DTNoTreatment());
        scnXml::DTRandom r;
        r.getOutcome().push_back(oT);
        r.getOutcome().push_back(oN);
        scnXml::DecisionTree dt;
        dt.setRandom(r);

        const int N = 1000;
        int nTreated = 0;
        for (int i = 0; i < N; ++i)
            nTreated += execOnce(dt).treated ? 1 : 0;
        TS_ASSERT_EQUALS(nTreated, 0); // golden (pre-rewrite)
    }

    void testScreenedSwallowed() {
        // CMDTMultiple quirk: treated ORs out of children, screened does
        // not. Diagnostic forces screened=true within its own subtree.
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back(treat1);

        scnXml::DTDiagnostic diag(simpleTreat, simpleTreat, "microscopy");
        scnXml::DecisionTree diagTree;
        diagTree.setDiagnostic(diag);

        scnXml::DTMultiple m;
        m.getDiagnostic().push_back(diag);
        scnXml::DecisionTree multi;
        multi.setMultiple(m);

        pgState = static_cast<Episode::State>(STATE_MALARIA);
        // Diagnostic alone: screened visible.
        CMDTOut out = execOnce(diagTree);
        TS_ASSERT_EQUALS(out.treated, true);
        TS_ASSERT_EQUALS(out.screened, true);
        // Inside Multiple: screened swallowed at the Multiple boundary,
        // treated survives (from either child).
        out = execOnce(multi);
        TS_ASSERT_EQUALS(out.treated, true);
        TS_ASSERT_EQUALS(out.screened, false); // golden (pre-rewrite quirk)
    }

    void testMockStateIsPerInstance() {
        auto second = UnitTest::makeWHMock(UnittestUtil::ctx());
        WHMockState& secondState = UnitTest::whMockState(*second);

        human->withinHostModel->treatment(*human, TreatmentId(0));

        TS_ASSERT_EQUALS(whm->nTreatments, 1);
        TS_ASSERT_EQUALS(secondState.nTreatments, 0);
    }

    void testDosing() {
        scnXml::DTTreatPKPD treat1("sched1", "dosage1");
        scnXml::DecisionTree dt1;
        dt1.getTreatPKPD().push_back(treat1);

        // Test our dosing table. Set with a multiplier of 1 below 5 and 5 from 5.
        TS_ASSERT_DELTA(runAndGetMgPrescribed(dt1, 0), 6, 1e-8);
        TS_ASSERT_DELTA(runAndGetMgPrescribed(dt1, 4.9), 6, 1e-8);
        TS_ASSERT_DELTA(runAndGetMgPrescribed(dt1, 5), 30, 1e-8);
        TS_ASSERT_DELTA(runAndGetMgPrescribed(dt1, 99), 30, 1e-8);

        scnXml::DTTreatPKPD treat2("sched2", "dosage1");
        scnXml::DecisionTree dt2;
        dt2.getTreatPKPD().push_back(treat2);

        // Test our dosing table. Set with a multiplier of 1 below 5 and 5 from 5.
        TS_ASSERT_DELTA(runAndGetMgPrescribed(dt2, 0), 7, 1e-8);
        TS_ASSERT_DELTA(runAndGetMgPrescribed(dt2, 99), 35, 1e-8);
    }

  private:
    unique_ptr<Human> human;
    WHInterface* wh;
    WHMockState* whm;
    double ageYears;
    Episode::State pgState;
};
