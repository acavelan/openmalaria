#pragma once
#include <cxxtest/TestSuite.h>
#include "PopulationAgeStructure.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ClinicalModel.h"
#include <future>
#include "mon/Continuous.h"

class PopulationStateSuite : public CxxTest::TestSuite {
    static std::vector<OM::DemographyAgeGroup> groups() {
        return {{1, 3.47}, {5, 12.76}, {10, 14.52}, {15, 12.75},
            {20, 10.84}, {25, 8.39}, {30, 7.00}, {35, 5.80},
            {40, 5.10}, {45, 4.18}, {50, 3.34}, {55, 2.99},
            {60, 2.56}, {65, 2.33}, {70, 1.77}, {75, 1.01},
            {80, 0.74}, {85, 0.27}, {90, 0.16}};
    }
public:
    void testContinuousReportingRegistriesHaveIndependentLifetimes() {
        using namespace OM;
        mon::Continuous::State second;
        int a = 2, b = 9;
        mon::Continuous::registerCallback(second, "count", "second",
            [&](ostream& out) { out << b++; });
        Population population;
        std::ostringstream output;
        {
            mon::Continuous::State first;
            mon::Continuous::registerCallback(first, "count", "first",
                [&](ostream& out) { out << a++; });
            mon::Continuous::registerCallback(first, "population", "size",
                [](Population& p, ostream& out) { out << p.targetSize; });
            mon::Continuous::registerCallback(first, "humans", "humans",
                [](const vector<Host::Human>& humans, ostream& out) { out << humans.size(); });
            first.registered.at("count").write(population, output);
            first.registered.at("population").write(population, output);
            first.registered.at("humans").write(population, output);
            TS_ASSERT_EQUALS(second.registered.size(), 1u);
            TS_ASSERT_EQUALS(b, 9);
        }
        second.registered.at("count").write(population, output);
        TS_ASSERT_EQUALS(output.str(), "2009");
        TS_ASSERT_EQUALS(a, 3);
        TS_ASSERT_EQUALS(b, 10);
    }

    void testInfantReportingHasIndependentCountersAndParameters() {
        OM::Clinical::InfantMortality first{{0, 0}, {0, 0}, 4.0};
        OM::Clinical::InfantMortality second{{0, 0, 0}, {0, 0, 0}, 9.0};
        first.reportRisk(0, true);
        second.reportRisk(0, false);
        first.reportRisk(0, false);
        second.reportRisk(1, false);
        first.reportRisk(1, false);
        second.reportRisk(2, false);
        TS_ASSERT_EQUALS(first.allCause(), 504.0);
        TS_ASSERT_EQUALS(second.allCause(), 9.0);
        first.reset();
        TS_ASSERT_EQUALS(first.deaths.size(), 2u);
        TS_ASSERT_EQUALS(first.deaths[0], 0);
        TS_ASSERT_EQUALS(first.intervalsAtRisk[0], 0);
        TS_ASSERT_EQUALS(first.nonMalariaMortality, 4.0);
        TS_ASSERT(std::isnan(first.allCause())); // Retain zero-exposure behavior.
        TS_ASSERT_EQUALS(second.allCause(), 9.0);
    }

    void testConcurrentFitsHaveIndependentStateAndTimeInputs() {
        auto a = groups(), b = groups();
        b[0].percentage = 6.0;
        b[10].percentage = 1.0;
        auto referenceA = OM::makeAgeStructure(a, 80 * 365, 5);
        auto referenceB = OM::makeAgeStructure(b, 70 * 365, 1);
        auto first = std::async(std::launch::async, [&] { return OM::makeAgeStructure(a, 80 * 365, 5); });
        auto second = std::async(std::launch::async, [&] { return OM::makeAgeStructure(b, 70 * 365, 1); });
        TS_ASSERT(first.get() == referenceA);
        TS_ASSERT(second.get() == referenceB);
        TS_ASSERT_EQUALS(referenceA.size(), 5841u);
        TS_ASSERT_EQUALS(referenceB.size(), 25551u);
        TS_ASSERT_EQUALS(referenceA.front(), 0.0);
        TS_ASSERT_EQUALS(referenceA.back(), 1.0);
        TS_ASSERT(std::is_sorted(referenceA.begin(), referenceA.end()));
    }

    void testPopulationTargetsKeepRoundingAndAgeOrder() {
        std::vector<double> proportions = {0.0, 0.125, 0.5, 1.0};
        TS_ASSERT_EQUALS(OM::targetCumPop(proportions, 0, 100), 100);
        TS_ASSERT_EQUALS(OM::targetCumPop(proportions, 1, 100), 50);
        TS_ASSERT_EQUALS(OM::targetCumPop(proportions, 2, 100), 13);
        TS_ASSERT_EQUALS(OM::targetCumPop(proportions, 3, 100), 0);
        TS_ASSERT_THROWS(OM::makeAgeStructure({}, 365, 5), const OM::util::xml_scenario_error&);
        TS_ASSERT_THROWS(OM::makeAgeStructure(groups(), 365, 0), const OM::util::xml_scenario_error&);
    }

    void testMaternalHistoryIsOwnedAndRetainsFiveMonths() {
        using namespace OM::Host;
        auto daily = makeNeonatalMortality(1), fiveDay = makeNeonatalMortality(5);
        TS_ASSERT_EQUALS(daily.prevByGestationalAge.size(), 150u);
        TS_ASSERT_EQUALS(fiveDay.prevByGestationalAge.size(), 30u);
        updateNeonatalRisk(daily, 0.8, 0);
        TS_ASSERT_EQUALS(fiveDay.riskFromMaternalInfection, 0.0);
        const double highRisk = daily.riskFromMaternalInfection;
        for (size_t i = 1; i < daily.prevByGestationalAge.size(); ++i)
            updateNeonatalRisk(daily, 0.1, i);
        TS_ASSERT_EQUALS(daily.riskFromMaternalInfection, highRisk);
        updateNeonatalRisk(daily, 0.1, 0);
        TS_ASSERT(daily.riskFromMaternalInfection < highRisk);
        double expected = (0.3 * 0.011) * (1.0 - exp((0.1 / (0.19 + 0.1)) * (-1.0 / 0.25)));
        TS_ASSERT_EQUALS(daily.riskFromMaternalInfection, expected);
        TS_ASSERT_THROWS(makeNeonatalMortality(0), const OM::util::xml_scenario_error&);
        TS_ASSERT_THROWS(makeNeonatalMortality(151), const OM::util::xml_scenario_error&);
    }

    void testMaternalHistoryResumesExactly() {
        using namespace OM::Host;
        auto original = makeNeonatalMortality(5);
        updateNeonatalRisk(original, 0.7, 0);
        std::ostringstream saved;
        original & saved;
        NeonatalMortality restored;
        std::istringstream input(saved.str());
        restored & input;
        for (int step = 1; step < 50; ++step) {
            updateNeonatalRisk(original, (step % 5) * 0.1, step % 30);
            updateNeonatalRisk(restored, (step % 5) * 0.1, step % 30);
        }
        std::ostringstream a, b;
        original & a;
        restored & b;
        TS_ASSERT_EQUALS(a.str(), b.str());
    }
};
