#pragma once

#include <cxxtest/TestSuite.h>
#include "interventions/Interventions.h"
#include "Transmission/NonVectorModel.h"
#include "Population.h"
#include "UnittestUtil.h"

class InterventionsSuite : public CxxTest::TestSuite {
public:
    void testEmptyAndExhaustedSchedules() {
        using namespace OM;
        UnittestUtil::initTime(5);
        sim::s_interv = 0;
        sim::s_end = sim::startDate() + 100;
        Transmission::HostTransmissionParameters parameters;
        parameters.availabilityByAge.set({util::AgeInterpolation::Constant, {{0, 1.0}}}, 90, "test");
        Transmission::NonVectorModel transmission(vector<double>(sim::stepsPerYear(), 0.0),
            Transmission::dynamicEIR, vector<double>(365, 1.0), 10, parameters);
        transmission.simulationMode = Transmission::dynamicEIR;
        Population population;
        Clinical::ClinicalParameters library;
        interventions::Interventions empty;
        empty.init(scnXml::Interventions("empty"), population, transmission, library);
        empty.loadFromCheckpoint(population.humans, transmission);
        empty.deploy(population.humans, transmission);
        TS_ASSERT_EQUALS(empty.nextTimed, 0u);
        TS_ASSERT(empty.timed.empty());

        scnXml::UninfectVectors clear;
        clear.getTimedDeployment().push_back(scnXml::TimedDeployment2("1t"));
        scnXml::Interventions input("test");
        input.setUninfectVectors(clear);
        interventions::Interventions scheduled;
        scheduled.init(input, population, transmission, library);
        transmission.laggedKappa[0] = 0.5;
        scheduled.deploy(population.humans, transmission);
        TS_ASSERT_EQUALS(transmission.laggedKappa[0], 0.5);
        sim::s_interv = 5;
        scheduled.deploy(population.humans, transmission);
        TS_ASSERT_EQUALS(scheduled.nextTimed, 1u);
        TS_ASSERT_EQUALS(transmission.laggedKappa[0], 0.0);
        transmission.laggedKappa[0] = 0.75;
        scheduled.deploy(population.humans, transmission);
        TS_ASSERT_EQUALS(transmission.laggedKappa[0], 0.75);

        interventions::Interventions resumed;
        resumed.init(input, population, transmission, library);
        sim::s_interv = 10;
        resumed.loadFromCheckpoint(population.humans, transmission);
        resumed.deploy(population.humans, transmission);
        TS_ASSERT_EQUALS(resumed.nextTimed, 1u);
        TS_ASSERT_EQUALS(transmission.laggedKappa[0], 0.75);
    }

    void testDefinitionsAndRemovalRulesBelongToTheirModel() {
        using namespace OM::interventions;
        Interventions first, second;
        first.identifierMap.emplace("net", ComponentId(0));
        second.identifierMap.emplace("net", ComponentId(7));
        first.removalRules[SubPopRemove::ON_FIRST_BOUT].push_back(ComponentId(0));

        TS_ASSERT_EQUALS(first.getComponentId("net").id, 0u);
        TS_ASSERT_EQUALS(second.getComponentId("net").id, 7u);
        TS_ASSERT(second.removalRules[SubPopRemove::ON_FIRST_BOUT].empty());
        TS_ASSERT_THROWS(second.getComponentId("missing"), const OM::util::xml_scenario_error&);
    }
};
