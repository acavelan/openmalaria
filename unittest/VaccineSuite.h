#pragma once
#include <cxxtest/TestSuite.h>
#include "interventions/Vaccine.h"

class VaccineSuite : public CxxTest::TestSuite {
    using Effect = OM::interventions::VaccineEffect;
    using Component = OM::interventions::VaccineComponent;
    using Rng = OM::util::LocalRng;

    static std::string state(std::vector<Effect>& effects, Rng& rng) {
        std::ostringstream stream;
        OM::interventions::checkpointVaccines(effects, stream);
        rng.checkpoint(stream);
        return stream.str();
    }
public:
    void testInterleavedDefinitionsWithSameIdStayIndependent() {
        using namespace OM::interventions;
        Component first(ComponentId(0), Vaccine::PEV, true);
        Component second(ComponentId(0), Vaccine::PEV, true);
        first.perGenotypeInitialMeanEfficacy = {{0.3, 0.4}};
        second.perGenotypeInitialMeanEfficacy = {{0.7, 0.8}};
        first.decayFunc = OM::util::DecayFunction(0.2);
        second.decayFunc = OM::util::DecayFunction(0.4);
        std::vector<Effect> a, b, isolatedA, isolatedB;
        Rng rngA(17, 0), rngB(23, 0), referenceA(17, 0), referenceB(23, 0);
        for (int dose = 0; dose < 3; ++dose) {
            vaccinate(a, rngA, first, {}, dose * 10, 1);
            vaccinate(b, rngB, second, {}, dose * 20, 1);
        }
        for (int dose = 0; dose < 3; ++dose) vaccinate(isolatedA, referenceA, first, {}, dose * 10, 1);
        for (int dose = 0; dose < 3; ++dose) vaccinate(isolatedB, referenceB, second, {}, dose * 20, 1);
        TS_ASSERT_EQUALS(state(a, rngA), state(isolatedA, referenceA));
        TS_ASSERT_EQUALS(state(b, rngB), state(isolatedB, referenceB));
        TS_ASSERT_EQUALS(a[0].parameters, &first);
        TS_ASSERT_EQUALS(b[0].parameters, &second);
        TS_ASSERT(first.reportDeployments && second.reportDeployments);
    }

    void testDoseLimitsAndRedeploymentKeepRandomSequence() {
        using namespace OM::interventions;
        Component parameters(ComponentId(0), Vaccine::PEV);
        parameters.perGenotypeInitialMeanEfficacy = {{1.0}};
        parameters.decayFunc = OM::util::DecayFunction(0.3);
        std::vector<Effect> effects;
        Rng rng(17, 0), expected(17, 0);
        VaccineLimits limits;
        limits.minPrevDoses = 1;
        TS_ASSERT(!vaccinate(effects, rng, parameters, limits, 10, 1));
        TS_ASSERT(effects.empty());
        limits.minPrevDoses = 0;
        limits.maxCumDoses = 1;
        TS_ASSERT(vaccinate(effects, rng, parameters, limits, 10, 1));
        double factor = parameters.decayFunc.heterogeneity.sample(expected);
        TS_ASSERT_EQUALS(*effects[0].hetSample, factor);
        TS_ASSERT(!vaccinate(effects, rng, parameters, limits, 20, 1));
        limits.maxCumDoses = 2;
        TS_ASSERT(vaccinate(effects, rng, parameters, limits, 30, 1));
        TS_ASSERT_EQUALS(*effects[0].hetSample, factor);
        TS_ASSERT_EQUALS(effects[0].numDosesAdministered, 2u);
        TS_ASSERT_EQUALS(effects[0].timeLastDeployment, 30);
        std::ostringstream actualState, expectedState;
        rng.checkpoint(actualState);
        expected.checkpoint(expectedState);
        TS_ASSERT_EQUALS(actualState.str(), expectedState.str());
    }

    void testCheckpointRebindsToSuppliedDefinitions() {
        using namespace OM::interventions;
        Component original(ComponentId(0), Vaccine::PEV);
        original.perGenotypeInitialMeanEfficacy = {{1.0}};
        std::vector<Effect> effects;
        Rng rng(17, 0);
        vaccinate(effects, rng, original, {}, 10, 1);
        std::ostringstream saved;
        checkpointVaccines(effects, saved);
        std::vector<std::unique_ptr<HumanInterventionComponent>> definitions;
        auto restoredDefinition = std::make_unique<Component>(ComponentId(0), Vaccine::PEV);
        restoredDefinition->decayFunc.kind = OM::util::DecayKind::Step;
        restoredDefinition->decayFunc.rate = 0.1;
        definitions.push_back(std::move(restoredDefinition));
        std::vector<Effect> restored;
        std::istringstream input(saved.str());
        checkpointVaccines(restored, input, definitions);
        TS_ASSERT_EQUALS(restored[0].parameters, definitions[0].get());
        TS_ASSERT_EQUALS(vaccineFactor(effects, Vaccine::PEV, 20), 0.0);
        TS_ASSERT_EQUALS(vaccineFactor(restored, Vaccine::PEV, 20), 1.0);
        std::ostringstream resaved;
        checkpointVaccines(restored, resaved);
        TS_ASSERT_EQUALS(saved.str(), resaved.str());
        definitions.clear();
        std::istringstream missing(saved.str());
        TS_ASSERT_THROWS(checkpointVaccines(restored, missing, definitions), const OM::util::checkpoint_error&);
        definitions.resize(1);
        std::istringstream nullDefinition(saved.str());
        TS_ASSERT_THROWS(checkpointVaccines(restored, nullDefinition, definitions), const OM::util::checkpoint_error&);
    }

    void testFactorUsesExplicitDateTypeAndGenotype() {
        using namespace OM::interventions;
        Component first(ComponentId(0), Vaccine::PEV), second(ComponentId(1), Vaccine::TBV);
        first.decayFunc.kind = OM::util::DecayKind::Linear;
        first.decayFunc.rate = 0.01;
        Effect a, b;
        a.parameters = &first; a.timeLastDeployment = 10; a.hetSample = 1.0;
        a.perGenotypeInitialEfficacy = {1.0, 0.5};
        b.parameters = &second; b.timeLastDeployment = 10; b.hetSample = 1.0;
        b.perGenotypeInitialEfficacy = {0.25, 0.5};
        std::vector<Effect> effects = {a, b};
        TS_ASSERT_EQUALS(vaccineFactor(effects, Vaccine::PEV, 60, 0), 0.5);
        TS_ASSERT_EQUALS(vaccineFactor(effects, Vaccine::PEV, 60, 1), 0.75);
        TS_ASSERT_EQUALS(vaccineFactor(effects, Vaccine::PEV, 110, 0), 1.0);
        TS_ASSERT_EQUALS(vaccineFactor(effects, Vaccine::TBV, 60, 0), 0.75);
        TS_ASSERT_EQUALS(vaccineFactor(effects, Vaccine::BSV, 60, 0), 1.0);
    }
};
