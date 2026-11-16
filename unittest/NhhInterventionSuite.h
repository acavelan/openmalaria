#pragma once
#include <cxxtest/TestSuite.h>
#include "Transmission/Anopheles/AnophelesModel.h"
#include "ScenarioXml.h"
#include "schema/interventions.h"
#include <sstream>

class NhhInterventionSuite : public CxxTest::TestSuite {
public:
    void testAnimalDeploymentUsesSpeciesOrderAndValidatesDescriptions() {
        using scnXml::DoubleValue;
        scnXml::Description3 description;
        description.getAnopheles().push_back(scnXml::NonHumanHostsVectorSpecies(
            DoubleValue(2.0), DoubleValue(0.3), DoubleValue(0.4), DoubleValue(0.5), DoubleValue(0.6), "second"));
        description.getAnopheles().push_back(scnXml::NonHumanHostsVectorSpecies(
            DoubleValue(3.0), DoubleValue(0.7), DoubleValue(0.8), DoubleValue(0.9), DoubleValue(1.0), "first"));
        scnXml::NonHumanHosts2 input(description, "cattle");
        const std::map<std::string, size_t> indices{{"first", 0}, {"second", 1}};
        auto values = OM::xml::readNhhDeployment(input, indices);
        TS_ASSERT_EQUALS(values.size(), 2u);
        TS_ASSERT_EQUALS(values[0].relativeAvailabilityHuman, 3.0);
        TS_ASSERT_EQUALS(values[0].biting, 0.7);
        TS_ASSERT_EQUALS(values[0].findingRestSite, 0.8);
        TS_ASSERT_EQUALS(values[0].resting, 0.9);
        TS_ASSERT_EQUALS(values[0].fecundity, 1.0);
        TS_ASSERT_EQUALS(values[1].relativeAvailabilityHuman, 2.0);
        input.getDescription().getAnopheles().pop_back();
        TS_ASSERT_THROWS(OM::xml::readNhhDeployment(input, indices), const OM::util::xml_scenario_error&);
        input.getDescription().getAnopheles().push_back(description.getAnopheles()[0]);
        TS_ASSERT_THROWS(OM::xml::readNhhDeployment(input, indices), const OM::util::xml_scenario_error&);
        input.getDescription().getAnopheles().back().setMosquito("unknown");
        TS_ASSERT_THROWS(OM::xml::readNhhDeployment(input, indices), const OM::util::xml_scenario_error&);
    }

    void testOptionalEffectsAndXmlValidation() {
        scnXml::NonHumanHostsSpeciesIntervention input("mosquito");
        scnXml::DecayFunction decay("constant");
        auto effect = OM::xml::readNhhIntervention(input, decay);
        TS_ASSERT(!effect.availability.decay);
        TS_ASSERT(!effect.biting.decay);
        input.setAvailabilityReduction(scnXml::AvailabilityReduction(-0.5));
        effect = OM::xml::readNhhIntervention(input, decay);
        TS_ASSERT_EQUALS(effect.availability.initial, -0.5);
        input.setAvailabilityReduction(scnXml::AvailabilityReduction(1.1));
        TS_ASSERT_THROWS(OM::xml::readNhhIntervention(input, decay), const OM::util::xml_scenario_error&);
        input.setAvailabilityReduction(scnXml::AvailabilityReduction(0.2));
        input.setPreprandialKillingEffect(scnXml::PreprandialKillingEffect(-0.1));
        TS_ASSERT_THROWS(OM::xml::readNhhIntervention(input, decay), const OM::util::xml_scenario_error&);
    }

    void testDeploymentRetainsSamplingOrderAndRedeploys() {
        using namespace OM::Transmission::Anopheles;
        NhhIntervention effect;
        OM::util::SimpleDecayingValue* values[] = {&effect.availability, &effect.biting,
            &effect.restingSite, &effect.resting, &effect.fecundity};
        for (auto* value : values) value->set(0.3, OM::util::DecayFunction(0.4));
        OM::util::LocalRng actual(17, 0), expected(17, 0);
        for (int time : {10, 20}) {
            deploy(effect, actual, time);
            for (auto* value : values) {
                TS_ASSERT_EQUALS(*value->het, value->decay->heterogeneity.sample(expected));
                TS_ASSERT_EQUALS(value->deploy_t, time);
            }
        }
        std::ostringstream a, b;
        actual.checkpoint(a); expected.checkpoint(b);
        TS_ASSERT_EQUALS(a.str(), b.str());
    }

    void testEffectsKeepBaselineAndCompoundInOrder() {
        using namespace OM::Transmission::Anopheles;
        Nhh baseline{2.0, 0.7, 0.8, 0.9, 1.1, 50};
        NhhIntervention first, second;
        first.availability.set(0.2, OM::util::DecayFunction());
        first.biting.set(0.3, OM::util::DecayFunction());
        second.availability.set(0.4, OM::util::DecayFunction());
        OM::util::LocalRng rng(17, 0);
        auto unchanged = applyIntervention(baseline, first, 0);
        TS_ASSERT_EQUALS(unchanged.avail_i, baseline.avail_i);
        deploy(first, rng, 0); deploy(second, rng, 0);
        auto modified = applyIntervention(applyIntervention(baseline, first, 10), second, 10);
        TS_ASSERT_EQUALS(modified.avail_i, (2.0 * (1.0 - 0.2)) * (1.0 - 0.4));
        TS_ASSERT_EQUALS(modified.P_B_I, 0.7 * (1.0 - 0.3));
        TS_ASSERT_EQUALS(modified.P_C_I, baseline.P_C_I);
        TS_ASSERT_EQUALS(modified.P_D_I, baseline.P_D_I);
        TS_ASSERT_EQUALS(modified.rel_fecundity, baseline.rel_fecundity);
        TS_ASSERT_EQUALS(modified.expiry, baseline.expiry);
        TS_ASSERT_EQUALS(baseline.avail_i, 2.0);
    }
};
