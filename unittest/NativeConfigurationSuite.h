#pragma once

#include <cxxtest/TestSuite.h>
#include "Host/InfectionIncidenceModel.h"
#include "Host/HumanTransmission.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "UnittestUtil.h"
#include "mon/Continuous.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"

class NativeConfigurationSuite : public CxxTest::TestSuite {
public:
    void testPathogenesisBranchesAndLegacyCheckpointFields() {
        using namespace OM::WithinHost::Pathogenesis;
        using namespace OM::util::checkpoint;
        for (auto kind : {PathogenesisParameters::Kind::Pyrogenic,
                          PathogenesisParameters::Kind::Mueller,
                          PathogenesisParameters::Kind::Predetermined}) {
            PathogenesisParameters config;
            config.kind = kind;
            config.initPyroThres = 10.0;
            config.Ystar1_26 = 100.0;
            config.Ystar2_13 = 20.0;
            config.a = config.b = 0.0;
            config.rateMultiplier_31 = 0.2;
            config.densityExponent_32 = 1.0;
            PathogenesisModel state(1.25, config);
            const double probability = state.getPEpisode(5.0, 2.0);
            if (kind == PathogenesisParameters::Kind::Mueller) {
                TS_ASSERT_DELTA(probability, 0.3296799539643607, 1e-15);
            } else if (kind == PathogenesisParameters::Kind::Pyrogenic) {
                TS_ASSERT_EQUALS(probability, 1.0 / 3.0);
            } else {
                TS_ASSERT_EQUALS(probability, 0.0);
                TS_ASSERT_EQUALS(state.getPEpisode(10.0, 2.0), 0.0);
                TS_ASSERT_EQUALS(state.getPEpisode(20.0, 2.0), 1.0);
            }
            state._pyrogenThres = 12.5;
            std::stringstream actual, expected;
            state & static_cast<std::ostream&>(actual);
            double comorbidity = 1.25, threshold = 12.5;
            comorbidity & static_cast<std::ostream&>(expected);
            if (kind != PathogenesisParameters::Kind::Mueller)
                threshold & static_cast<std::ostream&>(expected);
            TS_ASSERT_EQUALS(actual.str(), expected.str());
            PathogenesisModel restored(0.0, config);
            restored & static_cast<std::istream&>(actual);
            TS_ASSERT_EQUALS(restored._comorbidityFactor, 1.25);
            TS_ASSERT_EQUALS(restored._pyrogenThres,
                kind == PathogenesisParameters::Kind::Mueller ? 10.0 : 12.5);
        }
    }

    void testNonMalariaFeverParametersAreIndependent() {
        using namespace OM::WithinHost::Pathogenesis;
        UnittestUtil::initTime(5);
        PathogenesisParameters none, certain;
        certain.pg_NMF_incidence.set({OM::util::AgeInterpolation::Constant, {{0, 1.0}}}, 90, "test");
        OM::util::LocalRng rng(19, 0), untouched(19, 0);
        TS_ASSERT_EQUALS(sampleNMF(none, rng, 20.0), NONE);
        TS_ASSERT_EQUALS(rng.uniform_01(), untouched.uniform_01());
        TS_ASSERT_EQUALS(sampleNMF(certain, rng, 20.0), STATE_NMF);
        TS_ASSERT_EQUALS(sampleNMF(none, rng, 20.0), NONE);
    }

    void testContinuousOutputCanBeConfiguredWithoutXmlOrCommandLine() {
        UnittestUtil::initTime(5);
        OM::mon::Continuous::State state;
        OM::mon::Continuous::init(state, "", false);
        TS_ASSERT(!state.ctsOStream.is_open());
        state.ctsPeriod = OM::sim::oneTS();
        TS_ASSERT_THROWS(OM::mon::Continuous::init(state, "", true),
            const OM::util::checkpoint_error&);
    }

    void testNativeDailyEirMeanPreservesOrderAndValidation() {
        using OM::Transmission::averageEIR;
        TS_ASSERT_EQUALS(averageEIR({1.0, 2.0, 3.0}), 2.0);
        TS_ASSERT_EQUALS(averageEIR({1e16, 1.0, -1e16}), 0.0);
        TS_ASSERT_THROWS(averageEIR({}), const OM::util::xml_scenario_error&);
    }

    void testNativeVectorEirAggregationAndValidation() {
        UnittestUtil::initTime(5);
        OM::Transmission::Anopheles::AnophelesModel model;
        std::vector<double> daily(365, 2.0);
        model.changeEIRIntervention(daily);
        TS_ASSERT_EQUALS(model.interventionEIR.size(), 73u);
        for (double value : model.interventionEIR) TS_ASSERT_EQUALS(value, 10.0);
        daily.back() = 7.0;
        daily.front() = 0.0;
        model.changeEIRIntervention(daily);
        TS_ASSERT_EQUALS(model.interventionEIR.back(), 15.0);
        TS_ASSERT(model.interventionEIR.front() > 8.0);
        TS_ASSERT(model.interventionEIR.front() < 8.1);
        daily.resize(364);
        TS_ASSERT_THROWS(model.changeEIRIntervention(daily), const OM::util::xml_scenario_error&);
    }

    void testIncidenceOptionsAndValidation() {
        using namespace OM::Host;
        IncidenceConfiguration input;
        input.baselineAvailabilityShape = 0.2;
        input.gammaP = 2.0;
        input.negativeLogOneMinusSinf = 0.0;
        input.Simm = 0.4;
        input.EStar = 4.0;
        input.XStarP = 8.0;
        input.noPreErythrocytic = true;
        input.vaccineGenotype = true;
        for (auto kind : {IncidenceKind::Default, IncidenceKind::Heterogeneity,
                          IncidenceKind::NegativeBinomial, IncidenceKind::Lognormal}) {
            input.kind = kind;
            auto parameters = makeIncidenceParameters(input);
            TS_ASSERT_EQUALS(parameters.kind, kind);
            TS_ASSERT_EQUALS(parameters.Sinf, 0.0);
            TS_ASSERT_EQUALS(parameters.EstarInv, 0.25);
            TS_ASSERT_EQUALS(parameters.Xstar_pInv, 0.125);
            TS_ASSERT(parameters.noPreErythrocytic);
            TS_ASSERT(parameters.vaccineGenotype);
            TS_ASSERT(std::isfinite(parameters.inf_rate_shape_param));
        }
        input.baselineAvailabilityShape = 1.0;
        TS_ASSERT_THROWS(makeIncidenceParameters(input), const OM::util::xml_scenario_error&);
    }

    void testConcreteSamplerValuesAreIndependentAndPreserveDraws() {
        using namespace OM::util;
        auto original = LognormalSampler::fromMeanCV(2.0, 0.4);
        auto copy = original;
        copy.scaleMean(3.0);
        TS_ASSERT_DELTA(original.mean(), 2.0, 1e-14);
        TS_ASSERT_DELTA(copy.mean(), 6.0, 1e-14);
        LocalRng actual(17, 0), expected(17, 0);
        for (double cv : {0.0, 0.4}) {
            DistributionParameters input;
            input.kind = DistributionKind::Lognormal;
            input.cv = cv;
            input.truncate = 2.5;
            auto value = createSampler<LognormalSampler>(2.0, input);
            auto dynamic = createSampler(2.0, input);
            for (int i = 0; i < 8; ++i)
                TS_ASSERT_EQUALS(value.sample(actual), dynamic->sample(expected));
            input.kind = DistributionKind::Gamma;
            auto gamma = createSampler<GammaSampler>(2.0, input);
            dynamic = createSampler(2.0, input);
            for (int i = 0; i < 8; ++i)
                TS_ASSERT_EQUALS(gamma.sample(actual), dynamic->sample(expected));
        }
        std::ostringstream a, b;
        actual.checkpoint(a); expected.checkpoint(b);
        TS_ASSERT_EQUALS(a.str(), b.str());
    }

    void testSamplerInputs() {
        using namespace OM::util;
        LocalRng rng(17, 0);
        DistributionParameters input;
        TS_ASSERT_DELTA(createSampler(3.0, input)->sample(rng), 3.0, 1e-14);
        input.cv = 1.0;
        TS_ASSERT_THROWS(createSampler(3.0, input), const xml_scenario_error&);
        for (auto kind : {DistributionKind::Lognormal, DistributionKind::Gamma}) {
            input.kind = kind;
            TS_ASSERT_DELTA(createSampler(3.0, input)->mean(), 3.0, 1e-14);
            input.variance = 2.0;
            TS_ASSERT_THROWS(createSampler(3.0, input), const xml_scenario_error&);
            input.cv.reset();
            TS_ASSERT_DELTA(createSampler(3.0, input)->mean(), 3.0, 1e-14);
            input.variance.reset();
            TS_ASSERT_THROWS(createSampler(3.0, input), const xml_scenario_error&);
            input.cv = 1.0;
        }
    }

    void testHumanTransmissionFromNativeValues() {
        using namespace OM;
        Transmission::HostSpeciesConfiguration input{{}, {0.75, 0.0}, {0.5, 0.0}, {0.25, 0.0}};
        Transmission::HostTransmissionParameters parameters;
        parameters.species.push_back(Transmission::makeHostSpeciesParameters(input));
        parameters.species[0].entoAvailabilityFactor = 2.0;
        parameters.availabilityByAge.set({util::AgeInterpolation::Constant, {{0, 0.5}, {20, 1.0}}}, 90, "test");
        util::LocalRng rng(17, 0), untouched(17, 0);
        Transmission::calculateAvailabilityPercentiles(parameters);
        TS_ASSERT_EQUALS(parameters.availabilityPercentile(0), 0.0);
        TS_ASSERT_EQUALS(parameters.availabilityPercentile(50), 2.0);
        Host::HumanTransmission host;
        Host::initialize(host, rng, 3.0, parameters);
        TS_ASSERT_EQUALS(host.species[0].availability, 6.0);
        TS_ASSERT_EQUALS(host.species[0].biteProbability, 0.75);
        TS_ASSERT_EQUALS(host.species[0].restingProbability, 0.125);
        TS_ASSERT_EQUALS(parameters.availabilityByAge.eval(10), 0.5);
        std::ostringstream actual, expected;
        rng.checkpoint(actual);
        untouched.checkpoint(expected);
        TS_ASSERT_EQUALS(actual.str(), expected.str());
    }
};
