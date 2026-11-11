#pragma once

#include <cxxtest/TestSuite.h>
#include "Host/InfectionIncidenceModel.h"
#include <sstream>

class IncidenceSuite : public CxxTest::TestSuite {
    using State = OM::Host::IncidenceState;
    using Parameters = OM::Host::IncidenceParameters;
    using Rng = OM::util::LocalRng;

    static void step(State& state, const Parameters& parameters, Rng& rng) {
        const double expected = OM::Host::expectedInfections(
            state, parameters, rng, 2.5, 0.2, 0.8, 5);
        OM::Host::sampleNewInfections(rng, expected);
    }

    static std::string checkpoint(State& state, Rng& rng) {
        std::ostringstream out;
        state & out;
        rng.checkpoint(out);
        return out.str();
    }

public:
    void testInterleavedModelsMatchIsolatedModels() {
        Parameters gamma;
        gamma.kind = OM::Host::IncidenceKind::NegativeBinomial;
        gamma.noPreErythrocytic = true;
        gamma.inf_rate_shape_param = 2.0;
        Parameters lognormal = gamma;
        lognormal.kind = OM::Host::IncidenceKind::Lognormal;
        lognormal.inf_rate_shape_param = 0.3;
        lognormal.inf_rate_offset = 0.045;

        State a, b, isolatedA, isolatedB;
        Rng rngA(17, 0), rngB(23, 0), isolatedRngA(17, 0), isolatedRngB(23, 0);
        for (int i = 0; i < 20; ++i) {
            step(a, gamma, rngA);
            step(b, lognormal, rngB);
        }
        for (int i = 0; i < 20; ++i) step(isolatedA, gamma, isolatedRngA);
        for (int i = 0; i < 20; ++i) step(isolatedB, lognormal, isolatedRngB);
        TS_ASSERT_EQUALS(checkpoint(a, rngA), checkpoint(isolatedA, isolatedRngA));
        TS_ASSERT_EQUALS(checkpoint(b, rngB), checkpoint(isolatedB, isolatedRngB));
    }

    void testCheckpointResumesTheSameRandomSequence() {
        Parameters parameters;
        parameters.kind = OM::Host::IncidenceKind::Lognormal;
        parameters.noPreErythrocytic = true;
        parameters.inf_rate_shape_param = 0.3;
        parameters.inf_rate_offset = 0.045;
        State original, restored;
        Rng originalRng(17, 0), restoredRng(99, 0);
        step(original, parameters, originalRng);
        std::istringstream in(checkpoint(original, originalRng));
        restored & in;
        restoredRng.checkpoint(in);
        for (int i = 0; i < 20; ++i) {
            step(original, parameters, originalRng);
            step(restored, parameters, restoredRng);
        }
        TS_ASSERT_EQUALS(checkpoint(original, originalRng), checkpoint(restored, restoredRng));
    }
};
