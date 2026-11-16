#pragma once
#include <cxxtest/TestSuite.h>
#include "Transmission/Anopheles/SugarFeeding.h"
#include <cmath>

class SugarFeedingSuite : public CxxTest::TestSuite {
public:
    void testAvailabilityReproducesRequestedFeedingProbability() {
        using OM::Transmission::Anopheles::sugarFeedingAvailability;
        for (double human : {0.1, 1.0, 10.0})
        for (double animal : {0.0, 0.3})
        for (double mortality : {0.02, 0.5})
        for (double duration : {0.2, 1.0})
        for (double target : {0.01, 0.1, 0.5, 0.9}) {
            const double availability = sugarFeedingAvailability(
                {human, animal, mortality, duration, target});
            const double total = availability + human + animal + mortality;
            const double probability = (1.0 - std::exp(-total * duration)) * availability / total;
            TS_ASSERT(availability > 0.0);
            // Legacy convergence tests changes in availability, not this probability residual.
            TS_ASSERT_DELTA(probability, target, 1e-7);
        }
    }
};
