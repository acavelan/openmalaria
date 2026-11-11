#pragma once

#include <cxxtest/TestSuite.h>
#include "Host/HumanTransmission.h"
#include "util/checkpoint_containers.h"
#include <sstream>

class HumanTransmissionSuite : public CxxTest::TestSuite {
public:
    void testSpeciesRecordsPreserveLegacyCheckpointBytes() {
        using namespace OM;
        using namespace OM::util::checkpoint;
        Host::HumanTransmission original;
        original.species = {{0.2, 0.3, 0.4, 0.5}, {0.6, 0.7, 0.8, 0.9}};
        original.relativeAvailability = 1.25;
        original.outsideTransmission = true;

        std::ostringstream legacy;
        for (auto field : {std::vector<double>{0.2, 0.6}, {0.3, 0.7}, {0.4, 0.8}, {0.5, 0.9}})
            field & legacy;
        original.relativeAvailability & legacy;
        original.outsideTransmission & legacy;

        std::ostringstream written;
        Host::checkpoint(original, written);
        TS_ASSERT_EQUALS(written.str(), legacy.str());

        Host::HumanTransmission restored;
        std::istringstream input(legacy.str());
        Host::checkpoint(restored, input);
        TS_ASSERT_EQUALS(restored.species.size(), 2u);
        TS_ASSERT_EQUALS(restored.species[1].sampledAvailability, 0.6);
        TS_ASSERT_EQUALS(restored.species[0].availability, 0.3);
        TS_ASSERT_EQUALS(restored.species[1].biteProbability, 0.8);
        TS_ASSERT_EQUALS(restored.species[0].restingProbability, 0.5);
        std::ostringstream rewritten;
        Host::checkpoint(restored, rewritten);
        TS_ASSERT_EQUALS(rewritten.str(), legacy.str());
    }

    void testCheckpointRejectsInconsistentSpeciesCounts() {
        using namespace OM;
        using namespace OM::util::checkpoint;
        std::ostringstream corrupt;
        std::vector<double> first{0.2, 0.6}, second{0.3};
        first & corrupt;
        second & corrupt;
        std::istringstream input(corrupt.str());
        Host::HumanTransmission state;
        TS_ASSERT_THROWS(Host::checkpoint(state, input), const util::checkpoint_error&);
    }
};
