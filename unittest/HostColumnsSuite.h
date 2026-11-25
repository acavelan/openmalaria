#pragma once

#include <cxxtest/TestSuite.h>
#include "Host/HostColumns.h"
#include <sstream>

using namespace OpenMalaria;

class HostColumnsSuite : public CxxTest::TestSuite {
  public:
    // A leading survivor, overlapping runs, gaps, all-dead and later births.
    // RNG bytes must follow their human, independent of the species stride.
    void testStableCompaction() {
        for (size_t species : {0u, 1u, 3u}) {
            HumanColumns cols;
            MasterRng master(0, 1);
            std::vector<std::string> rng;
            for (int i = 0; i < 8; ++i) {
                cols.append_row(fromDays(i), master, species);
                std::ostringstream out;
                cols.rng.back().checkpoint(out);
                rng.push_back(out.str());
                for (size_t s = 0; s < species; ++s) {
                    double x = 10 * i + s;
                    cols.anoph[i * species + s] = {x, x + 1, x + 2, x + 3};
                }
            }
            const std::vector<uint32_t> keep{0, 2, 3, 5, 6, 7};
            cols.compact(keep);
            TS_ASSERT_EQUALS(cols.anoph.size(), keep.size() * species);
            for (size_t i = 0; i < keep.size(); ++i) {
                TS_ASSERT_EQUALS(cols.dob[i], fromDays(keep[i]));
                std::ostringstream out;
                cols.rng[i].checkpoint(out);
                TS_ASSERT_EQUALS(out.str(), rng[keep[i]]);
                for (size_t s = 0; s < species; ++s) {
                    const auto& p = cols.anoph[i * species + s];
                    double x = 10 * keep[i] + s;
                    TS_ASSERT_EQUALS(p.ento_raw, x);
                    TS_ASSERT_EQUALS(p.ento_avail, x + 1);
                    TS_ASSERT_EQUALS(p.prob_biting, x + 2);
                    TS_ASSERT_EQUALS(p.prob_resting, x + 3);
                }
            }
            cols.compact({});
            TS_ASSERT(cols.anoph.empty());
            TS_ASSERT_EQUALS(cols.append_row(zero(), master, species), 0u);
            TS_ASSERT_EQUALS(cols.anoph.size(), species);
        }
    }
};
