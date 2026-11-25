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

#pragma once

/* Revived 2026-11-19. The previous version of this file was a test-less
 * skeleton (only setUp/tearDown) and no longer compiled: it called
 * DescriptiveInfection::initParameters() (now init(const Parameters&)),
 * UnittestUtil::Infection_init() (removed), the default ctor
 * DescriptiveInfection() (now requires rng/genotype/origin) and the removed
 * gsl::setUp/tearDown (old RNG bootstrap). determine_densities_descriptive() remains
 * untested (many inputs, 2 outputs — needs full Human scaffolding; see
 * T2FixtureManifest.md).
 * 
 * What is pinned here is the small pure function DescriptiveInfection::sample_duration()
 * (DescriptiveInfection.cpp:128): log-normal draw over hard-coded params plus
 * integer time-step rounding — a T2 fixture over the RNG primitive layer,
 * goldens generated with the same adapter as rng_dist_v50.csv.
 */

#include <cxxtest/TestSuite.h>
#include "Host/WithinHost/Infection/Infection.h"
#include "util/random.h"
#include "UnittestUtil.h"
#include <sstream>
#include <string>

using namespace OM::WithinHost;
using util::LocalRng;

class DescriptiveInfectionSuite : public CxxTest::TestSuite {
  public:
    void setUp() { UnittestUtil::initTime(5); }

    /// DescriptiveInfection::sample_duration() over a fixed stream: per draw, the returned
    /// infection duration and the 4 post-draw Xoshiro256P state words.
    /// Values generated 2026-11-19 with golden-gen/RngGoldenAdapter.h (GCC
    /// 16.2.1, GSL 2.8, -O2 -ffp-contract=off), mirroring the hard-coded
    /// params in DescriptiveInfection.cpp:128-134.
    void testInfectionDuration() {
        // days = int(1 + lognormal) before the fromDays(...)-oneTS() rounding;
        // SimTime is days (sim.h:74), oneTS() = 5 after initTime(5).
        static const struct {
            int days;
            uint64_t s0, s1, s2, s3;
        } GOLDEN[] = {
            {300, 0x0000626c286e5305ULL, 0x01c42b2629420000ULL, 0x0e7989b244225680ULL, 0x817480389561d829ULL},
            {170, 0xf408399cc420592dULL, 0xd738b0d0bd0f8b2cULL, 0x5975fa34f30b8ea9ULL, 0xba402e971e1dca0dULL},
            {225, 0xab11f90cfedeb0c6ULL, 0x2fe96d14cc659320ULL, 0xb35d0a24ff11cf88ULL, 0x9519464c2bdaf33dULL},
            {251, 0x8a0d1aa3dc85c9e4ULL, 0xba7b7c7c629620a9ULL, 0x7fa91141e641fcc5ULL, 0xb202f8d40cdc3bfbULL},
            {95, 0x8ec1a89ada08caf7ULL, 0xc0f7275b910bf21fULL, 0x6014e061e249e797ULL, 0xe3082196a6d22d18ULL},
            {155, 0xb8fc2588a4426c74ULL, 0x237eab1b98e4e7efULL, 0xd21d962e8f503890ULL, 0x6f3082b8517be932ULL},
            {136, 0x5cf6bd1f1b1fa891ULL, 0x81fb96f93ae78546ULL, 0xf91fe9a9ef07364dULL, 0x5947150896268ed8ULL},
            {492, 0x812b27b6e95f8d0dULL, 0x2843ddd887b52649ULL, 0x88ce82697e7e3dd3ULL, 0x25c040ac232b0bd0ULL},
            {418, 0x68bde3752efba290ULL, 0x1f5a687086e686ddULL, 0xcef2315ebb5c104aULL, 0x40409c82ab36ec67ULL},
            {339, 0x83e5e63f4ef4b997ULL, 0xf81c723a88774ef7ULL, 0x35bfeee3f3387af0ULL, 0xee37b6885e3709bbULL},
            {79, 0xd361daaed5992983ULL, 0x89a740598d2db02cULL, 0x322d51c8a3b63dbcULL, 0xbaeb08d5ff0469a5ULL},
            {189, 0xb3f77f74d5a3da5cULL, 0xe939c8204ac54026ULL, 0x17adf51af9e1e435ULL, 0x654aca7b5daace42ULL},
        };

        LocalRng rng(3, 4); // Xoshiro256P(3, 4, default, default) — engine seed defaults in util/random.h
        for (size_t i = 0; i < sizeof(GOLDEN) / sizeof(GOLDEN[0]); ++i) {
            SimTime dur = DescriptiveInfection::sample_duration(UnittestUtil::ctx(), rng);
            TS_ASSERT_EQUALS(dur, sim::fromDays(GOLDEN[i].days) - UnittestUtil::clock().one_TS());

            // state words after this draw, via the model's own checkpoint API
            std::ostringstream os;
            rng.checkpoint(os);
            std::string bytes = os.str();
            uint64_t st[4];
            TS_ASSERT_EQUALS(bytes.size(), 4 * sizeof(uint64_t));
            if (bytes.size() == 4 * sizeof(uint64_t)) {
                memcpy(st, bytes.data(), sizeof st);
                TS_ASSERT_EQUALS(st[0], GOLDEN[i].s0);
                TS_ASSERT_EQUALS(st[1], GOLDEN[i].s1);
                TS_ASSERT_EQUALS(st[2], GOLDEN[i].s2);
                TS_ASSERT_EQUALS(st[3], GOLDEN[i].s3);
            }
        }
    }
};
