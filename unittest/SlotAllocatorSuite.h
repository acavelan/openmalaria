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

/* T2 fixture: pins the per-host infection-slot clamp (phase-6 SlotAllocator,
 * census fact 2 / latent-bug #4) BEFORE the refactor and forever after.
 *
 * The clamp has two RNG-visible modes which must never be conflated:
 *  A. import mode (importInfection_common/_descriptive): at the cap the whole
 *     body is skipped — the genotype draw is NOT consumed (0 draws).
 *  B. update mode (update_common/_descriptive): the caller's Poisson draw was
 *     already consumed, so counts are clamped (draw-then-drop): _l first,
 *     then _i against what remains.
 * Mode A is pinned bit-exactly via Xoshiro256P state words (checkpoint API,
 * same idiom as DescriptiveInfectionSuite); mode B is a pure function pinned
 * by a full table. T1 (scenario12/MolPairwise/Molineaux) covers the rest.
 */

#include <cxxtest/TestSuite.h>
#include "Host/WithinHost/WHInterface.h" // SlotAllocator (merged)
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Genotypes.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "UnittestUtil.h"
#include <sstream>
#include <string>

using namespace OpenMalaria;

class SlotAllocatorSuite : public CxxTest::TestSuite {
  public:
    void setUp() {
        UnittestUtil::initTime(5);
        UnittestUtil::resetModelOptions();
    }

    /// Mode A: importInfection returns true iff a genotype draw may follow.
    void testImportModeTable() {
        TS_ASSERT(SlotAllocator::try_import(0));
        TS_ASSERT(SlotAllocator::try_import(20));
        TS_ASSERT_EQUALS(SlotAllocator::try_import(21), false);
        TS_ASSERT_EQUALS(SlotAllocator::try_import(22), false);
        TS_ASSERT_EQUALS(SlotAllocator::MAX_INFECTIONS, 21);
    }

    /// Mode B: clamp of to-create counts, _l before _i (original line order:
    /// CommonWithinHost.cpp:131-132, DescriptiveWithinHost.cpp:104-105).
    void testUpdateClampTable() {
        struct {
            int n_i, n_l, have, e_i, e_l, e_free;
        } T[] = {
            {5, 5, 0, 5, 5, 21}, {5, 5, 19, 0, 2, 2},  {5, 5, 20, 0, 1, 1},
            {5, 5, 21, 0, 0, 0}, {4, 1, 18, 2, 1, 3},                          // ignored = 2+1-3 = 0
            {3, 0, 20, 1, 0, 1}, {0, 0, 10, 0, 0, 11}, {30, 30, 0, 0, 21, 21}, // _l fills all 21, _i gets none
        };
        for (size_t i = 0; i < sizeof(T) / sizeof(T[0]); ++i) {
            int n_i = T[i].n_i, n_l = T[i].n_l;
            SlotAllocator::clamp_new(n_i, n_l, T[i].have);
            TS_ASSERT_EQUALS(n_i, T[i].e_i);
            TS_ASSERT_EQUALS(n_l, T[i].e_l);
            TS_ASSERT_EQUALS(SlotAllocator::free_slots(T[i].have), T[i].e_free);
        }
        // stand-alone count clamp (InfectionIncidenceModel::numNewInfections):
        TS_ASSERT_EQUALS(SlotAllocator::clamp_count(25), 21);
        TS_ASSERT_EQUALS(SlotAllocator::clamp_count(21), 21);
        TS_ASSERT_EQUALS(SlotAllocator::clamp_count(0), 0);
    }

    static void stateWords(const LocalRng& rng, uint64_t st[4]) {
        std::ostringstream os;
        const_cast<LocalRng&>(rng).checkpoint(os);
        std::string bytes = os.str();
        TS_ASSERT_EQUALS(bytes.size(), 4 * sizeof(uint64_t));
        if (bytes.size() == 4 * sizeof(uint64_t))
            memcpy(st, bytes.data(), sizeof(uint64_t) * 4);
    }

    /// Mode A bit-exact: Descriptive-kind importInfection consumes exactly
    /// zero draws at the cap, and exactly the reference draw sequence below
    /// it (genotype draw + infection duration draw). State-word goldens
    /// captured 2026-11-19 from p5-green code in the FULL-suite run context
    /// (registration order is frozen by unittest/CMakeLists.txt), Xoshiro(3, 4).
    void testImportDrawContract() {
        UnittestUtil::ctx().genotypes.initSingle();
        UnittestUtil::Infection_init_latentP_and_NaN();

        LocalRng rng(3, 4);
        WHDescriptive wh(UnittestUtil::ctx(), rng, 1.0);

        // ctor consumption pinned too (any change here shifts every word below)
        static const uint64_t GOLD_CTOR[4] = {0x0000626c286e5305ULL, 0x01c42b2629420000ULL, 0x0e7989b244225680ULL,
                                              0x817480389561d829ULL};
        uint64_t st[4];
        stateWords(rng, st);
        for (int i = 0; i < 4; ++i)
            TS_ASSERT_EQUALS(st[i], GOLD_CTOR[i]);

        // below the cap: one import consumes the pinned draw sequence
        wh.numInfs = 20;
        wh.importInfection(rng);
        stateWords(rng, st);
        static const uint64_t GOLD_IMPORT[4] = {0xf408399cc420592dULL, 0xd738b0d0bd0f8b2cULL, 0x5975fa34f30b8ea9ULL,
                                                0xba402e971e1dca0dULL};
        for (int i = 0; i < 4; ++i)
            TS_ASSERT_EQUALS(st[i], GOLD_IMPORT[i]);
        TS_ASSERT_EQUALS(wh.numInfs, 21);

        // at the cap: nothing at all — state words identical, no increment
        // (order-independent: relative to the state we just measured)
        wh.numInfs = 21;
        uint64_t before[4];
        stateWords(rng, before);
        wh.importInfection(rng);
        stateWords(rng, st);
        for (int i = 0; i < 4; ++i)
            TS_ASSERT_EQUALS(st[i], before[i]);
        TS_ASSERT_EQUALS(wh.numInfs, 21);
    }
};
