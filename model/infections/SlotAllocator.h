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

#ifndef Hmod_infections_SlotAllocator
#define Hmod_infections_SlotAllocator

#include <algorithm>

namespace OM {
namespace WithinHost {

/** Sole owner of the per-host infection-slot cap (plan §5: "MAX_INFECTIONS=21
 * lives in one place"). The only real reason for the limit is to prevent
 * bad input from making the number of infections balloon stupidly.
 *
 * The cap has TWO RNG-visible modes which must never be conflated
 * (census fact 2, latent-bugs.md #4; pinned by unittest/SlotAllocatorSuite.h):
 *  A. import mode  — try_import(): at the cap the caller skips the entire
 *     create, including the genotype draw (0 draws consumed).
 *     // parity: CommonWithinHost.cpp:111, DescriptiveWithinHost.cpp:81
 *  B. update mode  — clamp_new(): the caller's Poisson draw was already
 *     consumed, so counts are clamped draw-then-drop; _l clamps first, then
 *     _i against what remains.
 *     // parity: CommonWithinHost.cpp:131-132, DescriptiveWithinHost.cpp:104-105
 */
struct SlotAllocator {
    /// Exact constraint is: _MOI <= MAX_INFECTIONS.
    static constexpr int MAX_INFECTIONS = 21;

    /// Free slots for a host with num_infs infections (0 at the cap).
    static constexpr int free_slots(int num_infs) { return MAX_INFECTIONS - num_infs; }

    /// Mode A: true iff an import may proceed (genotype draw included).
    static constexpr bool try_import(int num_infs) { return num_infs < MAX_INFECTIONS; }

    /// Mode B: clamp to-create counts in place (_l first — order is load-bearing).
    static inline void clamp_new(int& n_i, int& n_l, int num_infs) {
        n_l = std::min(n_l, free_slots(num_infs));
        n_i = std::min(n_i, free_slots(num_infs) - n_l);
    }

    /// Clamp a stand-alone drawn count (post-draw; no rng effect).
    /// parity: InfectionIncidenceModel.cpp:257-258
    static constexpr int clamp_count(int want) {
        return want > MAX_INFECTIONS ? MAX_INFECTIONS : want;
    }
};

} }
#endif
