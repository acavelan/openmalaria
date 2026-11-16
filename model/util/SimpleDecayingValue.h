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
#include "util/DecayFunction.h"

namespace OM::util {
struct SimpleDecayingValue {
    std::optional<DecayFunction> decay;
    double initial = 0.0;
    DecaySample het;
    SimTime deploy_t = sim::never();

    void set(double value, DecayFunction definition) {
        decay = std::move(definition);
        initial = value;
    }
    void set(double value) { initial = value; }

    void deploy(LocalRng& rng, SimTime time) {
        if (!decay) return;
        deploy_t = time;
        het = decay->heterogeneity.sample(rng);
    }

    double current_value(SimTime time) const {
        if (!decay || !het) return 0.0;
        return initial * eval(*decay, time - deploy_t, *het);
    }

    template<class S> void operator&(S& stream) {
        if (decay) checkpointDecaySample(het, stream);
        deploy_t & stream;
    }
};
}
