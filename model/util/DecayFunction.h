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

#include "util/sampler.h"
#include <optional>

namespace OM::util {
enum class DecayKind {
    Constant, Step, Linear, Exponential, Weibull, Hill, SmoothCompact,
    Plus, Minus, Divides, Multiplies, Emax
};

struct DecayFunction {
    explicit DecayFunction(double cv = 0.0)
        : heterogeneity(LognormalSampler::fromMeanCV(1.0, cv)) {}

    DecayKind kind = DecayKind::Constant;
    bool increasing = false;
    double initialEfficacy = 1.0;
    double rate = 0.0; // inverse duration, including the distribution's scale factor
    double shape = 1.0;
    double emax = 1.0, ic50Power = 1.0, slope = 1.0, initialConcentration = 0.0;
    LognormalSampler heterogeneity;
    vector<DecayFunction> children;
};

// A person stores only the sampled factor, not a copy of the configured curve.
using DecaySample = std::optional<double>;
double eval(const DecayFunction&, double ageDays, double factor);
SimTime sampleAgeOfDecay(const DecayFunction&, LocalRng&, SimTime interval);
void checkpointDecaySample(const DecaySample&, ostream&);
void checkpointDecaySample(DecaySample&, istream&);
}
