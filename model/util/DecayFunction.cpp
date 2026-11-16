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

#include "util/DecayFunction.h"

namespace OM::util {
static double compute(const DecayFunction& f, double age, double factor) {
    switch (f.kind) {
    case DecayKind::Constant:
        return age * factor == numeric_limits<double>::infinity() ? 0.0 : 1.0;
    case DecayKind::Step:
        return age * f.rate * factor < 1.0 ? 1.0 : 0.0;
    case DecayKind::Linear:
        return age * f.rate * factor < 1.0 ? 1.0 - age * f.rate * factor : 0.0;
    case DecayKind::Exponential:
        return exp(-age * f.rate * factor);
    case DecayKind::Weibull: {
        double p = -pow(age * f.rate * factor, f.shape);
        return p < -700.0 ? 0.0 : exp(p);
    }
    case DecayKind::Hill:
        return 1.0 / (1.0 + pow(age * f.rate * factor, f.shape));
    case DecayKind::SmoothCompact:
        return age * f.rate * factor < 1.0
            ? exp(f.shape - f.shape / (1.0 - pow(age * f.rate * factor, 2.0))) : 0.0;
    case DecayKind::Emax: {
        const double concentrationPow = pow(f.initialConcentration * eval(f.children.at(0), age, factor), f.slope);
        const double efficacy = f.emax * concentrationPow / (concentrationPow + f.ic50Power);
        return max(min(efficacy, 1.0), 0.0);
    }
    case DecayKind::Plus:
        return max(min(std::plus<double>{}(eval(f.children.at(0), age, factor),
            eval(f.children.at(1), age, factor)), 1.0), 0.0);
    case DecayKind::Minus:
        return max(min(std::minus<double>{}(eval(f.children.at(0), age, factor),
            eval(f.children.at(1), age, factor)), 1.0), 0.0);
    case DecayKind::Divides:
        return max(min(std::divides<double>{}(eval(f.children.at(0), age, factor),
            eval(f.children.at(1), age, factor)), 1.0), 0.0);
    case DecayKind::Multiplies:
        return max(min(std::multiplies<double>{}(eval(f.children.at(0), age, factor),
            eval(f.children.at(1), age, factor)), 1.0), 0.0);
    }
    throw TRACED_EXCEPTION_DEFAULT("invalid decay function");
}

double eval(const DecayFunction& f, double ageDays, double factor) {
    if (f.increasing) return 1.0 - compute(f, ageDays, factor) * f.initialEfficacy;
    return compute(f, ageDays, factor) * f.initialEfficacy;
}

SimTime sampleAgeOfDecay(const DecayFunction& f, LocalRng& rng, SimTime interval) {
    const auto round = [interval](double days) {
        return interval * floorToInt(days / interval + 0.5);
    };
    switch (f.kind) {
    case DecayKind::Constant: return sim::future();
    case DecayKind::Step: return round(1.0 / f.rate);
    case DecayKind::Linear: return round(rng.uniform_01() / f.rate);
    case DecayKind::Exponential: return round(-log(rng.uniform_01()) / f.rate);
    case DecayKind::Weibull: return round(pow(-log(rng.uniform_01()), 1.0/f.shape) / f.rate);
    case DecayKind::Hill: return round(pow(1.0 / rng.uniform_01() - 1.0, 1.0/f.shape) / f.rate);
    case DecayKind::SmoothCompact:
        return round(sqrt(1.0 - f.shape / (f.shape - log(rng.uniform_01()))) / f.rate);
    case DecayKind::Emax:
        return round(sampleAgeOfDecay(f.children.at(0), rng, interval));
    case DecayKind::Plus:
    case DecayKind::Minus:
    case DecayKind::Divides:
    case DecayKind::Multiplies:
        // Keep the original argument order: both calls consume the same RNG.
        return round(max(sampleAgeOfDecay(f.children.at(0), rng, interval),
                         sampleAgeOfDecay(f.children.at(1), rng, interval)));
    }
    throw TRACED_EXCEPTION_DEFAULT("invalid decay function");
}

void checkpointDecaySample(const DecaySample& sample, ostream& stream) {
    bool present = sample.has_value();
    present & stream;
    if (present) *sample & stream;
}

void checkpointDecaySample(DecaySample& sample, istream& stream) {
    bool present;
    present & stream;
    if (present) {
        double factor;
        factor & stream;
        sample = factor;
    } else sample.reset();
}
}
