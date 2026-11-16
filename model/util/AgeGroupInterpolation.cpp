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

#include "util/AgeGroupInterpolation.h"
#include "util/errors.h"
#include <fstream>

namespace OM::util {
void AgeGroupInterpolator::set(const AgeGroupValues& input, double maxAgeYears,
        const char* name, bool outputSamples) {
    reset();
    if (input.groups.empty())
        throw xml_scenario_error(string(name) + ": at least one age group required");
    map<double, double> data;
    auto pos = data.begin();
    if (input.interpolation == AgeInterpolation::Constant) {
        double greatestLbound = -1.0;
        for (const auto& group : input.groups) {
            double lbound = group.lowerBound;
            if (lbound >= greatestLbound) {
                greatestLbound = lbound;
                pos = data.insert(pos, make_pair(lbound, group.value));
            } else
                throw xml_scenario_error(string(name) + ": lower bound " + to_string(lbound) + " not greater than previous " + to_string(greatestLbound));
        }
        if (data.begin()->first != 0.0)
            throw xml_scenario_error(string(name) + ": first lower-bound must be 0");
    } else {
        if (input.groups.front().lowerBound != 0.0)
            throw xml_scenario_error(string(name) + ": first lower-bound must be 0");
        double greatestLbound = 0.0;
        double lastValue = input.groups.front().value;
        for (const auto& group : input.groups) {
            double lbound = group.lowerBound;
            if (lbound >= greatestLbound) {
                double insBound = 0.5*(greatestLbound + lbound);
                pos = data.insert(pos, make_pair(insBound, lastValue));
                greatestLbound = lbound;
                lastValue = group.value;
            } else
                throw xml_scenario_error(string(name) + ": lower bound " + to_string(lbound) + " less than previous " + to_string(greatestLbound));
        }
        data[0.5*(greatestLbound + maxAgeYears)] = lastValue;
        data[0.0] = data.begin()->second;
        data[numeric_limits<double>::infinity()] = data.rbegin()->second;
    }
    interpolation = input.interpolation;
    points = std::move(data);
    if (outputSamples) {
        ofstream stream(string(name) + ".csv");
        for (double age = 0.0; age < maxAgeYears; age += 0.1)
            stream << age << "," << eval(age) << endl;
    }
}

double AgeGroupInterpolator::eval(double ageYears) const {
    if (!isSet()) throw TRACED_EXCEPTION_DEFAULT("AgeGroupDummy");
    assert(ageYears >= 0.0);
    auto it = points.upper_bound(ageYears);
    if (interpolation == AgeInterpolation::Constant) return std::prev(it)->second;
    assert(it != points.end());
    double a1 = it->first;
    double f1 = it->second;
    --it;
    double a0 = it->first;
    double f0 = it->second;
    return (ageYears - a0) / (a1 - a0) * (f1 - f0) + f0;
}

void AgeGroupInterpolator::scale(double factor) {
    if (!isSet()) throw TRACED_EXCEPTION_DEFAULT("AgeGroupDummy");
    for (auto& point : points) point.second *= factor;
}

double AgeGroupInterpolator::firstGlobalMaximum() const {
    if (!isSet()) throw TRACED_EXCEPTION_DEFAULT("AgeGroupDummy");
    double age = 0.0, max = numeric_limits<double>::min();
    for (const auto& point : points) {
        if (point.second > max) {
            max = point.second;
            age = point.first;
        }
    }
    return age;
}
}
