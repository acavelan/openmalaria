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

#include "Transmission/Anopheles/SugarFeeding.h"
#include <cmath>
#include <memory>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

namespace OM::Transmission::Anopheles {
namespace {

double feedingProbabilityError(double x, void* params) {
    const auto& p = *static_cast<SugarFeedingInputs*>(params);
    double total = x + p.ah + p.an + p.mu_v_interv;
    return (1.0 - exp(-total * p.theta_d)) * x / total - p.PAt_wanted;
}

double feedingProbabilityDerivative(double x, void* params) {
    const auto& p = *static_cast<SugarFeedingInputs*>(params);
    double total = x + p.ah + p.an + p.mu_v_interv;
    double a = (1.0 - exp(-total * p.theta_d)) * x / (total * total);
    double b = (1.0 - exp(-total * p.theta_d)) / total;
    double c = exp(-total * p.theta_d) * x * p.theta_d / total;
    return -a + b + c;
}

// Keep the historical combined derivative, including its different exponential terms.
void feedingProbabilityAndDerivative(double x, void* params, double* y, double* dy) {
    const auto& p = *static_cast<SugarFeedingInputs*>(params);
    double total = x + p.ah + p.an + p.mu_v_interv;
    double a = (1.0 - exp(-total)) * x / (total * total);
    double b = (1.0 - exp(-total)) / total;
    double c = exp(-total * p.theta_d) * x * p.theta_d / total;
    *y = (1.0 - exp(-total * p.theta_d)) * x / total - p.PAt_wanted;
    *dy = -a + b + c;
}
}

double sugarFeedingAvailability(SugarFeedingInputs inputs) {
    std::unique_ptr<gsl_root_fdfsolver, decltype(&gsl_root_fdfsolver_free)> solver(
        gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_secant), gsl_root_fdfsolver_free);
    gsl_function_fdf function{feedingProbabilityError, feedingProbabilityDerivative,
        feedingProbabilityAndDerivative, &inputs};
    double availability = 0.5;
    gsl_root_fdfsolver_set(solver.get(), &function, availability);
    int status, iteration = 0;
    do {
        ++iteration;
        gsl_root_fdfsolver_iterate(solver.get());
        const double previous = availability;
        availability = gsl_root_fdfsolver_root(solver.get());
        status = gsl_root_test_delta(availability, previous, 0, 1e-4);
    } while (status == GSL_CONTINUE && iteration < 20);
    return availability;
}

}
