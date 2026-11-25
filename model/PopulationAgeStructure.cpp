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

// p8-files: PopulationAgeStructure.h rides the RunContext.h cone; demography
// spline/parse bodies stay here (wide-cone rule).
#include "PopulationAgeStructure.h"
#include "Clock.h"
#include "RunContext.h"
#include "util/errors.h"

#include <cmath>
#include <gsl_vector_double.h>
#include <gsl/gsl_multimin.h>
#include <fstream>

namespace OpenMalaria {

size_t AgeStructure::getMaxTStepsPerLife(const RunContext& ctx) { return ctx.age.cumAgeProp.size(); }

void AgeStructure::init(RunContext& ctx, const scnXml::Demography& demography) {
    auto& cumAgeProp = ctx.age.cumAgeProp;
    // this number of cells are needed:
    cumAgeProp.resize(ctx.clock.inSteps(ctx.clock.max_human_age) + 1);

    estimateRemovalRates(ctx, demography);
    calcCumAgeProp(ctx);
}

int AgeStructure::targetCumPop(RunContext& ctx, size_t ageTSteps, int targetPop) {
    const auto& cumAgeProp = ctx.age.cumAgeProp;
    size_t index = cumAgeProp.size() - 1 - ageTSteps;
    return (int)floor(cumAgeProp[index] * targetPop + 0.5);
}

/** Fit objective: set removal-rate parameters from the trial values
 * (param1, param2) and return residual sum-of-squares vs the field-data age
 * groups; splits the infant age groups as a side effect. Called back by the
 * GSL minimiser (wCalcRSS); formerly a public static AgeStructure member,
 * file-local since its only caller is in this TU. Body verbatim. */
static double setDemoParameters(RunContext& ctx, double param1, double param2) {
    auto& ageGroupBounds = ctx.age.ageGroupBounds;
    auto& ageGroupPercent = ctx.age.ageGroupPercent;
    double& mu0 = ctx.age.mu0;
    double& mu1 = ctx.age.mu1;
    double& alpha0 = ctx.age.alpha0;
    double& alpha1 = ctx.age.alpha1;
    constexpr double rho = AgeStructureState::rho;
    const double IMR = 0.1;
    double M_inf = -log(1 - IMR);

    mu1 = exp(param1) / 100;
    alpha1 = exp(param2) / 100;
    alpha0 = 4.0;
    mu0 = (M_inf - mu1 * (exp(alpha1 * 0.5) - 1) * alpha0) / (alpha1 * (1 - exp(-alpha0 * 0.5)));

    size_t ngroups = ageGroupPercent.size();
    double sumpred = 0.0;
    vector<double> M(ngroups, 0);
    vector<double> pred(ngroups, 0);

    for (size_t i = 0; i < ngroups - 1; i++) {
        double midpt = (ageGroupBounds[i + 1] + ageGroupBounds[i]) * 0.5;
        double M1 = mu0 * (1.0 - exp(-alpha0 * midpt)) / alpha0;
        double M2 = mu1 * (exp(alpha1 * midpt) - 1.0) / alpha1;
        M[i] = M1 + M2;
        pred[i] = (ageGroupBounds[i + 1] - ageGroupBounds[i]) * exp(-rho * midpt - M[i]);
        sumpred += pred[i];
    }
    for (size_t i = 0; i < ngroups - 1; i++) {
        pred[i] = pred[i] / sumpred * 100.0;
    }
    double L_inf = exp(-rho * 0.5 - M[1]);
    double M_nn = -log(1.0 - 0.4 * (1 - exp(-M[1])));
    double L1 = 1.0 / 12.0 * exp(-rho / 24.0 - M_nn);
    double perc_inf = ageGroupPercent[0] + ageGroupPercent[1];
    ageGroupPercent[0] = perc_inf * L1 / L_inf;
    ageGroupPercent[1] = perc_inf - ageGroupPercent[0];

    double result = 0.0;
    for (size_t i = 0; i < ngroups - 1; i++) {
        double residual = log(pred[i]) - log(ageGroupPercent[i]);
        result += residual * residual;
    }
    return result;
}

static double wCalcRSS(const gsl_vector* v, void* params) {
    double x, y;
    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);

    // params: the run (see minimizeCalc_rss)
    return setDemoParameters(*static_cast<RunContext*>(params), x, y);
}
/** Presumably from AR. Searches iteratively for a minimum to function func.
 *
 * @param func A function pointer to the function used for ...
 * @param param1 Initial guess of first param taken by func
 * @param param2 Initial guess of second
 */
void minimizeCalc_rss(RunContext& ctx, double param1, double param2) {
    gsl_vector* stepSize = gsl_vector_alloc(2);
    gsl_vector_set_all(stepSize, 0.1);

    gsl_vector* initialValues = gsl_vector_alloc(2);
    gsl_vector_set(initialValues, 0, param1);
    gsl_vector_set(initialValues, 1, param2);

    gsl_multimin_function minex_func;
    minex_func.f = &wCalcRSS;
    minex_func.params = &ctx;
    minex_func.n = 2;

    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(T, 2);
    gsl_multimin_fminimizer_set(minimizer, &minex_func, initialValues, stepSize);

    for (size_t iter = 0; iter < 100; ++iter) {
        if (gsl_multimin_fminimizer_iterate(minimizer))
            throw TRACED_EXCEPTION("gsl_multimin_fminimizer_iterate failed", Error::GSL);

        double size = gsl_multimin_fminimizer_size(minimizer);
        int status = gsl_multimin_test_size(size, 1e-2);
        if (status == GSL_SUCCESS) break;
    }
    // Call again to set final value.
    wCalcRSS(minimizer->x, &ctx);

    gsl_vector_free(initialValues);
    gsl_vector_free(stepSize);
    gsl_multimin_fminimizer_free(minimizer);
}

void AgeStructure::calcCumAgeProp(RunContext& ctx) {
    auto& cumAgeProp = ctx.age.cumAgeProp;
    const double mu0 = ctx.age.mu0, mu1 = ctx.age.mu1;
    const double alpha0 = ctx.age.alpha0, alpha1 = ctx.age.alpha1;
    constexpr double rho = AgeStructureState::rho;
    cumAgeProp[0] = 0.0;
    for (size_t j = 1; j < cumAgeProp.size(); ++j) {
        SimTime age = ctx.clock.fromTS(cumAgeProp.size() - j - 1);
        double ageYears = inYears(age);
        double M1s = (mu0 * (1.0 - exp(-alpha0 * ageYears)) / alpha0);
        double M2s = (mu1 * (exp(alpha1 * ageYears) - 1.0) / alpha1);
        double Ms = M1s + M2s;
        double predperc = exp(-rho * ageYears - Ms);
        if (age >= ctx.clock.max_human_age) predperc = 0.0;
        cumAgeProp[j] = cumAgeProp[j - 1] + predperc;
    }
    double totalCumPC = cumAgeProp[cumAgeProp.size() - 1];
    for (size_t j = 1; j < cumAgeProp.size(); ++j) {
        //Scale using the total cumAgeProp
        cumAgeProp[j] = cumAgeProp[j] / totalCumPC;
    }
}

}
