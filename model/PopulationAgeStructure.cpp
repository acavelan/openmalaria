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

#include "PopulationAgeStructure.h"
#include "util/errors.h"
#include <gsl/gsl_multimin.h>

namespace OM {
namespace {
constexpr double rho = 0.0;
struct DemographyFit {
    vector<double> ageGroupBounds, ageGroupPercent;
    double mu0, mu1, alpha0, alpha1;
    double residual(double param1, double param2);
};

double DemographyFit::residual (double param1, double param2)
{
    const double IMR = 0.1;
    double M_inf = -log (1 - IMR);

    mu1 = exp (param1) / 100;
    alpha1 = exp (param2) / 100;
    alpha0 = 4.0;
    mu0 = (M_inf - mu1 * (exp (alpha1 * 0.5) - 1) * alpha0) /
    (alpha1 * (1 - exp (-alpha0 * 0.5)));

    size_t ngroups = ageGroupPercent.size();
    double sumpred = 0.0;
    vector<double> M(ngroups,0);
    vector<double> pred(ngroups,0);

    for(size_t i = 0; i < ngroups - 1; i++) {
	double midpt = (ageGroupBounds[i+1] + ageGroupBounds[i]) * 0.5;
	double M1 = mu0 * (1.0 - exp (-alpha0 * midpt)) / alpha0;
	double M2 = mu1 * (exp (alpha1 * midpt) - 1.0) / alpha1;
	M[i]  = M1 + M2;
	pred[i] = (ageGroupBounds[i+1] - ageGroupBounds[i])
	* exp (-rho * midpt - M[i]);
	sumpred += pred[i];
    }
    for(size_t i = 0; i < ngroups - 1; i++) {
	pred[i] = pred[i] / sumpred * 100.0;
    }
    double L_inf = exp (-rho * 0.5 - M[1]);
    double M_nn  = -log (1.0 - 0.4 * (1 - exp (-M[1])));
    double L1    = 1.0 / 12.0 * exp (-rho / 24.0 - M_nn);
    double perc_inf = ageGroupPercent[0] + ageGroupPercent[1];
    ageGroupPercent[0] = perc_inf * L1 / L_inf;
    ageGroupPercent[1] = perc_inf - ageGroupPercent[0];

    double result = 0.0;
    for(size_t i = 0; i < ngroups - 1; i++) {
	double residual = log (pred[i]) - log (ageGroupPercent[i]);
	result += residual * residual;
    }
    return result;
}
static vector<double> cumulativeProportions(const DemographyFit& fit, SimTime maxAge, SimTime interval)
{
    vector<double> cumAgeProp(maxAge / interval + 1);
    cumAgeProp[0] = 0.0;
    for(size_t j=1;j < cumAgeProp.size(); ++j) {
	SimTime age = interval * int(cumAgeProp.size() - j - 1);
	double ageYears = sim::inYears(age);
	double M1s = (fit.mu0 * (1.0 - exp (-fit.alpha0 * ageYears)) / fit.alpha0);
	double M2s = (fit.mu1 * (exp (fit.alpha1 * ageYears) - 1.0) / fit.alpha1);
	double Ms = M1s + M2s;
	double predperc = exp (-rho * ageYears - Ms);
	if (age >= maxAge) {
	    predperc = 0.0;
	}
	cumAgeProp[j] = cumAgeProp[j-1] + predperc;
    }
    double totalCumPC = cumAgeProp[cumAgeProp.size()-1];
    for(size_t j=1;j < cumAgeProp.size(); ++j) {
	//Scale using the total cumAgeProp
	cumAgeProp[j] = cumAgeProp[j] / totalCumPC;
    }
    return cumAgeProp;
}

static double evaluateFit(const gsl_vector* values, void* context) {
    return static_cast<DemographyFit*>(context)->residual(
        gsl_vector_get(values, 0), gsl_vector_get(values, 1));
}
}

vector<double> makeAgeStructure(const vector<DemographyAgeGroup>& groups, SimTime maxAge, SimTime interval) {
    if (groups.empty() || interval <= 0 || maxAge < interval)
        throw util::xml_scenario_error("age structure requires age groups, a positive interval, and maximum age >= interval");
    DemographyFit fit;
    double sumperc = 0.0;

    size_t ngroups = groups.size() + 1;
    fit.ageGroupBounds.resize(ngroups+1,0);
    fit.ageGroupPercent.resize(ngroups,0);

    //Add age group for first month of life
    fit.ageGroupBounds[0] = 0.0;
    fit.ageGroupBounds[1] = 1.0 / 12.0;
    fit.ageGroupPercent[0] = 0.0;
    for(size_t i = 1;i < ngroups; i++) {
	fit.ageGroupBounds[i+1] = groups[i-1].upperBound;
	fit.ageGroupPercent[i] = groups[i-1].percentage;
	sumperc += fit.ageGroupPercent[i];
    }
    sumperc = 100.0 / sumperc; // multiplier to get percentages
    for(size_t i = 0;i < ngroups; i++) {
	fit.ageGroupPercent[i]  = fit.ageGroupPercent[i] * sumperc;
    }

    double stepValues[] = {0.1, 0.1};
    double initialValues[] = {0.371626412, 0.841209593};
    auto steps = gsl_vector_view_array(stepValues, 2);
    auto initial = gsl_vector_view_array(initialValues, 2);
    gsl_multimin_function objective;
    objective.f = &evaluateFit;
    objective.params = &fit;
    objective.n = 2;
    unique_ptr<gsl_multimin_fminimizer, decltype(&gsl_multimin_fminimizer_free)> minimizer(
        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, 2),
        &gsl_multimin_fminimizer_free);
    gsl_multimin_fminimizer_set(minimizer.get(), &objective, &initial.vector, &steps.vector);
    for (size_t iter = 0; iter < 100; ++iter) {
        if (gsl_multimin_fminimizer_iterate(minimizer.get()))
            throw TRACED_EXCEPTION("gsl_multimin_fminimizer_iterate failed", util::Error::GSL);
        double size = gsl_multimin_fminimizer_size(minimizer.get());
        if (gsl_multimin_test_size(size, 1e-2) == GSL_SUCCESS) break;
    }
    // The final evaluation also sets the fitted parameters used below.
    evaluateFit(minimizer->x, &fit);
    return cumulativeProportions(fit, maxAge, interval);
}

int targetCumPop(const vector<double>& proportions, size_t ageSteps, int targetPopulation) {
    size_t index = proportions.size() - 1 - ageSteps;
    return (int) floor(proportions[index] * targetPopulation + 0.5);
}
}
