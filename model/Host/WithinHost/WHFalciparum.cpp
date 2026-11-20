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

#include "RunContext.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/Genotypes.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/UnitParse.h"

#include <cmath>
#include <gsl/gsl_cdf.h>

namespace OM {
namespace WithinHost {

using namespace OM::util;

// -----  private model parameters  -----

//Standard dev innate immunity for densities
double sigma_i;

// contribution of parasite densities to acquired immunity in the presence of fever
double immPenalty_22;

// Remaining immunity against asexual parasites(after time step, each of 2 components y and h)
// This variable decays the effectors imm_h/imm_Y columns in a way that their
// effects on densities (1-Dh and 1-Dy) decay exponentially.
double asexImmRemain;

// Remaining immunity against asexual parasites(after each time step, each of 2 components y and h)
// This variable decays the imm_h/imm_Y columns exponentially.
double immEffectorRemain;

/// Critical value for immunity trigger (cumulative densities)
double invCumulativeYstar;
/// Critical value for immunity trigger (cumulative inoculations)
double invCumulativeHstar;

/// Maternal protection at birth
double alpha_m;

/// More or less (up to 0.693) inverse quantity of alphaMStar (AJTM p. 9 eq. 12),
/// decay rate of maternal protection in years^(-1).
double decayM;

SimTime Infection::s_latentP = sim::never();
int WHInterface::y_lag_len = 0;

// -----  static functions  -----

void WHInterface::setParams(double cumYStar, double cumHStar, double aM, double dM) {
    invCumulativeYstar = cumYStar;
    invCumulativeHstar = cumHStar;
    alpha_m = aM;
    decayM = dM;
}

// -----  Non-static  -----

void WHInterface::ctorFalciparum(LocalRng& rng, double comorbidityFactor) {
    pathogenesisModel.reset(Pathogenesis::PathogenesisModel::createPathogenesisModel(comorbidityFactor));

    // NOTE: negating a Gaussian sample with mean 0 is pointless — except that
    // the individual samples change. In any case the overhead is negligible.
    //FIXME: Should this be allowed to be greater than 1?
    // Oldest code on GoogleCode: _innateImmunity=(double)(W_GAUSS((0), (sigma_i)));
    // parity: order — local RNG stream: draw call order fixes stream state
    _innateImmSurvFact = exp(-rng.gauss(0.0, sigma_i));

    m_y_lag_i.resize(y_lag_len * ctx.genotypes.N());
    m_y_lag_l.resize(y_lag_len * ctx.genotypes.N());
}

double WHInterface::immunitySurvivalFactor(double ageInYears, double cumulativeExposureJ) {
    // Row reads hoisted to locals: one cols->deref per call instead of four
    // (nothing mutates them in here, so the values are bit-identical).
    const double imm_h = imm_h_ref(), imm_Y = imm_Y_ref();
    if (std::isnan(ageInYears) || std::isnan(cumulativeExposureJ) || std::isnan(imm_h) || std::isnan(imm_Y)) {
        throw base_exception("nan in immunitySurvivalFactor");
    }

    // Documentation: AJTMH pp22-23

    // Effect of cumulative Parasite density (named Dy in AJTM)
    double dY = 1.0;
    // Effect of number of infections experienced since birth (named Dh in AJTM)
    double dH = 1.0;
    if (imm_h > 1.0) {
        dH = 1.0 / (1.0 + (imm_h - 1.0) * invCumulativeHstar);
        dY = 1.0 / (1.0 + (imm_Y - cumulativeExposureJ) * invCumulativeYstar);
    }

    // Effect of age-dependent maternal immunity (named Dm in AJTM)
    double dA = 1.0 - alpha_m * exp(-decayM * ageInYears);

    return util::streamValidate(std::min(dY * dH * dA, 1.0));
}

WHInterface::ImmSurvTick WHInterface::immSurvivalTick(double ageInYears) {
    // p7c-C3: hoisted host part of immunitySurvivalFactor (see WHInterface.h).
    // Each value below equals the per-infection call's operand bit-for-bit.
    const double imm_h = imm_h_ref(), imm_Y = imm_Y_ref();
    if (std::isnan(ageInYears) || std::isnan(imm_h) || std::isnan(imm_Y))
        throw base_exception("nan in immunitySurvivalFactor");
    ImmSurvTick t;
    t.active = imm_h > 1.0;
    t.immY = imm_Y;
    t.invY = invCumulativeYstar;
    t.dH = t.active ? 1.0 / (1.0 + (imm_h - 1.0) * invCumulativeHstar) : 1.0;
    t.dA = 1.0 - alpha_m * exp(-decayM * ageInYears);
    return t;
}

double WHInterface::immunitySurvivalFactor(const ImmSurvTick& t, double cumulativeExposureJ) {
    if (std::isnan(cumulativeExposureJ))
        throw base_exception("nan in immunitySurvivalFactor");
    double dY = 1.0;
    if (t.active)
        dY = 1.0 / (1.0 + (t.immY - cumulativeExposureJ) * t.invY);
    return util::streamValidate(std::min(dY * t.dH * t.dA, 1.0));
}

// Infectiousness parameters: see AJTMH p.33; tau=1/sigmag**2
const double PTM_beta1 = 1.0;
const double PTM_beta2 = 0.46;
const double PTM_beta3 = 0.17;
const double PTM_tau = 0.066;
const double PTM_tau_prime = 1.0 / sqrt(1.0 / PTM_tau);
const double PTM_mu = -8.1;

inline double compute_y_lag(const SimClock& clock, const Genotypes& genotypes, const vector<double>& y_lag,
                            int y_lag_len, vector<double>& y_lag_g) {
    const size_t n = genotypes.N();

    // Take weighted sum of total asexual blood stage density 10, 15 and 20 days
    // before. Add y_lag_len to index to ensure positive.
    size_t d10 = mod_nn(y_lag_len + clock.inSteps(clock.ts1() - sim::fromDays(10)), y_lag_len);
    size_t d15 = mod_nn(y_lag_len + clock.inSteps(clock.ts1() - sim::fromDays(15)), y_lag_len);
    size_t d20 = mod_nn(y_lag_len + clock.inSteps(clock.ts1() - sim::fromDays(20)), y_lag_len);

    double y10_sum = 0.0, y15_sum = 0.0, y20_sum = 0.0;
    for (size_t g = 0; g < n; ++g) {
        double y10 = y_lag[d10 * n + g];
        double y15 = y_lag[d15 * n + g];
        double y20 = y_lag[d20 * n + g];
        y10_sum += y10;
        y15_sum += y15;
        y20_sum += y20;
        y_lag_g[g] = PTM_beta1 * y10 + PTM_beta2 * y15 + PTM_beta3 * y20;
    }

    // Weighted sum:
    return PTM_beta1 * y10_sum + PTM_beta2 * y15_sum + PTM_beta3 * y20_sum;
}

double WHInterface::probTransmissionToMosquito_falciparum(vector<double>& probTransGenotype_i,
                                                          vector<double>& probTransGenotype_l) const {
    // This model (often referred to as the gametocyte model) was designed for
    // 5-day time steps. We use the same model (sampling 10, 15 and 20 days
    // ago) for 1-day time steps to avoid having to design and analyse a new
    // model. Description: AJTMH pp.32-33 and p9.

    // Note: we don't allow for gametocydal treatments (e.g. Primaquine).
    const size_t n = ctx.genotypes.N();

    // Sum lagged densities across genotypes:
    // ponytail: reused scratch (83M/2 allocs-per-host-tick mallocs eliminated). Value-identical:
    // compute_y_lag overwrites all n entries before any read; nothing escapes; n fixed per run
    // (genotypes.N()) and TwoRunsSuite (two runs, one process) proves nothing leaks across runs.
    static std::vector<double> y_lag_g_i, y_lag_g_l;
    y_lag_g_i.resize(n);
    y_lag_g_l.resize(n);
    const double y_lag_sum_i = compute_y_lag(ctx.clock, ctx.genotypes, m_y_lag_i, y_lag_len, y_lag_g_i);
    const double y_lag_sum_l = compute_y_lag(ctx.clock, ctx.genotypes, m_y_lag_l, y_lag_len, y_lag_g_l);
    const double y_lag_sum = y_lag_sum_i + y_lag_sum_l;

    if (y_lag_sum < 0.001)
        return 0.0; // cut off for uninfectious humans

    // Get a zval, convert to equivalent Normal sample:
    const double zval = (log(y_lag_sum) + PTM_mu) * PTM_tau_prime;
    // parity: order — pure GSL CDF/quantile call (no stream consumption); keep the identical GSL function (F2)
    const double pone = gsl_cdf_ugaussian_P(zval);
    double pTransmit = pone * pone;

    // pTransmit has to be between 0 and 1:
    pTransmit = std::max(pTransmit, 0.0);
    pTransmit = std::min(pTransmit, 1.0);

    if (pTransmit <= 0.0)
        return pTransmit;

    for (size_t g = 0; g < n; ++g) {
        probTransGenotype_i[g] = pTransmit * y_lag_g_i[g] / y_lag_sum;
        probTransGenotype_l[g] = pTransmit * y_lag_g_l[g] / y_lag_sum;
    }

    // Include here the effect of transmission-blocking vaccination:
    util::streamValidate(pTransmit);
    return pTransmit;
}

bool WHInterface::diagnosticResult_falciparum(LocalRng& rng, const Diagnostic& diagnostic) const {
    return diagnostic.isPositive(rng, totalDensity, hrp2Density);
}

void WHInterface::treatment_falciparum(Host::Human& human, TreatmentId treatId) {
    const Treatments& treat = Treatments::select(ctx, treatId);
    treatSimple(human, treat.timeLiver, treat.timeBlood);

    // triggered intervention deployments:
    treat.deploy(ctx, human, mon::Deploy::TREAT, interventions::VaccineLimits(/*default initialise: no limits*/));
}
bool WHInterface::treatSimple_falciparum(Host::Human& human, SimTime timeLiver, SimTime timeBlood) {
    if (timeLiver != sim::zero()) {
        if (timeLiver < sim::zero())
            clearInfections(Treatments::LIVER);
        else
            treatExpiryLiver = max(int(treatExpiryLiver), int(ctx.clock.nowOrTs1()) + timeLiver);
        mon::recordEvent(mon::measure("nLiverStageTreatments"), human);
    }
    if (timeBlood != sim::zero()) {
        if (timeBlood < sim::zero())
            clearInfections(Treatments::BLOOD);
        else
            treatExpiryBlood = max(int(treatExpiryBlood), int(ctx.clock.nowOrTs1()) + timeBlood);
        return true; // blood stage treatment
    }
    return false; // no blood stage treatment
}

Pathogenesis::StatePair WHInterface::determineMorbidity_falciparum(Host::Human& human, double ageYears, bool isDoomed) {
    Pathogenesis::StatePair result =
        pathogenesisModel->determineState(human, ageYears, timeStepMaxDensity, totalDensity, isDoomed);

    /* Note: this model can easily be re-enabled, but is not used and not considered to be a good model.
    if( (result.state & Pathogenesis::MALARIA) && util::ModelOptions::option( util::PENALISATION_EPISODES ) ){
        // This does immunity penalisation:
        imm_Y_ref() = imm_Y_lag_ref() - immPenalty_22*(imm_Y_ref()-imm_Y_lag_ref());
        if (imm_Y_ref() < 0) {
            imm_Y_ref()=0.0;
        }
    }*/

    return result;
}

// -----  immunity  -----

void WHInterface::updateImmuneStatus() {
    // Same sequential updates on row-locals (bit-identical), one deref each side.
    double h = imm_h_ref(), Y = imm_Y_ref();
    if (immEffectorRemain < 1) {
        h *= immEffectorRemain;
        Y *= immEffectorRemain;
    }
    if (asexImmRemain < 1) {
        h *= asexImmRemain / (1 + (h * (1 - asexImmRemain) * invCumulativeHstar));
        Y *= asexImmRemain / (1 + (Y * (1 - asexImmRemain) * invCumulativeYstar));
    }
    imm_h_ref() = h;
    imm_Y_ref() = Y;
    imm_Y_lag_ref() = Y;
}

// -----  Checkpointing  -----

}
}
