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

/* T2/T3 sampler battery: for every distribution in model/util/random.h's API
 * (plus the direct-GSL binomial/exponential the RNG layer pins), draw 50k
 * samples from a seeded stream and check, against GSL's own CDFs / exact
 * moments:
 *   (a) a one-sample Kolmogorov-Smirnov test with p > 0.005;
 *   (b) sample mean and std within a 5-sigma z-test of theory (kurtosis-aware
 *       SE for the std, so heavy-tailed regimes aren't false-failed).
 * Regime boundaries are covered: gamma shape 0.5/2.0, Poisson mu 3.2/25.0,
 * binomial near p=0 and p=1, beta with both (a,b) <=1 and both >1.
 *
 * The seeds are fixed, so every assertion is deterministic; the battery locks
 * GSL sampler-call semantics and must stay green when phase 9 ports samplers.
 * scipy is not available, hence the in-file KS p-value (Numerical Recipes
 * ksccdf series). Runtime is well under the 20 s budget.
 */

#include <cxxtest/TestSuite.h>
#include "util/random.h"
#include "golden-gen/RngGoldenAdapter.h" // for binomial/exponential (no random.h wrapper)
#include <gsl/gsl_cdf.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace OM;
using util::LocalRng;

namespace {

constexpr size_t N = 50000;
const double KS_P_FLOOR = 0.005; // fail if p <= this
const double Z_LIMIT = 5.0;      // fail if |z| >= this

/// One-sample KS p-value for sorted samples vs a CDF (continuous or discrete).
/// Scans unique sample values v_j: sup|F̂−F| is attained at the steps, where
/// both functions are evaluated at the same x — the textbook d-minus from
/// continuous tables (F(x_i) − i/n) is NOT valid for step CDFs (it pairs F at
/// v_j with F̂ below v_j, different x).
template <class F> double ks_pvalue(const std::vector<double>& sorted, F cdf, bool discrete = false) {
    const double n = static_cast<double>(sorted.size());
    const double sqn = std::sqrt(n);
    double d = 0.0;
    double prevAvg = 0.0; // F̂ just below the current unique value
    size_t i = 0;
    while (i < sorted.size()) {
        size_t j = i;
        while (j < sorted.size() && sorted[j] == sorted[i])
            ++j;
        const double Fv = cdf(sorted[i]);
        const double avg = static_cast<double>(j) / n; // F̂ at v_j (inclusive)
        d = std::max(d, std::fabs(avg - Fv));
        // continuous: inside [v_{j-1}, v_j) the sup |F̂−F| is approached at v_j
        if (!discrete)
            d = std::max(d, std::fabs(prevAvg - Fv));
        prevAvg = avg;
        i = j;
    }
    const double lambda = (sqn + 0.12 + 0.11 / sqn) * d;
    if (lambda <= 0.0)
        return 1.0;
    // Q_KS(lambda) = 2 sum_{j>=1} (-1)^{j-1} exp(-2 j^2 lambda^2)
    double sum = 0.0, fac = -2.0 * lambda * lambda;
    for (int j = 1; j <= 100; ++j) {
        const double term = std::exp(fac * static_cast<double>(j) * j);
        sum += (j & 1) ? term : -term;
        if (term < 1e-12 * std::abs(sum) && j > 1)
            break;
    }
    double p = 2.0 * sum;
    return p < 0.0 ? 0.0 : (p > 1.0 ? 1.0 : p);
}

/// Mean/std z-tests. sigmaTheory = sqrt(varTheory); kappa is the
/// distribution's raw excess kurtosis + 3 (SE of s: sqrt((kappa-1)/(4(n-1)))·sigma;
/// equals the usual sigma/sqrt(2(n-1)) at kappa=3).
void check_mean_std(const std::vector<double>& x, double mu, double varTheory, double kappa, const char* what) {
    const double n = static_cast<double>(x.size());
    double m = 0.0;
    for (double v : x)
        m += v;
    m /= n;
    double m2 = 0.0;
    for (double v : x)
        m2 += (v - m) * (v - m);
    const double s = std::sqrt(m2 / (n - 1.0)); // sample std
    const double sigma = std::sqrt(varTheory);
    const double zMean = (m - mu) / std::sqrt(varTheory / n);
    const double zStd = (s - sigma) * 2.0 * std::sqrt(n - 1.0) / (sigma * std::sqrt(kappa - 1.0));
    char msg[256];
    std::snprintf(msg, sizeof(msg), "%s: |z_mean|=%.2f >= %.1f (mean=%.6g theory=%.6g)", what, std::fabs(zMean),
                  Z_LIMIT, m, mu);
    if (!(std::fabs(zMean) < Z_LIMIT))
        TS_FAIL(msg);
    std::snprintf(msg, sizeof(msg), "%s: |z_std|=%.2f >= %.1f (std=%.6g theory=%.6g)", what, std::fabs(zStd), Z_LIMIT,
                  s, sigma);
    if (!(std::fabs(zStd) < Z_LIMIT))
        TS_FAIL(msg);
}

/// Weibull excess kurtosis+3 from raw moments (r_j = Γ(1+j/k)λ^j).
double weibull_kappa(double lam, double k) {
    double r[5];
    for (int j = 0; j <= 4; ++j)
        r[j] = std::exp(std::lgamma(1.0 + j / k)) * std::pow(lam, j);
    const double m = r[1];
    const double c2 = r[2] - m * m;
    const double c4 = r[4] - 4 * m * r[3] + 6 * m * m * r[2] - 3 * std::pow(m, 4.0);
    return c4 / (c2 * c2);
}

/// Assert KS p > floor for sorted samples vs cdf.
template <class F> void check_ks(const std::vector<double>& sorted, F cdf, const char* what, bool discrete = false) {
    const double p = ks_pvalue(sorted, cdf, discrete);
    char msg[256];
    std::snprintf(msg, sizeof(msg), "%s: KS p=%.5g <= %.3g (D too large)", what, p, KS_P_FLOOR);
    if (!(p > KS_P_FLOOR))
        TS_FAIL(msg);
}

} // anonymous namespace

class SamplerBatterySuite : public CxxTest::TestSuite {
  public:
    void testUniform01() {
        LocalRng rng(11, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.uniform_01());
        std::sort(x.begin(), x.end());
        check_ks(x, [](double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); }, "uniform_01");
        check_mean_std(x, 0.5, 1.0 / 12.0, 1.8, "uniform_01");
    }

    void testGauss() {
        const double mu = -1.5, sig = 2.5;
        LocalRng rng(12, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.gauss(mu, sig));
        std::sort(x.begin(), x.end());
        check_ks(x, [&](double v) { return gsl_cdf_ugaussian_P((v - mu) / sig); }, "gauss");
        check_mean_std(x, mu, sig * sig, 3.0, "gauss");
    }

    void testGammaSmallShape() {
        const double a = 0.5, b = 2.0; // GSL small-shape branch
        LocalRng rng(13, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.gamma(a, b));
        std::sort(x.begin(), x.end());
        check_ks(x, [&](double v) { return gsl_cdf_gamma_P(v, a, b); }, "gamma(0.5,2.0)");
        check_mean_std(x, a * b, a * b * b, 3.0 + 6.0 / a, "gamma(0.5,2.0)");
    }

    void testGammaLargeShape() {
        const double a = 2.0, b = 1.0; // Marsaglia-Tsang branch
        LocalRng rng(14, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.gamma(a, b));
        std::sort(x.begin(), x.end());
        check_ks(x, [&](double v) { return gsl_cdf_gamma_P(v, a, b); }, "gamma(2.0,1.0)");
        check_mean_std(x, a * b, a * b * b, 3.0 + 6.0 / a, "gamma(2.0,1.0)");
    }

    void testBetaBothBelow1() {
        const double a = 0.5, b = 0.7; // Cheng BA branch
        LocalRng rng(15, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.beta(a, b));
        std::sort(x.begin(), x.end());
        check_ks(x, [&](double v) { return gsl_cdf_beta_P(v, a, b); }, "beta(0.5,0.7)");
        check_mean_std(x, a / (a + b), a * b / ((a + b) * (a + b) * (a + b + 1.0)),
                       3.0 + 6.0 * (((a - b) * (a - b) * (a + b + 1.0) - a * b * (a + b + 2.0)) /
                                    (a * b * (a + b + 2.0) * (a + b + 3.0))),
                       "beta(0.5,0.7)");
    }

    void testBetaBothAbove1() {
        const double a = 2.0, b = 3.0; // Cheng BB/BC branch
        LocalRng rng(16, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.beta(a, b));
        std::sort(x.begin(), x.end());
        check_ks(x, [&](double v) { return gsl_cdf_beta_P(v, a, b); }, "beta(2.0,3.0)");
        check_mean_std(x, a / (a + b), a * b / ((a + b) * (a + b) * (a + b + 1.0)),
                       3.0 + 6.0 * (((a - b) * (a - b) * (a + b + 1.0) - a * b * (a + b + 2.0)) /
                                    (a * b * (a + b + 2.0) * (a + b + 3.0))),
                       "beta(2.0,3.0)");
    }

    void testLogNormal() {
        const double mlog = 0.4, slog = 1.1;
        const double e2 = std::exp(slog * slog);
        LocalRng rng(17, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.log_normal(mlog, slog));
        std::sort(x.begin(), x.end());
        check_ks(x, [&](double v) { return gsl_cdf_lognormal_P(v, mlog, slog); }, "log_normal(0.4,1.1)");
        check_mean_std(x, std::exp(mlog + slog * slog / 2.0), (e2 - 1.0) * std::exp(2.0 * mlog + slog * slog),
                       e2 * e2 * e2 * e2 + 2.0 * e2 * e2 * e2 + 3.0 * e2 * e2 - 3.0, // e^{4s²}+2e^{3s²}+3e^{2s²}−3
                       "log_normal(0.4,1.1)");
    }

    void testPoissonSmallMu() {
        const double mu = 3.2; // direct method branch (mu < 10)
        LocalRng rng(18, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.poisson(mu));
        std::sort(x.begin(), x.end());
        check_ks(
            x, [&](double v) { return gsl_cdf_poisson_P(static_cast<unsigned int>(v), mu); }, "poisson(3.2)", true);
        check_mean_std(x, mu, mu, 3.0 + 1.0 / mu, "poisson(3.2)");
    }

    void testPoissonLargeMu() {
        const double mu = 25.0; // rejection branch
        LocalRng rng(19, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.poisson(mu));
        std::sort(x.begin(), x.end());
        check_ks(
            x, [&](double v) { return gsl_cdf_poisson_P(static_cast<unsigned int>(v), mu); }, "poisson(25.0)", true);
        check_mean_std(x, mu, mu, 3.0 + 1.0 / mu, "poisson(25.0)");
    }

    void testWeibull() {
        const double lam = 2.0, k = 1.5;
        const double g1 = std::tgamma(1.0 + 1.0 / k), g2 = std::tgamma(1.0 + 2.0 / k);
        LocalRng rng(20, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.weibull(lam, k));
        std::sort(x.begin(), x.end());
        check_ks(x, [&](double v) { return gsl_cdf_weibull_P(v, lam, k); }, "weibull(2.0,1.5)");
        check_mean_std(x, lam * g1, lam * lam * (g2 - g1 * g1), weibull_kappa(lam, k), "weibull(2.0,1.5)");
    }

    // binomial/exponential have no random.h wrapper: drawn through the
    // GoldenRng adapter (same Xoshiro256P primitives, direct GSL calls), which
    // is what rng_dist_v50.csv pins.
    void testBinomialNearZero() {
        const unsigned int n = 100;
        const double p = 0.01;
        GoldenRng::Rng rng(1, 2, 3, 4);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.binomial(n, p));
        std::sort(x.begin(), x.end());
        check_ks(
            x, [&](double v) { return gsl_cdf_binomial_P(static_cast<unsigned int>(v), p, n); }, "binomial(100,0.01)",
            true);
        check_mean_std(x, n * p, n * p * (1 - p), 3.0 + (1.0 - 6.0 * p * (1 - p)) / (n * p * (1 - p)),
                       "binomial(100,0.01)");
    }

    void testBinomialNearOne() {
        const unsigned int n = 100;
        const double p = 0.99;
        GoldenRng::Rng rng(5, 6, 7, 8);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.binomial(n, p));
        std::sort(x.begin(), x.end());
        check_ks(
            x, [&](double v) { return gsl_cdf_binomial_P(static_cast<unsigned int>(v), p, n); }, "binomial(100,0.99)",
            true);
        check_mean_std(x, n * p, n * p * (1 - p), 3.0 + (1.0 - 6.0 * p * (1 - p)) / (n * p * (1 - p)),
                       "binomial(100,0.99)");
    }

    void testExponential() {
        const double scale = 1.7;
        GoldenRng::Rng rng(9, 10, 11, 12);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.exponential(scale));
        std::sort(x.begin(), x.end());
        check_ks(x, [&](double v) { return gsl_cdf_exponential_P(v, scale); }, "exponential(1.7)");
        check_mean_std(x, scale, scale * scale, 9.0, "exponential(1.7)");
    }

    void testMaxMultiLogNormal() {
        // start=0 → max() never binds; the sample is exp(mlog + slog Φ⁻¹(U^{1/n})),
        // i.e. the max of n iid lognormal draws: CDF = Φ(z)^n. KS only (no
        // closed-form moments for the max).
        const int n = 5;
        const double mlog = 0.4, slog = 1.1;
        LocalRng rng(21, 20261119);
        std::vector<double> x;
        x.reserve(N);
        for (size_t i = 0; i < N; ++i)
            x.push_back(rng.max_multi_log_normal(0.0, n, mlog, slog));
        std::sort(x.begin(), x.end());
        check_ks(
            x,
            [&](double v) {
                double z = (std::log(v) - mlog) / slog;
                double f = gsl_cdf_ugaussian_P(z);
                return std::pow(f, static_cast<double>(n));
            },
            "max_multi_log_normal(0,5,0.4,1.1)");
    }
};
