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

/* Minimal adapter around the model's own RNG primitives, used by
 * golden-gen/rng_goldens_gen.cpp (CSV generator) and RngGoldenSuite.h (replay
 * gate). It mirrors model/util/random.h's GSL plumbing and distribution
 * wrappers verbatim (each mirrored block names its random.h counterpart) around a
 * directly 4-word-seedable Xoshiro256P — RNG<> (util/random.h) only exposes
 * 2-word seeding, and its state is private, while Xoshiro256P's public
 * binary_checkpoint (RNG state API in util/random.h) is our state-word reader.
 *
 * It is deliberately its own namespace (random.h is included for the engine only, never RNG<>) so it can be
 * compiled standalone (model + GSL headers only, no model library).
 */

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>

// The Xoshiro256P engine (util/random.h) uses unqualified ostream/istream/runtime_error;
// the model gets these via Clock.h. Provide
// the same before including it:
using std::istream;
using std::ostream;
using std::runtime_error;

#include "util/random.h"  // for the Xoshiro256P engine only
#include "util/errors.h"
// random.h opens namespace OM::util, so its bare Error:: resolves there;
// alias the same namespace here so the mirrored poisson wrapper
// stays character-identical:
namespace Error = ::OM::util::Error;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

namespace GoldenRng {

// --- mirrored from model/util/random.h (support functions) ---
template<class T>
long unsigned int sample_ulong (void *ptr) {
    T *rng = reinterpret_cast<T*>(ptr);
    return rng->gen_u32();
}
template<class T>
double sample_double01 (void *ptr) {
    T *rng = reinterpret_cast<T*>(ptr);
    return rng->gen_double();
}

// --- mirrored from model/util/(make_gsl_rng_type) — random.h ---
template<class T>
gsl_rng_type make_gsl_rng_type(T& rng) {
    return {
        "OM_RNG",		// name
        std::numeric_limits<uint32_t>::max(),
        std::numeric_limits<uint32_t>::min(),
        0,				// size of state; not used here
        nullptr,			// re-seed function; don't use
        &sample_ulong<T>,
        &sample_double01<T>
    };
}

/// Same construction as OM::util::RNG<Xoshiro256P>, but with
/// direct 4-word seeding and a state-word reader.
struct Rng {
    Xoshiro256P m_rng;
    
    // Hooks for GSL distributions
    gsl_rng_type m_gsl_type;
    gsl_rng m_gsl_gen;
    
    /// ctor body mirrored from RNG's seeding ctor
    Rng(uint64_t a, uint64_t b, uint64_t c, uint64_t d): m_rng(a, b, c, d) {
        m_gsl_type = make_gsl_rng_type(m_rng);
        m_gsl_gen.type = &m_gsl_type;
        m_gsl_gen.state = reinterpret_cast<void*>(&m_rng);
    }
    
    /// read the 4 Xoshiro256P state words via the model's own public
    /// checkpoint API (Xoshiro256P::binary_checkpoint, util/random.h)
    void state(uint64_t out[4]) {
        std::ostringstream o;
        m_rng.binary_checkpoint(o);
        std::string s = o.str();
        std::memcpy(out, s.data(), 4 * sizeof(uint64_t));
    }
    
    // --- distribution wrappers, mirrored from random.h (line cites below) ---
    
    /** Generate a random number in the range [0,1). */
    // mirrored from RNG's uniform_01
    inline double uniform_01 () {
        return m_rng.gen_double();
    }
    
    /** Gaussian with mean mean and standard deviation std (polar Box-Muller). */
    // mirrored from model/util/random.h
    double gauss (double mean, double std){
        return gsl_ran_gaussian(&m_gsl_gen,std)+mean;
    }
    
    /** Gamma distribution variate. */
    // mirrored from model/util/random.h
    double gamma (double a, double b){
        return gsl_ran_gamma(&m_gsl_gen, a, b);
    }
    
    /** Lognormal distribution variate. */
    // mirrored from model/util/random.h
    double log_normal (double meanlog, double stdlog){
        return gsl_ran_lognormal (&m_gsl_gen, meanlog, stdlog);
    }
    
    /** Maximum over multiple log-normal samples. */
    // mirrored from model/util/random.h (body)
    double max_multi_log_normal (double start, int n, double meanlog, double stdlog){
        double normp = pow( uniform_01(), 1.0 / n );
        double zval = gsl_cdf_ugaussian_Pinv (normp);
        double multi_sample = exp(meanlog + stdlog * zval);
        return std::max(start, multi_sample);
    }
    
    /** Beta distribution variate. */
    // mirrored from model/util/random.h
    double beta(double a, double b){
        return gsl_ran_beta (&m_gsl_gen,a,b);
    }
    
    /** Poisson distribution variate with mean lambda. */
    // mirrored from model/util/random.h
    int poisson(double lambda){
        if( !(std::isfinite)(lambda) ){
            //This would lead to an inifinite loop
            throw TRACED_EXCEPTION( "lambda is inf", Error::InfLambda );
        }
        return gsl_ran_poisson (&m_gsl_gen, lambda);
    }
    
    /** Weibull distribution variate. */
    // mirrored from model/util/random.h
    double weibull( double lambda, double k ){
        return gsl_ran_weibull( &m_gsl_gen, lambda, k );
    }
    
    // --- pinned here, but NOT part of random.h's API (no wrapper exists): ---
    
    /// Exponential: no random.h wrapper; direct GSL call.
    double exponential(double scale){
        return gsl_ran_exponential(&m_gsl_gen, scale);
    }
    /// Binomial: no random.h wrapper; direct GSL call (GSL internally picks
    /// Knuth vs BTPE from n and p).
    unsigned int binomial(unsigned int n, double p){
        return gsl_ran_binomial(&m_gsl_gen, p, n);
    }
};

/// The golden seed used everywhere (CSV header repeats it).
constexpr uint64_t SEED[4] = {1, 2, 3, 4};

/// One replayable draw series. params is an informational CSV column; p0..p3
/// are the actual arguments (interpreted per dist).
struct DrawSpec {
    const char* dist;
    const char* params;
    double p0, p1, p2, p3;
    int n;
};

/// Schedule shared verbatim between the CSV generator and the replay suite.
inline const DrawSpec* schedule(size_t& nOut) {
    static const DrawSpec S[] = {
        // dist, params, p0, p1, p2, p3, draws
        {"uniform_01",          "-",                               0.0,  0.0, 0.0, 0.0, 1000},
        {"gauss",               "mean=1.5;std=2.0",                1.5,  2.0, 0.0, 0.0, 1000},
        {"log_normal",          "meanlog=0.5;stdlog=1.2",          0.5,  1.2, 0.0, 0.0, 1000},
        // gamma: p0<1 → GSL small-shape branch; p0>=1 → Marsaglia-Tsang
        {"gamma",               "a=0.5;b=2.0",                     0.5,  2.0, 0.0, 0.0, 1000},
        {"gamma",               "a=2.5;b=0.8",                     2.5,  0.8, 0.0, 0.0, 1000},
        // beta: both-<=1 vs both->1 sampler branches (random.h → GSL)
        {"beta",                "a=0.5;b=0.7",                     0.5,  0.7, 0.0, 0.0, 1000},
        {"beta",                "a=2.0;b=3.0",                     2.0,  3.0, 0.0, 0.0, 1000},
        // poisson: mu<10 → direct method; mu>=10 → rejection branch
        {"poisson",             "mu=3.2",                          3.2,  0.0, 0.0, 0.0, 1000},
        {"poisson",             "mu=25.0",                        25.0,  0.0, 0.0, 0.0, 1000},
        {"weibull",             "lambda=2.0;k=1.5",                2.0,  1.5, 0.0, 0.0, 1000},
        // binomial/exponential: direct GSL (random.h has no wrapper for these)
        {"binomial",            "n=10;p=0.3",                     10.0,  0.3, 0.0, 0.0, 1000},
        {"binomial",            "n=60;p=0.05",                    60.0, 0.05, 0.0, 0.0, 1000},
        {"exponential",         "scale=1.7",                       1.7,  0.0, 0.0, 0.0, 1000},
        {"max_multi_lognormal", "start=0.3;n=5;meanlog=0.4;stdlog=1.1", 0.3, 5.0, 0.4, 1.1, 1000},
    };
    nOut = sizeof(S) / sizeof(S[0]);
    return S;
}

/// Deterministic gsl_cdf_ugaussian_Pinv spot values (used at random.h and
/// interventions/Deployments.h:195). Draw-independent: no RNG state consumed.
inline const double* pinv_spots(size_t& nOut) {
    static const double U[] = {
        1e-6, 1e-4, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 0.999999
    };
    nOut = sizeof(U) / sizeof(U[0]);
    return U;
}

/// Draw the next value of series s from rng.
inline double golden_draw(Rng& rng, const DrawSpec& s) {
    using std::strcmp;
    if (!strcmp(s.dist, "uniform_01"))          return rng.uniform_01();
    if (!strcmp(s.dist, "gauss"))               return rng.gauss(s.p0, s.p1);
    if (!strcmp(s.dist, "log_normal"))          return rng.log_normal(s.p0, s.p1);
    if (!strcmp(s.dist, "gamma"))               return rng.gamma(s.p0, s.p1);
    if (!strcmp(s.dist, "beta"))                return rng.beta(s.p0, s.p1);
    if (!strcmp(s.dist, "poisson"))             return static_cast<double>(rng.poisson(s.p0));
    if (!strcmp(s.dist, "weibull"))             return rng.weibull(s.p0, s.p1);
    if (!strcmp(s.dist, "binomial"))            return static_cast<double>(rng.binomial(static_cast<unsigned int>(s.p0), s.p1));
    if (!strcmp(s.dist, "exponential"))         return rng.exponential(s.p0);
    if (!strcmp(s.dist, "max_multi_lognormal")) return rng.max_multi_log_normal(s.p0, static_cast<int>(s.p1), s.p2, s.p3);
    return -1.0;  // unknown dist: the replay suite will fail loudly
}

}
