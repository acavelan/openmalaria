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

#ifndef UTIL_SAMPLER
#define UTIL_SAMPLER

#include <limits>
#include <optional>

#include "Global.h"
#include "util/random.h"

namespace scnXml {
    class SampledValueN; class SampledValueCV; class SampledValueLN;
    class BetaMeanSample; class WeibullSample;
}

namespace OM { namespace util {
    
        /** A normal sample, which can be turned into various log-normal samples.
     * 
     * Allows generation of correlated log-normal samples with different sigma.
     */
    class NormalSample {
    public:
        NormalSample() : x(numeric_limits<double>::signaling_NaN()) {}
        
        /// convert sample to N(mu,sigma)
        double asNormal( double mu, double sigma ) const;

        /// convert sample to lnN(mu,sigma)
        double asLognormal( double mu, double sigma ) const;
        
        static NormalSample generate(LocalRng& rng);
        
        /// Generate a sample correlated with base.
        /// 
        /// @param base sample correlated against
        /// @param correlation expected correlation between these parameters
        ///     (note: if this is used to sample a log-normal, then the
        ///     correlation is on the log-scale)
        /// @param factor should equal sqrt(1 - correlation^2); may be cached
        static NormalSample generate_correlated(NormalSample base, double correlation, double factor, LocalRng& rng);
        
    private:
        NormalSample( double variate ) : x(variate) {}
        
        double x;   // variate (sampled from N(0,1))
    };
    
    /** Sampler for normal values */
    class NormalSampler {
    public:
        NormalSampler() :
            mu( numeric_limits<double>::signaling_NaN() ),
            sigma( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters such that samples taken are:
         * X ~ N( m, s² )
         * 
         * @param m Mean of sampled variates
         * @param s Square-root of variance of sampled variates
         */
        void setParams( double m, double s );

        /** As above, using an XML element. */
        void setParams( const scnXml::SampledValueN& elt );
        
        /** Sample a value. */
        double sample(LocalRng& rng) const;

        /** Create a log-normal sample from an existing normal sample. */
        inline double sample(NormalSample sample) const{
            return sample.asNormal( mu, sigma );
        }
        
        /// Return mu / mean of distribution
        inline double getMu() const{
            return mu;
        }

        /// Return sigma / standard deviation of distribution
        inline double getSigma() const{
            return sigma;
        }
        
    private:
        double mu, sigma;
    };
    
    /** Sampler for log-normal or gamma values: one final class, distribution
     * chosen by Kind (was the LognormalSampler / GammaSampler subclasses).
     *
     * fromMeanCV docs (lognormal kind): for CV > 0 computes μ, σ such that
     * E[X] = exp(μ + σ²/2) = mean, CV = √(exp(σ²) − 1); CV == 0 collapses to a
     * point mass (μ = ln(mean), or −∞ if mean == 0). Gamma kind: k = 1/CV²,
     * θ = mean/k; CV == 0 is a point mass at the mean. fromMeanVariance:
     * lognormal uses CV = variance/mean; gamma uses k = mean²/variance.
     * Both throw util::xml_scenario_error on invalid mean/CV/variance
     * (lognormal allows mean == 0 iff CV == 0; gamma requires mean > 0). */
    class Sampler final {
    public:
        enum class Kind { Lognormal, Gamma };

        Sampler() = default;
        Sampler& operator=(const Sampler&) = default;

        std::unique_ptr<Sampler> clone() const { return std::make_unique<Sampler>(*this); }

        static unique_ptr<Sampler> fromMeanCV( Kind kind, double mean, double CV, std::optional<double> truncate = std::nullopt);
        static unique_ptr<Sampler> fromMeanVariance( Kind kind, double mean, double variance, std::optional<double> truncate = std::nullopt);

        /** Scale the mean (lognormal kind: mu is log-space, so sigma specified
         * via CV is unaffected). */
        void scaleMean( double scalar );

        /// Get the mean of the distribution.
        double mean() const;

        /** Sample a value. */
        double sample(LocalRng& rng) const;

        /** Create a log-normal sample from an existing normal sample. */
        inline double sample(NormalSample s) const{
            if (kind != Kind::Lognormal)
                throw std::runtime_error("sample(NormalSample) not implemented for this distribution");
            return s.asLognormal( mu, sigma );
        }

        /// P(X ≤ x) for this distribution.
        double cdf(double x) const;

    private:
        Kind kind = Kind::Lognormal;
        // Lognormal kind: log-space mu, sigma. Gamma kind: mean mu, k, theta.
        // CV / variance are the configuring inputs; truncate caps samples.
        double mu       = std::numeric_limits<double>::signaling_NaN();
        double sigma    = std::numeric_limits<double>::signaling_NaN();
        double k        = std::numeric_limits<double>::signaling_NaN();
        double theta    = std::numeric_limits<double>::signaling_NaN();
        double variance = std::numeric_limits<double>::signaling_NaN();
        double CV       = std::numeric_limits<double>::signaling_NaN();
        std::optional<double> truncate = std::nullopt;
    };

    /** Create a sampler from a SampledValueCV xml snippet (defined in
     * model/xml/util.cpp: the model/xml/ layer owns all schema access;
     * headers only forward-declare). */
    unique_ptr<Sampler> createSampler(double mean, const scnXml::SampledValueCV& elt);

    /** Create a lognormal Sampler from a SampledValueLN xml snippet. */
    unique_ptr<Sampler> createSampler(const scnXml::SampledValueLN& elt);

    /** Sampler for the Beta distribution.
     *
     * Input may be alpha and beta or mean and variance. As a special case,
     * variance zero is also supported, meaning return the mean. */
    class BetaSampler {
    public:
        BetaSampler() :
            a( numeric_limits<double>::signaling_NaN() ),
            b( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters: alpha, beta.
         * 
         * As a special case, when beta=0, alpha is interpreted as the mean and
         * is returned without sampling. */
        inline void setParams( double alpha, double beta ){
            a=alpha;
            b=beta;
        }
        
        /** Set parameters: alpha and beta are calculated such that the mean
         * and variance are as given. Variance may be zero, in which case the
         * mean is returned directly without sampling. */
        void setParamsMV( double mean, double variance );
        
        /** Set parameters from an XML element (defined in model/xml/util.cpp). */
        void setParams( const scnXml::BetaMeanSample& elt );
        
        /** Sample a value. */
        double sample(LocalRng& rng) const;
        
    private:
        //Note: if b is 0, then alpha is mean. Otherwise a and b are
        //the expected (α,β) parameters to the beta distribution.
        double a, b;
    };
    
    /** Sampler for the Weibull distribution. */
    class WeibullSampler {
    public:
        WeibullSampler() :
            scale( numeric_limits<double>::signaling_NaN() ),
            shape( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters from scale, shape. */
        inline void setScaleShape( double lambda, double k ){
            scale = lambda;
            shape = k;
        }
        
        /** Set parameters from an XML element (defined in model/xml/util.cpp). */
        void setParams( const scnXml::WeibullSample& elt );
        
        /** Sample a value. */
        inline double sample(LocalRng& rng) const{
            // parity: order — local RNG stream: draw call order fixes stream state
            return rng.weibull( scale, shape );
        }
        
    private:
        double scale, shape;    // λ, k
    };
    
} }
#endif
