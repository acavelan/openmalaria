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

#include "Clock.h"
#include "Context.h"
#include "Host/HostColumns.h"          // leaf: SoA columns owned by Population (IvRecord)
#include "interventions/ComponentId.h" // leaf: Component::Type + ComponentId only
#include "util/AgeGroupInterpolation.h"
#include "util/sampler.h"

#include "util/random.h"

namespace scnXml {
class Mosq;
}

namespace OM {
namespace Transmission {

using util::AgeGroupInterpolator;
using util::LocalRng;

/** Stores vector model data applicable between a category of host and a
 * mosquito species: intervention descriptions and model parameters.
 *
 * Parameters are read from XML, and the availability rate is adjusted. */
class PerHostAnophParams {
  public:
    PerHostAnophParams(const scnXml::Mosq& mosq);

    PerHostAnophParams(PerHostAnophParams&&) noexcept = default;
    PerHostAnophParams& operator=(PerHostAnophParams&&) noexcept = default;

    PerHostAnophParams(const PerHostAnophParams&) = delete;
    PerHostAnophParams& operator=(const PerHostAnophParams&) = delete;

    /// Register one species' parameters (registry: (*ctx.anophParams)).
    static void init(Context& ctx, const scnXml::Mosq& mosq);

    /// Get the number of vector species
    static size_t numSpecies(Context& ctx);

    /// Get parameters for the given vector species
    static const PerHostAnophParams& get(Context& ctx, size_t species);

    /** entoAvailability is calculated externally, then set after other
     * parameters have been initialised.
     * 
     * This function doesn't need to exist, but helps make this fact obvious.
     * 
     * It should be called exactly once.
      */
    static void setEntoAvailabilityFactor(Context& ctx, size_t species, double entoAvailability);

    /** Sample the population availability distribution and store its
     * percentiles in ctx.entoAvailabilityPercentiles. */
    static void calcAvailabilityPercentiles(Context& ctx);

    static double getEntoAvailabilityPercentile(Context& ctx, int p);

    /** @brief Probabilities of finding a host and surviving a feeding cycle
     * 
     * These parameters describe the mean and heterogeneity of α_i, P_B_i,
     * P_C_i and P_D_i across the human population. */
    //@{
    /** Availability rate (α_i) */
    double entoAvailabilityFactor = 1.0;

    /** Availability heterogeneity factor
     * This multiplies the entoAvailabilityFactor. */
    unique_ptr<util::Sampler> entoAvailability;

    /** Probability of mosquito successfully biting host (P_B_i) */
    util::BetaSampler probMosqBiting;

    /** Probability of mosquito escaping human and finding a resting site without
     * dying, after biting the human (P_C_i). */
    util::BetaSampler probMosqFindRestSite;

    /** Probability of mosquito successfully resting after finding a resting site
     * (P_D_i). */
    util::BetaSampler probMosqSurvivalResting;
    //@}

    // (no per-type data here: the registry lives in Context::anophParams)
};

/** Contains Model parameters which need to be stored per host.
 *
 * Currently many members are public and directly accessed. */
class PerHost {
  public:
    /// @brief Static member functions
    //@{
    /** Static initialisation. */
    static void init(Context& ctx, const scnXml::AgeGroupValues& availabilityToMosquitoes);
    //@}

    ///@brief Initialisation / checkpionting
    //@{
    Context* ctx = nullptr; ///< set by the owning Human ctor; used on checkpoint restore
    PerHost();
    ~PerHost() = default;
    PerHost(PerHost&&) = default;
    PerHost& operator=(PerHost&&) = default;
    void initialise(Context& ctx, LocalRng& rng, double availabilityFactor);
    //@}

    /// Call once per time step. Updates net holes (interventions/IvState.cpp).
    void update(Host::Human& human);

    /** Calculates the adjustment for body size in exposure to mosquitoes,
     * relative to an average adult.
     * 
     * The bites are assumed proportional to average surface area for hosts of
     * the given age. Linear interpolation is used to calculate this from the
     * input array of surface areas. 
     * 
     * @param ageYears Age of host
     * @return the ratio of bites received by the host to the average for an adult 
     *
     * This is the age factor of availiability; mean output should be
     * mean population availability (that is, 1.0/invMeanPopAvail).
     * 
     * Also has a switch to put individuals entirely outside transmission. */
    inline double relativeAvailabilityAge(double ageYears) const {
        return outside_trans() ? 0.0 : relAvailAge.eval(ageYears);
    }

    /** Get the availability of this host to mosquitoes relative to an average
     * adult (including heterogeneity and age effects).
     *
     * Used to drive a simulation from an input EIR.
     * Is rel_avail_het()*relativeAvailabilityAge(ageYears).
     * 
     * Mean output is less than 1.0 (roughly 1.0/invMeanPopAvail).
     */
    inline double relativeAvailabilityHetAge(double ageYears) const {
        return rel_avail_het() * relativeAvailabilityAge(ageYears);
    }

    /** Availability rate of human to mosquitoes (α_i). Equals 
     * entoAvailabilityHetVecItv()*getRelativeAvailability().
     *
     * To be clear, this includes effects from HetVecItc (het, interv, availability
     * rate) as well as age (avail. relative to an adult). It does not divide
     * by the average availability of the population, which was incorrectly done
     * in the past. */
    inline double entoAvailabilityFull(size_t species, double ageYears) const {
        return entoAvailabilityHetVecItv(species) * relativeAvailabilityAge(ageYears);
    }
    //@}

    /** Availability of host to mosquitoes (α_i) excluding age factor
     * 
     * (Includes heterogeneity, intervention, and human-to-vector availability
     * rate factors.)
     * 
     * Assume mean is human-to-vector availability rate factor. */
    double entoAvailabilityHetVecItv(size_t species) const;

    ///@brief Get effects of interventions pre/post biting
    //@{
    /** Probability of a mosquito succesfully biting a host (P_B_i). */
    double probMosqBiting(size_t species) const;
    /** Probability of a mosquito succesfully finding a resting
     * place after biting and then resting (P_C_i * P_D_i). */
    double probMosqResting(size_t species) const;
    /** Multiplicative factor for the number of fertile eggs laid by mosquitoes
     * after feeding on this host. Should be 1 normally, less than 1 to reduce
     * fertility, greater than 1 to increase. */
    double relMosqFecundity(size_t species) const;
    //@}

    ///@brief Convenience wrappers around several functions
    //@{
    /// entoAvailabilityHetVecItv * probMosqBiting
    inline double availBite(size_t species) const {
        return entoAvailabilityHetVecItv(species) * probMosqBiting(species);
    }
    //@}

    ///@brief Miscellaneous
    //@{
    /** Get the age at which individuals are considered adults (i.e. where
     * availability to mosquitoes reaches its maximum). */
    static inline SimTime adultAge(const SimClock& clock) {
        return clock.from_years_d(relAvailAge.firstGlobalMaximum());
    }

    /** Get whether the user has any active deployments of interventions of
     * the given type, where type is one of ITN/IRS/GVI (in other cases this
     * will always return false). */
    bool hasActiveInterv(interventions::Component::Type type) const;

    //@}

    /** Species availability rate of human to mosquitoes: per-species
     * vector at this row (old anoph* members; refs keep the old syntax). */
    std::vector<double>& anoph_ento_raw() const { return cols->anoph_ento_raw[row]; }
    std::vector<double>& anoph_ento_avail() const { return cols->anoph_ento_avail[row]; }
    std::vector<double>& anoph_prob_biting() const { return cols->anoph_prob_biting[row]; }
    std::vector<double>& anoph_prob_resting() const { return cols->anoph_prob_resting[row]; }
    /// Heterogeneity factor in availability (old relativeAvailabilityHet;
    /// already multiplied into the entoAvailability HostMosquitoInteraction uses).
    double& rel_avail_het() const { return cols->rel_avail_het[row]; }
    /// Determines whether human is outside transmission (old bool member;
    /// u8 column, same 1-byte wire shape, truthiness identical).
    uint8_t& outside_trans() const { return cols->outside_trans[row]; }

    /** Bind the host human's column row (Human ctor before initialise;
     * Human::set_row rebinds — stale mirror = wrong host after renumber). */
    void bind_columns(HumanColumns& c, size_t r) {
        cols = &c;
        row = r;
    }

    /// The host human's SoA columns/row (see bind_columns).
    HumanColumns* cols = nullptr;
    size_t row = 0;

  private:
    static AgeGroupInterpolator relAvailAge;
};

}
}
