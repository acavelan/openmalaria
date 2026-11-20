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

// p8-files: PerHost.h is read by interventions + WH + checkpoint paths; the
// bit-vector deploy/restore bodies stay here (wide-cone rule).
#include "Transmission/PerHost.h"
#include "Transmission/Model.h"
#include "Host/Human.h"
#include "interventions/IvState.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
using namespace OM::util;

// -----  PerHost static  -----

AgeGroupInterpolator PerHost::relAvailAge;

// -----  PerHost non-static -----

PerHost::PerHost() {
    // Transmission state lives in HumanColumns (bind_columns); old member
    // inits (outsideTransmission=false, relativeAvailabilityHet=signaling_NaN)
    // are now the column defaults in append_row.
}

void PerHost::initialise(Context& ctx, LocalRng& rng, double availabilityFactor) {
    rel_avail_het() = availabilityFactor;
    anoph_ento_raw().resize(PerHostAnophParams::numSpecies(ctx));
    anoph_ento_avail().resize(PerHostAnophParams::numSpecies(ctx));
    anoph_prob_biting().resize(PerHostAnophParams::numSpecies(ctx));
    anoph_prob_resting().resize(PerHostAnophParams::numSpecies(ctx));

    for (size_t i = 0; i < PerHostAnophParams::numSpecies(ctx); ++i) {
        const PerHostAnophParams& base = PerHostAnophParams::get(ctx, i);
        anoph_ento_raw()[i] = base.entoAvailability->sample(rng);
        anoph_ento_avail()[i] = base.entoAvailabilityFactor * anoph_ento_raw()[i] * availabilityFactor;
        anoph_prob_biting()[i] = base.probMosqBiting.sample(rng);
        auto pRest1 = base.probMosqFindRestSite.sample(rng);
        auto pRest2 = base.probMosqSurvivalResting.sample(rng);
        anoph_prob_resting()[i] = pRest1 * pRest2;
    }
}

void PerHost::update(Host::Human& human) { interventions::iv_update(*ctx, human); }

// Note: deployment goes straight to interventions::iv_deploy (record list on
// HumanColumns.iv); update and the effect accessors below multiply over
// those records in list order.

// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over the time span now - sim::never()
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
// parity: c — history-ordered iv-record FP product — multiply in vector order, not reassociated
double PerHost::entoAvailabilityHetVecItv(size_t species) const {
    double alpha_i = anoph_ento_avail()[species];
    for (const IvRecord& r : cols->iv[row])
        alpha_i *= interventions::iv_rel_attractiveness(*ctx, r, species);
    return alpha_i;
}
double PerHost::probMosqBiting(size_t species) const {
    double P_B_i = anoph_prob_biting()[species];
    for (const IvRecord& r : cols->iv[row])
        P_B_i *= interventions::iv_preprandial_survival(*ctx, r, species);
    return P_B_i;
}
double PerHost::probMosqResting(size_t species) const {
    double pRest = anoph_prob_resting()[species];
    for (const IvRecord& r : cols->iv[row])
        pRest *= interventions::iv_postprandial_survival(*ctx, r, species);
    return pRest;
}

double PerHost::relMosqFecundity(size_t species) const {
    double relFecundity = 1.0;
    for (const IvRecord& r : cols->iv[row])
        relFecundity *= interventions::iv_rel_fecundity(*ctx, r, species);
    return relFecundity;
}

bool PerHost::hasActiveInterv(interventions::Component::Type type) const {
    // kind is recorded on the record: no component-registry lookup needed
    for (const IvRecord& r : cols->iv[row])
        if (r.kind == type && r.deploy != sim::never())
            return true;
    return false;
}

}
}

// ———  PerHostAnophParams registry ((*ctx.anophParams))  ———

namespace OM {
namespace Transmission {

void PerHostAnophParams::init(Context& ctx, const scnXml::Mosq& mosq) { ctx.anophParams->emplace_back(mosq); }

size_t PerHostAnophParams::numSpecies(Context& ctx) { return ctx.anophParams->size(); }

const PerHostAnophParams& PerHostAnophParams::get(Context& ctx, size_t species) { return (*ctx.anophParams)[species]; }

void PerHostAnophParams::setEntoAvailabilityFactor(Context& ctx, size_t species, double entoAvailability) {
    (*ctx.anophParams)[species].entoAvailabilityFactor = entoAvailability;
}

void PerHostAnophParams::calcAvailabilityPercentiles(Context& ctx) {
    auto& params = (*ctx.anophParams);
    auto& entoAvailabilityPercentiles = ctx.entoAvailabilityPercentiles;
    int nSamples = 100000;
    vector<double> samples(nSamples, 0.0);

    // Using a different seed to make sure old results will not change
    // Does it really matter?
    // parity: order — MasterRng(0,0) availability-percentile seed is load-bearing for goldens — do not “clean up”
    util::MasterRng seed(0, 0);
    seed.seed(0, 0);

    util::LocalRng rng(seed);

    // Sample
    for (int i = 0; i < nSamples; i++) {
        double avail = 0.0;
        for (size_t s = 0; s < params.size(); ++s)
            avail += params[s].entoAvailabilityFactor * params[s].entoAvailability->sample(rng);
        samples[i] = avail;
    }

    // Sort samples
    sort(samples.begin(), samples.end());

    entoAvailabilityPercentiles.resize(100);

    // Calc percetiles threshold values
    for (int i = 0; i < 100; i++)
        entoAvailabilityPercentiles[i] = samples[int(i * nSamples / 100)];
    entoAvailabilityPercentiles[0] = 0.0;
}

double PerHostAnophParams::getEntoAvailabilityPercentile(Context& ctx, int p) {
    auto& entoAvailabilityPercentiles = ctx.entoAvailabilityPercentiles;
    if (p >= 0 && p < 100)
        return entoAvailabilityPercentiles[p];
    else if (p == 100)
        return std::numeric_limits<double>::max();
    else
        throw util::xml_error(
            "Availability percentiles in intervention specifications can only be in [0 and 100]");
}

}
}
