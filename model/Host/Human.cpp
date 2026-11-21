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
#include "Host/Human.h"
#include "Host/InfectionIncidenceModel.h"
#include "Clinical/ClinicalModel.h"
#include "Host/WithinHost/WHInterface.h"

#include "Transmission/TransmissionModel.h"
#include "util/ModelOptions.h"
#include "util/vectors.h"
#include "util/StreamValidator.h"
#include "Population.h"
#include "interventions/InterventionManager.h"
#include "interventions/Vaccine.h"
#include "mon/Monitoring.h"

namespace OpenMalaria {

struct HumanHet {
    /* Human heterogeneity; affects:
     * comorbidityFactor (stored in PathogenesisModel)
     * treatmentSeekingFactor (stored in CaseManagementModel)
     * availabilityFactor (stored in PerHost) */
    double comorbidityFactor;
    double treatmentSeekingFactor;
    double availabilityFactor;

    HumanHet() : comorbidityFactor(1.0), treatmentSeekingFactor(1.0), availabilityFactor(1.0) {}
};

HumanHet hetSample(LocalRng& rng)
// §9 draw contract (checked phase 6.4): max 4 draws — the three leading `if`s
// draw at most one each, then exactly one draw in the else-if chain.
// parity: Deployments.h:245-306 (het consumed by deploys in that order)
{
    HumanHet het;
    if (ModelOptions::option(TRANS_HET)) {
        het.availabilityFactor = 0.2;
        // parity: order — local RNG stream: draw call order fixes stream state
        if (rng.bernoulli(0.5)) {
            het.availabilityFactor = 1.8;
        }
    }
    if (ModelOptions::option(COMORB_HET)) {
        het.comorbidityFactor = 0.2;
        // parity: order — local RNG stream: draw call order fixes stream state
        if (rng.bernoulli(0.5)) {
            het.comorbidityFactor = 1.8;
        }
    }
    if (ModelOptions::option(TREAT_HET)) {
        het.treatmentSeekingFactor = 0.2;
        // parity: order — local RNG stream: draw call order fixes stream state
        if (rng.bernoulli(0.5)) {
            het.treatmentSeekingFactor = 1.8;
        }
    }
    if (ModelOptions::option(TRANS_TREAT_HET)) {
        het.treatmentSeekingFactor = 0.2;
        het.availabilityFactor = 1.8;
        // parity: order — local RNG stream: draw call order fixes stream state
        if (rng.bernoulli(0.5)) {
            het.treatmentSeekingFactor = 1.8;
            het.availabilityFactor = 0.2;
        }
    } else if (ModelOptions::option(COMORB_TREAT_HET)) {
        // parity: order — local RNG stream: draw call order fixes stream state
        if (rng.bernoulli(0.5)) {
            het.comorbidityFactor = 1.8;
            het.treatmentSeekingFactor = 0.2;
        } else {
            het.comorbidityFactor = 0.2;
            het.treatmentSeekingFactor = 1.8;
        }
    } else if (ModelOptions::option(COMORB_TRANS_HET)) {
        het.availabilityFactor = 1.8;
        het.comorbidityFactor = 1.8;
        // parity: order — local RNG stream: draw call order fixes stream state
        if (rng.bernoulli(0.5)) {
            het.availabilityFactor = 0.2;
            het.comorbidityFactor = 0.2;
        }
    } else if (ModelOptions::option(TRIPLE_HET)) {
        het.availabilityFactor = 1.8;
        het.comorbidityFactor = 1.8;
        het.treatmentSeekingFactor = 0.2;
        // parity: order — local RNG stream: draw call order fixes stream state
        if (rng.bernoulli(0.5)) {
            het.availabilityFactor = 0.2;
            het.comorbidityFactor = 0.2;
            het.treatmentSeekingFactor = 1.8;
        }
    }
    return het;
}

Human::Human(RunContext& ctx, SimTime dateOfBirth, HumanColumns& columns)
    : ctx(&ctx), cols(&columns), infIncidence(InfectionIncidenceModel::createModel(ctx)) {
    // parity: order — each Human copies+seeds its LocalRng from the master at
    // construction (append_row builds the rng column row; same ctor point as
    // the old member init): human-creation order == master consumption order;
    // never reorder Human ctor.
    row = columns.append_row(dateOfBirth, ctx.master_RNG);

    // Initial humans are created at time 0 and may have dateOfBirth in past. Otherwise dateOfBirth must be now.
    assert(dateOfBirth == ctx.clock.nowOrTs1() || (ctx.clock.now() == zero() && dateOfBirth < ctx.clock.now()));

    HumanHet het = hetSample(rng());
    withinHostModel = WHInterface::createWithinHostModel(ctx, rng(), het.comorbidityFactor);
    withinHostModel->bind_columns(columns, row);
    double iiFactor = infIncidence->getAvailabilityFactor(rng(), 1.0);
    perHostTransmission.ctx = &ctx;
    perHostTransmission.bind_columns(columns, row);
    perHostTransmission.initialise(ctx, rng(), het.availabilityFactor * iiFactor);
    clinicalModel = ClinicalModel::createClinicalModel(ctx, het.treatmentSeekingFactor);
    clinicalModel->bind_columns(columns, row);
}

void Human::set_row(size_t r) {
    row = r;
    // Rebind every submodel caching (cols,row): a stale mirror reads/writes
    // another host's columns after compaction (WHInterface immunity, Clinical).
    if (withinHostModel)
        withinHostModel->bind_columns(*cols, r);
    if (clinicalModel)
        clinicalModel->bind_columns(*cols, r);
    perHostTransmission.bind_columns(*cols, r);
}

double Human::vaccine_factor(Vaccine::Types type, uint32_t genotype) const {
    return OpenMalaria::vaccine_factor(*ctx, cols->vac[row], type, genotype);
}

namespace {
/** Insertion point of `id` in the row's spans (sorted ascending by id, so
 * iteration order matches the old std::map). */
vector<pair<ComponentId, SimTime>>::iterator find_cohort(vector<pair<ComponentId, SimTime>>& spans, ComponentId id) {
    return lower_bound(spans.begin(), spans.end(), id,
                       [](pair<ComponentId, SimTime> const& s, ComponentId i) { return s.first < i; });
}
}

void Human::addToCohort(ComponentId id, SimTime duration) {
    if (duration <= zero())
        return; // nothing to do
    auto& spans = cols->cohorts[row];
    auto it = find_cohort(spans, id);
    if (it != spans.end() && it->first == id)
        it->second = ctx->clock.nowOrTs1() + duration;
    else
        spans.insert(it, make_pair(id, ctx->clock.nowOrTs1() + duration));
    cols->cohort_set[row] = OpenMalaria::updateCohortSet(*ctx, cols->cohort_set[row], id, true);
}

void Human::removeFromCohort(ComponentId id) {
    auto& spans = cols->cohorts[row];
    auto it = find_cohort(spans, id);
    if (it != spans.end() && it->first == id)
        spans.erase(it);
}

void Human::removeFirstEvent(SubPopRemove::RemoveAtCode code) {
    const vector<ComponentId>& removeAtList = ctx->removeAtIds[code];
    for (auto it = removeAtList.begin(), end = removeAtList.end(); it != end; ++it) {
        auto& spans = cols->cohorts[row];
        auto expIt = find_cohort(spans, *it);
        if (expIt != spans.end() && expIt->first == *it) {
            if (expIt->second >= ctx->clock.nowOrTs0()) {
                // removeFirstEvent() is used for onFirstBout, onFirstTreatment
                // and onFirstInfection cohort options. Health system memory must
                // be reset for this to work properly; in theory the memory should
                // be independent for each cohort, but this is a usable approximation.
                clinicalModel->flushReports(); // reset HS memory

                // report removal due to first infection/bout/treatment
                recordEvent(measure("nSubPopRemovalFirstEvent"), *this);
            }
            cols->cohort_set[row] = OpenMalaria::updateCohortSet(*ctx, cols->cohort_set[row], expIt->first, false);
            // remove (affects reporting, restrictToSubPop and cumulative deployment):
            spans.erase(expIt);
        }
    }
}

void Human::updateCohortSet() {
    // check sub-pop expiry
    auto& spans = cols->cohorts[row];
    // NOTE: unlike the old std::map, vector erase() invalidates the end()
    // sentinel, so iterate by index (same order, same event sequence).
    for (size_t k = 0; k < spans.size();) {
        if (!(spans[k].second >= ctx->clock.ts0())) { // membership expired
            // don't flush reports
            // report removal due to expiry
            recordEvent(measure("nSubPopRemovalTooOld"), *this);
            cols->cohort_set[row] = OpenMalaria::updateCohortSet(*ctx, cols->cohort_set[row], spans[k].first, false);
            // erase element, but continue iteration
            spans.erase(spans.begin() + k);
        } else {
            ++k;
        }
    }
}

// -----  Non-static functions: per-time-step update  -----
void summarize(RunContext& ctx, Human& human, bool surveyOnlyNewEp) {
    if (surveyOnlyNewEp && human.clinicalModel->isExistingCase()) {
        // This modifies the denominator to treat the health-system-memory
        // period immediately after a bout as 'not at risk'.
        return;
    }

    recordStat(measure("nHost"), human);
    recordStat(measure("sumAge"), human, inYears(human.age(ctx.clock.now())));
    bool patent = human.withinHostModel->summarize(human);
    human.infIncidence->summarize(human);

    if (patent && isReported(ctx)) {
        // this should happen after all other reporting!
        human.removeFirstEvent(SubPopRemove::ON_FIRST_INFECTION);
    }
}

// —— The 7 §9 host-update stages (bodies verbatim from the old update;
// statements, order and `// parity:` comments untouched) ——

bool host_is_dead(HostView& v) {
    if (v.human.clinicalModel->isDead(v.scratch.age0)) {
        v.human.kill();
        return true;
    }
    return false;
}

void host_update_mon_cohorts(HostView& v) {
    streamValidate(v.scratch.age0);

    // monitoringAgeGroup is the group for the start of the time step.
    updateAgeGroup(v.ctx, v.human.monitoringAgeGroup(), v.scratch.age0);

    v.human.updateCohortSet();
}

void host_calculate_eir(HostView& v) {
    // Age at  the end of the update period. In most cases
    // the difference between this and age at the start is not especially
    // important in the model design, but since we parameterised with
    // age1 we should stick with it.
    v.scratch.age1 = inYears(v.human.age(v.ctx.clock.ts1()));

    // age1 used only in PerHost::relativeAvailabilityAge(); difference to age0 should be minor
    v.transmission.getEIR(v.human, v.scratch.age0, v.scratch.age1, v.scratch.EIR_per_genotype_i,
                          v.scratch.EIR_per_genotype_l);

    v.scratch.eir_i = vectors::sum(v.scratch.EIR_per_genotype_i);
    v.scratch.eir_l = vectors::sum(v.scratch.EIR_per_genotype_l);
    v.scratch.eir = v.scratch.eir_i + v.scratch.eir_l;
}

void host_sample_new_infections(HostView& v) {
    double x = 0.0;

    if (v.scratch.eir > 0)
        x = v.scratch.eir_i / v.scratch.eir;

    double expectedNNewInfs = v.human.infIncidence->expectedNumNewInfections(v.human, v.scratch.eir);
    double expectedNNewInfs_i = x * expectedNNewInfs;
    double expectedNNewInfs_l = (1.0 - x) * expectedNNewInfs;

    // One more rng call if imported infection
    v.scratch.nNewInfs_l = v.human.infIncidence->numNewInfections(v.human, expectedNNewInfs_l);
    v.scratch.nNewInfs_i = v.human.infIncidence->numNewInfections(v.human, expectedNNewInfs_i);
}

void host_within_host_update(HostView& v) {
    // age1 used when medicating drugs (small effect) and in immunity model (which was parameterised for it)
    v.human.withinHostModel->update(v.human, v.human.rng(), v.scratch.nNewInfs_i, v.scratch.nNewInfs_l,
                                    v.scratch.EIR_per_genotype_i, v.scratch.EIR_per_genotype_l, v.scratch.age1);
    v.human.infIncidence->reportNumNewInfections(v.human, v.scratch.nNewInfs_i, v.scratch.nNewInfs_l);
}

void host_clinical_update(HostView& v) {
    // age1 used to get case fatality and sequelae probabilities, determine pathogenesis
    v.human.clinicalModel->update(v.human, v.scratch.age1, v.scratch.age0 == zero());
}

void host_update_infant_deaths(HostView& v) { v.human.clinicalModel->updateInfantDeaths(v.scratch.age0); }

void update(RunContext& ctx, Human& human, TransmissionModel& transmission, HostScratch& scratch) {
    HostView v{ctx, human, transmission, scratch};
    // For integer age checks we use age0 to e.g. get 73 steps comparing less than 1 year old
    scratch.age0 = human.age(ctx.clock.ts0());
    if (host_is_dead(v))
        return;

    host_update_mon_cohorts(v);
    host_calculate_eir(v);
    host_sample_new_infections(v);
    host_within_host_update(v);
    host_clinical_update(v);
    host_update_infant_deaths(v);
}

}
