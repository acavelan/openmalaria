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

// Wire section: the living human population.
//
// Restore reuses the real Human constructor (row append + copy-seed +
// kind submodels, so the RNG RULE construction moments are reproduced),
// then overwrites every field from the wire; the per-human RNG is read
// last, after the subobject creation drew from those fresh streams.
//
// NOT on the wire (rewrite contract): iv/vac/cohorts/copula columns —
// per-host intervention state; see HANDOFF.md "Round-trip limits".

#include "sections.h"

#include "Population.h"
#include "Host/Human.h"
#include "Host/InfectionIncidenceModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/Episode.h"
#include "PkPd/LSTMDrug.h"
#include "RunContext.h"

namespace OM {

namespace Host {

void checkpoint(checkpoint::Archive& a, InfectionIncidenceModel& m) {
    a.value(m.m_pInfected);
    a.value(m.m_cumulativeEIRa);
}

} // namespace Host

namespace WithinHost {

void checkpoint(checkpoint::Archive& a, VivaxBrood& b) {
    a.vector(b.releaseDates);
    a.value(b.bloodStageClearDate);
    a.boolean(b.primaryHasStarted);
    a.boolean(b.relapseHasStarted);
    a.boolean(b.relapsebHasStarted);
    a.boolean(b.hadEvent);
    a.boolean(b.hadRelapse);
    a.value(b.origin);
}

namespace Pathogenesis {
void checkpoint(checkpoint::Archive& a, PathogenesisModel& pm) {
    a.value(pm._comorbidityFactor);
    a.value(pm._pyrogenThres);
}
} // namespace Pathogenesis

void checkpoint(checkpoint::Archive& a, Infection& inf) {
    a.value(inf.m_startDate);
    a.value(inf.m_density);
    a.value(inf.m_cumulativeExposureJ);
    a.value(inf.m_genotype);
    a.value(inf.m_origin);
    switch (inf.kind) {
    case Infection::Kind::Descriptive:
        a.value(inf.descriptive.duration);
        a.boolean(inf.descriptive.notPrintedMDWarning);
        break;
    case Infection::Kind::Empirical:
        for (double& d : inf.empirical->_laggedLogDensities)
            a.value(d);
        a.value(inf.empirical->_patentGrowthRateMultiplier);
        break;
    case Infection::Kind::Molineaux:
        for (float& f : inf.molineaux->mi)
            a.value(f); // mi: plain float array
        a.value(inf.molineaux->Sm_summation);
        for (float& f : inf.molineaux->lagged_Pc)
            a.value(f);
        a.value(inf.molineaux->Pc_star);
        a.value(inf.molineaux->Pm_star);
        { // variants grow with infection age (clone per relapse), so the
            // count goes on the wire; Variant has a user default ctor (not
            // trivially copyable) hence field-by-field.
            uint64_t nv = inf.molineaux->variants.size();
            a.value(nv);
            if (a.reading()) {
                if (nv > uint64_t(a.remaining()) / 64u)
                    throw OM::checkpoint::checkpoint_error("variant count too large");
                inf.molineaux->variants.resize(size_t(nv));
            }
            for (auto& v : inf.molineaux->variants) {
                a.value(v.Pi1);
                a.value(v.Pi2);
                a.value(v.Si_summation);
                for (float& f : v.lagged_Pi)
                    a.value(f);
            }
        }
        break;
    case Infection::Kind::Penny:
        for (double& d : inf.penny->cirDensities)
            a.value(d);
        for (double& d : inf.penny->seqDensities)
            a.value(d);
        a.value(inf.penny->threshold_N);
        a.value(inf.penny->threshold_V);
        a.value(inf.penny->threshold_C);
        a.value(inf.penny->variantSpecificSummation);
        a.value(inf.penny->clonalSummation);
        break;
    }
}

} // namespace WithinHost

namespace PkPd {

void checkpoint(checkpoint::Archive& a, LSTMDrug& drug) {
    a.enumerated(drug.kind);
    { // pair isn't trivially copyable: element-wise
        uint64_t n = drug.doses.size();
        a.value(n);
        if (a.reading()) {
            drug.doses.clear();
            for (uint64_t i = 0; i < n; ++i) {
                LSTMDrug::DoseVec::value_type d;
                a.value(d.first);
                a.value(d.second);
                drug.doses.push_back(d);
            }
        } else {
            for (auto& d : drug.doses) {
                a.value(d.first);
                a.value(d.second);
            }
        }
    }
    a.value(drug.vol_dist);
    switch (drug.kind) {
    case LSTMDrug::Kind::OneComp:
        a.value(drug.concentration);
        a.value(drug.neg_elim_sample);
        break;
    case LSTMDrug::Kind::ThreeComp:
        a.value(drug.concA);
        a.value(drug.concB);
        a.value(drug.concC);
        a.value(drug.concABC);
        a.value(drug.elim_sample);
        a.value(drug.k12);
        a.value(drug.k21);
        a.value(drug.k13);
        a.value(drug.k31);
        a.value(drug.nka);
        a.value(drug.last_bm);
        a.value(drug.na);
        a.value(drug.nb);
        a.value(drug.ng);
        a.value(drug.A);
        a.value(drug.B);
        a.value(drug.C);
        break;
    case LSTMDrug::Kind::Conversion:
        a.value(drug.qtyG);
        a.value(drug.qtyP);
        a.value(drug.qtyM);
        a.value(drug.nkP_sample);
        a.value(drug.nconv_sample);
        a.value(drug.nkM_sample);
        a.value(drug.vol_dist_metabolite);
        a.value(drug.nka);
        a.value(drug.last_bm);
        break;
    }
}

namespace {
void ckLSTM(checkpoint::Archive& a, RunContext& ctx, LSTMModel& m, util::LocalRng& hostRng) {
    uint64_t n = m.m_drugs.size();
    a.value(n);
    if (a.reading()) {
        for (uint64_t i = 0; i < n; ++i) {
            size_t typeIndex;
            a.value(typeIndex);
            // Fresh instance from the drug type config (kind + metabolite
            // structure); its sampled parameters come from the wire next.
            m.m_drugs.push_back(ctx.drugs.createInstance(hostRng, typeIndex));
            checkpoint(a, *m.m_drugs.back());
        }
    } else {
        for (auto& drug : m.m_drugs) {
            size_t typeIndex = drug->getIndex();
            a.value(typeIndex);
            checkpoint(a, *drug);
        }
    }
    // Pending medications (queue order matters: medicate() consumes front).
    uint64_t q = m.medicateQueue.size();
    a.value(q);
    if (a.reading()) {
        m.medicateQueue.clear();
        for (uint64_t i = 0; i < q; ++i) {
            MedicateData md;
            a.value(md.drug);
            a.value(md.qty);
            a.value(md.time);
            m.medicateQueue.push_back(md);
        }
    } else {
        for (auto& md : m.medicateQueue) {
            a.value(md.drug);
            a.value(md.qty);
            a.value(md.time);
        }
    }
}
} // namespace

} // namespace PkPd

namespace Clinical {

namespace {
void ckEpisode(checkpoint::Archive& a, Episode& ep) {
    a.value(ep.time);
    if (ep.time != sim::never()) {
        a.value(ep.surveyPeriod);
        a.value(ep.ageGroup);
        a.value(ep.cohortSet);
        int32_t state = ep.state;
        a.value(state);
        int32_t infType = ep.infectionType;
        a.value(infType);
        if (a.reading()) {
            ep.state = Episode::State(state);
            ep.infectionType = WithinHost::InfectionOrigin(infType);
        }
    }
}
} // namespace

void checkpoint(checkpoint::Archive& a, ClinicalModel& cm) {
    a.value(cm.m_treatmentSeekingFactor);
    a.value(cm.doomed);
    ckEpisode(a, cm.latestReport);
    {
        int32_t s = cm.latestState;
        a.value(s);
        cm.latestState = Episode::State(s);
        s = cm.esPgState;
        a.value(s);
        cm.esPgState = Episode::State(s);
    }
}

} // namespace Clinical

namespace checkpoint {

namespace {

using WithinHost::WHInterface;

void ckWH(Archive& a, RunContext& ctx, WHInterface& wh, util::LocalRng& rng) {
    // Shared P. falciparum block (written for every kind; Vivax leaves most
    // at their defaults — same bytes, simplest symmetric code).
    a.value(wh.numInfs);
    a.value(wh._innateImmSurvFact);
    a.value(wh.totalDensity);
    a.value(wh.hrp2Density);
    a.value(wh.timeStepMaxDensity);
    a.vector(wh.m_y_lag_i);
    a.vector(wh.m_y_lag_l);
    a.value(wh.treatExpiryLiver);
    a.value(wh.treatExpiryBlood);
    {
        bool pm = wh.pathogenesisModel != nullptr;
        a.boolean(pm);
        if (a.reading() && pm && !wh.pathogenesisModel)
            throw checkpoint_error("pathogenesis model kind mismatch");
        if (pm)
            WithinHost::Pathogenesis::checkpoint(a, *wh.pathogenesisModel);
    }
    switch (wh.m_kind) {
    case WHInterface::Kind::Common:
    case WHInterface::Kind::Descriptive: {
        a.value(wh.hetMassMultiplier);
        PkPd::ckLSTM(a, ctx, wh.pkpdModel, rng);
        uint64_t n = wh.infections.size();
        a.value(n);
        if (a.reading()) {
            // Kind-correct creation: Descriptive builds directly (no factory
            // set); the other Common kinds use the factory from WH init.
            // Constructor draws are overwritten, like the Human ctor's.
            const bool desc = wh.m_kind == WHInterface::Kind::Descriptive;
            if (!desc && ctx.createInfection == nullptr)
                throw checkpoint_error("within-host infection factory unset");
            for (uint64_t i = 0; i < n; ++i) {
                WithinHost::Infection* inf =
                    desc ? WithinHost::infection_new(ctx, rng, 0, 0, WithinHost::Infection::Kind::Descriptive)
                         : ctx.createInfection(ctx, rng, 0, 0);
                WithinHost::checkpoint(a, *inf);
                wh.infections.push_back(inf);
            }
            if (uint64_t(wh.infections.size()) != n)
                throw checkpoint_error("internal: infection restore count");
        } else {
            for (WithinHost::Infection* inf : wh.infections)
                WithinHost::checkpoint(a, *inf);
        }
        break;
    }
    case WHInterface::Kind::Vivax: {
        a.boolean(wh.noPQ);
        {
            int32_t m = wh.morbidity;
            a.value(m);
            if (a.reading())
                wh.morbidity = WithinHost::Pathogenesis::State(m);
        }
        a.value(wh.treatExpiryLiverVivax);
        a.value(wh.treatExpiryBloodVivax);
        a.value(wh.pEvent);
        a.value(wh.pFirstRelapseEvent);
        a.value(wh.pSecondRelapseEvent);
        a.value(wh.pSevere);
        uint64_t n = wh.vivaxInfections.size();
        a.value(n);
        if (a.reading()) {
            for (uint64_t i = 0; i < n; ++i) {
                // Ctor samples a fresh brood (draws overwritten below).
                wh.vivaxInfections.emplace_back(ctx, rng, 0, &wh);
                WithinHost::checkpoint(a, wh.vivaxInfections.back());
            }
        } else {
            for (auto& brood : wh.vivaxInfections)
                WithinHost::checkpoint(a, brood);
        }
        break;
    }
    case WHInterface::Kind::Mock:
        throw checkpoint_error("unit-test mock within-host model on wire");
    }
}

void ckHuman(Archive& a, RunContext& ctx, Host::Human& h) {
    auto& c = *h.cols;
    const size_t r = h.row;
    a.value(c.dob[r]);
    a.value(c.cohort_set[r]);
    a.value(c.mon_age_group[r]);
    a.value(c.next_cts_dist[r]);
    a.value(c.imm_h[r]);
    a.value(c.imm_Y[r]);
    a.value(c.imm_Y_lag[r]);
    a.value(c.cum_prim_inf[r]);
    a.value(c.t_last_treatment[r]);
    a.value(c.es_case_start[r]);
    a.value(c.es_time_recovery[r]);
    a.value(c.es_time_treatment[r]);
    a.value(c.es_prev_density[r]);
    a.vector(c.anoph_ento_raw[r]);
    a.vector(c.anoph_ento_avail[r]);
    a.vector(c.anoph_prob_biting[r]);
    a.vector(c.anoph_prob_resting[r]);
    a.value(c.rel_avail_het[r]);
    a.value(c.outside_trans[r]);
    Host::checkpoint(a, *h.infIncidence);
    ckWH(a, ctx, *h.withinHostModel, h.rng());
    Clinical::checkpoint(a, *h.clinicalModel);
    // parity: order — this human's RNG last: everything above may create
    // objects drawing from the fresh copy-seeded stream; the wire restores
    // the real state after those draws.
    a.rng(h.rng());
}

} // namespace

void ckHumans(Archive& a, RunContext& ctx, Population& population) {
    uint64_t n = population.humans.size();
    a.value(n);
    a.value(population.recentBirths);
    if (a.reading()) {
        if (n > uint64_t(a.remaining())) // each human costs > 0 bytes
            throw checkpoint_error("human count exceeds file size");
        for (uint64_t i = 0; i < n; ++i) {
            // Real ctor at the real moment (row append + copy-seed from
            // master); doB is overwritten from the wire immediately after.
            population.humans.emplace_back(ctx, ctx.clock.nowOrTs1(), population.cols);
            ckHuman(a, ctx, population.humans.back());
        }
    } else {
        for (auto& h : population.humans)
            ckHuman(a, ctx, h);
    }
}

} // namespace checkpoint
} // namespace OM
