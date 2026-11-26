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

/* Closes two hard gaps from T2FixtureManifest.md (percentile table, CFR):
 * both functions are reachable through public headers, so the goldens are
 * generated in-process rather than from the binary.
 *
 * 1. PerHostAnophParams::calcAvailabilityPercentiles() — the load-bearing
 *    MasterRng(0,0)-seeded 100-entry ento-availability percentile table
 *    (PerHost.h:82-111; main.cpp calls it before the sim, and interventions
 *    interpolate on it). Pinned bit-exactly for one species with
 *    availability distr=gamma CV=0.3, factor 1.0. Bit-exact against both
 *    the same-platform GSL reference and the stored table.
 * 2. ClinicalModel::getCommunityCFR() (ClinicalModel.cpp:115-118) — the
 *    community-vs-hospital CFR map x·θ/(1−cfr+x). It reads the file-global
 *    oddsRatioThreshold (ClinicalModel.cpp:45, external
 *    linkage), set directly here instead of via Parameters → no model
 *    changes needed. Goldens evaluated in IEEE-754 double, same op order.
 */

#include <cxxtest/TestSuite.h>
#include "Transmission/PerHost.h"
#include "Clinical/ClinicalModel.h"
#include "UnittestUtil.h"
#include "util/ModelOptions.h"
#include "golden-gen/RngGoldenAdapter.h"
#include <cstring>
#include <limits>

// file-global defined by ClinicalModel.cpp:45; no header declares it
namespace OpenMalaria {
extern double oddsRatioThreshold;
}

// getCommunityCFR is a public static member since the hierarchy collapse
// (phase 4b) — no probe needed.

class PercentileCFRSuite : public CxxTest::TestSuite {
  public:
    void testNewClinicalModelState() {
        auto clinical = std::make_unique<CMImmediateOutcomes>(UnittestUtil::ctx(), 1.0);
        TS_ASSERT_EQUALS(clinical->latestState, NONE);
    }

    void testEntoAvailabilityPercentiles() {
        UnittestUtil::initTime(5);
        UnittestUtil::DescriptiveInfection_init(); // = ModelOptions::reset()
        if (PerHostAnophParams::numSpecies(UnittestUtil::ctx()) == 0) {
            // minimal in-memory scnXml::Mosq (schema classes are constructible
            // programmatically — pattern used by UnittestUtil for PKPD types)
            scnXml::IntValue rest(3), eip(10);
            scnXml::DoubleValue laidEggs(0.95), seeking(0.15), survFeed(0.62);
            scnXml::SampledValueCV avail;
            avail.setDistr("gamma");
            avail.setCV(0.3);
            scnXml::BetaMeanSample biting(0.11, 0.0), restSite(0.88, 0.0), resting(0.89, 0.0);
            scnXml::DoubleValue oviposit(0.0), bloodIndex(0.9);
            scnXml::Mosq mosq(rest, eip, laidEggs, seeking, survFeed, avail, biting, restSite, resting, oviposit,
                              bloodIndex);
            UnittestUtil::ctx().anophParams->emplace_back(mosq);
            PerHostAnophParams::setEntoAvailabilityFactor(UnittestUtil::ctx(), 0, 1.0);
        }
        PerHostAnophParams::calcAvailabilityPercentiles(UnittestUtil::ctx());

        // Frozen original-model capture for the explicitly selected platform.
        // MasterRng(0,0), gamma CV=0.3, factor 1.0, 100k samples.
        static const double GOLDEN[100] = {
#include "configured/availability_percentiles.inc"
        };
        // Independent reference for the fixed gamma fixture: baseline seed,
        // k=1/CV², theta=mean/k, 100k draws, sorted order statistics.
        // Local GSL/libm rounding is shared; model sampling and percentile
        // selection must agree bit-for-bit, even off the CSV's platform.
        MasterRng master(0, 0);
        LocalRng initial(master);
        const auto& s = initial.m_rng.s;
        GoldenRng::Rng reference(s[0], s[1], s[2], s[3]);
        const double k = 1.0 / (0.3 * 0.3);
        std::vector<double> samples(100000);
        for (double& x : samples) x = reference.gamma(k, 1.0 / k);
        std::sort(samples.begin(), samples.end());
        for (int p = 0; p < 100; ++p) {
            double v = PerHostAnophParams::getEntoAvailabilityPercentile(UnittestUtil::ctx(), p);
            double expected = p == 0 ? 0.0 : samples[p * 1000];
            TSM_ASSERT_SAME_DATA("same-platform percentile mismatch", &v, &expected, sizeof(double));
            TSM_ASSERT_SAME_DATA("stored percentile mismatch", &v, &GOLDEN[p], sizeof(double));
        }
        TS_ASSERT_EQUALS(PerHostAnophParams::getEntoAvailabilityPercentile(UnittestUtil::ctx(), 100),
                         std::numeric_limits<double>::max());
        TS_ASSERT_THROWS_ANYTHING(PerHostAnophParams::getEntoAvailabilityPercentile(UnittestUtil::ctx(), -1));
        TS_ASSERT_THROWS_ANYTHING(PerHostAnophParams::getEntoAvailabilityPercentile(UnittestUtil::ctx(), 101));
    }

    void testCommunityCFR() {
        oddsRatioThreshold = 1.5;
        // golden[cfr] = x/((1-cfr)+x), x=cfr*1.5 — evaluated 2026-11-19 in
        // IEEE-754 double with the op order of ClinicalModel.cpp:116-117
        static const double TAB[][2] = {
            {0.0, 0.0},
            {0.01, 0.014925373134328359},
            {0.1, 0.14285714285714288},
            {0.3, 0.39130434782608697},
            {0.5, 0.59999999999999998},
            {0.9, 0.93103448275862066},
            {1.0, 1.0},
        };
        for (const auto& row : TAB) {
            double got = ClinicalModel::getCommunityCFR(row[0]);
            if (!(std::memcmp(&got, &row[1], sizeof(double)) == 0))
                TS_FAIL("CFR mismatch");
        }
    }
};
