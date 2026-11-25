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

#ifndef Hmod_PercentileCFRSuite
#define Hmod_PercentileCFRSuite

/* Closes two hard gaps from T2FixtureManifest.md (percentile table, CFR):
 * both functions are reachable through public headers, so the goldens are
 * generated in-process rather than from the binary.
 *
 * 1. PerHostAnophParams::calcAvailabilityPercentiles() — the load-bearing
 *    MasterRng(0,0)-seeded 100-entry ento-availability percentile table
 *    (PerHost.h:82-111; main.cpp calls it before the sim, and interventions
 *    interpolate on it). Pinned bit-exactly for one species with
 *    availability distr=gamma CV=0.3, factor 1.0.
 * 2. ClinicalModel::getCommunityCFR() (ClinicalModel.cpp:115-118) — the
 *    community-vs-hospital CFR map x·θ/(1−cfr+x). It reads the file-global
 *    OM::Clinical::oddsRatioThreshold (ClinicalModel.cpp:45, external
 *    linkage), set directly here instead of via Parameters → no model
 *    changes needed. Goldens evaluated in IEEE-754 double, same op order.
 */

#include <cxxtest/TestSuite.h>
#include "Transmission/PerHost.h"
#include "Clinical/ClinicalModel.h"
#include "UnittestUtil.h"
#include "util/ModelOptions.h"
#include <cstring>
#include <limits>

// file-global defined by ClinicalModel.cpp:45; no header declares it
namespace OM { namespace Clinical { extern double oddsRatioThreshold; } }

// getCommunityCFR is a public static member since the hierarchy collapse
// (phase 4b) — no probe needed.

class PercentileCFRSuite : public CxxTest::TestSuite
{
public:
    void testEntoAvailabilityPercentiles () {
        UnittestUtil::initTime( 5 );
        UnittestUtil::DescriptiveInfection_init();  // = ModelOptions::reset()
        if (Transmission::PerHostAnophParams::numSpecies(UnittestUtil::ctx()) == 0) {
            // minimal in-memory scnXml::Mosq (schema classes are constructible
            // programmatically — pattern used by UnittestUtil for PKPD types)
            scnXml::IntValue rest(3), eip(10);
            scnXml::DoubleValue laidEggs(0.95), seeking(0.15), survFeed(0.62);
            scnXml::SampledValueCV avail;
            avail.setDistr("gamma");
            avail.setCV(0.3);
            scnXml::BetaMeanSample biting(0.11, 0.0), restSite(0.88, 0.0), resting(0.89, 0.0);
            scnXml::DoubleValue oviposit(0.0), bloodIndex(0.9);
            scnXml::Mosq mosq(rest, eip, laidEggs, seeking, survFeed, avail,
                              biting, restSite, resting, oviposit, bloodIndex);
            Transmission::PerHostAnophParams::init(UnittestUtil::ctx(), mosq);
            Transmission::PerHostAnophParams::setEntoAvailabilityFactor(UnittestUtil::ctx(), 0, 1.0);
        }
        Transmission::PerHostAnophParams::calcAvailabilityPercentiles(UnittestUtil::ctx());
        
        // blessed 2026-11-19: MasterRng(0,0) seed, gamma CV=0.3, factor 1.0,
        // 100k samples (GCC 16.2.1, GSL 2.8, -O3 -ffp-contract=off)
        static const double GOLDEN[100] = {
            0, 0.43394887354156364, 0.48155911570088317, 0.51478085081433911, 0.53853386245530854, 0.56093440861090682, 0.57921041632433945, 0.59666331377567738, 0.61096509664026266, 0.62474109576452841,
            0.6373535913417997, 0.64904645249680437, 0.66114385891010929, 0.67212191709311575, 0.68279476608169154, 0.69342067948658759, 0.70316609067292923, 0.7135134741258925, 0.72315806794309312, 0.73218780882108303,
            0.74120973436097926, 0.75014966789447446, 0.75807492438341817, 0.76616752141553079, 0.77464506183506188, 0.78321587827137351, 0.7907223125274212, 0.79858169244512112, 0.80678193246302587, 0.81458861851279285,
            0.82226371207261306, 0.83017429638193441, 0.83799836028067898, 0.84537534727493768, 0.85275169796791461, 0.86000881817830888, 0.86733054996520254, 0.87470711300917425, 0.88179413519243777, 0.88891082675910915,
            0.89619313563513947, 0.90351286463701652, 0.91080274449808851, 0.91800670409113916, 0.92518862945701108, 0.93260841381189308, 0.94021661817552749, 0.94750477186489634, 0.95493516487396268, 0.96229540206068265,
            0.96951641686188383, 0.97639183407398722, 0.98387845816752917, 0.99135903899092703, 0.9987680804319804, 1.0062276318798744, 1.0136798781084249, 1.0215430848026046, 1.0293700640930896, 1.037243887303626,
            1.0451810385817, 1.0537038749681469, 1.062128568587084, 1.0701611982564256, 1.0783178476480444, 1.0867152272751792, 1.095549583548473, 1.1043837209249041, 1.1130157917014556, 1.1219423143187583,
            1.1312605388592833, 1.1408823800578267, 1.1509509251825902, 1.1610677988743474, 1.1715310949914011, 1.1827326358798205, 1.193752817352254, 1.2049620838730206, 1.2158540972819174, 1.2278329202391507,
            1.2400827622396284, 1.2525094016897724, 1.2660278802152989, 1.2797776176399633, 1.2941682694203034, 1.3094082731539498, 1.3244725773136092, 1.3407280966903099, 1.3581208009020067, 1.3771194925603183,
            1.3975358792579338, 1.4195050438190029, 1.4445177393128859, 1.4720539546419755, 1.5033293220923645, 1.5380241756892483, 1.5808360459617221, 1.6340908655830597, 1.7063221886130431, 1.8280769186712187,
        };
        for (int p = 0; p < 100; ++p) {
            double v = Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(UnittestUtil::ctx(),  p );
            if (!(std::memcmp(&v, &GOLDEN[p], sizeof(double)) == 0)) TS_FAIL("percentile mismatch");
        }
        TS_ASSERT_EQUALS( Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(UnittestUtil::ctx(), 100),
                          std::numeric_limits<double>::max() );
        TS_ASSERT_THROWS_ANYTHING( Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(UnittestUtil::ctx(), -1) );
        TS_ASSERT_THROWS_ANYTHING( Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(UnittestUtil::ctx(), 101) );
    }
    
    void testCommunityCFR () {
        Clinical::oddsRatioThreshold = 1.5;
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
            double got = OM::Clinical::ClinicalModel::getCommunityCFR( row[0] );
            if (!(std::memcmp(&got, &row[1], sizeof(double)) == 0)) TS_FAIL("CFR mismatch");
        }
    }
};

#endif
