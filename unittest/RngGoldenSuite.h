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

#ifndef Hmod_RngGoldenSuite
#define Hmod_RngGoldenSuite

/* T2 fixture: replays unittest/goldens/rng_dist_v50.csv through the same calls
 * that produced it (golden-gen/RngGoldenAdapter.h — a verbatim mirror of
 * model/util/random.h's GSL plumbing and distribution wrappers over
 * Xoshiro256P), and asserts bit-exact equality of every drawn value AND of the
 * 4 post-draw xoshiro state words.
 *
 * This locks the primitive layer (seeding, state advance, draw counts per
 * sampler call). When phase 9 replaces the GSL samplers with inline ports,
 * this suite is the gate: any change in per-call consumption or output bits
 * fails here, per distribution and parameter set.
 *
 * Regenerate goldens only via unittest/golden-gen/build-goldens.sh (never by
 * editing the CSV). Comparison is bitwise (memcmp), never tolerance-based.
 */

#include <cxxtest/TestSuite.h>
#include "golden-gen/RngGoldenAdapter.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

class RngGoldenSuite : public CxxTest::TestSuite
{
public:
    static const char* CSV_FILE;
    
    /// next non-comment, non-empty line
    static bool nextDataRow(std::ifstream& f, std::string& out) {
        std::string l;
        while (std::getline( f, l )) {
            if (l.empty() || l[0] == '#') continue;
            out = l;
            return true;
        }
        return false;
    }
    
    void testGoldenFilePresent() {
        std::ifstream f( CSV_FILE );
        if (!(f.is_open())) TS_FAIL("rng_dist_v50.csv not found next to the unittest binary (see CMakeLists configure_file)");
        // expected data-row count: all schedule series + Pinv spots
        size_t nSeries;
        const GoldenRng::DrawSpec* S = GoldenRng::schedule( nSeries );
        size_t want = 0;
        for (size_t i = 0; i < nSeries; ++i) want += S[i].n;
        size_t nSpots;
        GoldenRng::pinv_spots( nSpots );
        want += nSpots;
        
        size_t got = 0;
        std::string line;
        while (std::getline( f, line )) {
            if (!line.empty() && line[0] != '#') ++got;  // '#' = header comment
        }
        if (!(got == want)) TS_FAIL("CSV data-row count differs from the replay schedule");
    }
    
    void testBitExactReplay() {
        std::ifstream f( CSV_FILE );
        if (!(f.is_open())) TS_FAIL("rng_dist_v50.csv not found");
        
        std::string line;
        size_t nSeries;
        const GoldenRng::DrawSpec* S = GoldenRng::schedule( nSeries );
        for (size_t i = 0; i < nSeries; ++i) {
            const GoldenRng::DrawSpec& spec = S[i];
            GoldenRng::Rng rng( GoldenRng::SEED[0], GoldenRng::SEED[1],
                                 GoldenRng::SEED[2], GoldenRng::SEED[3] );
            for (int d = 0; d < spec.n; ++d) {
                if (!nextDataRow( f, line )) TS_FAIL("CSV ended early");
                char dist[64], params[128], hex[4][17];
                int idx;
                char value[64];
                int nFld = std::sscanf( line.c_str(), "%63[^,],%127[^,],%d,%63[^,],%16[0-9a-f],%16[0-9a-f],%16[0-9a-f],%16[0-9a-f]",
                    dist, params, &idx, value, hex[0], hex[1], hex[2], hex[3] );
                if (!(nFld == 8)) TS_FAIL("malformed CSV row");
                if (nFld != 8) continue;
                if (!(std::strcmp( dist, spec.dist ) == 0)) TS_FAIL("CSV row out of order: dist mismatch");
                if (!(std::strcmp( params, spec.params ) == 0)) TS_FAIL("CSV row out of order: params mismatch");
                if (!(idx == d)) TS_FAIL("CSV row out of order: draw_index mismatch");
                
                double want = std::strtod( value, nullptr );  // %.17g round-trips exactly
                uint64_t wantState[4];
                for (int w = 0; w < 4; ++w) wantState[w] = std::strtoull( hex[w], nullptr, 16 );
                
                double v = GoldenRng::golden_draw( rng, spec );
                uint64_t st[4];
                rng.state( st );
                
                char msg[256];
                std::snprintf( msg, sizeof(msg), "%s(%s) draw %d: value bits differ", spec.dist, spec.params, d );
                if (!(std::memcmp( &v, &want, sizeof(double) ) == 0)) TS_FAIL(msg);
                for (int w = 0; w < 4; ++w) {
                    std::snprintf( msg, sizeof(msg), "%s(%s) draw %d: state word %d differs", spec.dist, spec.params, d, w );
                    if (!(st[w] == wantState[w])) TS_FAIL(msg);
                }
            }
        }
        
        // Pinv spots: deterministic quantile evaluations, no state consumed.
        size_t nSpots;
        const double* U = GoldenRng::pinv_spots( nSpots );
        for (size_t i = 0; i < nSpots; ++i) {
            if (!nextDataRow( f, line )) TS_FAIL("CSV ended early (Pinv block)");
            char dist[64], params[128], hex[4][17];
            int idx;
            char value[64];
            int nFld = std::sscanf( line.c_str(), "%63[^,],%127[^,],%d,%63[^,],%16[0-9a-f],%16[0-9a-f],%16[0-9a-f],%16[0-9a-f]",
                dist, params, &idx, value, hex[0], hex[1], hex[2], hex[3] );
            if (!(nFld == 8)) TS_FAIL("malformed Pinv CSV row");
            if (nFld != 8) continue;
            if (!(std::strcmp( dist, "ugaussian_Pinv" ) == 0)) TS_FAIL("Pinv row out of order");
            double want = std::strtod( value, nullptr );
            double v = gsl_cdf_ugaussian_Pinv( U[i] );
            char msg[128];
            std::snprintf( msg, sizeof(msg), "ugaussian_Pinv spot %zu (%s): value bits differ", i, params );
            if (!(std::memcmp( &v, &want, sizeof(double) ) == 0)) TS_FAIL(msg);
        }
    }
};

const char* RngGoldenSuite::CSV_FILE = "rng_dist_v50.csv";

#endif
