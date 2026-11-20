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

// Phase-3 acceptance test (plan §7 row 3): run the full engine twice, in one
// process, through the same code path main() uses (OM::run_simulation).
// Run 2 must be bit-identical to run 1 — any per-run state that survives run
// 1 (i.e. any static that should live in Context) makes this red.
//
// Capture: the monitoring MeasureStore report values (exact IEEE-754 bit
// comparison, NaN included) plus the survey-output bytes written by
// mon::writeSurveyData(). scenario1.xml: pop 100, timed interventions,
// dynamic EIR — full run ≈ 0.3 s.

#include <cxxtest/TestSuite.h>

#include "configured/TestPaths.h"

#include "Context.h"
#include "util/CommandLine.h"
#include "mon/Monitoring.h"

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace {

/// Raw IEEE-754 bit pattern of every monitoring report value, in store order
/// (contract (g): store order is deterministic). NaN-safe: compares bits.
void captureMonStores(OM::Context& ctx, std::vector<uint64_t>& out) {
    for (const OM::mon::internal::MeasureStore& store : (*ctx.monitoring).stores) {
        out.push_back(uint64_t(INT64_MAX) + 1 + static_cast<uint64_t>(store.output.measure));
        out.push_back(store.size());
        std::visit(
            [&out](const auto& reports) {
                for (const auto v : reports) {
                    double d = static_cast<double>(v);
                    uint64_t bits;
                    std::memcpy(&bits, &d, sizeof bits);
                    out.push_back(bits);
                }
            },
            store.reports);
    }
}

struct RunCapture {
    std::vector<uint64_t> monValues;
    std::vector<char> surveyBytes;
};

RunCapture captureRun(Context& ctx, const std::string& surveyFile) {
    RunCapture c;
    captureMonStores(ctx, c.monValues);
    std::ifstream in(surveyFile, std::ios::in | std::ios::binary);
    c.surveyBytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    return c;
}

} // namespace

class TwoRunsSuite : public CxxTest::TestSuite {
  public:
    void testTwoRunsInOneProcessAreBitIdentical() {
        const std::string src(UnittestSourceDir); // <repo>/unittest/
        const std::string scenario = src + "../test/scenario1.xml";
        const std::string xsd = src + "../schema/scenario_50.xsd";
        const std::string resources = src + "../test";
        TS_ASSERT(std::filesystem::exists(scenario));
        TS_ASSERT(std::filesystem::exists(xsd));

        // isolated cwd: survey/ctsout files land here; XSD must sit next to cwd
        std::error_code ec;
        const std::filesystem::path oldCwd = std::filesystem::current_path();
        std::filesystem::path work = std::filesystem::temp_directory_path() / "om-tworuns";
        std::filesystem::remove_all(work, ec);
        std::filesystem::create_directories(work);
        std::filesystem::copy_file(xsd, work / "scenario_current.xsd");
        std::filesystem::current_path(work);

        const char* argv[] = {"unittest",        "--scenario", scenario.c_str(), "--resource-path", resources.c_str(),
                              "--output-format", "bin",        "--output",       "tworuns.bin"};
        std::string scenarioFile;
        try {
            scenarioFile = OM::util::CommandLine::parse(5, const_cast<char**>(argv));
        } catch (const std::exception& e) {
            std::filesystem::current_path(oldCwd);
            TS_FAIL(std::string("CommandLine::parse: ") + e.what());
            return;
        }

        RunCapture run1, run2;
        bool ok = true;
        std::string err;
        for (int run = 1; run <= 2 && ok; ++run) {
            try {
                Context ctx; // each run gets a fresh context: runs must not share state
                OM::run_simulation(ctx, scenarioFile);
                (run == 1 ? run1 : run2) = captureRun(ctx, "tworuns.bin");
            } catch (const std::exception& e) {
                ok = false;
                err = "run " + std::to_string(run) + ": " + e.what();
            }
        }
        std::filesystem::current_path(oldCwd);

        if (!ok) {
            TS_FAIL(err);
            return;
        }

        TS_ASSERT_EQUALS(run1.monValues.size(), run2.monValues.size());
        size_t firstDiff = run1.monValues.size();
        for (size_t i = 0; i < std::min(run1.monValues.size(), run2.monValues.size()); ++i) {
            if (run1.monValues[i] != run2.monValues[i]) {
                firstDiff = i;
                break;
            }
        }
        if (firstDiff < run1.monValues.size()) {
            const auto hex = [](uint64_t v) {
                char b[19];
                std::snprintf(b, sizeof b, "0x%016llx", static_cast<unsigned long long>(v));
                return std::string(b);
            };
            TS_FAIL("monitoring store value differs at #" + std::to_string(firstDiff) +
                    ": run1=" + hex(run1.monValues[firstDiff]) + " run2=" + hex(run2.monValues[firstDiff]));
        }
        TS_ASSERT_EQUALS(run1.surveyBytes.size(), run2.surveyBytes.size());
        if (!run1.surveyBytes.empty()) { // sanity: run 1 produced output at all
            TS_ASSERT(run1.surveyBytes == run2.surveyBytes);
            // in-memory value capture must be non-empty too, else the test proves nothing
            TS_ASSERT(!run1.monValues.empty());
        }
    }
};
