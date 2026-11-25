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

// Wire sections: statics (run-modified globals) and monitoring (continuous
// output position + survey/report accumulators).
//
// Free wire helpers for the subsystem statics deliberately live here (in
// their home namespaces) so all wire knowledge stays in model/checkpoint/:
// the model headers declare nothing checkpoint-related.

#include "sections.h"

#include "Clock.h"
#include "RunContext.h"
#include "Host/NeonatalMortality.h"
#include "Host/InfectionIncidenceModel.h"
#include "Host/ImportedInfections.h"
#include "Host/WithinHost/Genotypes.h"
#include "Clinical/ClinicalModel.h"
#include "mon/Monitoring.h"

namespace OpenMalaria {

// Run-modified globals, defined in their subsystems' TUs (namespace-scope,
// external linkage; the checkpoint module is the only code outside those TUs
// that touches them):
extern std::vector<int> infantDeaths;            // ClinicalModel.cpp
extern std::vector<int> infantIntervalsAtRisk;   // ditto
extern double riskFromMaternalInfection;         // NeonatalMortality.cpp
extern std::vector<double> prevByGestationalAge; // ditto

/// NeonatalMortality run state (was NeonatalMortality::staticCheckpoint).
void neonatalStatics(Archive& a) {
    a.value(riskFromMaternalInfection);
    a.vector(prevByGestationalAge);
}
/// InfectionIncidenceModel run counter (was InfectionIncidenceModel::staticCheckpoint).
void infectionIncidenceStatics(Archive& a) { a.value(InfectionIncidenceModel::ctsNewInfections); }
/// ImportedInfections rate-table cursor (was ImportedInfections::checkpoint).
void checkpoint(Archive& a, ImportedInfections& ii) { a.value(ii.lastIndex); }

/// InfantMortality death/risk summaries (was InfantMortality::staticCheckpoint).
void infantMortalityStatics(Archive& a) {
    a.vector(infantDeaths);
    a.vector(infantIntervalsAtRisk);
}

void ckStatics(Archive& a, RunContext& ctx) {
    neonatalStatics(a);
    infectionIncidenceStatics(a);
    infantMortalityStatics(a);
    int32_t mode = ctx.genotypes.current_mode; // int on wire (enum domain)
    a.value(mode);
    if (a.reading()) ctx.genotypes.current_mode = Genotypes::SampleMode(mode);
    checkpoint(a, *ctx.importedInfections);
}

namespace {
/// A vector whose length is config-fixed: write contents only, checking the
/// count matches the already-initialised size on read.
template <class T> void fixedVector(Archive& a, std::vector<T>& v) {
    uint64_t n = v.size();
    a.value(n);
    if (a.reading() && n != v.size())
        throw checkpoint_error("monitoring vector: checkpoint has " + std::to_string(n) + " entries, scenario has " +
                               std::to_string(v.size()));
    for (auto& e : v)
        a.value(e);
}
} // namespace

void ckMonitoring(Archive& a, RunContext& ctx) {
    // Continuous output: resume position in the cts file (header written by
    // Continuous::init; rows are flushed per record, so this is exact).
    a.value(ctx.ctsStreamOffset);
    if (a.reading() && ctx.ctsPeriod != zero()) {
        ctx.ctsOStream.seekp(ctx.ctsStreamOffset, std::ios::beg);
        if (ctx.ctsOStream.fail()) throw checkpoint_error("continuous output: can't seek to saved position");
    }

    auto& rt = *ctx.monitoring;
    a.boolean(rt.isInit);
    a.value(rt.surveyIndex);
    a.value(rt.survNumEvent);
    fixedVector(a, rt.cohortSubPopNumbers);
    { // surveyDates: config-sized vector; only .num changes during a run
        uint64_t n = rt.surveyDates.size();
        a.value(n);
        if (a.reading() && n != rt.surveyDates.size()) throw checkpoint_error("survey date count mismatch");
        for (auto& sd : rt.surveyDates)
            a.value(sd.num);
    }
    { // conditions: config-sized vector; only .value changes during a run
        uint64_t n = rt.conditions.size();
        a.value(n);
        if (a.reading() && n != rt.conditions.size()) throw checkpoint_error("condition count mismatch");
        for (auto& c : rt.conditions)
            a.boolean(c.value);
    }
    {
        uint64_t n = rt.stores.size();
        a.value(n);
        if (a.reading() && n != rt.stores.size()) throw checkpoint_error("measure store count mismatch");
        for (auto& store : rt.stores) {
            // variant index: 0 = int reports, 1 = double reports (both runs
            // agree by config; checked anyway)
            int32_t idx = std::holds_alternative<std::vector<int>>(store.reports) ? 0 : 1;
            a.value(idx);
            if (a.reading() && idx != 0 && idx != 1) throw checkpoint_error("bad report store kind");
            std::visit(
                [&a](auto& v) {
                    uint64_t n = v.size();
                    a.value(n);
                    if (a.reading() && n != v.size()) throw checkpoint_error("report store size mismatch");
                    for (auto& e : v)
                        a.value(e);
                },
                store.reports);
        }
    }
}

} // namespace OpenMalaria
