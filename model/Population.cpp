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
#include "Population.h"
#include "mon/Continuous.h"

#include "Host/Human.h"
#include "Host/NeonatalMortality.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Clinical/ClinicalModel.h"
#include "Transmission/TransmissionModel.h"

#include "util/errors.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"

#include <cmath>

namespace OpenMalaria {
using namespace std::placeholders;

/** Variables for continuous reporting */
const vector<double> ctsDemogAgeGroups = {1.0, 5.0, 10.0, 15.0, 25.0};

// -----  non-static methods: creation/destruction, checkpointing  -----

Population::Population(RunContext& ctx, size_t size) : ctx(ctx), size(size) {}

void Population::createInitialHumans() {
    /* We create a whole population here, regardless of whether humans can
    survive until start of vector init (vector model needs a whole population
    structure in any case). However, we don't update humans known not to survive
    until vector init, which saves computation and memory (no infections). */

    int cumulativePop = 0;
    for (size_t iage_prev = AgeStructure::getMaxTStepsPerLife(ctx), iage = iage_prev - 1; iage_prev > 0;
         iage_prev = iage, iage -= 1) {
        int targetPop = AgeStructure::targetCumPop(ctx, iage, size);
        while (cumulativePop < targetPop) {
            SimTime dob = zero() - ctx.clock.fromTS(iage);
            streamValidate(dob);
            humans.push_back(Human(ctx, dob, cols));
            ++cumulativePop;
        }
    }

    // Vector setup dependant on human population structure (we *want* to
    // include all humans, whether they'll survive to vector init phase or not).
    assert(ctx.clock.now() == zero()); // assumed below

    // Huge-page hint for the heap segments just grown by the appends above
    // (best-effort, see Host/HostColumns.h).
    advise_heap_huge();
}

void Population::update() {
    ++g_compact_stats.pop_updates; // TEMP instrumentation
    //NOTE: other parts of code are not set up to handle changing population size. Also
    // size is assumed to be the _actual and exact_ population size by other code.
    int cumPop = 0;

    // NOTE: this pass only decides (reads per-host state through POSITION,
    // cols[i]); Human::row goes stale for rows the compaction below moves and
    // is fixed by the renumber pass, so don't use Human accessors here.
    std::vector<uint32_t> keep;
    keep.reserve(humans.size());
    for (size_t i = 0; i < humans.size(); ++i) {
        bool isDead = cols.death[i] != never();

        // if (Actual number of people so far > target population size for this age)
        // "outmigrate" some to maintain population shape
        //NOTE: better to use age(ctx.clock.ts0())? Possibly, but the difference will not be very significant.
        // Also see targetPop = ... comment above
        SimTime age1 = ctx.clock.ts1() - cols.dob[i];
        bool outmigrate = cumPop >= AgeStructure::targetCumPop(ctx, ctx.clock.inSteps(age1), size);

        if (!(isDead || outmigrate)) {
            keep.push_back(static_cast<uint32_t>(i));
            ++cumPop;
        }
    } // end of per-human updates

    if (keep.size() != humans.size()) {
        // TEMP instrumentation
        ++g_compact_stats.steps_with_removals;
        g_compact_stats.removals_total += humans.size() - keep.size();
        // One stable batched compaction == the old per-death erase_at +
        // humans.erase (same survivor order), but each column shifts once.
        // The w != r guard avoids self-move (empties nested containers).
        cols.compact(keep);
        size_t w = 0;
        for (const uint32_t r : keep) {
            if (w != r) humans[w] = std::move(humans[r]);
            ++w;
        }
        humans.erase(humans.begin() + w, humans.end()); // no default ctor

        // Compaction shifted rows left; renumber so Human::row matches
        // position again (newborns below append at the correct index
        // themselves).
        for (size_t i = 0; i < humans.size(); ++i)
            humans[i].set_row(i);
    }

    // increase population size to targetPop
    recentBirths += (size - cumPop);
    while (cumPop < (int)size) {
        // humans born at end of this time step = beginning of next, hence ts1
        humans.push_back(Human(ctx, ctx.clock.ts1(), cols));
        ++cumPop;
    }

    // Yearly re-scan for heap segments grown since the last pass (arena
    // extensions, monitoring buffers); the madvise hint is idempotent and
    // costs ~1 ms, so doing it every week isn't worth the tracking (see
    // Host/HostColumns.h).
    if (mod_nn(ctx.clock.now(), 52) == 0) advise_heap_huge();
}

void registerContinousPopulationCallbacks(RunContext& ctx) {
    ostringstream ctsDemogTitle;
    for (double ubound : ctsDemogAgeGroups)
        ctsDemogTitle << "\thost % ≤ " << ubound;

    Continuous::registerCallback(ctx, "hosts", "\thosts", &ctsHosts);
    Continuous::registerCallback(ctx, "host demography", ctsDemogTitle.str(),
                                 [&ctx](Population& pop, ostream& s) { ctsHostDemography(ctx, pop, s); });
    Continuous::registerCallback(ctx, "recent births", "\trecent births", &ctsRecentBirths);
    Continuous::registerCallback(ctx, "patent hosts", "\tpatent hosts",
                                 [&ctx](Population& pop, ostream& s) { ctsPatentHosts(ctx, pop, s); });
    Continuous::registerCallback(ctx, "immunity h", "\timmunity h", &ctsImmunityh);
    Continuous::registerCallback(ctx, "immunity Y", "\timmunity Y", &ctsImmunityY);
    Continuous::registerCallback(ctx, "median immunity Y", "\tmedian immunity Y", &ctsMedianImmunityY);
    Continuous::registerCallback(ctx, "human age availability", "\thuman age availability",
                                 [&ctx](Population& pop, ostream& s) { ctsMeanAgeAvailEffect(ctx, pop, s); });
    Continuous::registerCallback(ctx, "ITN coverage", "\tITN coverage", &ctsITNCoverage);
    Continuous::registerCallback(ctx, "IRS coverage", "\tIRS coverage", &ctsIRSCoverage);
    Continuous::registerCallback(ctx, "GVI coverage", "\tGVI coverage", &ctsGVICoverage);
}

void ctsHosts(Population& population, ostream& stream) {
    // this option is intended for debugging human initialization; normally this should equal size.
    stream << '\t' << population.size;
}

void ctsHostDemography(RunContext& ctx, Population& population, ostream& stream) {
    auto iter = population.humans.crbegin();
    int cumCount = 0;
    for (double ubound : ctsDemogAgeGroups) {
        while (iter != population.humans.crend() && inYears(iter->age(ctx.clock.now())) < ubound) {
            ++cumCount;
            ++iter;
        }
        stream << '\t' << cumCount;
    }
}

void ctsRecentBirths(Population& population, ostream& stream) {
    stream << '\t' << population.recentBirths;
    population.recentBirths = 0;
}

void ctsPatentHosts(RunContext& ctx, Population& population, ostream& stream) {
    int patent = 0;
    for (Human& human : population.humans) {
        auto diag = diagnostics::monitoringDiagnostic(ctx);
        if (human.withinHostModel->diagnosticResult(human.rng(), diag))
            ++patent;
    }
    stream << '\t' << patent;
}

void ctsImmunityh(Population& population, ostream& stream) {
    double x = 0.0;
    for (const Human& human : population.humans) {
        x += human.withinHostModel->getCumulative_h();
    }
    x /= population.size;
    stream << '\t' << x;
}

void ctsImmunityY(Population& population, ostream& stream) {
    double x = 0.0;
    for (const Human& human : population.humans) {
        x += human.withinHostModel->getCumulative_Y();
    }
    x /= population.size;
    stream << '\t' << x;
}

void ctsMedianImmunityY(Population& population, ostream& stream) {
    vector<double> list;
    list.reserve(population.size);
    for (const Human& human : population.humans) {
        list.push_back(human.withinHostModel->getCumulative_Y());
    }
    sort(list.begin(), list.end());
    double x;
    if (mod_nn(population.size, 2) == 0) {
        size_t i = population.size / 2;
        x = (list[i - 1] + list[i]) / 2.0;
    } else {
        x = list[population.size / 2];
    }
    stream << '\t' << x;
}

void ctsMeanAgeAvailEffect(RunContext& ctx, Population& population, ostream& stream) {
    int nHumans = 0;
    double avail = 0.0;
    for (Human& human : population.humans) {
        if (!human.perHostTransmission.outside_trans()) {
            ++nHumans;
            avail += human.perHostTransmission.relativeAvailabilityAge(inYears(human.age(ctx.clock.now())));
        }
    }
    stream << '\t' << avail / nHumans;
}

void ctsITNCoverage(Population& population, ostream& stream) {
    int nActive = 0;
    for (const Human& human : population.humans) {
        nActive += human.perHostTransmission.hasActiveInterv(Component::ITN);
    }
    double coverage = static_cast<double>(nActive) / population.size;
    stream << '\t' << coverage;
}

void ctsIRSCoverage(Population& population, ostream& stream) {
    int nActive = 0;
    for (const Human& human : population.humans) {
        nActive += human.perHostTransmission.hasActiveInterv(Component::IRS);
    }
    double coverage = static_cast<double>(nActive) / population.size;
    stream << '\t' << coverage;
}

void ctsGVICoverage(Population& population, ostream& stream) {
    int nActive = 0;
    for (const Human& human : population.humans) {
        nActive += human.perHostTransmission.hasActiveInterv(Component::GVI);
    }
    double coverage = static_cast<double>(nActive) / population.size;
    stream << '\t' << coverage;
}

}
