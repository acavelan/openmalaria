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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "WHMock.h"
#include "Host/WithinHost/Diagnostic.h"
#include "util/errors.h"

using namespace std;

namespace OM {
namespace UnitTest {

using util::LocalRng;

namespace {
WHMockState* state = nullptr;

void treatment(void*) { state->nTreatments += 1; }

bool treatSimple(void*, SimTime timeLiver, SimTime timeBlood) {
    state->nTreatments += 1;
    state->lastTimeLiver = timeLiver;
    state->lastTimeBlood = timeBlood;
    return timeBlood != sim::zero();
}

void treatPkPd(void*, RunContext& ctx, size_t schedule, size_t dosages, double age, double delay_d) {
    state->nTreatments += 1;
    state->pkpd.prescribe(ctx, schedule, dosages, age, numeric_limits<double>::quiet_NaN(), delay_d);
}

double getTotalDensity(const void*) { return state->totalDensity; }

bool diagnosticResult(const void*, LocalRng& rng, const WithinHost::Diagnostic& diagnostic) {
    return diagnostic.isPositive(rng, state->totalDensity, numeric_limits<double>::quiet_NaN());
}

RunContext::WHMockHooks hooks = {nullptr, &treatment, &treatSimple, &treatPkPd, &getTotalDensity, &diagnosticResult};

/** Within-host model for unit tests (see WHMockState): treatment/treatSimple/
 * treatPkPd/density/diagnostic run through the hook table on RunContext
 * (installed by makeWHMock); everything else behaves exactly as the old
 * Kind::Mock dispatch branch did — throws unimplemented_exception. No
 * infections, no draws, no state beyond the hook table. */
class WHMock final : public WithinHost::WHInterface {
  public:
    explicit WHMock(RunContext& ctx) : WHInterface(ctx, Kind::Mock) {}

    void treatment(Host::Human&, WithinHost::TreatmentId) override {
        const RunContext::WHMockHooks* h = ctx.whMock;
        h->treatment(h->state);
    }
    bool treatSimple(Host::Human&, SimTime timeLiver, SimTime timeBlood) override {
        const RunContext::WHMockHooks* h = ctx.whMock;
        return h->treatSimple(h->state, timeLiver, timeBlood);
    }
    void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d) override {
        const RunContext::WHMockHooks* h = ctx.whMock;
        h->treatPkPd(h->state, ctx, schedule, dosages, age, delay_d);
    }
    double getTotalDensity() const override {
        const RunContext::WHMockHooks* h = ctx.whMock;
        return h->getTotalDensity(h->state);
    }
    bool diagnosticResult(LocalRng& rng, const WithinHost::Diagnostic& diagnostic) const override {
        const RunContext::WHMockHooks* h = ctx.whMock;
        return h->diagnosticResult(h->state, rng, diagnostic);
    }

    double probTransmissionToMosquito(vector<double>&, vector<double>&) const override {
        throw util::unimplemented_exception("not needed in unit test");
    }
    bool summarize(Host::Human&) const override { throw util::unimplemented_exception("not needed in unit test"); }
    void importInfection(LocalRng&) override { throw util::unimplemented_exception("not needed in unit test"); }
    void optionalPqTreatment(Host::Human&) override { throw util::unimplemented_exception("not needed in unit test"); }
    void update(Host::Human&, LocalRng&, int&, int&, vector<double>&, vector<double>&, double) override {
        throw util::unimplemented_exception("not needed in unit test");
    }
    WithinHost::Pathogenesis::StatePair determineMorbidity(Host::Human&, double, bool) override {
        throw util::unimplemented_exception("not needed in unit test");
    }
    void clearImmunity() override { throw util::unimplemented_exception("not needed in unit test"); }
    double getCumulative_h() const override { throw util::unimplemented_exception("not needed in unit test"); }
    double getCumulative_Y() const override { throw util::unimplemented_exception("not needed in unit test"); }
};
}

unique_ptr<WithinHost::WHInterface> makeWHMock(RunContext& ctx) {
    delete state;
    state = new WHMockState();
    hooks.state = state;
    ctx.whMock = &hooks;
    return unique_ptr<WithinHost::WHInterface>(new WHMock(ctx));
}

}
}
