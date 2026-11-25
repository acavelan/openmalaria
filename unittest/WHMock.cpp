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

namespace OpenMalaria {
namespace UnitTest {

namespace {
/** Within-host model for unit tests; unneeded operations throw just as the
 * old Kind::Mock dispatch branch did. Each instance owns its test state. */
class WHMock final : public WHInterface {
  public:
    explicit WHMock(RunContext& ctx) : WHInterface(ctx, Kind::Mock) {}

    WHMockState state;

    void treatment(Human&, TreatmentId) override { state.nTreatments += 1; }
    bool treatSimple(Human&, SimTime timeLiver, SimTime timeBlood) override {
        state.nTreatments += 1;
        state.lastTimeLiver = timeLiver;
        state.lastTimeBlood = timeBlood;
        return timeBlood != zero();
    }
    void treatPkPd(size_t schedule, size_t dosages, double age, double delay_d) override {
        state.nTreatments += 1;
        state.pkpd.prescribe(ctx, schedule, dosages, age, numeric_limits<double>::quiet_NaN(), delay_d);
    }
    double getTotalDensity() const override { return state.totalDensity; }
    bool diagnosticResult(LocalRng& rng, const Diagnostic& diagnostic) const override {
        return diagnostic.isPositive(rng, state.totalDensity, numeric_limits<double>::quiet_NaN());
    }

    double probTransmissionToMosquito(vector<double>&, vector<double>&) const override {
        throw unimplemented_exception("not needed in unit test");
    }
    bool summarize(Human&) const override { throw unimplemented_exception("not needed in unit test"); }
    void importInfection(LocalRng&) override { throw unimplemented_exception("not needed in unit test"); }
    void optionalPqTreatment(Human&) override { throw unimplemented_exception("not needed in unit test"); }
    void update(Human&, LocalRng&, int&, int&, vector<double>&, vector<double>&, double) override {
        throw unimplemented_exception("not needed in unit test");
    }
    StatePair determineMorbidity(Human&, double, bool) override {
        throw unimplemented_exception("not needed in unit test");
    }
    void clearImmunity() override { throw unimplemented_exception("not needed in unit test"); }
    double getCumulative_h() const override { throw unimplemented_exception("not needed in unit test"); }
    double getCumulative_Y() const override { throw unimplemented_exception("not needed in unit test"); }
};
}

unique_ptr<WHInterface> makeWHMock(RunContext& ctx) { return std::make_unique<WHMock>(ctx); }

WHMockState& whMockState(WHInterface& mock) { return static_cast<WHMock&>(mock).state; }

}
}
