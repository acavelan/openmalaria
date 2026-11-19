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

#include "WHMock.h"
#include "core/RunExtensions.h"
#include "Host/WithinHost/Diagnostic.h"

using namespace std;

namespace OM {
namespace UnitTest {

using util::LocalRng;

namespace {
WHMockState* state = nullptr;

void treatment(void*){
    state->nTreatments += 1;
}

bool treatSimple(void*, SimTime timeLiver, SimTime timeBlood){
    state->nTreatments += 1;
    state->lastTimeLiver = timeLiver;
    state->lastTimeBlood = timeBlood;
    return timeBlood != sim::zero();
}

void treatPkPd(void*, RunContext& ctx, size_t schedule, size_t dosages, double age, double delay_d){
    state->nTreatments += 1;
    state->pkpd.prescribe( ctx, schedule, dosages, age, numeric_limits<double>::quiet_NaN(), delay_d );
}

double getTotalDensity(const void*){
    return state->totalDensity;
}

bool diagnosticResult(const void*, LocalRng& rng, const WithinHost::Diagnostic& diagnostic){
    return diagnostic.isPositive( rng, state->totalDensity, numeric_limits<double>::quiet_NaN() );
}

RunExtensions::WHMockHooks hooks = { nullptr, &treatment, &treatSimple, &treatPkPd, &getTotalDensity, &diagnosticResult };
}

unique_ptr<WithinHost::WHInterface> makeWHMock(RunContext& ctx){
    delete state;
    state = new WHMockState();
    hooks.state = state;
    ctx.ext->whMock = &hooks;
    LocalRng rng(0, 0);   // unused: the Mock kind's constructor draws nothing
    return unique_ptr<WithinHost::WHInterface>(
        new WithinHost::WHInterface( ctx, rng, 1.0, WithinHost::WHInterface::Kind::Mock ));
}

}
}
