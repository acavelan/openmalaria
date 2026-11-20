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

// Checkpoint wire format — table of contents. Read top to bottom:
//
//   header      magic "OMCk" + format version
//   phase       Warmup | Main (what run.cpp skips on resume)
//   clock       t0, t1, interv, endTime, estEndTime
//   statics     (statics.cpp) run-modified globals + imported-infections pos
//   monitoring  (statics.cpp) cts stream offset + survey/report state
//   transmission(transmission.cpp) model kind state, vectors or non-vector
//   humans      (humans.cpp)    living humans' full per-host state
//   master RNG  last: sections above may create objects (Human ctor
//               copy-seeds from master); the wire then overwrites every
//               stream touched.
//
// Everything config/XML-derived is rebuilt by init, never read back.

#include "checkpoint.h"
#include "sections.h"

#include "Context.h"

#include <fstream>

namespace OM {
namespace checkpoint {

namespace {
constexpr uint32_t MAGIC = 0x4B634D4F; // "OMCk" little-endian bytes
constexpr uint32_t FORMAT_VERSION = 1; // bump on any incompatible change

void ckHeader(Archive& a) {
    uint32_t magic = MAGIC, version = FORMAT_VERSION;
    a.value(magic);
    a.value(version);
    if (a.reading()) {
        if (magic != MAGIC)
            throw checkpoint_error("not an OpenMalaria checkpoint file");
        if (version != FORMAT_VERSION)
            throw checkpoint_error("checkpoint format version " + std::to_string(version) + ", this build reads " +
                                   std::to_string(FORMAT_VERSION));
    }
}

void ckClock(Archive& a, Context& ctx, Phase& phase, SimTime& endTime, SimTime& estEndTime) {
    uint8_t p = uint8_t(phase);
    a.value(p);
    if (a.reading()) {
        if (p != uint8_t(Phase::Warmup) && p != uint8_t(Phase::Main))
            throw checkpoint_error("bad checkpoint phase");
        phase = Phase(p);
    }
    a.value(ctx.clock.t0);
    a.value(ctx.clock.t1);
    a.value(ctx.clock.interv);
    a.value(endTime);
    a.value(estEndTime);
}
} // namespace

void writeCheckpoint(Context& ctx, Population& population, Transmission::TransmissionModel& transmission, Phase phase,
                     SimTime endTime, SimTime estEndTime, const std::string& path) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out.is_open())
        throw checkpoint_error("cannot open \"" + path + "\" for writing");
    Archive a(out);
    Phase p = phase;
    ckHeader(a);
    ckClock(a, ctx, p, endTime, estEndTime);
    ckStatics(a, ctx);
    ckMonitoring(a, ctx);
    ckTransmission(a, ctx, transmission);
    ckHumans(a, ctx, population);
    // parity: order — master RNG last (restore above may draw from it; the
    // wire restores the saved state after those draws).
    a.rng(ctx.master_RNG);
    a.checkExhausted(); // no-op on write; keeps sections symmetric
    out.flush();
    if (out.fail())
        throw checkpoint_error("error writing \"" + path + "\"");
}

void readCheckpoint(Context& ctx, Population& population, Transmission::TransmissionModel& transmission, Phase& phase,
                    SimTime& endTime, SimTime& estEndTime, const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open())
        throw checkpoint_error("cannot open checkpoint \"" + path + "\"");
    Archive a(in);
    ckHeader(a);
    ckClock(a, ctx, phase, endTime, estEndTime);
    ckStatics(a, ctx);
    ckMonitoring(a, ctx);
    ckTransmission(a, ctx, transmission);
    ckHumans(a, ctx, population);
    a.rng(ctx.master_RNG);
    a.checkExhausted();
}

} // namespace checkpoint
} // namespace OM
