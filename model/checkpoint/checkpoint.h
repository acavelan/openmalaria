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

#ifndef Hmod_checkpoint_h
#define Hmod_checkpoint_h

// Run-state checkpointing (plan v3.3 phase 8-CK, redone xml/-style).
//
// writeCheckpoint() serialises the complete running-model state to one
// plain binary file; readCheckpoint() restores it into a freshly
// XML-initialised run, so a resumed run is bit-identical to an
// uninterrupted one (save at day D + resume == continuous run to the end).
//
// Config rebuilt from the scenario XML never goes on the wire (kind
// switches, drug types, decay curves, component registries...): only true
// run state. Per-host intervention/vaccine/cohort state is deliberately
// NOT serialised (rewrite contract): deploys are deterministic replays from
// the scenario; saves mid-intervention-period are exact only for the
// no-deploy-since-scenario-init cases — see HANDOFF.md limits.
//
// File format: magic + version, then the sections below in this exact
// order (checkpoint.cpp is the table of contents):
//   clock, statics, monitoring, transmission, humans, master RNG (last:
//   every restore step may draw from any stream; the wire overwrites).
// There is no compatibility with the removed v2 (boost-archive/gzstream)
// format — old files simply fail the magic check.

#include "Clock.h"

#include <cstdint>
#include <string>

namespace OM {

class RunContext;
class Population;
namespace Transmission { class TransmissionModel; }

namespace checkpoint {

/// Which simulation phase the snapshot was taken in; drives what run.cpp
/// skips on resume (WARMUP: continue warmup/calibration; MAIN: skip the
/// whole init/warmup/calibration/reset block).
enum class Phase : uint8_t { Warmup = 1, Main = 2 };

/// Write the checkpoint (phase + loop-end times included). Throws
/// checkpoint_error on I/O trouble.
void writeCheckpoint (RunContext& ctx, Population& population,
    Transmission::TransmissionModel& transmission, Phase phase,
    SimTime endTime, SimTime estEndTime, const std::string& path);

/// Restore into the current (XML-initialised) run. On return: clock, all
/// RNG streams, statics, monitoring, transmission and humans match the
/// save point; run.cpp then calls calcAvailabilityPercentiles and either
/// continues the warmup loop (Phase::Warmup) or goes straight to the main
/// loop (Phase::Main). Throws checkpoint_error on a bad file.
/// Outputs: save-time phase and loop-end times.
void readCheckpoint (RunContext& ctx, Population& population,
    Transmission::TransmissionModel& transmission, Phase& phase,
    SimTime& endTime, SimTime& estEndTime, const std::string& path);

} // namespace checkpoint
} // namespace OM
#endif
