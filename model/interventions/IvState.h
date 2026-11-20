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

// Phase 7 (7.2a): per-host vector-intervention state = IvRecord columns
// (HostColumns.h), no per-host classes. Bodies here moved verbatim from the
// deleted HumanITN/HumanIRS/HumanGVI (PerHostInterventionData) classes.
// Leaf header: IvRecord + free functions only (implementation pulls in the
// ITN/IRS/GVI parameter registries; 7.2b made the parameters plain public
// structs, so no friend declarations are needed).

#ifndef Hmod_interventions_IvState
#define Hmod_interventions_IvState

#include "Global.h"
#include "humans/HostColumns.h"  // IvRecord

namespace OM {
class RunContext;
namespace Host { class Human; }
namespace interventions {
struct ComponentConfig;  // interventions/Components.h

/** Deploy a vector intervention on a human: the component's existing record
 * is reset in place (redeploy draws), else a record is appended (first-use
 * draws). Draw order pinned (§9); list order is load-bearing (parity: c). */
void iv_deploy( RunContext& ctx, Host::Human& human, const ComponentConfig& params );

/// Call once per time step (old PerHost::update loop; only ITNs mutate).
void iv_update( RunContext& ctx, Host::Human& human );

/// Per-record vector effects, 1.0 if not active (the history-ordered
/// products of these are the PerHost accessors — parity: c).
double iv_rel_attractiveness( RunContext& ctx, const IvRecord& r, size_t species );
double iv_preprandial_survival( RunContext& ctx, const IvRecord& r, size_t species );
double iv_postprandial_survival( RunContext& ctx, const IvRecord& r, size_t species );
double iv_rel_fecundity( RunContext& ctx, const IvRecord& r, size_t species );

/// Record-list checkpointing. Wire NOT frozen (user ruling, phase 7): one
/// uniform record layout replaces the old per-class streams; decay pointer
/// is recovered from the component registry on read.
void iv_checkpoint( RunContext& ctx, std::vector<IvRecord>& iv, ostream& stream );
void iv_checkpoint( RunContext& ctx, std::vector<IvRecord>& iv, istream& stream );

} }
#endif
