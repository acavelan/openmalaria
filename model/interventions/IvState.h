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

#pragma once

#include "Clock.h"
#include "Host/HostColumns.h" // IvRecord

namespace OpenMalaria {
class RunContext;
class Human;
struct ComponentConfig; // interventions/Components.h

/** Deploy a vector intervention on a human: the component's existing record
 * is reset in place (redeploy draws), else a record is appended (first-use
 * draws). Draw order pinned (§9); list order is load-bearing (parity: c). */
void iv_deploy(RunContext& ctx, Human& human, const ComponentConfig& params);

/// Call once per time step (old PerHost::update loop; only ITNs mutate).
void iv_update(RunContext& ctx, Human& human);

/// Per-record vector effects, 1.0 if not active.
double iv_rel_attractiveness(RunContext& ctx, const IvRecord& r, size_t species);
double iv_preprandial_survival(RunContext& ctx, const IvRecord& r, size_t species);
double iv_postprandial_survival(RunContext& ctx, const IvRecord& r, size_t species);
double iv_rel_fecundity(RunContext& ctx, const IvRecord& r, size_t species);

// parity: c — history-ordered scalar left folds; strict record iteration and
// *= order are load-bearing. The initial value is supplied by the caller.
inline double iv_rel_attractiveness(RunContext& ctx, const std::vector<IvRecord>& iv, double alpha_i,
                                    size_t species) {
    for (const IvRecord& r : iv)
        alpha_i *= iv_rel_attractiveness(ctx, r, species);
    return alpha_i;
}
inline double iv_preprandial_survival(RunContext& ctx, const std::vector<IvRecord>& iv, double P_B_i,
                                      size_t species) {
    for (const IvRecord& r : iv)
        P_B_i *= iv_preprandial_survival(ctx, r, species);
    return P_B_i;
}
inline double iv_postprandial_survival(RunContext& ctx, const std::vector<IvRecord>& iv, double pRest,
                                       size_t species) {
    for (const IvRecord& r : iv)
        pRest *= iv_postprandial_survival(ctx, r, species);
    return pRest;
}
inline double iv_rel_fecundity(RunContext& ctx, const std::vector<IvRecord>& iv, double relFecundity,
                               size_t species) {
    for (const IvRecord& r : iv)
        relFecundity *= iv_rel_fecundity(ctx, r, species);
    return relFecundity;
}

}
