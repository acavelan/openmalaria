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

#include "util/random.h"
#include <iosfwd>
#include <vector>

namespace OpenMalaria {
class Human;
class Parameters;
class RunContext;

/** Translate per-host EIR into new infections. Per-host state lives directly
 * in Human; model options select one set of equations for the whole run. */
namespace InfectionIncidenceModel {
void init(RunContext& ctx, const Parameters& parameters);
void ctsReportNewInfections(const std::vector<Human>&, std::ostream& stream);
extern int ctsNewInfections;
inline void preMainSimInit() { ctsNewInfections = 0; }

double getAvailabilityFactor(LocalRng& rng, double baseAvailability = 1.0);
void summarize(const Human& human);
double getModelExpectedInfections(Human& human, double effectiveEIR);
double susceptibility(const Human& human);
double expectedNumNewInfections(Human& human, double effectiveEIR);
int numNewInfections(Human& human, double expectedNumNewInfections);
void reportNumNewInfections(Human& human, int newNumInfections_i, int newNumInfections_l);
}

}
