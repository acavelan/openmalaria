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

// Module-private section functions (one per wire section; the Archive
// carries the direction — write each field exactly once, in this order:
// see checkpoint.cpp for the full table of contents).

#include "Stream.h"

namespace OpenMalaria {
class RunContext;
class Population;
class TransmissionModel;

void ckStatics(Archive&, RunContext&);
void ckMonitoring(Archive&, RunContext&);
void ckTransmission(Archive&, RunContext&, TransmissionModel&);
void ckHumans(Archive&, RunContext&, Population&);

} // namespace OpenMalaria
