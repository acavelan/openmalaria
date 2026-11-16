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

namespace OM::Transmission::Anopheles {

struct SugarFeedingInputs
{
    double ah; // The availability rate of hosts h. Time-1
    double an; // The availability rate of hosts n. Time-1
    double mu_v_interv; // Mosquito mortality rate while host-seeking. Time-1
    double theta_d; // Proportion of the night that the mosquito spends host-seeking. Time.
    double PAt_wanted; // PAt is the daily probability of a mosquito feeding on an ATSB
};

// Solve for the availability rate giving the requested sugar-feeding probability.
double sugarFeedingAvailability(SugarFeedingInputs inputs);

}
