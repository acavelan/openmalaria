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
#include "util/AgeGroupInterpolation.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"

namespace OM::WithinHost {
struct WithinHostParameters {
    enum class Kind { Descriptive, Dummy, Empirical, Molineaux, Penny, Vivax };
    Kind kind = Kind::Descriptive;
    Pathogenesis::PathogenesisParameters pathogenesis;
    bool reportInfectionsByGenotype = false;
    bool reportInfectedOrPatentInfected = false;
    bool reportPatentInfected = false;
    bool opt_vaccine_genotype = false;

    double hetMassMultStdDev = numeric_limits<double>::signaling_NaN();
    double minHetMassMult = numeric_limits<double>::signaling_NaN();
    util::AgeGroupInterpolator massByAge;

    // Innate immunity standard deviation and fever-associated acquired-immunity penalty.
    double sigma_i = 0.0;
    double immPenalty_22 = 0.0;
    // Remaining effect on densities and remaining immune effectors per step, respectively.
    double asexImmRemain = 0.0;
    double immEffectorRemain = 0.0;
    // Inverse cumulative-density and inoculation thresholds for acquired immunity.
    double invCumulativeYstar = 0.0;
    double invCumulativeHstar = 0.0;
    double alpha_m = 0.0; // maternal protection at birth
    double decayM = 0.0; // maternal-protection decay rate per year
    int y_lag_len = 0; // density history for transmission to mosquitoes
};
}
