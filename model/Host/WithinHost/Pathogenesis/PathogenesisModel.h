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
#include "Global.h"
#include "Host/WithinHost/Pathogenesis/State.h"
#include "util/AgeGroupInterpolation.h"
#include "util/random.h"

namespace OM::Host { struct Human; }
namespace OM::WithinHost::Pathogenesis {
struct PathogenesisParameters {
    enum class Kind { Pyrogenic, Mueller, Predetermined };
    Kind kind = Kind::Pyrogenic;
    double pg_indirRiskCoFactor = 0.0; // indirect-mortality comorbidity at birth
    double pg_severeMalThreshold = 0.0; // critical density for severe malaria
    double pg_comorbIntercept = 0.0; // severe-malaria comorbidity at birth
    double pg_inv_critAgeComorb = 0.0; // inverse critical age for comorbidity
    util::AgeGroupInterpolator pg_NMF_incidence;
    double rateMultiplier_31 = 0.0;
    double densityExponent_32 = 0.0;
    double initPyroThres = 0.0;
    double Ystar2_13 = 0.0;
    double Ystar1_26 = 0.0;
    double a = numeric_limits<double>::signaling_NaN();
    double b = numeric_limits<double>::signaling_NaN();
};

State sampleNMF(const PathogenesisParameters&, util::LocalRng&, double ageYears);

struct PathogenesisModel {
    PathogenesisModel(double comorbidityFactor, const PathogenesisParameters& config);
    const PathogenesisParameters& config;
    double _comorbidityFactor;
    double _pyrogenThres;

    StatePair determineState(Host::Human&, double ageYears,
        double timeStepMaxDensity, double endDensity, bool isDoomed);
    double getPEpisode(double timeStepMaxDensity, double totalDensity);
    void updatePyrogenThres(double totalDensity);
    void summarize(const Host::Human&) const;

    template<class S> void operator&(S& stream) {
        _comorbidityFactor & stream;
        // Preserve the historical model-specific checkpoint layout.
        if (config.kind != PathogenesisParameters::Kind::Mueller) _pyrogenThres & stream;
    }
};
}
