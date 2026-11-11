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

#ifndef Hmod_InfectionIncidence
#define Hmod_InfectionIncidence

#include "Global.h"
#include "util/random.h"

namespace OM::Host {

enum class IncidenceKind { Default, Heterogeneity, NegativeBinomial, Lognormal };

struct IncidenceParameters {
    IncidenceKind kind = IncidenceKind::Default;
    double baseline_avail_shape_param = 0.0, baseline_avail_offset = 0.0;
    double inf_rate_shape_param = 0.0, inf_rate_offset = 0.0;
    double gamma_p = 0.0, Sinf = 0.0, Simm = 0.0;
    double Xstar_pInv = 0.0, EstarInv = 0.0;
    bool noPreErythrocytic = false;
    bool vaccineGenotype = false;
};

struct IncidenceState {
    double probabilityInfected = 0.0;
    double cumulativeEIR = 0.0;

    template<class S> void operator&(S& stream) {
        probabilityInfected & stream;
        cumulativeEIR & stream;
    }
};

double sampleAvailability(const IncidenceParameters&, util::LocalRng&);
double expectedInfections(IncidenceState&, const IncidenceParameters&,
    util::LocalRng&, double effectiveEIR, double relativeAvailability,
    double vaccineFactor, SimTime interval);
int sampleNewInfections(util::LocalRng&, double expectedNumInfections);

}
#endif
