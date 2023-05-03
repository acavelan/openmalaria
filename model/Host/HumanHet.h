/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

// This is a simple header, for direct inclusion in Human.cpp only.

#include "util/random.h"
#include "util/ModelOptions.h"

namespace OM {
namespace Host {

/// (Some) heterogeneity parameters of humans
struct HumanHet {
    static bool opt_trans_het, opt_comorb_het, opt_treat_het,
    opt_trans_treat_het, opt_comorb_treat_het,
    opt_comorb_trans_het, opt_triple_het;

    /* Human heterogeneity; affects:
     * comorbidityFactor (stored in PathogenesisModel)
     * treatmentSeekingFactor (stored in CaseManagementModel)
     * availabilityFactor (stored in Transmission::PerHost) */
    double comorbidityFactor;
    double treatmentSeekingFactor;
    double availabilityFactor;
    HumanHet() : comorbidityFactor(1.0), treatmentSeekingFactor(1.0),
        availabilityFactor(1.0) {}
    
    static void init();

    static HumanHet sample(util::LocalRng& rng);
};

}
}
