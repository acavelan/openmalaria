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

#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "Host/WithinHost/DescriptiveWithinHost.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "Host/WithinHost/Treatments.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
//using namespace std;

#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

bool reportInfectionsByGenotype = false;

bool opt_vivax_simple = false,
        opt_dummy_whm = false, opt_empirical_whm = false,
        opt_molineaux_whm = false, opt_penny_whm = false,
        opt_common_whm = false;

// -----  static functions  -----


TreatmentId WHInterface::addTreatment(RunContext& ctx, const scnXml::TreatmentOption& desc){
    return Treatments::addTreatment( ctx, desc );
}

unique_ptr<WHInterface> WHInterface::createWithinHostModel(RunContext& ctx, LocalRng& rng, double comorbidityFactor) {
    if( opt_vivax_simple ) {
        return unique_ptr<WHInterface>(new WHVivax( ctx, rng, comorbidityFactor ));
    } else if( opt_common_whm ) {
        return unique_ptr<WHInterface>(new CommonWithinHost( ctx, rng, comorbidityFactor ));
    } else {
        return unique_ptr<WHInterface>(new DescriptiveWithinHostModel( ctx, rng, comorbidityFactor ));
    }
}


// -----  Non-static  -----


void WHInterface::checkpoint (istream& stream) {
    numInfs & stream;

    if (numInfs > MAX_INFECTIONS)
        throw util::checkpoint_error("numInfs: " + to_string(numInfs));
}
void WHInterface::checkpoint (ostream& stream) {
    numInfs & stream;
}

}
}
