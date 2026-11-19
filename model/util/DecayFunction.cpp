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

#include <algorithm>
#include <limits>
#include <memory>

#include "Global.h"
#include "util/sampler.h"
#include "util/UnitParse.h"
#include "util/DecayFunction.h"

namespace OM {
namespace util {
    
// -----  interface / static functions  -----
void DecayFunction::checkpointSample(
        const unique_ptr<DecayFunction>& sample, ostream& stream
) const {
    bool present = bool(sample);
    present & stream;
    if (present) sample->hetFactor & stream;
}

void DecayFunction::checkpointSample(
        unique_ptr<DecayFunction>& sample, istream& stream
) const {
    bool present;
    present & stream;
    if (present) {
        double factor;
        factor & stream;
        sample = hetSample(factor);
    } else {
        sample.reset();
    }
}


} }
