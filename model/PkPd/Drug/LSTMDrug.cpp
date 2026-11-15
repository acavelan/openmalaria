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

#include "PkPd/Drug/LSTMDrug.h"
#include "PkPd/Drug/LSTMDrugType.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <algorithm>

using namespace std;

namespace OM {
namespace PkPd {

LSTMDrug::LSTMDrug(Kind kind, const LSTMDrugType& typeData, double vol_dist)
    : kind(kind), typeData(typeData), vol_dist(vol_dist) {}

size_t LSTMDrug::getIndex() const { return typeData.getIndex(); }

bool comp(pair<double,double> lhs, pair<double,double> rhs){
    return lhs.first < rhs.first;
}

// Create our list of doses. Optimise for the case where we only have 1 or less
// per day, but be able to handle more.
// If two doses coincide, they can be combined but doing so is not essential.

void LSTMDrug::medicate(double time, double qty){
    // Insert in the right position to maintain sorting:
    auto elt = make_pair (time, qty);
    auto pos = lower_bound(doses.begin(), doses.end(), elt, comp);
    doses.insert(pos, std::move(elt));
    assert(is_sorted(doses.begin(), doses.end(), comp));
}


double LSTMDrug::getConcentration(size_t index) const {
    switch (kind) {
    case ONE_COMP:   return lstmdrug::oneComp::getConcentration(*this, index);
    case THREE_COMP: return lstmdrug::threeComp::getConcentration(*this, index);
    default:         return lstmdrug::conversion::getConcentration(*this, index);
    }
}

double LSTMDrug::calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) const {
    switch (kind) {
    case ONE_COMP:   return lstmdrug::oneComp::drugFactor(*this, rng, inf, body_mass);
    case THREE_COMP: return lstmdrug::threeComp::drugFactor(*this, rng, inf, body_mass);
    default:         return lstmdrug::conversion::drugFactor(*this, rng, inf, body_mass);
    }
}

void LSTMDrug::updateConcentration(double body_mass) {
    switch (kind) {
    case ONE_COMP:   lstmdrug::oneComp::update(*this, body_mass); break;
    case THREE_COMP: lstmdrug::threeComp::update(*this, body_mass); break;
    default:         lstmdrug::conversion::update(*this, body_mass); break;
    }
}

void LSTMDrug::checkpoint(ostream& stream) {
    switch (kind) {
    case ONE_COMP:   lstmdrug::oneComp::checkpoint(*this, stream); break;
    case THREE_COMP: lstmdrug::threeComp::checkpoint(*this, stream); break;
    default:         lstmdrug::conversion::checkpoint(*this, stream); break;
    }
}
void LSTMDrug::checkpoint(istream& stream) {
    switch (kind) {
    case ONE_COMP:   lstmdrug::oneComp::read(*this, stream); break;
    case THREE_COMP: lstmdrug::threeComp::read(*this, stream); break;
    default:         lstmdrug::conversion::read(*this, stream); break;
    }
}


}
}
