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
#include <limits>

using namespace std;

namespace OM {
namespace PkPd {

LSTMDrug::LSTMDrug (Kind kind, const LSTMDrugType& type,
                    const LSTMDrugType* metabolite, LocalRng& rng) :
    kind(kind), typeData(&type), metaboliteData(metabolite)
{
    switch (kind) {
    case Kind::OneComp:
        vol_dist = type.sample_Vd(rng);
        concentration = 0.0;
        neg_elim_sample = -type.sample_elim_rate(rng);
        break;
    case Kind::ThreeComp:
        vol_dist = type.sample_Vd(rng);
        concA = 0.0; concB = 0.0; concC = 0.0; concABC = 0.0;
        elim_sample = type.sample_elim_rate(rng);
        k12 = type.sample_k12(rng);
        k21 = type.sample_k21(rng);
        k13 = type.sample_k13(rng);
        k31 = type.sample_k31(rng);
        nka = -type.sample_ka(rng);
        last_bm = numeric_limits<double>::quiet_NaN();
        na = numeric_limits<double>::quiet_NaN();
        nb = numeric_limits<double>::quiet_NaN();
        ng = numeric_limits<double>::quiet_NaN();
        A = numeric_limits<double>::quiet_NaN();
        B = numeric_limits<double>::quiet_NaN();
        C = numeric_limits<double>::quiet_NaN();
        break;
    case Kind::Conversion:
        vol_dist = type.sample_Vd(rng);
        qtyG = 0.0; qtyP = 0.0; qtyM = 0.0;
        nka = -type.sample_ka(rng);
        nkP_sample = -type.sample_elim_rate(rng);
        nconv_sample = -type.sample_conv_rate(rng);
        nkM_sample = -metaboliteData->sample_elim_rate(rng);
        vol_dist_metabolite = metaboliteData->sample_Vd(rng);
        break;
    }
}

size_t LSTMDrug::getIndex() const {
    return typeData->getIndex();   // Conversion: parent drug's index
}

double LSTMDrug::getConcentration(size_t index) const {
    switch (kind) {
    case Kind::OneComp:
        return index == typeData->getIndex() ? concentration : 0.0;
    case Kind::ThreeComp:
        return index == typeData->getIndex() ? conc() : 0.0;
    case Kind::Conversion:
        if( index == typeData->getIndex() ) return getParentConcentration();
        else if( index == metaboliteData->getIndex() ) return getMetaboliteConcentration();
        else return 0.0;
    }
    return 0.0;
}

double LSTMDrug::calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) const {
    switch (kind) {
    case Kind::OneComp: return drugFactor_oneComp (rng, inf, body_mass);
    case Kind::ThreeComp: return drugFactor_threeComp (rng, inf, body_mass);
    default: return drugFactor_conversion (rng, inf, body_mass);
    }
}

void LSTMDrug::updateConcentration (double body_mass) {
    switch (kind) {
    case Kind::OneComp: updateConcentration_oneComp (body_mass); break;
    case Kind::ThreeComp: updateConcentration_threeComp (body_mass); break;
    case Kind::Conversion: updateConcentration_conversion (body_mass); break;
    }
}

void LSTMDrug::checkpoint (istream& stream) {
    switch (kind) {
    case Kind::OneComp: checkpoint_oneComp (stream); break;
    case Kind::ThreeComp: checkpoint_threeComp (stream); break;
    case Kind::Conversion: checkpoint_conversion (stream); break;
    }
}
void LSTMDrug::checkpoint (ostream& stream) {
    switch (kind) {
    case Kind::OneComp: checkpoint_oneComp (stream); break;
    case Kind::ThreeComp: checkpoint_threeComp (stream); break;
    case Kind::Conversion: checkpoint_conversion (stream); break;
    }
}

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

}
}
