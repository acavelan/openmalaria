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

// p8-files: LSTMDrugType.h rides the Context.h cone; the type-registry XML
// init + parameter sampling stay here (wide-cone rule).
#include "PkPd/LSTMDrugType.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "PkPd/LSTMDrug.h"

#include <cmath>
#include <iomanip>

using namespace std;

namespace OM {
namespace PkPd {
using namespace OM::util;
using WithinHost::Genotypes;

// -----  Static variables and functions  -----

double LSTMDrugPD::calcFactor(double Kn, double neg_elim_rate, double* conc, double duration) const {
    const double C0 = *conc;
    const double C1 = C0 * exp(neg_elim_rate * duration);

    // From Hastings & Winter 2011 paper
    // Note: these look a little different from original equations because Kn
    // is calculated when parameters are read from the scenario document instead of now.
    const double numerator = Kn + pow(C1, n);
    const double denominator = Kn + pow(C0, n);
    const double power = V / (-neg_elim_rate * n);

    *conc = C1;                                 // conc is an in/out parameter
    return pow(numerator / denominator, power); // unitless
}

double LSTMDrugPD::IC50_pow_slope(LocalRng& rng, size_t index, WithinHost::Infection* inf) const {
    double Kn; // gets sampled once per infection
    if (inf->Kn().count(index) != 0) {
        Kn = inf->Kn().at(index);
    } else {
        Kn = pow(IC50->sample(rng), n);
        inf->Kn()[index] = Kn;
    }
    return Kn;
}

void DrugRegistry::clear() {
    drugTypes.clear();
    drugTypeNames.clear();
}

// Add index to drugsInUse if not already present. Not fast but doesn't need to be.
static void drugIsUsed(vector<size_t>& drugsInUse, size_t index) {
    for (size_t i : drugsInUse) {
        if (i == index)
            return; // already in list
    }
    drugsInUse.push_back(index);
}
size_t DrugRegistry::findDrug(string _abbreviation) {
    auto iter = drugTypeNames.find(_abbreviation);
    if (iter == drugTypeNames.end())
        throw util::xml_error(string("attempt to use drug without description: ").append(_abbreviation));
    size_t index = iter->second;

    // We assume that drugs are used when and only when findDrug returns their
    // index or they are a metabolite of a drug returned here.
    drugIsUsed(drugsInUse, index);
    if (drugTypes[index].conversion_rate /*conversion model used*/) {
        drugIsUsed(drugsInUse, drugTypes[index].metabolite);
    }

    return index;
}

LSTMDrugType& DrugRegistry::get(size_t index) { return drugTypes.at(index); }

const vector<size_t>& DrugRegistry::getDrugsInUse() const { return drugsInUse; }

unique_ptr<LSTMDrug> DrugRegistry::createInstance(LocalRng& rng, size_t index) {
    // One create-line per variant; the drug type's XML structure selects it
    // (same selection the old Kind switch made at construction).
    LSTMDrugType& typeData = drugTypes[index];
    if (typeData.conversion_rate) {
        return unique_ptr<LSTMDrug>(new LSTMDrugConversion(typeData, &drugTypes[typeData.metabolite], rng));
    } else if (typeData.k12) {
        // k21 is set when k12 is set; k13 and k31 may be set
        return unique_ptr<LSTMDrug>(new LSTMDrugThreeComp(typeData, rng));
    } else {
        // none of k12/k21/k13/k31 should be set in this case
        return unique_ptr<LSTMDrug>(new LSTMDrugOneComp(typeData, rng));
    }
}

// -----  Non-static LSTMDrugType functions  -----

LSTMDrugType::~LSTMDrugType() {}

const LSTMDrugPD& LSTMDrugType::getPD(uint32_t genotype) const { return PD[genotype_mapping[genotype]]; }

}
}
