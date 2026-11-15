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

#ifndef Hmod_LSTMDrug
#define Hmod_LSTMDrug

#include "Global.h"
#include "util/checkpoint_containers.h"
#include "util/random.h"

#include <limits>
#include <vector>

namespace OM {
namespace WithinHost {
    class CommonInfection;
}
namespace PkPd {

using util::LocalRng;

class LSTMDrugType;

/** PK/PD state of one drug type in one human: one Kind-tagged data holder
 * for all three PK models (was LSTMDrugOneComp/ThreeComp/Conversion). The
 * per-model logic is the lstmdrug free functions below, defined in the same
 * per-model .cpp files and dispatched by the methods here.
 *
 * Checkpoint streams keep the old per-kind field order (LSTMModel writes the
 * drug type index first, which selects the kind when loading). */
struct LSTMDrug {
    enum Kind { ONE_COMP, THREE_COMP, CONVERSION };

    /** Volume of distribution must be specified here (from sample or mean).
     * Remaining fields: zero, except the cached quantities (NaN). */
    LSTMDrug(Kind kind, const LSTMDrugType& typeData, double vol_dist);

    /// The drug type's index (the CONVERSION model reports the parent's).
    size_t getIndex() const;

    /** Indicate a new medication this time step: stores the (time, qty) pair
     * in doses, keeping the list ordered by time.
     *
     * @param time Time of administration, as the delay since the start of the
     *  current time step in units of days.
     * @param qty Amount of active ingredient, in mg (total) */
    void medicate (double time, double qty);

    /** Get the concentration of the given drug type contained in this model
     * (only compartments with active PD), in mg/l; 0 if the drug index
     * doesn't match that used. */
    double getConcentration(size_t index) const;

    /** Returns the total drug factor for this drug over one day (multiply the
     * factors of all drugs). Doesn't change any state: may be called several
     * times per time step (once per infection) or not at all.
     *
     * @param body_mass Weight of patient in kg */
    double calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) const;

    /** Update the state with today's doses and one day of decay, and clear
     * the day's doses.
     *
     * @param body_mass Weight of patient in kg */
    void updateConcentration (double body_mass);

    /// Checkpointing
    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);
    template<class S>
    void operator& (S& stream) {
        doses & stream;
        checkpoint(stream);
    }

    Kind kind;

    /// Always links a drug instance to its drug-type data. For CONVERSION
    /// this is the parent drug; see metaboliteType.
    const LSTMDrugType& typeData;

    /// First is time (days), second is additional concentration (mg / l; for
    /// one- and three-compartment models) or quantity (mg; for conversion
    /// model). List of each dose given today (and possibly tomorrow),
    /// ordered by time. Read in calculateDrugFactor, updated in
    /// updateConcentration().
    typedef std::vector<std::pair<double,double> > DoseVec;
    DoseVec doses;

    /// Volume of distribution (CONVERSION: parent drug), sampled on creation.
    double vol_dist;

    /// ONE_COMP: concentration in blood; units: mg / l.
    double concentration = 0.0;
    /// ONE_COMP: sampled elimination rate constant (negative).
    double neg_elim_sample = 0.0;

    /// THREE_COMP: concentrations in the blood and other compartments; see
    /// "Permutation of the three-compartment equation", Diggory Hardy, 2015.
    double concA = 0.0, concB = 0.0, concC = 0.0, concABC = 0.0;
    /// THREE_COMP: sampled rate constants (elim_sample positive).
    double elim_sample = 0.0, k12 = 0.0, k21 = 0.0, k13 = 0.0, k31 = 0.0;
    /// THREE_COMP and CONVERSION: minus the absorption rate constant.
    double nka = 0.0;

    /// CONVERSION: the metabolite's drug type.
    const LSTMDrugType* metaboliteType = nullptr;
    /// CONVERSION: amounts of parent drug in gut, parent drug in circulation
    /// and metabolite in circulation; units: mg (A, B, C in the paper).
    double qtyG = 0.0, qtyP = 0.0, qtyM = 0.0;
    /// CONVERSION: sampled elimination (parent, metabolite) and conversion
    /// rate constants (negative), and the metabolite's volume of
    /// distribution (the parent's is vol_dist).
    double nkP_sample = 0.0, nconv_sample = 0.0, nkM_sample = 0.0;
    double vol_dist_metabolite = 0.0;

    /// THREE_COMP and CONVERSION: body mass at the last update (NaN before).
    /// THREE_COMP also caches the derived decay constants na, nb, ng and the
    /// dose-response coefficients A, B, C, recomputed whenever the body mass
    /// changes. Mutable: calculateDrugFactor() is const.
    mutable double last_bm = std::numeric_limits<double>::quiet_NaN();
    mutable double na = std::numeric_limits<double>::quiet_NaN();
    mutable double nb = std::numeric_limits<double>::quiet_NaN();
    mutable double ng = std::numeric_limits<double>::quiet_NaN();
    mutable double A = std::numeric_limits<double>::quiet_NaN();
    mutable double B = std::numeric_limits<double>::quiet_NaN();
    mutable double C = std::numeric_limits<double>::quiet_NaN();

    /// THREE_COMP: total blood concentration.
    double conc() const { return concA + concB + concC - concABC; }
};

/** Per-PK-model implementations of the LSTMDrug methods: the old
 * LSTMDrugOneComp/ThreeComp/Conversion virtual-method bodies with member
 * accesses turned into d.field. Checkpoint/read keep the per-kind order. */
//@{
namespace lstmdrug {
    namespace oneComp {
        double getConcentration(const LSTMDrug&, size_t index);
        double drugFactor(const LSTMDrug&, LocalRng&, WithinHost::CommonInfection*, double body_mass);
        void update(LSTMDrug&, double body_mass);
        void checkpoint(LSTMDrug&, ostream&);
        void read(LSTMDrug&, istream&);
    }
    namespace threeComp {
        double getConcentration(const LSTMDrug&, size_t index);
        double drugFactor(const LSTMDrug&, LocalRng&, WithinHost::CommonInfection*, double body_mass);
        void update(LSTMDrug&, double body_mass);
        void checkpoint(LSTMDrug&, ostream&);
        void read(LSTMDrug&, istream&);
    }
    namespace conversion {
        double getConcentration(const LSTMDrug&, size_t index);
        double drugFactor(const LSTMDrug&, LocalRng&, WithinHost::CommonInfection*, double body_mass);
        void update(LSTMDrug&, double body_mass);
        void checkpoint(LSTMDrug&, ostream&);
        void read(LSTMDrug&, istream&);
    }
}
//@}

}
}
#endif
