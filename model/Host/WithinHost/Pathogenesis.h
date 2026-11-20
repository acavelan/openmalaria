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

#include "Clock.h"
#include "Parameters.h"
#include "util/random.h"

namespace scnXml {
class HSESNMF;
class Clinical;
}
namespace OM {
class Context;
}

namespace OM {
namespace Host {
class Human;
}
namespace WithinHost {
namespace Pathogenesis {

using util::LocalRng;

/** Types of sickness; used by case management.
   *
   * Most values are flags which can be combined in any form. A few
   * combinations set follow. */
enum State {
    /* Values here are written in hexadecimal: http://en.wikipedia.org/wiki/Hexadecimal
     * Many are designed to be "flags", so the value corresponds to a single bit:
     * http://en.wikipedia.org/wiki/Flag_byte
     * Max: 0x80000000
     * (note & | ^ are C++'s binary AND, OR and XOR operators). */
    NONE = 0, ///< Not sick

    // Flags for current state/worst state to report:
    SICK = 0x1,         ///< Sick (may or may not be from malaria)
    MALARIA = 0x2,      ///< Malaria sickness
    SEVERE = 0x8,       ///< Severe malaria case
    COINFECTION = 0x4,  ///< Malaria with a coinfection
    COMPLICATED = 0x20, ///< Flag used to indicate SEVERE and/or COINFECTION

    // Flags for outcome reporting:
    STATE_NMF = SICK, ///< Alias: SICK
    //STATE_NMF_NEED_AB    = SICK | NEED_ANTIBIOTIC,       ///< Combination: SICK, NEED_ANTIBIOTIC
    STATE_MALARIA = SICK | MALARIA,                      ///< Combination: SICK, MALARIA
    STATE_SEVERE = STATE_MALARIA | COMPLICATED | SEVERE, ///< Combination: SICK, MALARIA, COMPLICATED, SEVERE
    STATE_COINFECTION = STATE_MALARIA | COMPLICATED |
        COINFECTION, ///< Combination: SICK, MALARIA, COMPLICATED, COINFECTION
};
struct StatePair {
    StatePair() : state(NONE), indirectMortality(false) {}
    State state;
    bool indirectMortality; // implies death a little later on which may not be seen as a malaria death
};

/*! Pathogenesis model: one final class, behaviour chosen by DATA (Kind switch
 * — the old Mueller/Pyrogen/Predet subclass trio, incl. the Submodels.h merge,
 * differed only in getPEpisode + one extra state field + checkpoint bytes).
 *
 * Previously named MorbidityModel and PresentationModel. */
class Model {
  public:
    enum class Kind { Mueller, Pyrogen, Predet };

    /// Calls static init on correct Model.
    static void init(Context& ctx, const Parameters& parameters, const scnXml::Clinical& clinical, bool nmfOnly);

    /** Create an instance, Kind dependant on global options.
     *
     * @param cF = Comorbidity factor (currently set in Human). */
    static Model* createPathogenesisModel(double cF);

    /** Determines the health of the individual based on his/her parasitemia.
     *
     * May introduce severe or uncomplicated cases of malaria, as well as non-
     * malaria fevers. */
    StatePair determineState(Host::Human& human, double ageYears, double timeStepMaxDensity, double endDensity,
                             bool isDoomed);

    /** For Vivax: determine the chance of a NMF and sample, returning either
     * NONE or STATE_NMF. */
    static Pathogenesis::State sampleNMF(LocalRng& rng, double ageYears);

    /** Summarize Model details
     *
     * Only the pyrogen-threshold kinds (Pyrogen, Predet) have anything to add
     * to the summary (parity: old PyrogenPathogenesis::summarize). */
    void summarize(const Host::Human& human);

    /// Read parameters (parity: old MuellerPathogenesis::init /
    /// PyrogenPathogenesis::init; called from Model::init in xml).
    static void initMueller(Context& ctx, const Parameters& parameters);
    static void initPyrogen(Context& ctx, const Parameters& parameters);

  protected:
    /** Create a Model (_pyrogenThres = initPyroThres for the
     * pyrogenic kinds, per old PyrogenPathogenesis ctor). */
    Model(double cF, Kind kind);

    /// Episode probability: one case per Kind (were three virtual overrides).
    double getPEpisode(double timeStepMaxDensity, double totalDensity);

  public:
    //! comorbidity factor for heterogeneity
    double _comorbidityFactor;

  private:
    /// Determine the current pyrogenic threshold (Pyrogen and Predet kinds).
    void updatePyrogenThres(double totalDensity);

    Kind kind;

  public:
    /// Critical density for fever (clinical episodes); Pyrogen and Predet only
    double _pyrogenThres;
};

}
}
}
