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

#ifndef Hmod_PathogenesisModel
#define Hmod_PathogenesisModel

#include "Clock.h"
#include "Parameters.h"
#include "Host/WithinHost/State.h"
#include "util/random.h"

namespace scnXml{
    class HSESNMF;
    class Clinical;
}
namespace OM { class RunContext; }

namespace OM {
    namespace Host { class Human; }
    namespace checkpoint { class Archive; }
    namespace WithinHost { namespace Pathogenesis {

using util::LocalRng;

/*! Pathogenesis model: one final class, behaviour chosen by DATA (Kind switch
 * — the old Mueller/Pyrogen/Predet subclass trio, incl. the Submodels.h merge,
 * differed only in getPEpisode + one extra state field + checkpoint bytes).
 *
 * Previously named MorbidityModel and PresentationModel. */
class PathogenesisModel {
public:
    enum class Kind { Mueller, Pyrogen, Predet };

    /// Calls static init on correct PathogenesisModel.
    static void init( RunContext& ctx, const Parameters& parameters, const scnXml::Clinical& clinical, bool nmfOnly );

    /** Create an instance, Kind dependant on global options.
     *
     * @param cF = Comorbidity factor (currently set in Human). */
    static PathogenesisModel* createPathogenesisModel(double cF);

    /** Determines the health of the individual based on his/her parasitemia.
     *
     * May introduce severe or uncomplicated cases of malaria, as well as non-
     * malaria fevers. */
    StatePair determineState( Host::Human& human, double ageYears,
            double timeStepMaxDensity, double endDensity, bool isDoomed );
    
    /** For Vivax: determine the chance of a NMF and sample, returning either
     * NONE or STATE_NMF. */
    static Pathogenesis::State sampleNMF( LocalRng& rng, double ageYears );
    
    /** Summarize PathogenesisModel details
     *
     * Only the pyrogen-threshold kinds (Pyrogen, Predet) have anything to add
     * to the summary (parity: old PyrogenPathogenesis::summarize). */
    void summarize (const Host::Human& human);

    /// Read parameters (parity: old MuellerPathogenesis::init /
    /// PyrogenPathogenesis::init; called from PathogenesisModel::init in xml).
    static void initMueller( RunContext& ctx, const Parameters& parameters );
    static void initPyrogen( RunContext& ctx, const Parameters& parameters );

    /// Save/restore wire state (model/checkpoint/humans.cpp).
    void checkpoint (checkpoint::Archive& a);

protected:
    /** Create a PathogenesisModel (_pyrogenThres = initPyroThres for the
     * pyrogenic kinds, per old PyrogenPathogenesis ctor). */
    PathogenesisModel(double cF, Kind kind);

    /// Episode probability: one case per Kind (were three virtual overrides).
    double getPEpisode(double timeStepMaxDensity, double totalDensity);

    //! comorbidity factor for heterogeneity
    double _comorbidityFactor;

private:
    /// Determine the current pyrogenic threshold (Pyrogen and Predet kinds).
    void updatePyrogenThres(double totalDensity);

    Kind kind;
    /// Critical density for fever (clinical episodes); Pyrogen and Predet only
    double _pyrogenThres;
};

} } }
#endif
