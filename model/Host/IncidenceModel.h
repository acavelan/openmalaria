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
#include "Transmission/PerHost.h"
#include "util/random.h"

namespace OM {
class Parameters;
namespace Host {
class Human;

using util::LocalRng;

/** Models how a per-host EIR translates into new infections
 * (roughly when bites from infected mosquitos infect the host).
 *
 * One final class, behaviour chosen by DATA (Kind switch — the four old
 * subclasses differed only in two formulas and added no state; this was the
 * file's own TODO):
 * - Default: Smith et al, AJTMH 2006 75 Suppl 2
 * - HetWorkaround: emulates old, presumably unintended, behaviour — the EIR
 *   was not in one place adjusted by the availability factor used in
 *   transmission heterogeneity, where it possibly should have been
 * - NegBinom: NEGATIVE_BINOMIAL_MASS_ACTION
 * - LogNormal: LOGNORMAL_MASS_ACTION
 *
 * There are also two susceptibility models which should be compatible with all
 * of these (see susceptibility()). */
class IncidenceModel {
  public:
    enum class Kind { Default, HetWorkaround, NegBinom, LogNormal };

  protected:
    /// The run this model belongs to (clock accessed as ctx.clock).
    Context& ctx;

  public:
    ///@brief Static initialisation & constructors
    //@{
    /// Read in/initialise parameters
    static void init(Context& ctx, const Parameters& parameters);
    /// Create a new instance (Kind selected from the model options)
    static IncidenceModel* createModel(Context& ctx);

    /// Reset summary outputs at beginning of main simulation
    static inline void preMainSimInit() { ctsNewInfections = 0; }
    //@}

  protected:
    /// Create a new model
    explicit IncidenceModel(Context& ctx, Kind kind) : ctx(ctx), kind(kind), m_pInfected(0.0), m_cumulativeEIRa(0.0) {}

  public:
    /** Return an availability multiplier, dependant on the model (NegBinom
   * and LogNormal models use this). Ideally availability adjustments
   * should have nothing to do with the IncidenceModel though.
   * 
   * @param baseAvailability This was BaselineAvailabilityMean from Constant.h,
   *	and had the value 1.0. Whether it should be anything else I don't know.
   */
    double getAvailabilityFactor(LocalRng& rng, double baseAvailability = 1.0);

    /// Output _pinfected to the summary
    void summarize(const Host::Human& human);

    /** Calculate the expected number of new infections to introduce.
 * 
 * Firstly converts the EIR into an expected number of infections (what
 * getExpectedNumberOfInfections() used to do):
 * 1. Calculates h from the EIR measured on adults where h is the expected
 * number of epidemiological inoculations.
 * 2. Calculates the updated values of the pre-erythrocytic exposure.
 * 
 * Secondly calculates the number of new infections to introduce via a
 * stochastic process.
 * 
 * @param human   Reference to human
 * @param effectiveEIR    EIR adjusted for availability (age, het,
 *        outside transmission due to hospitalisation)
 */
    double expectedNumNewInfections(OM::Host::Human& human, double effectiveEIR);

    /** Calculate the number of new infections to introduce.
   * 
   * @param human   Reference to human
   * @param expectedNumNewInfections   Expected number of infections
   * to be introduced
   */
    int numNewInfections(OM::Host::Human& human, double expectedNumNewInfections);

    /** Report the number of new infections
   * 
   * This number might be lower than the numNewInfections from this model
   * due to vaccine factors being applied in the WithinHost model
   * 
   * @param int newNumInfections_i  Number of introduced infections to be reported
   * @param int newNumInfections_l  Number of local (indigenous) infections to be reported
   */
    void reportNumNewInfections(Human& human, int newNumInfections_i, int newNumInfections_l);

  public: // run state; wire in checkpoint/humans.cpp + statics.cpp
    /// Calculates the expected number of infections, excluding vaccine effects.
    /// Kind switch: one case per old-subclass formula.
    double getModelExpectedInfections(LocalRng& rng, double effectiveEIR, const Transmission::PerHost& phTrans);

    double susceptibility();

    static void ctsReportNewInfections(const vector<Host::Human>&, ostream& stream);

    /// Which formulas this instance runs (set by createModel from ModelOptions).
    Kind kind;

    /** Probability of infection (cumulative or reset to zero in massTreatment).
   *
   * Appears to be used only for calculating expected inoculations for the
   * analysis of pre-erythrocytic immunity. */
    double m_pInfected;

    //!Number of infective bites since birth
    double m_cumulativeEIRa; //TODO(memory opt): not needed by NegBinom and LogNormal

    /// Number of new infections introduced, per continuous reporting period
    static int ctsNewInfections;
};

}
}
