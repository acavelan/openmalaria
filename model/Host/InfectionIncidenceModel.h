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

#ifndef Hmod_InfectionIncidence
#define Hmod_InfectionIncidence

#include "Global.h"
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
 * There are four versions of this model, with different availability models:
 * - Plain (the default): Smith et al, AJTMH 2006 75 Suppl 2
 * - HeterogeneityWorkaround: emulates old, presumably unintended, behaviour
 * - NegBinomMAII: NEGATIVE_BINOMIAL_MASS_ACTION
 * - LogNormalMAII: LOGNORMAL_MASS_ACTION
 * 
 * There are also two susceptibility models which should be compatible with all
 * of these (see susceptibility()). */
struct InfectionIncidenceParams {
    enum Kind { Plain, HeterogeneityWorkaround, NegBinomMAII, LogNormalMAII };
    
    /** Read in/initialise parameters (kind selected from model options). */
    void init( const OM::Parameters& parameters );
    
    Kind kind = Plain;
    
    /// Shape constant of (Gamma) distribution of availability
    double baselineAvailShapeParam = 0.0, baselineAvailOffset = 0.0;
    
    //! Describes the shape of the Infectionrate distribution, related to the
    //! baseline availabilty distr.
    double infRateShapeParam = 0.0, infRateOffset = 0.0;
    
    ///@brief Parameters used in the "expected number of infections" model
    //@{
    /// Steepness of relationship between success of inoculation and Xp in Phase A model
    double gammaP = 0.0;
    /// Lower limit of success probability of inoculations at high exposure in Phase A model
    double Sinf = 0.0;
    /// Lower limit of success probability of inoculations in immune individuals in Phase A model
    double Simm = 0.0;
    /// 1 over the critical value of cumulative number of entomologic inoculations in Phase A model
    double XstarPInv = 0.0;
    /// 1 over the critical value of EIR in Phase A pre-erythrocytic model
    double EstarInv = 0.0;
    //@}
    
    /// Model options (NO_PRE_ERYTHROCYTIC, VACCINE_GENOTYPE):
    bool noPreErythrocytic = false, vaccineGenotype = false;
};

/** Per-host infection incidence state; the free functions below do the work,
 * dispatched by params.kind where the models differ. */
struct InfectionIncidenceModel {
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        m_pInfected & stream;
        m_cumulativeEIRa & stream;
    }
    
    /// Output _pinfected to the summary
    void summarize (const Host::Human& human) const;
    
    /** Probability of infection (cumulative or reset to zero in massTreatment).
     *
     * Appears to be used only for calculating expected inoculations for the
     * analysis of pre-erythrocytic immunity. */
    double m_pInfected = 0.0;
    
    //!Number of infective bites since birth
    double m_cumulativeEIRa = 0.0;//TODO(memory opt): not needed by NegBinomMAII and LogNormalMAII
};

/** Return an availability multiplier, dependant on the model kind (NegBinomMAII
 * and LogNormalMAII differ from the rest).
 * 
 * @param baseAvailability This was BaselineAvailabilityMean from Constant.h,
 *	and had the value 1.0. Whether it should be anything else I don't know.
 */
double getAvailabilityFactor (LocalRng& rng, double baseAvailability,
                              const InfectionIncidenceParams& params);

/** Calculates the expected number of infections, excluding vaccine effects */
double getModelExpectedInfections (const InfectionIncidenceModel& model,
                                   LocalRng& rng, double effectiveEIR,
                                   const Transmission::PerHost& phTrans,
                                   const InfectionIncidenceParams& params);

double susceptibility (const InfectionIncidenceModel& model,
                       const InfectionIncidenceParams& params);

/** Calculate the expected number of new infections to introduce: converts the
 * EIR into an expected number of infections, introduces the effect of
 * vaccination and updates the pre-erythrocytic immunity state.
 * 
 * @param human   Reference to human
 * @param effectiveEIR    EIR adjusted for availability (age, het,
 *        outside transmission due to hospitalisation) */
double expectedNumNewInfections (InfectionIncidenceModel& model, Human& human,
                                 double effectiveEIR,
                                 const InfectionIncidenceParams& params);

/** Calculate the number of new infections to introduce.
 * 
 * @param expectedNumNewInfections   Expected number of infections
 * to be introduced */
int numNewInfections (Human& human, double expectedNumNewInfections);

/** Report the number of new infections
 * 
 * This number might be lower than the numNewInfections from this model
 * due to vaccine factors being applied in the WithinHost model
 * 
 * @param ctsNewInfections  per-continuous-reporting-period counter (Model-owned)
 * @param int newNumInfections_i  Number of introduced infections to be reported
 * @param int newNumInfections_l  Number of local (indigenous) infections to be reported
 */
void reportNumNewInfections (Human& human, int newNumInfections_i,
                             int newNumInfections_l, int& ctsNewInfections);

} }
#endif
