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

// p8-files: InfectionIncidenceModel.h sits in the Human/WH chain; the per-Kind
// incidence formula bodies stay here (wide-cone rule).
#include "RunContext.h"
#include "Host/InfectionIncidenceModel.h"
#include "Host/Human.h"
#include "Transmission/PerHost.h"
#include "Host/WithinHost/WHInterface.h"
#include "Parameters.h"
#include "mon/Continuous.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/errors.h"

#include <stdexcept>
#include <cmath>

namespace OpenMalaria {

// ———  model constants set by init()  ———

/* Shape constant of (Gamma) distribution of availability
real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
double baseline_avail_shape_param, baseline_avail_offset;

//VARIABLES INCLUDED IN CORE GETs of number of infections
//! Describes the shape of the Infectionrate distribution, related to the baseline availabilty distr.
double inf_rate_shape_param, inf_rate_offset;

/** @brief Parameters used in the "expected number of infections" model */
//@{
//!Steepness of relationship between success of inoculation and Xp in Phase A model
double gamma_p;
//!Lower limit of success probability of inoculations at high exposure in Phase A model
double Sinf;
//!Lower limit of success probability of inoculations in immune individuals in Phase A model
double Simm;
//!1 over the critical value of cumulative number of entomologic inoculations in Phase A model
double Xstar_pInv;
//!1 over the critical value of EIR in Phase A pre-erythrocytic model
double EstarInv;
//@}

// model options:
bool opt_neg_bin_mass_action = false, opt_lognormal_mass_action = false, opt_no_pre_erythrocytic = false,
     opt_any_het = false, opt_vaccine_genotype = false;
enum class Kind { Default, HetWorkaround, NegBinom, LogNormal };
Kind kind = Kind::Default;

// ———  variables  ———
int InfectionIncidenceModel::ctsNewInfections = 0;

// -----  static initialisation  -----

void InfectionIncidenceModel::init(RunContext& ctx, const Parameters& parameters) {
    baseline_avail_shape_param = parameters[Parameter::BASELINE_AVAILABILITY_SHAPE];
    baseline_avail_offset = 0.5 * pow(baseline_avail_shape_param, 2);

    gamma_p = parameters[Parameter::GAMMA_P];
    Sinf = 1 - exp(-parameters[Parameter::NEG_LOG_ONE_MINUS_SINF]);
    Simm = parameters[Parameter::SIMM];
    EstarInv = 1.0 / parameters[Parameter::E_STAR];
    Xstar_pInv = 1.0 / parameters[Parameter::X_STAR_P];

    //! constant defining the constraint for the Gamma shape parameters
    /// Used for the case where availability is assumed gamma distributed
    double r_square_Gamma;
    /*
    //! Expected number of inoculations
    /// Product of measured EIR, susceptibility and length of time Global::interval
    double gsi = 1.0;
    
    r_square_Gamma=(totalInfectionrateVariance**2-gsi*BaselineAvailabilityMean)/(gsi*BaselineAvailabilityMean)**2
    r_square_Gamma must be greater than zero, so r_square_LogNormal is also. 
    */
    r_square_Gamma = 0.649; //such that r_square_LogNormal =0.5

    opt_no_pre_erythrocytic = ModelOptions::option(NO_PRE_ERYTHROCYTIC);
    opt_neg_bin_mass_action = ModelOptions::option(NEGATIVE_BINOMIAL_MASS_ACTION);
    opt_vaccine_genotype = ModelOptions::option(VACCINE_GENOTYPE);

    if (opt_neg_bin_mass_action) {
        inf_rate_shape_param = (baseline_avail_shape_param + 1.0) / (r_square_Gamma * baseline_avail_shape_param - 1.0);
        inf_rate_shape_param = std::max(inf_rate_shape_param, 0.0);
    } else if (ModelOptions::option(LOGNORMAL_MASS_ACTION)) {
        opt_lognormal_mass_action = true;

        //! constant defining the constraint for the log Normal variance
        /// Used for the case where availability is assumed log Normally distributed
        double r_square_LogNormal = log(1.0 + r_square_Gamma);

        inf_rate_shape_param = sqrt(r_square_LogNormal - 1.86 * pow(baseline_avail_shape_param, 2));
        inf_rate_shape_param = std::max(inf_rate_shape_param, 0.0);
        inf_rate_offset = 0.5 * pow(inf_rate_shape_param, 2);
        if ((std::isnan)(inf_rate_shape_param)) {
            throw xml_scenario_error("bad parameter 16 (BASELINE_AVAILABILITY_SHAPE)");
        }
    } else {
        if (ModelOptions::option(TRANS_HET) || ModelOptions::option(COMORB_TRANS_HET) ||
            ModelOptions::option(TRANS_TREAT_HET) || ModelOptions::option(TRIPLE_HET)) {
            opt_any_het = true;
            cerr << "Warning: will use heterogeneity workaround." << endl;
        }
    }

    kind = opt_neg_bin_mass_action   ? Kind::NegBinom
           : opt_lognormal_mass_action ? Kind::LogNormal
           : opt_any_het               ? Kind::HetWorkaround
                                      : Kind::Default;

    Continuous::registerCallback(ctx, "new infections", "\tnew infections",
                                 &InfectionIncidenceModel::ctsReportNewInfections);
}

// -----  other static methods  -----

void InfectionIncidenceModel::ctsReportNewInfections(const vector<Human>&, ostream& stream) {
    stream << '\t' << ctsNewInfections;
    ctsNewInfections = 0;
}

// -----  non-static methods  -----

double InfectionIncidenceModel::getAvailabilityFactor(LocalRng& rng, double baseAvailability) {
    switch (kind) {
    case Kind::NegBinom:
        // Gamma sample with k=BaselineAvailabilityShapeParam; mean is baseAvailability
        // parity: order — local RNG stream: draw call order fixes stream state
        return rng.gamma(baseline_avail_shape_param, baseAvailability / baseline_avail_shape_param);
    case Kind::LogNormal:
        // given BaselineAvailabilityShapeParam = sqrt (log (1 + variance/mean²))
        // and baseAvailability = mean, this is a draw from the log-normal distribution.
        if (baseAvailability != 1.0) {
            // NOTE: shouldn't the normal_mean parameter be adjusted when baseAvailability != 1.0?
            throw TRACED_EXCEPTION_DEFAULT("LogNormalMAII::getAvailabilityFactor");
        }
        return rng.log_normal(log(baseAvailability) - baseline_avail_offset, baseline_avail_shape_param);
    default: // Default and HetWorkaround have no availability adjustment
        return baseAvailability;
    }
}

void InfectionIncidenceModel::summarize(const Human& human) {
    recordStat(measure("nExpectd"), human, human.p_infected);
}

double InfectionIncidenceModel::getModelExpectedInfections(Human& human, double effectiveEIR) {
    LocalRng& rng = human.rng();
    RunContext& ctx = *human.ctx;
    const PerHost& phTrans = human.perHostTransmission;
    switch (kind) {
    case Kind::HetWorkaround: // parity: old HeterogeneityWorkaroundII formula
        return (Sinf + (1 - Sinf) / (1 + effectiveEIR / (ctx.clock.one_TS() * phTrans.rel_avail_het()) * EstarInv)) *
               susceptibility(human) * effectiveEIR;
    case Kind::NegBinom: // parity: old NegBinomMAII formula
        // Documentation: http://www.plosmedicine.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pmed.1001157.s009
        // parity: order — local RNG stream: draw call order fixes stream state
        return rng.gamma(inf_rate_shape_param, effectiveEIR * susceptibility(human) / inf_rate_shape_param);
    case Kind::LogNormal: { // parity: old LogNormalMAII formula
        // Documentation: http://www.plosmedicine.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pmed.1001157.s009
        if (effectiveEIR < 1e-10)
            return 0.0;
        double meanlog = log(effectiveEIR * susceptibility(human)) - inf_rate_offset;
        return rng.log_normal(meanlog, inf_rate_shape_param);
    }
    default:
        // First two lines are availability adjustment: S_1(i,t) from AJTMH 75 (suppl 2) p12 eqn. (5)
        // Note that NegBinom and LogNormal supercede this model; see above
        return (Sinf + (1 - Sinf) / (1 + effectiveEIR / ctx.clock.one_TS() * EstarInv)) * susceptibility(human) *
               effectiveEIR;
    }
}

double InfectionIncidenceModel::susceptibility(const Human& human) {
    if (opt_no_pre_erythrocytic) {
        //! The average proportion of bites from sporozoite positive mosquitoes resulting in infection.
        /*! 
    This is computed as 0.19 (the value S from a neg bin mass action model fitted 
    to Saradidi data, divided by 0.302 (the ratio of body surface area in a 
    0.5-6 year old child (as per Saradidi) to adult) 
    \sa getExpectedNumberOfInfections() 
    */
        return 0.702;
    } else {
        // S_2(i,t) from AJTMH 75 (suppl 2) p12 eqn. (7)
        return Simm + (1.0 - Simm) / (1.0 + pow(human.cumulative_eir_a * Xstar_pInv, gamma_p));
    }
}

double InfectionIncidenceModel::expectedNumNewInfections(Human& human, double effectiveEIR) {
    // error check (should be OK if kappa is checked, for nonVector model):
    if (!(std::isfinite)(effectiveEIR)) {
        ostringstream out;
        out << "effectiveEIR is not finite: " << effectiveEIR << endl;
        throw TRACED_EXCEPTION(out.str(), Error::EffectiveEIR);
    }

    double expectedNumInfections = getModelExpectedInfections(human, effectiveEIR);

    // Only to be consistent with old simulation runs when set to false
    // Setting this option to true will only affect reporting
    if (opt_vaccine_genotype == false)
        //Introduce the effect of vaccination. Note that this does not affect cumEIR.
        expectedNumInfections *= human.vaccine_factor(Vaccine::PEV);

    //Update pre-erythrocytic immunity
    human.cumulative_eir_a += effectiveEIR;

    human.p_infected = 1.0 - exp(-expectedNumInfections) * (1.0 - human.p_infected);
    if (human.p_infected < 0.0)
        human.p_infected = 0.0;
    else if (human.p_infected > 1.0)
        human.p_infected = 1.0;

    return expectedNumInfections;
}

int InfectionIncidenceModel::numNewInfections(Human& human, double expectedNumInfections) {
    if (expectedNumInfections > 0.0000001) {
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        int n = human.rng().poisson(expectedNumInfections);

        n = SlotAllocator::clamp_count(n); // post-draw clamp
        return n;
    }

    if ((std::isnan)(expectedNumInfections)) { // check for not-a-number
        // bad Params::BASELINE_AVAILABILITY_SHAPE ?
        throw TRACED_EXCEPTION("numNewInfections: NaN", Error::NumNewInfections);
    }
    return 0;
}

void InfectionIncidenceModel::reportNumNewInfections(Human& human, int newNumInfections_i, int newNumInfections_l) {
    recordEvent(measure("nNewInfections"), human, newNumInfections_i + newNumInfections_l);
    recordEvent(measure("nNewInfections_Introduced"), human, newNumInfections_i);
    recordEvent(measure("nNewInfections_Indigenous"), human, newNumInfections_l);
    ctsNewInfections += newNumInfections_i + newNumInfections_l;
}

}
