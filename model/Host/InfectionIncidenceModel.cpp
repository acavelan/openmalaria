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

#include "Host/InfectionIncidenceModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "util/errors.h"

namespace OM::Host {
using util::LocalRng;

IncidenceParameters makeIncidenceParameters(const IncidenceConfiguration& config) {
    IncidenceParameters result;
    result.baseline_avail_shape_param = config.baselineAvailabilityShape;
    result.baseline_avail_offset = 0.5 * pow(result.baseline_avail_shape_param, 2);

    result.gamma_p=config.gammaP;
    result.Sinf=1-exp(-config.negativeLogOneMinusSinf);
    result.Simm=config.Simm;
    result.EstarInv = 1.0/config.EStar;
    result.Xstar_pInv = 1.0/config.XStarP;

    // Constraint on Gamma availability variance; log(1 + 0.649) is approximately 0.5.
    const double r_square_Gamma = 0.649;

    result.noPreErythrocytic = config.noPreErythrocytic;
    const bool negativeBinomial = config.kind == IncidenceKind::NegativeBinomial;
    result.vaccineGenotype = config.vaccineGenotype;

    if (negativeBinomial) {
        result.kind = IncidenceKind::NegativeBinomial;
        result.inf_rate_shape_param = (result.baseline_avail_shape_param+1.0) / (r_square_Gamma*result.baseline_avail_shape_param - 1.0);
        result.inf_rate_shape_param=std::max(result.inf_rate_shape_param, 0.0);
    } else if (config.kind == IncidenceKind::Lognormal) {
        result.kind = IncidenceKind::Lognormal;

        //! constant defining the constraint for the log Normal variance
        /// Used for the case where availability is assumed log Normally distributed
        double r_square_LogNormal = log(1.0+r_square_Gamma);

        result.inf_rate_shape_param = sqrt(r_square_LogNormal - 1.86 * pow(result.baseline_avail_shape_param, 2));
        result.inf_rate_shape_param=std::max(result.inf_rate_shape_param, 0.0);
        result.inf_rate_offset = 0.5 * pow(result.inf_rate_shape_param, 2);
        if( (std::isnan)(result.inf_rate_shape_param) ){
            throw util::xml_scenario_error( "bad parameter 16 (BASELINE_AVAILABILITY_SHAPE)" );
        }
    }else{
        if( config.kind == IncidenceKind::Heterogeneity )
        {
            result.kind = IncidenceKind::Heterogeneity;
            cerr << "Warning: will use heterogeneity workaround." << endl;
        }
    }

    return result;
}


double sampleAvailability(const IncidenceParameters& p, LocalRng& rng) {
    switch (p.kind) {
    case IncidenceKind::NegativeBinomial:
        return rng.gamma(p.baseline_avail_shape_param, 1.0 / p.baseline_avail_shape_param);
    case IncidenceKind::Lognormal:
        return rng.log_normal(log(1.0) - p.baseline_avail_offset, p.baseline_avail_shape_param);
    default:
        return 1.0;
    }
}

static double susceptibility(const IncidenceState& state, const IncidenceParameters& p) {
    // Smith et al., AJTMH 75 (suppl 2), equation (7).
    if (p.noPreErythrocytic) return 0.702;
    return p.Simm + (1.0-p.Simm) /
        (1.0 + pow(state.cumulativeEIR*p.Xstar_pInv, p.gamma_p));
}

static double modelExpectedInfections(const IncidenceState& state,
        const IncidenceParameters& p, LocalRng& rng, double effectiveEIR,
        double relativeAvailability, SimTime interval) {
    switch (p.kind) {
    case IncidenceKind::Default:
        return (p.Sinf+(1-p.Sinf) /
            (1 + effectiveEIR/interval*p.EstarInv)) *
            susceptibility(state, p) * effectiveEIR;
    case IncidenceKind::Heterogeneity:
        // Preserve the historical availability adjustment.
        return (p.Sinf+(1-p.Sinf) /
            (1 + effectiveEIR/(interval*relativeAvailability)*p.EstarInv)) *
            susceptibility(state, p) * effectiveEIR;
    case IncidenceKind::NegativeBinomial:
        return rng.gamma(p.inf_rate_shape_param,
            effectiveEIR * susceptibility(state, p) / p.inf_rate_shape_param);
    case IncidenceKind::Lognormal:
        if (effectiveEIR < 1e-10) return 0.0;
        double meanlog = log(effectiveEIR * susceptibility(state, p)) - p.inf_rate_offset;
        return rng.log_normal(meanlog, p.inf_rate_shape_param);
    }
    throw TRACED_EXCEPTION_DEFAULT("invalid incidence model");
}

double expectedInfections(IncidenceState& state, const IncidenceParameters& p,
        LocalRng& rng, double effectiveEIR, double relativeAvailability,
        double vaccineFactor, SimTime interval)
{
    // error check (should be OK if kappa is checked, for nonVector model):
    if( !(std::isfinite)(effectiveEIR) ){
        ostringstream out;
        out << "effectiveEIR is not finite: " << effectiveEIR << endl;
        throw TRACED_EXCEPTION (out.str(), util::Error::EffectiveEIR);
    }

    double expectedNumInfections = modelExpectedInfections(state, p, rng, effectiveEIR, relativeAvailability, interval);

    // Only to be consistent with old simulation runs when set to false
    // Setting this option to true will only affect reporting
    if(p.vaccineGenotype == false)
        //Introduce the effect of vaccination. Note that this does not affect cumEIR.
        expectedNumInfections *= vaccineFactor;

    //Update pre-erythrocytic immunity
    state.cumulativeEIR += effectiveEIR;

    state.probabilityInfected = 1.0 - exp(-expectedNumInfections) * (1.0-state.probabilityInfected);
    if (state.probabilityInfected < 0.0)
        state.probabilityInfected = 0.0;
    else if (state.probabilityInfected > 1.0)
        state.probabilityInfected = 1.0;

    return expectedNumInfections;
}

int sampleNewInfections(LocalRng& rng, double expectedNumInfections)
{
  if (expectedNumInfections > 0.0000001){
    int n = rng.poisson(expectedNumInfections);

    if( n > WithinHost::WHInterface::MAX_INFECTIONS )
        n = WithinHost::WHInterface::MAX_INFECTIONS;
    return n;
  }

  if ( (std::isnan)(expectedNumInfections) ){	// check for not-a-number
      // bad Params::BASELINE_AVAILABILITY_SHAPE ?
      throw TRACED_EXCEPTION( "numNewInfections: NaN", util::Error::NumNewInfections );
  }
  return 0;
}

}
