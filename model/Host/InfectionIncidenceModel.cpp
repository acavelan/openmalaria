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
#include "Host/Human.h"
#include "Transmission/PerHost.h"
#include "Host/WithinHost/WHInterface.h"
#include "Parameters.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/errors.h"

#include <stdexcept>
#include <cmath>

namespace OM { namespace Host {
    using namespace OM::util;

// -----  initialisation  -----

void InfectionIncidenceParams::init ( const OM::Parameters& parameters ) {
    baselineAvailShapeParam = parameters[Parameter::BASELINE_AVAILABILITY_SHAPE];
    baselineAvailOffset = 0.5 * pow(baselineAvailShapeParam, 2);
    
    gammaP=parameters[Parameter::GAMMA_P];
    Sinf=1-exp(-parameters[Parameter::NEG_LOG_ONE_MINUS_SINF]);
    Simm=parameters[Parameter::SIMM];
    EstarInv = 1.0/parameters[Parameter::E_STAR];
    XstarPInv = 1.0/parameters[Parameter::X_STAR_P];
    
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
    r_square_Gamma=0.649; //such that r_square_LogNormal =0.5
    
    noPreErythrocytic = util::ModelOptions::option (util::NO_PRE_ERYTHROCYTIC);
    vaccineGenotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    if (util::ModelOptions::option (util::NEGATIVE_BINOMIAL_MASS_ACTION)) {
        kind = NegBinomMAII;
        infRateShapeParam = (baselineAvailShapeParam+1.0) / (r_square_Gamma*baselineAvailShapeParam - 1.0);
        infRateShapeParam=std::max(infRateShapeParam, 0.0);
    } else if (util::ModelOptions::option (util::LOGNORMAL_MASS_ACTION)) {
        kind = LogNormalMAII;
        
        //! constant defining the constraint for the log Normal variance
        /// Used for the case where availability is assumed log Normally distributed
        double r_square_LogNormal = log(1.0+r_square_Gamma);
        
        infRateShapeParam = sqrt(r_square_LogNormal - 1.86 * pow(baselineAvailShapeParam, 2));
        infRateShapeParam=std::max(infRateShapeParam, 0.0);
        infRateOffset = 0.5 * pow(infRateShapeParam, 2);
        if( (std::isnan)(infRateShapeParam) ){
            throw util::xml_scenario_error( "bad parameter 16 (BASELINE_AVAILABILITY_SHAPE)" );
        }
    }else{
        if( util::ModelOptions::option( util::TRANS_HET ) ||
            util::ModelOptions::option( util::COMORB_TRANS_HET ) ||
            util::ModelOptions::option( util::TRANS_TREAT_HET ) ||
            util::ModelOptions::option( util::TRIPLE_HET ) )
        {
            kind = HeterogeneityWorkaround;
            cerr << "Warning: will use heterogeneity workaround." << endl;
        }
    }
}

// -----  free functions  -----

void InfectionIncidenceModel::summarize (const Host::Human& human) const {
    mon::recordStat(mon::measure("nExpectd"), human, m_pInfected);
}

double getAvailabilityFactor (LocalRng& rng, double baseAvailability,
                              const InfectionIncidenceParams& params) {
    switch (params.kind) {
    case InfectionIncidenceParams::NegBinomMAII:
        // Gamma sample with k=BaselineAvailabilityShapeParam; mean is baseAvailability
        return rng.gamma(params.baselineAvailShapeParam,
                         baseAvailability/params.baselineAvailShapeParam);
    case InfectionIncidenceParams::LogNormalMAII:
        // given BaselineAvailabilityShapeParam = sqrt (log (1 + variance/mean²))
        // and baseAvailability = mean, this is a draw from the log-normal distribution.
        if( baseAvailability != 1.0 ){
            // NOTE: shouldn't the normal_mean parameter be adjusted when baseAvailability != 1.0?
            throw TRACED_EXCEPTION_DEFAULT("LogNormalMAII::getAvailabilityFactor");
        }
        return rng.log_normal (log(baseAvailability) - params.baselineAvailOffset,
                               params.baselineAvailShapeParam);
    default:
        return baseAvailability;
    }
}

double susceptibility (const InfectionIncidenceModel& model,
                       const InfectionIncidenceParams& params) {
    if (params.noPreErythrocytic) {
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
        return params.Simm + (1.0-params.Simm) /
            (1.0 + pow(model.m_cumulativeEIRa*params.XstarPInv, params.gammaP));
    }
}

double getModelExpectedInfections (const InfectionIncidenceModel& model,
                                   LocalRng& rng, double effectiveEIR,
                                   const Transmission::PerHost& phTrans,
                                   const InfectionIncidenceParams& params) {
    switch (params.kind) {
    case InfectionIncidenceParams::HeterogeneityWorkaround:
        return (params.Sinf+(1-params.Sinf) / 
            (1 + effectiveEIR/(sim::oneTS()*phTrans.relativeAvailabilityHet)*params.EstarInv)) *
            susceptibility (model, params) * effectiveEIR;
    case InfectionIncidenceParams::NegBinomMAII:
        // Documentation: http://www.plosmedicine.org/article/fetchSingleRepresentation.action?uri:info:doi:10.1371/journal.pmed.1001157.s009
        return rng.gamma(params.infRateShapeParam,
            effectiveEIR * susceptibility (model, params) / params.infRateShapeParam);
    case InfectionIncidenceParams::LogNormalMAII: {
        // Documentation: http://www.plosmedicine.org/article/fetchSingleRepresentation.action?uri:info:doi:10.1371/journal.pmed.1001157.s009
        if (effectiveEIR < 1e-10)
            return 0.0;
        double meanlog = log(effectiveEIR * susceptibility (model, params)) - params.infRateOffset;
        return rng.log_normal(meanlog, params.infRateShapeParam);
    }
    default:
        // First two lines are availability adjustment: S_1(i,t) from AJTMH 75 (suppl 2) p12 eqn. (5)
        // Note that NegBinomMAII and LogNormalMAII supercede this model; see above
        return (params.Sinf+(1-params.Sinf) / 
            (1 + effectiveEIR/sim::oneTS()*params.EstarInv)) *
            susceptibility (model, params) * effectiveEIR;
    }
}

double expectedNumNewInfections (InfectionIncidenceModel& model, Human& human,
                                 double effectiveEIR,
                                 const InfectionIncidenceParams& params, SimTime ts1)
{
    // error check (should be OK if kappa is checked, for nonVector model):
    if( !(std::isfinite)(effectiveEIR) ){
        ostringstream out;
        out << "effectiveEIR is not finite: " << effectiveEIR << endl;
        throw TRACED_EXCEPTION (out.str(), util::Error::EffectiveEIR);
    }

    double expectedNumInfections = getModelExpectedInfections (model, human.rng, effectiveEIR, human.perHostTransmission, params);    
    
    // Only to be consistent with old simulation runs when set to false
    // Setting this option to true will only affect reporting
    if(params.vaccineGenotype == false)
        //Introduce the effect of vaccination. Note that this does not affect cumEIR.
        expectedNumInfections *= human.vaccine.getFactor( interventions::Vaccine::PEV, ts1 );
  
    //Update pre-erythrocytic immunity
    model.m_cumulativeEIRa += effectiveEIR;

    model.m_pInfected = 1.0 - exp(-expectedNumInfections) * (1.0-model.m_pInfected);
    if (model.m_pInfected < 0.0)
        model.m_pInfected = 0.0;
    else if (model.m_pInfected > 1.0)
        model.m_pInfected = 1.0;
    
    return expectedNumInfections;
}

int numNewInfections (Human& human, double expectedNumInfections)
{ 
  if (expectedNumInfections > 0.0000001){
    int n = human.rng.poisson(expectedNumInfections);

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

void reportNumNewInfections (Human& human, int newNumInfections_i,
                             int newNumInfections_l, int& ctsNewInfections)
{
    mon::recordEvent(mon::measure("nNewInfections"), human, newNumInfections_i+newNumInfections_l);
    mon::recordEvent(mon::measure("nNewInfections_Introduced"), human, newNumInfections_i);
    mon::recordEvent(mon::measure("nNewInfections_Indigenous"), human, newNumInfections_l);
    ctsNewInfections += newNumInfections_i+newNumInfections_l;
}

} }
