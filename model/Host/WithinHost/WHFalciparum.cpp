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

#include "Host/WithinHost/WHFalciparum.h"
#include "Host/WithinHost/DescriptiveWithinHost.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Treatment.h"
#include "Host/WithinHost/Genotypes.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/checkpoint_containers.h"

#include <cmath>
#include <gsl/gsl_cdf.h>


namespace OM {
namespace WithinHost {

using namespace OM::util;

SimTime Infection::s_latentP = sim::never();

// -----  Non-static  -----

WHFalciparum::WHFalciparum( LocalRng& rng, double comorbidityFactor, const WithinHostParameters& config ):
    WHInterface(),
    m_cumulative_h(0.0), m_cumulative_Y(0.0), m_cumulative_Y_lag(0.0),
    totalDensity(0.0), hrp2Density(0.0), timeStepMaxDensity(0.0),
    pathogenesisModel(comorbidityFactor, config.pathogenesis),
    treatExpiryLiver(0), treatExpiryBlood(0), config(config)
{
    // NOTE: negating a Gaussian sample with mean 0 is pointless — except that
    // the individual samples change. In any case the overhead is negligible.
    //FIXME: Should this be allowed to be greater than 1?
    // Oldest code on GoogleCode: _innateImmunity=(double)(W_GAUSS((0), (sigma_i)));
    _innateImmSurvFact = exp(-rng.gauss(0.0, config.sigma_i));
    
    m_y_lag_i.resize(config.y_lag_len * Genotypes::N());
    m_y_lag_l.resize(config.y_lag_len * Genotypes::N());
}

WHFalciparum::~WHFalciparum()
{
}

double WHFalciparum::immunitySurvivalFactor (double ageInYears, double cumulativeExposureJ) {
    if (std::isnan(ageInYears) || std::isnan(cumulativeExposureJ) ||
        std::isnan(m_cumulative_h) || std::isnan(m_cumulative_Y)) {
        throw base_exception("nan in immunitySurvivalFactor");
    }
    
    // Documentation: AJTMH pp22-23
    
    // Effect of cumulative Parasite density (named Dy in AJTM)
    double dY = 1.0;
    // Effect of number of infections experienced since birth (named Dh in AJTM)
    double dH = 1.0;
    if (m_cumulative_h > 1.0) {
        dH = 1.0 / (1.0 + (m_cumulative_h - 1.0) * config.invCumulativeHstar);
        dY = 1.0 / (1.0 + (m_cumulative_Y - cumulativeExposureJ) * config.invCumulativeYstar);
    }
    
    // Effect of age-dependent maternal immunity (named Dm in AJTM)
    double dA = 1.0 - config.alpha_m * exp(-config.decayM * ageInYears);
    
    return util::streamValidate( std::min(dY*dH*dA, 1.0) );
}

// Infectiousness parameters: see AJTMH p.33; tau=1/sigmag**2 
const double PTM_beta1=1.0;
const double PTM_beta2=0.46;
const double PTM_beta3=0.17;
const double PTM_tau= 0.066;
const double PTM_tau_prime = 1.0 / sqrt(1.0 / PTM_tau);
const double PTM_mu= -8.1;

inline double compute_y_lag(const vector<double> &y_lag, int y_lag_len, vector<double> &y_lag_g)
{
    const size_t n = Genotypes::N();

    // Take weighted sum of total asexual blood stage density 10, 15 and 20 days
    // before. Add y_lag_len to index to ensure positive.
    size_t d10 = mod_nn(y_lag_len + sim::inSteps(sim::ts1() - sim::fromDays(10)), y_lag_len);
    size_t d15 = mod_nn(y_lag_len + sim::inSteps(sim::ts1() - sim::fromDays(15)), y_lag_len);
    size_t d20 = mod_nn(y_lag_len + sim::inSteps(sim::ts1() - sim::fromDays(20)), y_lag_len);

    double y10_sum = 0.0, y15_sum = 0.0, y20_sum = 0.0;
    for( size_t g = 0; g < n; ++g ){
        double y10 = y_lag[d10 * n + g];
        double y15 = y_lag[d15 * n + g];
        double y20 = y_lag[d20 * n + g];
        y10_sum += y10;
        y15_sum += y15;
        y20_sum += y20;
        y_lag_g[g] = PTM_beta1 * y10 + PTM_beta2 * y15 + PTM_beta3 * y20;
    }
    
    // Weighted sum:
    return PTM_beta1 * y10_sum + PTM_beta2 * y15_sum + PTM_beta3 * y20_sum;
}

double WHFalciparum::probTransmissionToMosquito(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l) const{
    // This model (often referred to as the gametocyte model) was designed for
    // 5-day time steps. We use the same model (sampling 10, 15 and 20 days
    // ago) for 1-day time steps to avoid having to design and analyse a new
    // model. Description: AJTMH pp.32-33 and p9.
    
    // Note: we don't allow for gametocydal treatments (e.g. Primaquine).
    const size_t n = Genotypes::N();

    // Sum lagged densities across genotypes:
    vector<double> y_lag_g_i(n), y_lag_g_l(n);
    const double y_lag_sum_i = compute_y_lag(m_y_lag_i, config.y_lag_len, y_lag_g_i);
    const double y_lag_sum_l = compute_y_lag(m_y_lag_l, config.y_lag_len, y_lag_g_l);
    const double y_lag_sum = y_lag_sum_i + y_lag_sum_l;

    if( y_lag_sum < 0.001 ) return 0.0; // cut off for uninfectious humans
    
    // Get a zval, convert to equivalent Normal sample:
    const double zval = (log(y_lag_sum) + PTM_mu) * PTM_tau_prime;
    const double pone = gsl_cdf_ugaussian_P(zval);
    double pTransmit = pone*pone;

    // pTransmit has to be between 0 and 1:
    pTransmit=std::max(pTransmit, 0.0);
    pTransmit=std::min(pTransmit, 1.0);
    
    if(pTransmit <= 0.0)
        return pTransmit;

    for( size_t g = 0; g < n; ++g )
    {
        probTransGenotype_i[g] = pTransmit * y_lag_g_i[g] / y_lag_sum;
        probTransGenotype_l[g] = pTransmit * y_lag_g_l[g] / y_lag_sum;
    }

    // Include here the effect of transmission-blocking vaccination:
    util::streamValidate( pTransmit );
    return pTransmit;
}

bool WHFalciparum::diagnosticResult( LocalRng& rng, const Diagnostic& diagnostic ) const{
    return diagnostic.isPositive( rng, totalDensity, hrp2Density );
}


bool WHFalciparum::treatSimple( Host::Human& human, SimTime timeLiver, SimTime timeBlood ){
    if( timeLiver != sim::zero() ){
        if( timeLiver < sim::zero() )
            clearInfections( Treatment::LIVER );
        else
            treatExpiryLiver = max( int(treatExpiryLiver), int(sim::nowOrTs1()) + timeLiver );
        mon::recordEvent(mon::measure("nLiverStageTreatments"), human);
    }
    if( timeBlood != sim::zero() ){
        if( timeBlood < sim::zero() )
            clearInfections( Treatment::BLOOD );
        else
            treatExpiryBlood = max( int(treatExpiryBlood), int(sim::nowOrTs1()) + timeBlood );
        return true;    // blood stage treatment
    }
    return false;    // no blood stage treatment
}

Pathogenesis::StatePair WHFalciparum::determineMorbidity( Host::Human& human,
        double ageYears, bool isDoomed )
{
    Pathogenesis::StatePair result = pathogenesisModel.determineState( human,
            ageYears, timeStepMaxDensity, totalDensity, isDoomed );
    
    /* Note: this model can easily be re-enabled, but is not used and not considered to be a good model.
    if( (result.state & Pathogenesis::MALARIA) && util::ModelOptions::option( util::PENALISATION_EPISODES ) ){
        // This does immunity penalisation:
        m_cumulative_Y = m_cumulative_Y_lag - immPenalty_22*(m_cumulative_Y-m_cumulative_Y_lag);
        if (m_cumulative_Y < 0) {
            m_cumulative_Y=0.0;
        }
    }*/
    
    return result;
}


// -----  immunity  -----

void WHFalciparum::updateImmuneStatus() {
    if (config.immEffectorRemain < 1) {
        m_cumulative_h *= config.immEffectorRemain;
        m_cumulative_Y *= config.immEffectorRemain;
    }
    if (config.asexImmRemain < 1) {
        m_cumulative_h *= config.asexImmRemain /
                      (1+(m_cumulative_h*(1-config.asexImmRemain) * config.invCumulativeHstar));
        m_cumulative_Y *= config.asexImmRemain /
                      (1+(m_cumulative_Y*(1-config.asexImmRemain) * config.invCumulativeYstar));
    }
    m_cumulative_Y_lag = m_cumulative_Y;
}


// -----  Checkpointing  -----

void WHFalciparum::checkpoint (istream& stream) {
    WHInterface::checkpoint( stream );
    _innateImmSurvFact & stream;
    m_cumulative_h & stream;
    m_cumulative_Y & stream;
    m_cumulative_Y_lag & stream;
    totalDensity & stream;
    hrp2Density & stream;
    timeStepMaxDensity & stream;
    m_y_lag_i & stream;
    m_y_lag_l & stream;
    pathogenesisModel & stream;
    treatExpiryLiver & stream;
    treatExpiryBlood & stream;
}
void WHFalciparum::checkpoint (ostream& stream) {
    WHInterface::checkpoint( stream );
    _innateImmSurvFact & stream;
    m_cumulative_h & stream;
    m_cumulative_Y & stream;
    m_cumulative_Y_lag & stream;
    totalDensity & stream;
    hrp2Density & stream;
    timeStepMaxDensity & stream;
    m_y_lag_i & stream;
    m_y_lag_l & stream;
    pathogenesisModel & stream;
    treatExpiryLiver & stream;
    treatExpiryBlood & stream;
}

}
}
