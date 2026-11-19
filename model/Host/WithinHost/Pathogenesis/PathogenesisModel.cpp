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

#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "Host/Human.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/AgeGroupInterpolation.h"
#include "Parameters.h"

#include <cmath>
using namespace std;

namespace OM { namespace WithinHost { namespace Pathogenesis {

using namespace OM::util;
using util::AgeGroupInterpolator;

//BEGIN static
/// Comorbidity prevalence at birth as a risk factor for indirect mortality
double pg_indirRiskCoFactor;
/// sevMal: critical density for severe malaria bout (Y*B1)
double pg_severeMalThreshold;
/// Comorbidity prevalence at birth as a risk factor for severe
double pg_comorbIntercept;
/// One over critical age for co-morbidity (for both severe and indirect)
double pg_inv_critAgeComorb;

/// Rate of Non-Malaria Fever incidence by age. Non-seasonal.
AgeGroupInterpolator pg_NMF_incidence;

bool opt_predetermined_episodes = false, opt_mueller_pres_model = false;

// ———  Müller model parameters (origin: Submodels.cpp)  ———
double rateMultiplier_31;
double densityExponent_32;

//Pyrogenic threshold at birth (Y*0)
double initPyroThres;
// Ystar2: critical value in determining increase in pyrogenic threshold
double Ystar2_13;
double Ystar1_26;

//Number of categories in the numerical approximation used in updatePyrogenThres()
const size_t n = 11;
// Derived parameters without good names:
double a = numeric_limits<double>::signaling_NaN(),
    b = numeric_limits<double>::signaling_NaN();
//END static


PathogenesisModel* PathogenesisModel::createPathogenesisModel(double cF) {
    // Kind from model options (same dispatch chain as the old subclass factory)
    if (opt_predetermined_episodes) return new PathogenesisModel(cF, Kind::Predet);
    if (opt_mueller_pres_model) return new PathogenesisModel(cF, Kind::Mueller);
    return new PathogenesisModel(cF, Kind::Pyrogen);
}

PathogenesisModel::PathogenesisModel(double cF, Kind kind) :
    _comorbidityFactor(cF), kind(kind),
    _pyrogenThres(kind == Kind::Mueller ? 0.0 : initPyroThres)
{}

void PathogenesisModel::initMueller( RunContext& ctx, const Parameters& parameters ){
    // parity: old MuellerPathogenesis::init (Submodels.cpp:41)
    rateMultiplier_31 = parameters[Parameter::MUELLER_RATE_MULTIPLIER]
        * ctx.clock.years_per_step;
    densityExponent_32 = parameters[Parameter::MUELLER_DENSITY_EXPONENT];
}

void PathogenesisModel::initPyrogen( RunContext& ctx, const Parameters& parameters ){
    // parity: old PyrogenPathogenesis::init (Submodels.cpp:68)
    initPyroThres = parameters[Parameter::Y_STAR_0];
    
    double delt = 1.0 / n;
    double smuY = -log(0.5) /
        (ctx.clock.steps_per_year * parameters[Parameter::Y_STAR_HALF_LIFE]);
    b = -smuY * delt;
    
    Ystar2_13 = parameters[Parameter::Y_STAR_SQ];
    
    //alpha: factor determining increase in pyrogenic threshold
    double alpha14 = parameters[Parameter::ALPHA];
    a = alpha14 * ctx.clock.one_TS() * delt;  
    
    //Ystar1: critical value of parasite density in determing increase in pyrog t
    Ystar1_26 = parameters[Parameter::Y_STAR_1];
}

double PathogenesisModel::getPEpisode(double timeStepMaxDensity, double totalDensity) {
    switch(kind) {
    case Kind::Mueller: {
        // parity: old MuellerPathogenesis::getPEpisode (Submodels.cpp:47)
        double incidenceDensity = rateMultiplier_31
                * (pow(totalDensity, densityExponent_32));
        return 1.0 - exp(-incidenceDensity);
    }
    case Kind::Predet:
        // parity: old PredetPathogenesis::getPEpisode (Submodels.cpp:125)
        updatePyrogenThres(totalDensity);
        return (timeStepMaxDensity > _pyrogenThres) ? 1 : 0;
    case Kind::Pyrogen:
        // parity: old PyrogenPathogenesis::getPEpisode (Submodels.cpp:91)
        updatePyrogenThres(totalDensity);
        return timeStepMaxDensity / (timeStepMaxDensity + _pyrogenThres);
    }
    return 0.0;  // unreachable: every Kind has a case
}

void PathogenesisModel::updatePyrogenThres(double totalDensity){
    // Note: this calculation is slow (something like 5% of runtime)
    
    //Numerical approximation to equation 2, AJTMH p.57
    for( size_t i = 1; i <= n; ++i ){
        _pyrogenThres += totalDensity * a /
            ( (Ystar1_26 + totalDensity) * (Ystar2_13 + _pyrogenThres) )
            + b * _pyrogenThres;
    }
}

void PathogenesisModel::summarize (const Host::Human& human) {
    if( kind == Kind::Mueller ) return;  // old base summarize was empty
    mon::recordStat(mon::measure("sumPyrogenThresh"), human, _pyrogenThres);
    mon::recordStat(mon::measure("sumLogPyrogenThres"), human, log(_pyrogenThres+1.0) );
}

Pathogenesis::StatePair PathogenesisModel::determineState(Host::Human& human,
        double ageYears, double timeStepMaxDensity, double endDensity, bool isDoomed)
{
    double pMalariaFever = getPEpisode(timeStepMaxDensity, endDensity);
    StatePair result;
    //TODO(performance): would using a single RNG sample and manipulating probabilities be faster?
    //Decide whether a clinical episode occurs and if so, which type
    // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
    if( human.rng.bernoulli( pMalariaFever ) ){
        const double prSevereEpisode = timeStepMaxDensity / (timeStepMaxDensity + pg_severeMalThreshold);
        const double comorb_factor = _comorbidityFactor / (1.0 + ageYears * pg_inv_critAgeComorb);
        const double pCoinfection = pg_comorbIntercept * comorb_factor;
        
        // Expectation of a severe bout:
        const double exSevere = prSevereEpisode + (1.0 - prSevereEpisode) * pCoinfection;
        mon::recordStat(mon::measure("expectedSevere"), human, exSevere);
        mon::recordStat(mon::measure("expectedSevereWithoutComorbidities"), human, prSevereEpisode);
        
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if( human.rng.bernoulli( prSevereEpisode ) )
            result.state = STATE_SEVERE;
        else {
            // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
            if( human.rng.bernoulli( pCoinfection ) )
                result.state = STATE_COINFECTION;
            else
                result.state = STATE_MALARIA;
        }

        // Indirect mortality:
        // IndirectRisk is the probability of dying from the indirect effects
        // of malaria conditional on not having an acute attack of malaria
        double indirectRisk = pg_indirRiskCoFactor * comorb_factor;
        if( !isDoomed ){
            mon::recordStat(mon::measure("expectedIndirectDeaths"), human, indirectRisk);
        }
        //NOTE: we could only evaluate indirectMortality if not already doomed,
        // except that (a) it affects random numbers, and (b) it affects
        // evaluation of uncomplicated cases with the 5-day HS when
        // indirectMortBugfix is not enabled.
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if( human.rng.bernoulli( indirectRisk ) )
            result.indirectMortality = true;        
    }else{
        result.state = sampleNMF( human.rng, ageYears );
    }
    return result;
}

Pathogenesis::State PathogenesisModel::sampleNMF( LocalRng& rng, double ageYears ){
    if ( pg_NMF_incidence.isSet() ) {
        double pNMF = pg_NMF_incidence.eval( ageYears );
        // parity: order — local RNG stream: draw call order fixes stream state
        if( rng.bernoulli( pNMF ) )
            return Pathogenesis::STATE_NMF;
    }
    return Pathogenesis::NONE;
}

} } }
