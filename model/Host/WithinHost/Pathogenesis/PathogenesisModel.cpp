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
#include <cmath>

namespace OM::WithinHost::Pathogenesis {
using util::LocalRng;

PathogenesisModel::PathogenesisModel(double cF, const PathogenesisParameters& config) :
    config(config), _comorbidityFactor(cF), _pyrogenThres(config.initPyroThres) {}

Pathogenesis::StatePair PathogenesisModel::determineState(Host::Human& human,
        double ageYears, double timeStepMaxDensity, double endDensity, bool isDoomed)
{
    double pMalariaFever = getPEpisode(timeStepMaxDensity, endDensity);
    StatePair result;
    //TODO(performance): would using a single RNG sample and manipulating probabilities be faster?
    //Decide whether a clinical episode occurs and if so, which type
    if( human.rng.bernoulli( pMalariaFever ) ){
        const double prSevereEpisode = timeStepMaxDensity / (timeStepMaxDensity + config.pg_severeMalThreshold);
        const double comorb_factor = _comorbidityFactor / (1.0 + ageYears * config.pg_inv_critAgeComorb);
        const double pCoinfection = config.pg_comorbIntercept * comorb_factor;

        // Expectation of a severe bout:
        const double exSevere = prSevereEpisode + (1.0 - prSevereEpisode) * pCoinfection;
        mon::recordStat(mon::measure("expectedSevere"), human, exSevere);
        mon::recordStat(mon::measure("expectedSevereWithoutComorbidities"), human, prSevereEpisode);

        if( human.rng.bernoulli( prSevereEpisode ) )
            result.state = STATE_SEVERE;
        else {
            if( human.rng.bernoulli( pCoinfection ) )
                result.state = STATE_COINFECTION;
            else
                result.state = STATE_MALARIA;
        }

        // Indirect mortality:
        // IndirectRisk is the probability of dying from the indirect effects
        // of malaria conditional on not having an acute attack of malaria
        double indirectRisk = config.pg_indirRiskCoFactor * comorb_factor;
        if( !isDoomed ){
            mon::recordStat(mon::measure("expectedIndirectDeaths"), human, indirectRisk);
        }
        //NOTE: we could only evaluate indirectMortality if not already doomed,
        // except that (a) it affects random numbers, and (b) it affects
        // evaluation of uncomplicated cases with the 5-day HS when
        // indirectMortBugfix is not enabled.
        if( human.rng.bernoulli( indirectRisk ) )
            result.indirectMortality = true;
    }else{
        result.state = sampleNMF( config, human.rng, ageYears );
    }
    return result;
}

Pathogenesis::State sampleNMF(const PathogenesisParameters& config, LocalRng& rng, double ageYears){
    if ( config.pg_NMF_incidence.isSet() ) {
        double pNMF = config.pg_NMF_incidence.eval( ageYears );
        if( rng.bernoulli( pNMF ) )
            return Pathogenesis::STATE_NMF;
    }
    return Pathogenesis::NONE;
}

double PathogenesisModel::getPEpisode(double timeStepMaxDensity, double totalDensity) {
    if (config.kind == PathogenesisParameters::Kind::Mueller) {
        double incidenceDensity = config.rateMultiplier_31
                * (pow(totalDensity, config.densityExponent_32));
        return 1.0 - exp(-incidenceDensity);
    }
    updatePyrogenThres(totalDensity);
    if (config.kind == PathogenesisParameters::Kind::Predetermined)
        return (timeStepMaxDensity > _pyrogenThres) ? 1 : 0;
    return timeStepMaxDensity / (timeStepMaxDensity + _pyrogenThres);
}

void PathogenesisModel::updatePyrogenThres(double totalDensity){
    // Note: this calculation is slow (something like 5% of runtime)

    //Numerical approximation to equation 2, AJTMH p.57
    for( size_t i = 1; i <= 11; ++i ){
        _pyrogenThres += totalDensity * config.a /
            ( (config.Ystar1_26 + totalDensity) * (config.Ystar2_13 + _pyrogenThres) )
            + config.b * _pyrogenThres;
    }
}

void PathogenesisModel::summarize (const Host::Human& human) const {
    if (config.kind == PathogenesisParameters::Kind::Mueller) return;
    mon::recordStat(mon::measure("sumPyrogenThresh"), human, _pyrogenThres);
    mon::recordStat(mon::measure("sumLogPyrogenThres"), human, log(_pyrogenThres+1.0) );
}
}
