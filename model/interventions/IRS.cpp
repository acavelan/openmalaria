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

#include "Global.h"
#include "core/RunExtensions.h"
#include "core/RunContext.h"
#include "interventions/IRS.h"
#include "Host/Human.h"
#include "util/errors.h"
#include "util/SpeciesIndexChecker.h"
#include "util/CommandLine.h"
#include "R_nmath/qnorm.h"
#include <cmath>

namespace OM { namespace interventions {

void IRSComponent::deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method, VaccineLimits )const{
    iv_deploy( ctx, human, *this );
    mon::recordDeploy(mon::measure("nMassIRS"), mon::measure("nCtsIRS"), human, method);
}

Component::Type IRSComponent::componentType()const{ return Component::IRS; }

void IRSComponent::print_details( std::ostream& out )const{
    out << id().id << "\tIRS";
}

IRSComponent::IRSAnopheles::RelativeAttractiveness::RelativeAttractiveness() :
    lPF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{}
IRSComponent::IRSAnopheles::SurvivalFactor::SurvivalFactor() :
    BF( numeric_limits< double >::signaling_NaN() ),
    PF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() ),
    invBaseSurvival( numeric_limits< double >::signaling_NaN() )
{}
void IRSComponent::IRSAnopheles::SurvivalFactor::init1(){
    BF = 0.0;
    PF = 0.0;
    insecticideScaling = 0.0;
    invBaseSurvival = 1.0;
}
double IRSComponent::IRSAnopheles::RelativeAttractiveness::relativeAttractiveness(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lPF*insecticideComponent );
    assert( relAvail>=0.0 );
    return relAvail;
}
double IRSComponent::IRSAnopheles::SurvivalFactor::survivalFactor(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double killingEffect = BF + PF*insecticideComponent;
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    assert( survivalFactor >= 0.0 );
    assert( survivalFactor <= 1.0 );
    return survivalFactor;
}

double IRSComponent::sampleInitialInsecticide(LocalRng& rng) const{
    double sample = initialInsecticide.sample(rng);
    if( sample < 0.0 )
        sample = 0.0;       // avoid negative samples
    if( sample > maxInsecticide )
        sample = maxInsecticide;
    return sample;
}


// ———  per-human data  ———
// (phase 7.2a: the HumanIRS class is gone; its draws live in iv_deploy and
// its effects in iv_* — interventions/IvState.cpp)

} }
