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

// model/xml/ conversion layer — vivax within-host model: <vivax> description
// and liver-stage drug parameters from scnXml. Moved verbatim from
// Host/WithinHost/WHVivax.cpp (the HypnozoiteReleaseDistribution struct moved
// to WHVivax.h, shared with the model code). State stays defined in
// WHVivax.cpp (extern re-declared below).

#include "Clock.h"
#include "RunContext.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "util/sampler.h"
#include "xml/latentP.h"
#include "schema/scenario.h"

#include <cassert>
#include <limits>
#include <memory>
#include <set>
#include <string>

using namespace std;

namespace OM {
namespace WithinHost {
using namespace OM::util;

// State defined in WHVivax.cpp (external linkage):
extern SimTime latentP, bloodStageProtectionLatency;
extern double probBloodStageInfectiousToMosq;
extern int maxNumberHypnozoites;
extern double baseNumberHypnozoites;
extern unique_ptr<HypnozoiteReleaseDistribution> latentRelapse1st, latentRelapse2nd;
extern double pSecondRelease;
extern WeibullSampler bloodStageLength;    // units: days
extern double pPrimaryA, pPrimaryB, pRelapseOneA, pRelapseOneB,
    pRelapseTwoA, pRelapseTwoB, pEventIsSevere;
extern string vivaxClinOption;
extern map<double,int> nHypnozoitesProbMap;
extern bool ignoreNoPQ;
extern double pHetNoPQ, pReceivePQ, effectivenessPQ;

std::unique_ptr<HypnozoiteReleaseDistribution> createHypnozoiteReleaseDistribution(RunContext& ctx, const scnXml::HypnozoiteReleaseDistribution& elt)
{
    return std::make_unique<HypnozoiteReleaseDistribution>(ctx, createSampler(elt),  elt.getLatentRelapse());
}

void initNHypnozoites(){
    assert(baseNumberHypnozoites <= 1 && baseNumberHypnozoites >= 0);
    double total = 0.0;
    for( int n = 0; n <= maxNumberHypnozoites; ++n )
        total += pow( baseNumberHypnozoites, n );
    
    double cumP = 0.0;
    for( int n = 0; n <= maxNumberHypnozoites; ++n ){
        cumP += pow( baseNumberHypnozoites, n ) / total;
        // pair n with the cumulative probability of sampling n:
        nHypnozoitesProbMap[cumP] = n;
    }
}

void WHInterface::initVivax( RunContext& ctx, const OM::Parameters& parameters, const scnXml::Model& model ){

    // Default value to use if no <parameters> element is present in XML.
    latentP = read_latentP( ctx.clock, model );

    if( !model.getVivax().present() )
        throw util::xml_scenario_error( "no vivax model description in scenario XML" );
    const scnXml::Vivax& elt = model.getVivax().get();
    probBloodStageInfectiousToMosq = elt.getProbBloodStageInfectiousToMosq().getValue();
    maxNumberHypnozoites = elt.getHypnozoiteRelease().getNumberHypnozoites().getMax();
    baseNumberHypnozoites = elt.getHypnozoiteRelease().getNumberHypnozoites().getBase();
    const scnXml::HypnozoiteRelease& hr = elt.getHypnozoiteRelease();
    latentRelapse1st = createHypnozoiteReleaseDistribution(ctx, hr.getFirstReleaseDays());
    if(hr.getSecondReleaseDays().present()){
        latentRelapse2nd = createHypnozoiteReleaseDistribution(ctx, hr.getSecondReleaseDays().get());
        pSecondRelease = hr.getPSecondRelease();
        assert( pSecondRelease >= 0 && pSecondRelease <= 1 );
    } // else pSecondRelease is NaN and other values don't get used
    bloodStageProtectionLatency = ctx.clock.roundToTSFromDays( elt.getBloodStageProtectionLatency().getValue() );
    bloodStageLength.setParams(elt.getBloodStageLengthDays());
    
    const scnXml::ClinicalEvents& ce = elt.getClinicalEvents();
    pPrimaryA = ce.getPPrimaryInfection().getA();
    pPrimaryB = ce.getPPrimaryInfection().getB();
    pRelapseOneA = ce.getPRelapseOne().getA();
    pRelapseOneB = ce.getPRelapseOne().getB();
    pRelapseTwoA = ce.getPRelapseTwoPlus().getA();
    pRelapseTwoB = ce.getPRelapseTwoPlus().getB();
    pEventIsSevere = ce.getPEventIsSevere().getValue();
    vivaxClinOption = ce.getVivaxClinOption();
    
	// Define accepted values for vivaxClinOption
    const std::set<std::string> acceptedOptions = {"A1j", "A2j", "B1j", "B2j"};

   // Check if vivaxClinOption is valid
    if (acceptedOptions.find(vivaxClinOption) == acceptedOptions.end()) {
        throw util::xml_scenario_error("Invalid vivaxClinOption: " + vivaxClinOption +
                                     ". Accepted values are: A1j, A2j, B1j, B2j.");
    }
		
    initNHypnozoites();
    Pathogenesis::PathogenesisModel::init( ctx, parameters, model.getClinical(), true );
}

void WHInterface::setHSParametersVivax(const scnXml::LiverStageDrug* elt){
    double oldPHetNoPQ = pHetNoPQ;
    if( elt == 0 ){
        ignoreNoPQ = false;
        pHetNoPQ = 0.0;
        pReceivePQ = 0.0;
        effectivenessPQ = 1.0;  // sensible default: does not affect simple liver stage treatment option
    }else{
        ignoreNoPQ = elt->getIgnoreCannotReceive().present() ?
            elt->getIgnoreCannotReceive().get().getValue() : false;
        pHetNoPQ = elt->getPHumanCannotReceive().getValue();
        pReceivePQ = elt->getPUseUncomplicated().present() ?
            elt->getPUseUncomplicated().get().getValue() :  0.0;
        if( pReceivePQ > 0.0 ){
            if( util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS) ){
                cerr << "Deprecation warning: pUseUncomplicated is deprecated; it is "
                    "suggested to use the liver stage simple treatment option instead." << endl;
            }
        }
        effectivenessPQ = elt->getEffectivenessOnUse().getValue();
    }
    
    if( !std::isnan(oldPHetNoPQ) && oldPHetNoPQ != pHetNoPQ ){
        throw util::xml_scenario_error( "changeHS cannot change pHumanCannotReceive value" );
    }
}
}
}
