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

#include "RunContext.h"
#include "util/checkpoint.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "Host/Human.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/checkpoint_containers.h"
#include "util/UnitParse.h"
#include "util/CommandLine.h"
#include "util/sampler.h"
#include <algorithm>
#include <limits>
#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;



// ———  parameters  ———

// Set from the parameters block, or a default value is assumed.
SimTime latentP = sim::never();       // attribute on <parameters> element.

// Set from <vivax .../> element:
double probBloodStageInfectiousToMosq = numeric_limits<double>::signaling_NaN();
int maxNumberHypnozoites = -1;
double baseNumberHypnozoites = numeric_limits<double>::signaling_NaN();
std::unique_ptr<HypnozoiteReleaseDistribution> latentRelapse1st, latentRelapse2nd;
double pSecondRelease = numeric_limits<double>::signaling_NaN();
SimTime bloodStageProtectionLatency = sim::never();
WeibullSampler bloodStageLength;    // units: days
double pPrimaryA = numeric_limits<double>::signaling_NaN(),
    pPrimaryB = numeric_limits<double>::signaling_NaN();
double pRelapseOneA = numeric_limits<double>::signaling_NaN(),
    pRelapseOneB = numeric_limits<double>::signaling_NaN();
double pRelapseTwoA = numeric_limits<double>::signaling_NaN(),
    pRelapseTwoB = numeric_limits<double>::signaling_NaN();
double pEventIsSevere = numeric_limits<double>::signaling_NaN();
std::string vivaxClinOption = " ";

// Set from healthSystem element:
bool ignoreNoPQ = false;
double pHetNoPQ = numeric_limits<double>::signaling_NaN();
double pReceivePQ = numeric_limits<double>::signaling_NaN();
double effectivenessPQ = numeric_limits<double>::signaling_NaN();

#ifdef WHVivaxSamples
WHInterface *sampleHost = 0;
VivaxBrood *sampleBrood = 0;
#endif


// ———  individual models  ———


// number of hypnozoites per brood:
map<double,int> nHypnozoitesProbMap;


int sampleNHypnozoites(LocalRng& rng){
    // parity: order — local RNG stream: draw call order fixes stream state
    double x = rng.uniform_01();
    // upper_bound finds the first key (cumulative probability) greater than x:
    auto it = nHypnozoitesProbMap.upper_bound( x );
    assert( it != nHypnozoitesProbMap.end() );  // i.e. that we did find a real key
    return it->second;  // corresponding n
}


// time to hypnozoite release after initial release:
SimTime sampleReleaseDelay(LocalRng& rng){
    bool isFirstRelease = true;
    if( pSecondRelease == 1.0 ||
        // parity: order — local RNG stream: draw call order fixes stream state
        (pSecondRelease > 0.0 && rng.bernoulli(pSecondRelease)) ){
        isFirstRelease = false;
    }
    
    if (isFirstRelease){
        return latentRelapse1st->sampleReleaseDelay(rng);
    } else {
        return latentRelapse2nd->sampleReleaseDelay(rng);
    }
}


// ———  per-brood code  ———

VivaxBrood::VivaxBrood( RunContext& ctx, LocalRng& rng, int origin, WHInterface *host ) :
        ctx(ctx),
        primaryHasStarted( false ),
        relapseHasStarted( false ),
		relapsebHasStarted( false),
        hadEvent( false ),
        hadRelapse( false ),
		origin(origin)
 {
    set<SimTime> releases;     // used to initialise releaseDates; a set is better to use now but a vector later
    
    // primary blood stage plus hypnozoites (relapses)
    releases.insert( ctx.clock.nowOrTs0() + latentP );
    int numberHypnozoites = sampleNHypnozoites(rng);
    
    for( int i = 0; i < numberHypnozoites; ){
        SimTime randomReleaseDelay = sampleReleaseDelay(rng);
        SimTime timeToRelease = ctx.clock.nowOrTs0() + latentP + randomReleaseDelay;
        bool inserted = releases.insert( timeToRelease ).second;
        if( inserted ) ++i;     // successful
        // else: sample clash with an existing release date, so resample
    }

    // Copy times to the vector, backwards (smallest last):
    releaseDates.insert( releaseDates.end(), releases.rbegin(), releases.rend() );
    
#ifdef WHVivaxSamples
    if( sampleHost == host && sampleBrood == 0 ){
        sampleBrood = this;
        cout << "New sample brood";
        for( auto it = releaseDates.begin(); it != releaseDates.end(); ++it )
            cout << '\t' << *it;
        cout << endl;
    }
#endif
}
VivaxBrood::~VivaxBrood(){
#ifdef WHVivaxSamples
    if( sampleBrood == this ){
        sampleBrood = 0;
        cout << "Brood terminated" << endl;
    }
#endif
}

void VivaxBrood::checkpoint( ostream& stream ){
    releaseDates & stream;
    bloodStageClearDate & stream;
    primaryHasStarted & stream;
    relapseHasStarted & stream;
	relapsebHasStarted & stream;
    hadEvent & stream;
    hadRelapse & stream;
	origin & stream;
}
VivaxBrood::VivaxBrood( RunContext& ctx, istream& stream ) : ctx(ctx) {
    (void)ctx;
    releaseDates & stream;
    bloodStageClearDate & stream;
    primaryHasStarted & stream;
    relapseHasStarted & stream;
	relapsebHasStarted & stream;
    hadEvent & stream;
    hadRelapse & stream;
	origin & stream;
}


VivaxBrood::UpdResult VivaxBrood::update(LocalRng& rng){
    if( bloodStageClearDate == ctx.clock.ts0() ){
        //NOTE: this effectively means that both asexual and sexual stage
        // parasites self-terminate. It also means the immune system can
        // protect against new blood-stage vivaxInfections for a short time.
    }
    
    UpdResult result;
    while( releaseDates.size() > 0 && releaseDates.back() == ctx.clock.ts0() ){
        releaseDates.pop_back();
        
#ifdef WHVivaxSamples
        if( sampleBrood == this ){
            cout << "Time\t" << ctx.clock.ts0();
            for( auto it = releaseDates.begin(); it != releaseDates.end(); ++it )
                cout << '\t' << *it;
            cout << endl;
        }
#endif
        
        // an existing or recently terminated blood stage from the same brood
        // protects against a newly released hypnozoite for relapse classification B 
        if( (vivaxClinOption=="B1j" || vivaxClinOption=="B2j") && (bloodStageClearDate + bloodStageProtectionLatency >= ctx.clock.ts0()) ) continue;
        		
        if( !relapsebHasStarted && relapseHasStarted ){
            relapsebHasStarted = true;
            result.newRelapsebBS = true;
        }		
        if( !relapseHasStarted && primaryHasStarted ){
            relapseHasStarted = true;
            result.newRelapseBS = true;
        }
        if( !primaryHasStarted ){
            primaryHasStarted = true;
            result.newPrimaryBS = true;
        }
        result.newBS = true;
        
        double lengthDays = bloodStageLength.sample(rng);
        bloodStageClearDate = ctx.clock.ts0() + ctx.clock.roundToTSFromDays( lengthDays );
        // Assume gametocytes emerge at the same time (they mature quickly and
        // we have little data, thus assume coincedence of start)
    }
    
    result.isFinished = (!isPatent()) && (releaseDates.size() == 0);
    return result;
}

void VivaxBrood::treatmentBS(){
    // Blood stage treatment: clear both asexual and sexual parasites from the
    // blood. NOTE: we assume vivaxInfections removed via treatment do not leave
    // protective immunity since the patient was unable to self-clear.
    bloodStageClearDate = sim::never();
}

void VivaxBrood::treatmentLS(){
    releaseDates.clear();       // 100% clearance
    
    /* partial clearance code, in case of need:
    vector<SimTime> survivingZoites;
    survivingZoites.reserve( releaseDates.size() );   // maximum size we need
    for( auto it = releaseDates.begin(); it != releaseDates.end(); ++it ){
        // parity: order — local RNG stream: draw call order fixes stream state
        if( !rng.bernoulli( pClearEachHypnozoite ) ){
            survivingZoites.push_back( *it );    // copy            
        }
    }
    swap( releaseDates, survivingZoites );
    */
}


// ———  per-host code  ———

void WHInterface::ctorVivax( LocalRng& rng, double comorbidityFactor )
{
    if( comorbidityFactor != 1.0 ){
        throw TRACED_EXCEPTION_DEFAULT( "This vivax model cannot be used with comorbidity heterogeneity" );
    }
#ifdef WHVivaxSamples
    if( sampleHost == 0 ){
        sampleHost = this;
        cout << "New host" << endl;
    }
#endif
    // parity: order — local RNG stream: draw call order fixes stream state
    noPQ = ( pHetNoPQ > 0.0 && rng.bernoulli(pHetNoPQ) );
}

void WHInterface::dtorVivax(){
#ifdef WHVivaxSamples
    if( this == sampleHost ){
        sampleHost = 0;
        cout << "Host terminates" << endl;
    }
#endif
}



double WHInterface::probTransmissionToMosquito_vivax(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l)const{
    assert( ctx.genotypes.N() == 1 );
    for(auto inf = vivaxInfections.begin();
         inf != vivaxInfections.end(); ++inf)
    {
        if( inf->isPatent() ){
            // we have gametocytes from at least one brood
            return probBloodStageInfectiousToMosq;
        }
    }
    return 0;   // no gametocytes
}


bool WHInterface::summarize_vivax(Host::Human& human) const{
    if( vivaxInfections.size() == 0 ) return false;  // no vivaxInfections: not patent, nothing to report
    mon::recordStat(mon::measure("nInfect"), human);
    bool patentHost = false;
    // (patent) vivaxInfections are reported by genotype, even though we don't have
    // genotype in this model
    mon::recordStat(mon::measure("totalInfs"), human, static_cast<int>(vivaxInfections.size()));
    for(auto inf = vivaxInfections.begin(); inf != vivaxInfections.end(); ++inf) {
        if (inf->isPatent()){
            mon::recordStat(mon::measure("totalPatentInf"), human);
            patentHost = true;
        }
    }
    if( patentHost ) mon::recordStat(mon::measure("nPatent"), human);
    return patentHost;
}

void WHInterface::importInfection_vivax(LocalRng& rng){
    // this means one new liver stage infection, which can result in multiple blood stages
    vivaxInfections.push_back( VivaxBrood( ctx, rng, WithinHost::InfectionOrigin::Imported, this ) );
}

void WHInterface::update_vivax(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears)
{
    pSevere = 0.0;
    
    // create new vivaxInfections, letting the constructor do the initialisation work:
    for( int i = 0; i < nNewInfs_i; ++i )
        vivaxInfections.push_back( VivaxBrood( ctx, rng, InfectionOrigin::Introduced, this ) );

    for( int i = 0; i < nNewInfs_l; ++i )
        vivaxInfections.push_back( VivaxBrood( ctx, rng, InfectionOrigin::Indigenous, this ) );

	
    
    // update vivaxInfections
    // NOTE: currently no BSV model
    morbidity = Pathogenesis::NONE;
    uint32_t oldCumInf = cum_prim_inf_ref();
    bool treatmentLiver = treatExpiryLiverVivax > ctx.clock.ts0();
    bool treatmentBlood = treatExpiryBloodVivax > ctx.clock.ts0();
    double matImmClin = 1 - (0.90 * exp(-2.53*ageInYears));
    //double matImmClin = 1 - (1 * exp(-(0.639*8*ageInYears)/2.53));
    auto inf = vivaxInfections.begin();
    while( inf != vivaxInfections.end() ){
        if( treatmentLiver ) inf->treatmentLS();
        if( treatmentBlood ) inf->treatmentBS();        // clearance due to treatment; no protection against reemergence
        VivaxBrood::UpdResult result = inf->update(rng);
        if( result.newPrimaryBS ) cum_prim_inf_ref() += 1;
        
        if( result.newBS ){
            // Sample for each new blood stage infection: the chance of some clinical event.
            // model variant: no illness from relapses possible unless there was illness from the primary infection
            
            bool clinicalEvent = false;
            if( result.newPrimaryBS ){
                // Blood stage is primary. oldCumInf wasn't updated yet.
                //double pPrimaryInfEvent = matImmClin * pPrimaryA * pPrimaryB / (pPrimaryB+oldCumInf);
				double pPrimaryInfEvent = matImmClin * pPrimaryA * exp(-pPrimaryB * oldCumInf);
                // parity: order — local RNG stream: draw call order fixes stream state
                clinicalEvent = rng.bernoulli( pPrimaryInfEvent );
                inf->hadEvent = clinicalEvent;
					
			} else if ( result.newRelapseBS ){
                    // Blood stage is a relapse. oldCumInf wasn't updated yet. Subtract 1 from oldCumInf not
                    // to count the current brood in the number of cumulative primary vivaxInfections 
					
                    double pFirstRelapseEvent = matImmClin * pRelapseOneA * exp(-pRelapseOneB * (oldCumInf-1));
                    // parity: order — local RNG stream: draw call order fixes stream state
                    clinicalEvent = rng.bernoulli( pFirstRelapseEvent );     
                    inf->hadRelapse = clinicalEvent;
      
            } else if ( result.newRelapsebBS ){
                             
				    if (vivaxClinOption=="A1j" || vivaxClinOption=="A2j"){			 
                       double pFirstRelapseEvent = matImmClin * pRelapseOneA * exp(-pRelapseOneB * (oldCumInf-1));
                       // parity: order — local RNG stream: draw call order fixes stream state
                       clinicalEvent = rng.bernoulli( pFirstRelapseEvent );     
                       inf->hadRelapse = clinicalEvent;
				    }
				    if (vivaxClinOption=="B1j" || vivaxClinOption=="B2j"){			 
                       double pSecondRelapseEvent = matImmClin * pRelapseTwoA * exp(-pRelapseTwoB * (oldCumInf-1));
                       // parity: order — local RNG stream: draw call order fixes stream state
                       clinicalEvent = rng.bernoulli( pSecondRelapseEvent );
					   inf->hadRelapse = clinicalEvent;
			        }
				                    
            } else {
                     double pSecondRelapseEvent = matImmClin * pRelapseTwoA * exp(-pRelapseTwoB * (oldCumInf-1));
                     // parity: order — local RNG stream: draw call order fixes stream state
                     clinicalEvent = rng.bernoulli( pSecondRelapseEvent );
                    }
             

            if( clinicalEvent ){
                pSevere = pSevere + (1.0 - pSevere) * pEventIsSevere;
                // parity: order — local RNG stream: draw call order fixes stream state
                if( rng.bernoulli( pEventIsSevere ) )
                    morbidity = static_cast<Pathogenesis::State>( morbidity | Pathogenesis::STATE_SEVERE );
                else
                    morbidity = static_cast<Pathogenesis::State>( morbidity | Pathogenesis::STATE_MALARIA );
            }
        }
        
        if( result.isFinished ) inf = vivaxInfections.erase( inf );
        else ++inf;
    }
    
    
    
    //NOTE: currently we don't model co-infection or indirect deaths
    if( morbidity == Pathogenesis::NONE ){
        morbidity = Pathogenesis::PathogenesisModel::sampleNMF( rng, ageInYears );
    }
}

bool WHInterface::diagnosticResult_vivax( LocalRng& rng, const Diagnostic& diagnostic ) const{
    //TODO(monitoring): this shouldn't ignore the diagnostic (especially since
    // it should always return true if diagnostic.density=0)
    for(auto inf = vivaxInfections.begin(); inf != vivaxInfections.end(); ++inf) {
        if (inf->isPatent())
            return true;        // at least one patent infection
    }
    return false;
}

Pathogenesis::StatePair WHInterface::determineMorbidity_vivax( Host::Human& human, double ageYears, bool ){
    mon::recordStat(mon::measure("expectedSevere"), human, pSevere);
    Pathogenesis::StatePair result;     // no indirect mortality in the vivax model
    result.state = morbidity;
    return result;
}

void WHInterface::clearImmunity_vivax(){
    throw TRACED_EXCEPTION_DEFAULT( "vivax model does not include immune suppression" );
}

void WHInterface::optionalPqTreatment_vivax( Host::Human& human ){
    // PQ clears liver stages. We don't worry about the effect of PQ on
    // gametocytes, because these are always cleared by blood-stage drugs with
    // Vivax, and PQ is not given without BS drugs. NOTE: this ignores drug failure.
    // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
    if (pReceivePQ > 0.0 && (ignoreNoPQ || !noPQ) && human.rng().bernoulli(pReceivePQ)){
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if( human.rng().bernoulli(effectivenessPQ) ){
            for( auto it = vivaxInfections.begin(); it != vivaxInfections.end(); ++it ){
                it->treatmentLS();
            }
        }
        mon::recordEvent(mon::measure("nLiverStageTreatments"), human);
    }
}
bool WHInterface::treatSimple_vivax( Host::Human& human, SimTime timeLiver, SimTime timeBlood ){
    //TODO: this should be implemented properly (allowing effects on next
    // update instead of now)
    
    // liver-stage treatment is only via "LiverStageDrug" option, if at all
    if( timeLiver != sim::zero() ){
        if( pReceivePQ > 0.0 ){
            // This is only really disallowed to prevent simultaneous treatment through both methods
            throw util::xml_scenario_error("simple treatment for vivax liver "
            "stages is incompatible with case-management pUseUncomplicated "
            "(liverStageDrug) option; it is suggested to use the former over the latter");
        }
        // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
        if( (ignoreNoPQ || !noPQ) && (effectivenessPQ == 1.0 || human.rng().bernoulli(effectivenessPQ)) ){
            if( timeLiver >= sim::zero() ){
                treatExpiryLiverVivax = max( int(treatExpiryLiverVivax), ctx.clock.nowOrTs1() + timeLiver );
            }else{
                for( auto it = vivaxInfections.begin(); it != vivaxInfections.end(); ++it ){
                    it->treatmentLS();
                }
            }
        }
        mon::recordEvent(mon::measure("nLiverStageTreatments"), human);
    }
    
    // there probably will be blood-stage treatment
    if( timeBlood != sim::zero() ){
        if( timeBlood < sim::zero() ){
            // legacy mode: retroactive clearance
            for( auto it = vivaxInfections.begin(); it != vivaxInfections.end(); ++it ){
                it->treatmentBS();
            }
        }else{
            treatExpiryBloodVivax = max( int(treatExpiryBloodVivax), ctx.clock.nowOrTs1() + timeBlood );
        }
        return true;    // blood stage treatment
    }
    return false;    // no blood stage treatment
}


// ———  boring stuff: checkpointing and set-up  ———

void WHInterface::checkpoint_vivax(istream& stream){
    size_t len;
    len & stream;
    for( size_t i = 0; i < len; ++i ){
        vivaxInfections.push_back( VivaxBrood( ctx, stream ) );
    }
    noPQ & stream;
    int morbidity_i;
    morbidity_i & stream;
    morbidity = static_cast<Pathogenesis::State>( morbidity_i );
    cum_prim_inf_ref() & stream;
    treatExpiryLiverVivax & stream;
    treatExpiryBloodVivax & stream;
    pEvent & stream;
    pFirstRelapseEvent & stream;
}
void WHInterface::checkpoint_vivax(ostream& stream){
    vivaxInfections.size() & stream;
    for( auto it = vivaxInfections.begin(); it != vivaxInfections.end(); ++it ){
        it->checkpoint( stream );
    }
    noPQ & stream;
    static_cast<int>( morbidity ) & stream;
    cum_prim_inf_ref() & stream;
    treatExpiryLiverVivax & stream;
    treatExpiryBloodVivax & stream;
    pEvent & stream;
    pFirstRelapseEvent & stream;
}

char const*const not_impl = "feature not available in Vivax model";
void WHInterface::treatPkPd_vivax(size_t schedule, size_t dosages, double age, double delay_d){
    throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHInterface::getTotalDensity_vivax() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHInterface::getCumulative_h_vivax() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHInterface::getCumulative_Y_vivax() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }


}
}
