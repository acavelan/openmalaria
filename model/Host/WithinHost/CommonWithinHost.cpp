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

#include "core/RunContext.h"
#include "util/checkpoint.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "core/RunExtensions.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/AgeGroupInterpolation.h"
#include "util/random.h"
#include "util/StreamValidator.h"

using namespace std;

namespace OM {
namespace WithinHost {


double hetMassMultStdDev = std::numeric_limits<double>::signaling_NaN();
double minHetMassMult = std::numeric_limits<double>::signaling_NaN();
util::AgeGroupInterpolator massByAge;

bool reportInfectedOrPatentInfected = false;

// -----  Initialization  -----


void WHInterface::ctorCommon( LocalRng& rng, double comorbidityFactor )
{
    ctorFalciparum( rng, comorbidityFactor );
    
    assert( ctx.clock.one_TS() == sim::fromDays(1) || ctx.clock.one_TS() == sim::fromDays(5) );
    
    opt_vaccine_genotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    // Sample a weight heterogeneity factor
#ifndef NDEBUG
    int counter = 0;
#endif
    do {
        // parity: order — local RNG stream: draw call order fixes stream state
        hetMassMultiplier = rng.gauss( 1.0, hetMassMultStdDev );
#ifndef NDEBUG
        assert( counter < 100 );        // too many resamples: resamples should rarely be needed...
        ++counter;
#endif
    } while( hetMassMultiplier < minHetMassMult );
}

void WHInterface::dtorCommon() {
    for( auto inf = infections.begin(); inf != infections.end(); ++inf ){
        delete *inf;
    }
    infections.clear();
}

// -----  Simple infection adders/removers  -----

void WHInterface::clearInfections( Treatments::Stages stage ){
    for(auto inf = infections.begin(); inf != infections.end();) {
        if( stage == Treatments::BOTH ||
            (stage == Treatments::LIVER && !(*inf)->bloodStage()) ||
            (stage == Treatments::BLOOD && (*inf)->bloodStage())
        ){
            delete *inf;
            inf = infections.erase( inf );
        }else{
            ++inf;
        }
    }
    numInfs = infections.size();
}

// -----  interventions -----

void WHInterface::treatPkPd_common(size_t schedule, size_t dosage, double age, double delay_d){
    double mass = massByAge.eval( age ) * hetMassMultiplier;
    pkpdModel.prescribe( ctx, schedule, dosage, age, mass, delay_d );
}
void WHInterface::clearImmunity_falciparum() {
    for(auto inf = infections.begin(); inf != infections.end(); ++inf) {
        (*inf)->clearImmunity();
    }
    imm_h_ref() = 0.0;
    imm_Y_lag_ref() = 0.0;
}
void WHInterface::importInfection_common(LocalRng& rng){
    if( SlotAllocator::try_import(numInfs) ){   // mode A: at the cap, no draw at all
        imm_h_ref() += 1;
        numInfs += 1;
        // This is a hook, used by interventions. The newly imported infections
        // should use initial frequencies to select genotypes.
        vector<double> weights( 0 );        // zero length: signal to use initial frequencies
        uint32_t genotype = ctx.genotypes.sampleGenotype(rng, weights);
        infections.push_back(ctx.ext->createInfection(ctx, rng, genotype, WithinHost::InfectionOrigin::Imported));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void WHInterface::update_common(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears)
{
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    SlotAllocator::clamp_new( nNewInfs_i, nNewInfs_l, numInfs );   // mode B

    int nNewInfsIgnored = nNewInfs_i + nNewInfs_l - SlotAllocator::free_slots(numInfs);

    assert( numInfs>=0 && numInfs<=SlotAllocator::MAX_INFECTIONS );

    int nNewInfsDiscarded = 0;
    for( int i=0; i<nNewInfs_i; ++i ) {
        uint32_t genotype = ctx.genotypes.sampleGenotype(rng, genotype_weights_i);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( opt_vaccine_genotype )
        {
            double vaccineFactor = human.vaccine_factor( interventions::Vaccine::PEV, genotype );
            // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
            if(vaccineFactor == 1.0 || human.rng().bernoulli(vaccineFactor))
                infections.push_back(ctx.ext->createInfection (ctx, rng, genotype, InfectionOrigin::Introduced));
            else
                nNewInfsDiscarded++;
        }
        else if (opt_vaccine_genotype == false)
            infections.push_back(ctx.ext->createInfection (ctx, rng, genotype, InfectionOrigin::Introduced));
    }
    // Update nNewInfs, this is the number that will be reported in Human
    nNewInfs_i -= nNewInfsDiscarded;
    numInfs += nNewInfs_i;

    nNewInfsDiscarded = 0;
    for( int i=0; i<nNewInfs_l; ++i ) {
        uint32_t genotype = ctx.genotypes.sampleGenotype(rng, genotype_weights_l);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( opt_vaccine_genotype )
        {
            double vaccineFactor = human.vaccine_factor( interventions::Vaccine::PEV, genotype );
            // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
            if(vaccineFactor == 1.0 || human.rng().bernoulli(vaccineFactor))
                infections.push_back(ctx.ext->createInfection (ctx, rng, genotype, InfectionOrigin::Indigenous));
            else
                nNewInfsDiscarded++;
        }
        else if (opt_vaccine_genotype == false)
            infections.push_back(ctx.ext->createInfection (ctx, rng, genotype, InfectionOrigin::Indigenous));
    }
    // Update nNewInfs, this is the number that will be reported in Human
    nNewInfs_l -= nNewInfsDiscarded;
    numInfs += nNewInfs_l;

    assert( numInfs == static_cast<int>(infections.size()) );
    
    updateImmuneStatus ();

    totalDensity = 0.0;
    hrp2Density = 0.0;
    timeStepMaxDensity = 0.0;
    
    bool treatmentLiver = treatExpiryLiver > ctx.clock.ts0();
    bool treatmentBlood = treatExpiryBlood > ctx.clock.ts0();
    
    double body_mass = massByAge.eval( ageInYears ) * hetMassMultiplier;
    
    // p7b B4 hoist: pure + empty⇒1.0, so per-infection calls collapse to one
    // check per host (vac column is frozen during the human update).
    const bool anyVac = human.has_vaccines();
    for( SimTime now = ctx.clock.ts0(), end = ctx.clock.ts0() + ctx.clock.one_TS(); now < end; now = now + sim::oneDay() ){
        // every day, medicate drugs, update each infection, then decay drugs
        pkpdModel.medicate(ctx, rng);
        
        double sumLogDens = 0.0;
        
        for(auto inf = infections.begin(); inf != infections.end();) {
            // Note: this is only one treatment model; there is also the PK/PD model
            bool expires = ((*inf)->bloodStage() ? treatmentBlood : treatmentLiver);
            
            if( !expires ){     /* no expiry due to simple treatment model; do update */
                const double drugFactor = pkpdModel.getDrugFactor(rng, *inf, body_mass);
                const double immFactor = immunitySurvivalFactor(ageInYears, (*inf)->cumulativeExposureJ());
                const double bsvFactor = anyVac ? human.vaccine_factor(interventions::Vaccine::BSV, opt_vaccine_genotype? (*inf)->genotype() : 0) : 1.0;
                const double survivalFactor = bsvFactor * _innateImmSurvFact * immFactor * drugFactor;
                // update, may result in termination of infection:
                expires = infection_update(**inf, rng, survivalFactor, now, body_mass);
            }
            
            if( expires ){
                delete *inf;
                inf = infections.erase(inf);        // inf points to next infection now so don't increment with ++inf
                --numInfs;
            } else {
                double density = (*inf)->getDensity();
                totalDensity += density;
                if( !(*inf)->isHrp2Deficient() ){
                    hrp2Density += density;
                }
                timeStepMaxDensity = max(timeStepMaxDensity, density);
                if( density > 0 ){
                    // Base 10 logarithms are usually used; +1 because it avoids negatives in output while having very little affect on high densities
                    sumLogDens += log10(1.0 + density);
                }
                ++inf;
            }
        }
        pkpdModel.decayDrugs (body_mass);
    }
    
    // As in AJTMH p22, cumulative_h (X_h + 1) doesn't include infections added
    // this time-step and cumulative_Y only includes past densities, thus we
    // increment these after the update.
    imm_h_ref() += nNewInfs_i + nNewInfs_l;
    imm_Y_ref() += totalDensity;
    
    util::streamValidate(totalDensity);
    util::streamValidate(hrp2Density);
    assert( (std::isfinite)(totalDensity) );        // inf probably wouldn't be a problem but NaN would be
    
    // Cache total density for infectiousness calculations
    int y_lag_i = ctx.clock.moduloSteps(ctx.clock.ts1(), y_lag_len);
    for( size_t g = 0; g < ctx.genotypes.N(); ++g )
    {
        m_y_lag_i[y_lag_i * ctx.genotypes.N() + g] = 0.0;
        m_y_lag_l[y_lag_i * ctx.genotypes.N() + g] = 0.0;
    }

    for( auto inf = infections.begin(); inf != infections.end(); ++inf )
    {
        if((*inf)->origin() == InfectionOrigin::Imported)
            m_y_lag_i[y_lag_i * ctx.genotypes.N() + (*inf)->genotype()] += (*inf)->getDensity();
        else
            m_y_lag_l[y_lag_i * ctx.genotypes.N() + (*inf)->genotype()] += (*inf)->getDensity();
    }

    // This is a bug, we keep it this way to be consistent with old simulations
    if(nNewInfsIgnored > 0)
        nNewInfs_l += nNewInfsIgnored;
}

// -----  Summarize  -----

// Used in summarizeInfs.
vector<Infection*> sortedInfs;
struct InfGenotypeSorter {
    bool operator() (Infection* i, Infection* j){
        return i->genotype() < j->genotype();
    }
} infGenotypeSorter;

bool WHInterface::summarize_common( Host::Human& human )const{
    pathogenesisModel->summarize( human );
    pkpdModel.summarize( ctx, human );
    
    InfectionOrigin infectionType = get_infection_origin(infections);
    
    // If the number of infections is 0 and parasite density is positive we default to Indigenous
    if( infections.size() > 0 ){
        mon::recordStat(mon::measure("nInfect"), human);
        if(infectionType == InfectionOrigin::Indigenous)
            mon::recordStat(mon::measure("nInfect_Indigenous"), human);
        else if(infectionType == InfectionOrigin::Introduced)
            mon::recordStat(mon::measure("nInfect_Introduced"), human);
        else
            mon::recordStat(mon::measure("nInfect_Imported"), human);

        if( reportInfectedOrPatentInfected ){
            for(auto inf = infections.begin(); inf != infections.end(); ++inf)
            {
                uint32_t genotype = (*inf)->genotype();
                mon::recordStat(mon::measure("totalInfs"), human, 1, 0, genotype, 0);
                if((*inf)->origin() == InfectionOrigin::Indigenous)
                    mon::recordStat(mon::measure("totalInfs_Indigenous"), human, 1, 0, genotype, 0);
                else if((*inf)->origin() == InfectionOrigin::Introduced)
                    mon::recordStat(mon::measure("totalInfs_Introduced"), human, 1, 0, genotype, 0);
                else
                    mon::recordStat(mon::measure("totalInfs_Imported"), human, 1, 0, genotype, 0);

                if( diagnostics::monitoringDiagnostic( ctx ).isPositive( human.rng(), (*inf)->getDensity(), std::numeric_limits<double>::quiet_NaN() ) ){
                    mon::recordStat(mon::measure("totalPatentInf"), human, 1, 0, genotype, 0);
                    if((*inf)->origin() == InfectionOrigin::Indigenous)
                        mon::recordStat(mon::measure("totalPatentInf_Indigenous"), human, 1, 0, genotype, 0);
                    else if((*inf)->origin() == InfectionOrigin::Introduced)
                        mon::recordStat(mon::measure("totalPatentInf_Introduced"), human, 1, 0, genotype, 0);
                    else
                        mon::recordStat(mon::measure("totalPatentInf_Imported"), human, 1, 0, genotype, 0);
                }
            }
        }
        
        if( reportInfectionsByGenotype ){
            // Instead of storing nInfs and total density by genotype we sort
            // infections by genotype and report each in sequence.
            // We don't sort in place since that would affect random number sampling
            // order when updating, and the monitoring system should not in my
            // opinion affect outputs (since it would make testing harder).
            sortedInfs.assign( infections.begin(), infections.end() );
            sort( sortedInfs.begin(), sortedInfs.end(), infGenotypeSorter );
            auto inf = sortedInfs.begin();
            while( inf != sortedInfs.end() ){
                uint32_t genotype = (*inf)->genotype();
                double dens = 0.0;
                do{     // at start: genotype is that of the current infection, dens is 0
                    dens += (*inf)->getDensity();
                    ++inf;
                }while( inf != sortedInfs.end() && (*inf)->genotype() == genotype );
                // we had at least one infection of this genotype
                mon::recordStat(mon::measure("nInfectByGenotype"), human, 1, 0, genotype, 0);
                if( diagnostics::monitoringDiagnostic( ctx ).isPositive(human.rng(), dens, std::numeric_limits<double>::quiet_NaN()) ){
                    mon::recordStat(mon::measure("nPatentByGenotype"), human, 1, 0, genotype, 0);
                    mon::recordStat(mon::measure("logDensByGenotype"), human, log(dens), 0, genotype, 0);
                }
            }
        }
    }
    
    // Some treatments (simpleTreat with steps=-1) clear infections immediately
    // (and are applied after update()), thus infections.size() may be 0 while
    // totalDensity > 0. Here we report the last calculated density.
    if( diagnostics::monitoringDiagnostic( ctx ).isPositive(human.rng(), totalDensity, std::numeric_limits<double>::quiet_NaN()) ){
        mon::recordStat(mon::measure("nPatent"), human);
        if(infectionType == InfectionOrigin::Imported)
            mon::recordStat(mon::measure("nPatent_Imported"), human);
        else if(infectionType == InfectionOrigin::Introduced)
            mon::recordStat(mon::measure("nPatent_Introduced"), human);
        else if(infectionType == InfectionOrigin::Indigenous)
            mon::recordStat(mon::measure("nPatent_Indigenous"), human);

        if(totalDensity == 0.0)
            mon::recordStat(mon::measure("sumlogDens"), human, 0.0);
        else
            mon::recordStat(mon::measure("sumlogDens"), human, log(totalDensity));
        return true;    // patent
    }
    return false;       // not patent
}


void WHInterface::checkpoint_common (istream& stream) {
    checkpoint_falciparum (stream);
    hetMassMultiplier & stream;
    pkpdModel.checkpoint(ctx, stream);
    for(int i = 0; i < numInfs; ++i) {
        infections.push_back (ctx.ext->checkpointedInfection (ctx, stream));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}

void WHInterface::checkpoint_common (ostream& stream) {
    checkpoint_falciparum (stream);
    hetMassMultiplier & stream;
    pkpdModel.checkpoint(stream);
    for(auto inf = infections.begin(); inf != infections.end(); ++inf) {
        (**inf) & stream;
    }
}
}
}
