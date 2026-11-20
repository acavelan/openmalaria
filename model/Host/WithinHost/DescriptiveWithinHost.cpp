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

#include "core/RunExtensions.h"
#include "Global.h"
#include "core/RunContext.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include "util/errors.h"

#include <cassert>

using namespace std;

namespace OM {
namespace WithinHost {

extern bool bugfix_max_dens;    // DescriptiveInfection.cpp
bool reportPatentInfected = false;
// -----  Initialization  -----

void WHInterface::initDescriptive(RunContext& ctx){
    reportPatentInfected = (
                                mon::isUsed(ctx, mon::measure("totalPatentInf"))
                                || mon::isUsed(ctx, mon::measure("totalPatentInf_Imported"))
                                || mon::isUsed(ctx, mon::measure("totalPatentInf_Introduced"))
                                || mon::isUsed(ctx, mon::measure("totalPatentInf_Indigenous"))
                            );
}

void WHInterface::ctorDescriptive( LocalRng& rng, double comorbidityFactor )
{
    ctorFalciparum( rng, comorbidityFactor );
    
    assert( ctx.clock.one_TS() == sim::fromDays(5) );
    opt_vaccine_genotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);
}

void WHInterface::dtorDescriptive() {
    for( auto inf = infections.begin(); inf != infections.end(); ++inf ){
        delete *inf;
    }
    infections.clear();
}


// -----  Simple infection adders/removers  -----

void WHInterface::loadInfection_descriptive(istream& stream) {
    infections.push_back(infection_new_checkpointed(ctx, stream, Infection::Kind::Descriptive));
}

// clearInfections() and clearImmunity_falciparum() are shared with the Common
// kind (identical code, defined in CommonWithinHost.cpp).

void WHInterface::importInfection_descriptive(LocalRng& rng){
    if( SlotAllocator::try_import(numInfs) ){   // mode A: at the cap, no draw at all
        imm_h_ref() += 1;
        numInfs += 1;
        // This is a hook, used by interventions. The newly imported infections
        // should use initial frequencies to select genotypes.
        vector<double> weights( 0 );        // zero length: signal to use initial frequencies
        uint32_t genotype = ctx.genotypes.sampleGenotype(rng, weights);
        infections.push_back(infection_new(ctx, rng, genotype, InfectionOrigin::Imported, Infection::Kind::Descriptive));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void WHInterface::update_descriptive(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears)
{
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    int nNewInfsToBeCreated_i = nNewInfs_i;
    int nNewInfsToBeCreated_l = nNewInfs_l;

    SlotAllocator::clamp_new( nNewInfs_i, nNewInfs_l, numInfs );   // mode B

    numInfs += nNewInfs_i;
    assert( numInfs>=0 && numInfs<=SlotAllocator::MAX_INFECTIONS );
    for( int i=0; i<nNewInfs_i; ++i ) {
        uint32_t genotype = ctx.genotypes.sampleGenotype(rng, genotype_weights_i);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( opt_vaccine_genotype )
        {
            double vaccineFactor = human.vaccine_factor( interventions::Vaccine::PEV, genotype );
            // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
            if(vaccineFactor == 1.0 || human.rng().bernoulli(vaccineFactor))
                infections.push_back(infection_new (ctx, rng, genotype, InfectionOrigin::Introduced, Infection::Kind::Descriptive));
        }
        else if (opt_vaccine_genotype == false)
            infections.push_back(infection_new (ctx, rng, genotype, InfectionOrigin::Introduced, Infection::Kind::Descriptive));
    }
    assert( numInfs == static_cast<int>(infections.size()) );

    numInfs += nNewInfs_l;
    assert( numInfs>=0 && numInfs<=SlotAllocator::MAX_INFECTIONS );
    for( int i=0; i<nNewInfs_l; ++i ) {
        uint32_t genotype = ctx.genotypes.sampleGenotype(rng, genotype_weights_l);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( opt_vaccine_genotype )
        {
            double vaccineFactor = human.vaccine_factor( interventions::Vaccine::PEV, genotype );
            // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
            if(vaccineFactor == 1.0 || human.rng().bernoulli(vaccineFactor))
                infections.push_back(infection_new (ctx, rng, genotype, InfectionOrigin::Indigenous, Infection::Kind::Descriptive));
        }
        else if (opt_vaccine_genotype == false)
            infections.push_back(infection_new (ctx, rng, genotype, InfectionOrigin::Indigenous, Infection::Kind::Descriptive));
    }
    assert( numInfs == static_cast<int>(infections.size()) );

    updateImmuneStatus ();

    totalDensity = 0.0;
    hrp2Density = 0.0;
    timeStepMaxDensity = 0.0;

    bool treatmentLiver = treatExpiryLiver > ctx.clock.ts0();
    bool treatmentBlood = treatExpiryBlood > ctx.clock.ts0();

    // p7b B4 hoist: vaccine_factor is pure over the host's (frozen during host
    // update) vac column and returns exactly 1.0 when it is empty — same value,
    // one check per host instead of one call per infection.
    const bool anyVac = human.has_vaccines();
    // p7c-C3 hoist: dH/dA/exp of immunitySurvivalFactor are host-tick constants
    // (updateImmuneStatus above is this tick's last mutation of imm_h/imm_Y;
    // ageInYears is fixed) — one computation per host tick, same IEEE values
    // the per-infection call used. Computed on first use so a NaN still throws
    // at exactly the first infection that reaches the call (original point).
    WHInterface::ImmSurvTick immTick = {};
    bool immTickReady = false;
    for(auto inf = infections.begin(); inf != infections.end();) {
        //NOTE: it would be nice to combine this code with that in
        // CommonWithinHost.cpp, but a few changes would be needed:
        // INNATE_MAX_DENS and MAX_DENS_CORRECTION would need to be required
        // (couldn't support old parameterisations using buggy versions of code
        // any more).
        // SP drug action and the PK/PD model would need to be abstracted
        // behind a common interface.
        if ( (*inf)->expired() /* infection has self-terminated */ ||
            ((*inf)->bloodStage() ? treatmentBlood : treatmentLiver) )
        {
            delete *inf;
            inf=infections.erase(inf);
            numInfs--;
            continue;
        }
        
        // Should be: infStepMaxDens = 0.0, but has some history.
        // See MAX_DENS_CORRECTION in DescriptiveInfection.cpp.
        double infStepMaxDens = timeStepMaxDensity;
        if( !immTickReady ) { immTick = immSurvivalTick(ageInYears); immTickReady = true; }
        double immSurvFact = immunitySurvivalFactor(immTick, (*inf)->cumulativeExposureJ());
        double bsvFactor = anyVac ? human.vaccine_factor(interventions::Vaccine::BSV, opt_vaccine_genotype? (*inf)->genotype() : 0) : 1.0;

        determine_densities_descriptive(**inf, rng, imm_h_ref(), infStepMaxDens, immSurvFact, _innateImmSurvFact, bsvFactor);

        if (bugfix_max_dens)
            infStepMaxDens = std::max(infStepMaxDens, timeStepMaxDensity);
        timeStepMaxDensity = infStepMaxDens;

        double density = (*inf)->getDensity();
        totalDensity += density;
        if( !(*inf)->isHrp2Deficient() ){
            hrp2Density += density;
        }

        ++inf;
    }
    
    // As in AJTMH p22, cumulative_h (X_h + 1) doesn't include infections added
    // this time-step and cumulative_Y only includes past densities.
    imm_h_ref() += nNewInfs_i + nNewInfs_l;
    imm_Y_ref() += ctx.clock.one_TS() * totalDensity;
    
    util::streamValidate( totalDensity );
    util::streamValidate( hrp2Density );
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
    if(opt_vaccine_genotype == false)
    {
        nNewInfs_i = nNewInfsToBeCreated_i;
        nNewInfs_l = nNewInfsToBeCreated_l;
    }
}

// -----  Summarize  -----

bool WHInterface::summarize_descriptive( Host::Human& human )const{
    pathogenesisModel->summarize( human );
    
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

        int nImported = 0, nIntroduced = 0, nIndigenous = 0;
        for( auto inf = infections.begin(); inf != infections.end(); ++inf )
        {
            if((*inf)->origin() == InfectionOrigin::Indigenous) nIndigenous++;
            else if((*inf)->origin() == InfectionOrigin::Introduced) nIntroduced++;
            else nImported++;
        }

        // (patent) infections are reported by genotype, even though we don't have
        // genotype in this model
        mon::recordStat(mon::measure("totalInfs"), human, static_cast<int>(infections.size()));
        mon::recordStat(mon::measure("totalInfs_Imported"), human, nImported);
        mon::recordStat(mon::measure("totalInfs_Introduced"), human, nIntroduced);
        mon::recordStat(mon::measure("totalInfs_Indigenous"), human, nIndigenous);

        if( reportPatentInfected ){
            for(auto inf = infections.begin(); inf != infections.end(); ++inf)
            {
                if( diagnostics::monitoringDiagnostic( ctx ).isPositive( human.rng(), (*inf)->getDensity(), std::numeric_limits<double>::quiet_NaN() ) )
                {
                    mon::recordStat(mon::measure("totalPatentInf"), human);
                    if((*inf)->origin() == InfectionOrigin::Indigenous)
                        mon::recordStat(mon::measure("totalPatentInf_Indigenous"), human);
                    else if((*inf)->origin() == InfectionOrigin::Introduced)
                        mon::recordStat(mon::measure("totalPatentInf_Introduced"), human);
                    else
                        mon::recordStat(mon::measure("totalPatentInf_Imported"), human);
                }
            }
        }
        if( reportInfectionsByGenotype ){
            // accumulate total density by genotype
            map<uint32_t, double> dens_by_gtype;
            for(auto inf = infections.begin(); inf != infections.end(); ++inf)
                dens_by_gtype[(*inf)->genotype()] += (*inf)->getDensity();
            
            for( auto gtype: dens_by_gtype ){
                // we had at least one infection of this genotype
                mon::recordStat(mon::measure("nInfectByGenotype"), human, 1, 0, gtype.first, 0);
                if( diagnostics::monitoringDiagnostic( ctx ).isPositive(human.rng(), gtype.second, std::numeric_limits<double>::quiet_NaN()) ){
                    mon::recordStat(mon::measure("nPatentByGenotype"), human, 1, 0, gtype.first, 0);
                    mon::recordStat(mon::measure("logDensByGenotype"), human, log(gtype.second), 0, gtype.first, 0);
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

        if(totalDensity > 1e-10)
            mon::recordStat(mon::measure("sumlogDens"), human, log(totalDensity));
        return true;    // patent
    }
    return false;       // not patent
}


// -----  Data checkpointing  -----

void WHInterface::checkpoint_descriptive (istream& stream) {
    checkpoint_falciparum (stream);
    for(int i=0; i<numInfs; ++i) {
        loadInfection_descriptive(stream);  // create infections from the stream
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}
void WHInterface::checkpoint_descriptive (ostream& stream) {
    checkpoint_falciparum (stream);
    for(Infection* inf : infections)
        *inf & stream;
}

char const*const not_impl = "feature not available with the \"descriptive\" within-host model";
void WHInterface::treatPkPd_descriptive(size_t schedule, size_t dosages, double age, double delay_d){
    throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }

}
}
