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

#include "Host/WithinHost/Genotypes.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/vectors.h"
#include "util/CommandLine.h"

#include <iomanip>

namespace OM {
namespace WithinHost {

namespace GT /*for genotype impl details*/{
// ———  Model constants (after init)  ———
// keys are cumulative probabilities; last entry should equal 1; values are genotype codes
map<double,uint32_t> cum_initial_freqs;

// we give each allele of each loci a unique code
map<string, map<string, uint32_t> > alleleCodes;
uint32_t nextAlleleCode = 0;


vector<Genotypes::Genotype> genotypes;
// ———  Model variables  ———
enum SampleMode{
    SAMPLE_FIRST,      // always choose first genotype (essentially the off switch)
    SAMPLE_INITIAL,    // sample from initial probabilities
    SAMPLE_TRACKING    // sample from tracked success at genotype level (no recombination)
};
// Mode to use now (until switched) and from the start of the intervention period.
SampleMode current_mode = SAMPLE_FIRST, interv_mode = SAMPLE_FIRST;
}
size_t Genotypes::N_genotypes = 1;

Genotypes::Genotype Genotypes::Genotype::cross( const Genotype& that )const{
    Genotype result(*that.alleles.begin(),
                    init_freq*that.init_freq,
                    fitness*that.fitness,
                    hrp2_deficient || that.hrp2_deficient);
    result.alleles.insert( alleles.begin(), alleles.end() );
    result.alleles.insert( that.alleles.begin(), that.alleles.end() );
    return result;
}

void Genotypes::initSingle()
{
    // no specification implies there is a single genotype
    GT::genotypes.assign( 1, Genotypes::Genotype(
        0 /*allele code*/, 1.0/*frequency*/, 1.0/*fitness*/, false /*hrp2 deficiency*/) );
    N_genotypes = 1;
}



void Genotypes::preMainSimInit(){
    GT::current_mode = GT::interv_mode;
}

uint32_t Genotypes::findAlleleCode(const string& locus, const string& allele){
    auto it = GT::alleleCodes.find(locus);
    if( it == GT::alleleCodes.end() ) return numeric_limits<uint32_t>::max();
    auto it2 = it->second.find( allele );
    if( it2 == it->second.end() ) return numeric_limits<uint32_t>::max();
    return it2->second;
}

const vector< Genotypes::Genotype >& Genotypes::getGenotypes(){
    return GT::genotypes;
}

uint32_t Genotypes::sampleGenotype( LocalRng& rng, vector<double>& genotype_weights ){
    if( GT::current_mode == GT::SAMPLE_FIRST ){
        return 0;       // always the first genotype code
    }else if( GT::current_mode == GT::SAMPLE_INITIAL
            || genotype_weights.size() == 0 )
    {
        // parity: order — local RNG stream: draw call order fixes stream state
        double sample = rng.uniform_01();
        auto it = GT::cum_initial_freqs.upper_bound( sample );
        assert( it != GT::cum_initial_freqs.end() );
        return it->second;
    }else{
        assert( GT::current_mode == GT::SAMPLE_TRACKING );
        assert( genotype_weights.size() == N_genotypes );
        double weight_sum = util::vectors::sum( genotype_weights );
        assert( weight_sum >= 0.0 && weight_sum < 1e5 );        // possible loss of precision or other error
        double sample = rng.uniform_01() * weight_sum;
        double cum = 0.0;
        for( size_t g = 0; g < N_genotypes; ++g ){
            cum += genotype_weights[g];
            if( sample < cum ) return g;
        }
        return 0;       // just to be safe (could happen if weight_sum == 0.0)
    }
}

double Genotypes::initialFreq( size_t genotype ){
    if( GT::genotypes.size() == 0 ){
        assert( genotype == 0 );
        assert( N_genotypes == 1 );
        return 1.0;     // only 1 genotype, thus initial frequency is 100%
    }else{
        return GT::genotypes[genotype].init_freq;
    }
}


// ———  checkpointing  ———

void Genotypes::staticCheckpoint( ostream& stream ){
    int t = GT::current_mode;
    t & stream;
}
void Genotypes::staticCheckpoint( istream& stream ){
    int t;
    t & stream;
    GT::current_mode = static_cast<GT::SampleMode>(t);
    assert( GT::current_mode == GT::SAMPLE_FIRST ||
        GT::current_mode == GT::SAMPLE_INITIAL ||
        GT::current_mode == GT::SAMPLE_TRACKING );
}

}
}
