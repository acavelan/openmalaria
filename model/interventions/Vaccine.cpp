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
#include "core/RunExtensions.h"
#include "interventions/Vaccine.h"
#include "interventions/HumanComponents.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include "util/checkpoint.h"
#include "util/checkpoint_containers.h"
#include "Host/WithinHost/Genotypes.h"

#include <limits>
#include <cmath>

namespace OM {
namespace interventions {
using namespace OM::util;


const VaccineParams& VaccineParams::getParams( RunContext& ctx, ComponentId component ){
    auto& params = ctx.ext->vaccineComponents;
    assert( component.id < params.size() && params[component.id] != 0 );
    return *params[component.id];
}

double VaccineParams::getInitialEfficacy (LocalRng& rng, size_t numPrevDoses, uint32_t genotype) const
{
    /* If initialMeanEfficacy.size or more doses have already been given, use
     * the last efficacy. */
    if (numPrevDoses >= perGenotypeInitialMeanEfficacy[genotype].size())
        numPrevDoses = perGenotypeInitialMeanEfficacy[genotype].size() - 1;
    double ime = perGenotypeInitialMeanEfficacy[genotype][numPrevDoses];
    //NOTE(validation): With extra valiadation in random, the first difference is noticed here:
    util::streamValidate(ime);
    util::streamValidate(efficacyB);
    if (ime == 0.0){
        return 0.0;
    } else if (ime < 1.0) {
        // parity: order — local RNG stream: draw call order fixes stream state
        double result = rng.betaWithMean (ime, efficacyB);
        //NOTE(validation):: Without extra validation in random, the first difference is noticed here:
        util::streamValidate(result);
        //TODO(validation):: Why the difference? Bug in/limitation of StreamValidatior?
        // Extra memory allocation due to the extra logging causes some bad
        // memory usage to manifest differently? Am I forgetting to initialise
        // some variable?
        return result;
    } else {
        return 1.0;
    }
}

#if 0
/* R_0: this used to use _the_ vaccine configuration; this is just a check.
 * If reimplementing, use a new vaccine instance; implement some hook to ensure
 * that TBV is given to everyone after enabling the intervention, etc. */
void Vaccine::verifyEnabledForR_0 (){
    if( !types[PEV].active || !types[TBV].active )
        throw util::xml_scenario_error("PEV and TBV vaccines must have a "
                "description to use the insertR_0Case intervention");
}
#endif

// ———  per-human vaccine state: VacRecord column (old PerHumanVaccine)  ———

double vaccine_factor( RunContext& ctx, const std::vector<VacRecord>& vac,
                       Vaccine::Types type, uint32_t genotype )
{
    // parity: c — records multiply in insertion (deployment) order
    double factor = 1.0;
    for( const VacRecord& r : vac ){
        if( VaccineParams::getParams(ctx, r.id).type == type ){
            SimTime age = ctx.clock.ts1() - r.time_last;  // implies age 1 TS on first use
            double decayFactor = r.decay->eval_with( age, r.decay_factor );
            factor *= 1.0 - r.per_genotype_eff[genotype] * decayFactor;
        }
    }
    return factor;
}

bool vaccinate( RunContext& ctx, Host::Human& human,
                ComponentId componentId, interventions::VaccineLimits vaccLimits )
{
    auto& vac = human.cols->vac[human.row];
    VacRecord* record = 0;
    for( VacRecord& r : vac ){
        if( r.id == componentId ){
            record = &r;
            break;
        }
    }
    
    uint32_t numDosesAdministered = (record == 0) ? 0 : record->num_doses;
    if( numDosesAdministered < vaccLimits.minPrevDoses ||
        numDosesAdministered >= vaccLimits.maxCumDoses )
        return false;   // no vaccination (this replaces the old schedule for continuous doses)
    
    const VaccineParams& params = VaccineParams::getParams( ctx, componentId );
    
    if( record == 0 ){
        // parity: draw order — first dose per component samples the decay het
        // BEFORE the per-genotype efficacy draws (old PerEffectPerHumanVaccine ctor)
        VacRecord r;
        r.id = componentId;
        r.decay = params.decayFunc.get();
        r.decay_factor = params.decayFunc->sample_het_factor( human.rng() );
        vac.push_back( r );
        record = &vac.back();
    }
    
    record->per_genotype_eff.resize(ctx.genotypes.N());

    for(size_t g=0; g<params.perGenotypeInitialMeanEfficacy.size(); g++)
        record->per_genotype_eff[g] = params.getInitialEfficacy(human.rng(), numDosesAdministered, g);

    util::streamValidate(record->per_genotype_eff);
    
    record->num_doses = numDosesAdministered + 1;
    record->time_last = ctx.clock.nowOrTs1();
    
    return true;
}

}
}
