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

#include "Host/WithinHost/Diagnostic.h"
#include "core/RunExtensions.h"
#include "Parameters.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include <limits>

namespace OM { namespace WithinHost {

// ———  Diagnostic (non-static)  ———


Diagnostic::Diagnostic(double minDens){
    specificity = numeric_limits<double>::quiet_NaN();
    dens_lim = minDens;
    uses_hrp2 = false;
}

bool Diagnostic::isPositive( LocalRng& rng, double dens, double densHRP2 ) const {
    if( uses_hrp2 ){
        assert( densHRP2 == densHRP2 ); // monitoring diagnostic passes NaN; use of HRP2 is not supported
        dens = densHRP2;
    }
    if( (std::isnan)(specificity) ){
        // use deterministic test
        return dens >= dens_lim;
    }else{
        // dens_lim is dens_50 in this case
        double pPositive = 1.0 + specificity * (dens / (dens + dens_lim) - 1.0);
//         double pPositive = (dens + dens_lim - dens_lim * specificity) / (dens + dens_lim);       // equivalent
        // parity: order — local RNG stream: draw call order fixes stream state
        return rng.bernoulli( pPositive );
    }
}

bool Diagnostic::allowsFalsePositives() const{
    if( (std::isnan)(specificity) ) return dens_lim <= 0.0;
    return specificity < 1.0;
}


// ———  diagnostics (static)  ———

void diagnostics::clear( RunContext& ctx ){
    ctx.ext->diagnostic_set.clear();
}


const Diagnostic& diagnostics::get( RunContext& ctx, const string& name ){
    auto& diagnostic_set = ctx.ext->diagnostic_set;
    auto it = diagnostic_set.find(name);
    if( it == diagnostic_set.end() ){
        throw util::xml_scenario_error( string("diagnostic not found: ").append(name) );
    }
    return it->second;
}

const Diagnostic& diagnostics::make_deterministic( RunContext& ctx, double minDens ){
    auto& diagnostic_set = ctx.ext->diagnostic_set;
    string name = "";   // anything not matching an existing name is fine
    assert( diagnostic_set.count(name) == 0 );  // we shouldn't need to call make_deterministic more than once, so I think this will do
    auto ret = diagnostic_set.insert( make_pair(std::move(name), Diagnostic( minDens )) );
    assert(ret.second);
    return ret.first->second;
}


const Diagnostic& diagnostics::monitoringDiagnostic( RunContext& ctx ){
    assert( ctx.ext->monitoring_diagnostic != 0 );
    return *ctx.ext->monitoring_diagnostic;
}
} }
