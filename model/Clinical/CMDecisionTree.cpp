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

#include "Clinical/CMDecisionTree.h"
#include "Clinical/ClinicalModel.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"
#include "util/errors.h"
#include "util/random.h"
#include <algorithm>
#include <cassert>

namespace OM {
namespace Clinical {
using interventions::ComponentId;

// ———  flat node-array interpreter (was: CMDecisionTree subclass exec())  ———

namespace {

/** First branch with key strictly > value (== map<double>::upper_bound on the
 * old branches map; keys are sorted unique doubles, no NaN). Returns the tree
 * size sentinel (nodes.size() is never a valid child of a non-leaf ... use
 * UINT32_MAX) if none. */
inline uint32_t branch_select( const FlatTree& t, const FlatTree::Node& nd,
    double value )
{
    const auto b = t.branches.begin();
    auto it = std::upper_bound( b + nd.aux, b + nd.end, value,
        [](double v, const DTBranch& e){ return v < e.key; } );
    return it == b + nd.end ? UINT32_MAX : it->child;
}

/** Walk from node n, accumulating outcomes into e. Iterative tail walk for
 * every node kind except DTK_Multiple (which sequences children and needs one
 * direct recursive call; nesting is compile-bounded and shallow). */
void walk( const FlatTree& t, DTExecCtx& e, uint32_t n )
{
    for( ;; ){
        const FlatTree::Node& nd = t.nodes[n];
        switch( nd.kind ){
        case DTK_Random: {
            // parity: d — decision-tree contract: one uniform_01 per Random
            // node, then first cumulative key strictly greater (was:
            // branches.upper_bound(uniform_01()) on map<double>).
            uint32_t c = branch_select( t, nd, e.rng.uniform_01() );
            assert( c != UINT32_MAX );  // last key >= 1.0 > uniform_01()
            n = c;
            continue;
        }
        case DTK_Age:
            // age is that of human at start of time step (i.e. may be 0)
            n = branch_select( t, nd, e.ageYears );
            if( n == UINT32_MAX )
                throw TRACED_EXCEPTION( "bad age-based decision tree switch",
                                        util::Error::PkPd );
            continue;
        case DTK_Multiple: {
            // parity: CMDTMultiple::exec quirk — treated ORs out of children
            // but screened does not leave the Multiple boundary: save/restore
            // around the children (nothing reads screened during descent, so
            // in-place accumulation with this restore is equivalent).
            bool sc = e.screened;
            for( uint32_t i = nd.aux; i != nd.end; ++i )
                walk( t, e, t.children[i] );
            e.screened = sc;
            return;
        }
        case DTK_CaseType:
            // Uses of this in complicated cases trigger an exception during
            // initialisation (compile).
            assert( (e.pgState & Episode::SICK)
                    && !(e.pgState & Episode::COMPLICATED) );
            n = (e.pgState & Episode::SECOND_CASE) ? nd.c1 : nd.c0;
            continue;
        case DTK_InfectionOrigin: {
            auto o = e.human.withinHostModel->getInfectionOrigin();
            n = o == WithinHost::InfectionOrigin::Imported ? nd.c0 :
                o == WithinHost::InfectionOrigin::Introduced ? nd.c1 : nd.c2;
            continue;
        }
        case DTK_Diagnostic: {
            // parity: the diagnostic's own draws happen here, before descent.
            const bool pos = e.human.withinHostModel->diagnosticResult(
                e.rng, *t.diagnostics[nd.aux] );
            walk( t, e, pos ? nd.c0 : nd.c1 );
            e.screened = true;    // set after the child ran (was: exec result)
            return;
        }
        case DTK_Uncomplicated: {
            const Clinical::Episode& latest = e.human.clinicalModel->getLatestReport();
            const SimTime memory = static_cast<SimTime>( nd.aux );
            const bool recent = ((e.pgState & Episode::SICK)
                                    && !(e.pgState & Episode::COMPLICATED))
                                || (e.pgState & Episode::MALARIA);
            n = (recent && latest.time + memory >= e.ctx.clock.nowOrTs0())
                ? nd.c0 : nd.c1;
            continue;
        }
        case DTK_Severe:
            n = (e.pgState & Episode::COMPLICATED) ? nd.c0 : nd.c1;
            continue;
        case DTK_Cohort:
            // Rely on the health system memory to not count the same episode twice
            n = e.human.isInSubPop( ComponentId( nd.aux ) ) ? nd.c0 : nd.c1;
            continue;
        case DTK_NoTreatment:
            return;
        case DTK_TreatFailure:
            e.treated = true;     // report treatment, no parasite effect
            return;
        case DTK_Report:
            for( uint32_t i = nd.aux; i != nd.end; ++i )
                mon::recordEvent( mon::measure("nCMDTReport"), e.human, 1,
                                  t.reports[i] );
            return;
        case DTK_TreatPKPD:
            for( uint32_t i = nd.aux; i != nd.end; ++i )
                e.human.withinHostModel->treatPkPd( t.pkpd[i].schedule,
                    t.pkpd[i].dosage, e.ageYears, 0.0 );
            e.treated = true;
            return;
        case DTK_TreatSimple: {
            bool bsTreatment = false;
            for( uint32_t i = nd.aux; i != nd.end; ++i )
                // last simple treatment's result is the outcome (parity:
                // CMDTTreatSimple::exec overwrote, did not OR)
                bsTreatment = e.human.withinHostModel->treatSimple( e.human,
                    t.simple[i].first, t.simple[i].second );
            e.treated = bsTreatment;
            return;
        }
        case DTK_Deploy:
            //NOTE: it's not intuitively obvious what treated should mean for
            // intervention deployment; deploying reports nothing (CMDTOut(false)).
            t.deploys[nd.aux].deploy( e.ctx, e.human, mon::Deploy::TREAT,
                                      interventions::VaccineLimits() );
            return;
        }
    }
}

} // unnamed namespace

CMDTOut exec_decision_tree( const FlatTree& tree, DTExecCtx& e )
{
    walk( tree, e, tree.root );
    CMDTOut out( e.treated );
    out.screened = e.screened;
    return out;
}

} }
