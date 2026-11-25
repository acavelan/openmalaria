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

#pragma once

#include "RunContext.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"
#include "util/errors.h"
#include <algorithm>
#include "Host/Human.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Clinical/Episode.h"
#include "interventions/Interfaces.h"
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace scnXml {
class DecisionTree;
}

namespace OpenMalaria {
class Diagnostic;

/** All output data from the decision tree. */
struct CMDTOut {
    CMDTOut() : treated(false), screened(false) {}
    explicit CMDTOut(bool t) : treated(t), screened(false) {}
    bool treated;  // true iff some blood-stage treatment was administered
    bool screened; // true iff some diagnostic was used
};

/** Node kinds of the compiled flat decision tree. See exec_decision_tree(). */
enum DTKind : uint8_t {
    // leaves (actions):
    DTK_NoTreatment,  // do nothing
    DTK_TreatFailure, // report treatment without affecting parasites
    DTK_Report,       // record monitoring events for each output id
    DTK_TreatPKPD,    // deliver one or more treatments via the PK/PD model
    DTK_TreatSimple,  // deliver one or more simple treatments
    DTK_Deploy,       // deploy one or more interventions
    // branching nodes:
    DTK_Random,          // one uniform_01, upper-bound of cumulative-p keys
    DTK_Age,             // upper-bound of age-category keys (years)
    DTK_Multiple,        // execute all children in order
    DTK_CaseType,        // first/second line (uncomplicated only)
    DTK_InfectionOrigin, // imported/introduced/indigenous
    DTK_Diagnostic,      // run a diagnostic; positive/negative child
    DTK_Uncomplicated,   // recent uncomplicated report within memory
    DTK_Severe,          // pgState has COMPLICATED
    DTK_Cohort,          // human is in the component's sub-population
};

/** One keyed branch of a Random/Age node: select first key strictly > value. */
struct DTBranch {
    double key;
    uint32_t child;
};

/** A decision tree compiled from XML into a flat node array plus side tables.
 * Built once per health-system load by compile_decision_tree() (xml layer);
 * owned by RunContext::decision_trees for the run's lifetime. Contains no
 * per-run host state, so it is never checkpointed (parity: was:
 * CMDecisionTree.cpp statics, likewise not checkpointed).
 *
 * Field use per kind (unused fields left 0):
 *  DTK_Random/DTK_Age:   branches[aux..end)   (sorted unique keys, folded)
 *  DTK_Multiple:         children[aux..end)
 *  DTK_Report:           reports[aux..end)
 *  DTK_TreatPKPD:        pkpd[aux..end)
 *  DTK_TreatSimple:      simple[aux..end)
 *  DTK_Deploy:           deploys[aux] (HumanIntervention list, pre-sorted)
 *  DTK_Diagnostic:       diagnostics[aux]
 *  DTK_Cohort:           aux = ComponentId number
 *  DTK_CaseType:         c0=firstLine, c1=secondLine
 *  DTK_InfectionOrigin:  c0=imported, c1=introduced, c2=indigenous
 *  DTK_Uncomplicated:    aux=memory (SimTime), c0=positive, c1=negative
 *  DTK_Severe/DTK_Cohort: c0=positive, c1=negative */
struct FlatTree {
    struct Node {
        uint32_t aux = 0, end = 0;       // side-table range [aux,end) or scalar
        uint32_t c0 = 0, c1 = 0, c2 = 0; // fixed children (node-array indices)
        uint8_t kind = 0;                // DTKind
    };
    struct PkPdTreat {
        size_t schedule;
        size_t dosage;
    };

    uint32_t root = 0;
    std::vector<Node> nodes;
    std::vector<DTBranch> branches;                  // Random: cumulative p; Age: upper bounds
    std::vector<uint32_t> children;                  // Multiple
    std::vector<int> reports;                        // Report: monitoring output numbers
    std::vector<PkPdTreat> pkpd;                     // TreatPKPD (delay_h is not applied here)
    std::vector<std::pair<SimTime, SimTime>> simple; // TreatSimple liver/blood
    /** THE clinical→interventions seam (phase 5, plan §5): compiled leaves
     * hold interventions by value (the only sanctioned direction; the
     * interventions module must not include Clinical/ — the last two
     * runtime edges live in the xml load layer: DecisionTreeComponent in
     * xml/HumanInterventionComponents.h and changeHS setHS in
     * xml/interventions.cpp; phase 7 owns moving those behind ctx.ext). */
    std::vector<HumanIntervention> deploys;
    std::vector<const Diagnostic*> diagnostics;
};

/** Per-execution state: one host, one walk, results accumulated in place.
 * Replaces the old CMHostData parameter-passing struct; the reference fields
 * must outlive the walk (they do: all call sites run it synchronously). */
struct DTExecCtx {
    DTExecCtx(RunContext& ctx_, Human& human_, double ageYears_, Episode::State pgState_)
        : ctx(ctx_), human(human_), rng(human_.rng()), ageYears(ageYears_), pgState(pgState_) {}
    RunContext& ctx;
    Human& human;
    LocalRng& rng;   // the host's stream; Random nodes draw from here
    double ageYears; // human age at start of the time step
    Episode::State pgState;
    bool treated = false; // accumulated outcome (see CMDTOut)
    bool screened = false;
};

/** Compile a decision tree from XML into the flat representation, resolving
 * drug/cohort/diagnostic names to ids, precomputing cumulative-p keys (in
 * map<double> numeric order, with tie-drop and identical-leaf-run folding)
 * and validating all constraints that the old node classes checked in their
 * constructors. Owned by (*ctx.decision_trees); identical content compiled
 * twice is stored twice (the old decision-library dedup only saved memory).
 * @param isUC see the old CMDecisionTree::create docs: false rejects caseType. */
const FlatTree& compile_decision_tree(RunContext& ctx, const scnXml::DecisionTree& node, bool isUC);

/** Execute the decision tree: one free-function interpreter over the flat
 * node array, no virtual dispatch during treatment decisions.
 *
 * Reporting: use of diagnostics is reported (CMDTOut::screened); treatment
 * is not, but CMDTOut::treated says whether any took place. */
CMDTOut exec_decision_tree(const FlatTree& tree, DTExecCtx& e);

// ———  flat node-array interpreter (was: CMDecisionTree subclass exec())  ———

namespace {

/** First branch with key strictly > value (== map<double>::upper_bound on the
 * old branches map; keys are sorted unique doubles, no NaN). Returns the tree
 * size sentinel (nodes.size() is never a valid child of a non-leaf ... use
 * UINT32_MAX) if none. */
inline uint32_t branch_select(const FlatTree& t, const FlatTree::Node& nd, double value) {
    const auto b = t.branches.begin();
    auto it = std::upper_bound(b + nd.aux, b + nd.end, value, [](double v, const DTBranch& e) { return v < e.key; });
    return it == b + nd.end ? UINT32_MAX : it->child;
}

/** Walk from node n, accumulating outcomes into e. Iterative tail walk for
 * every node kind except DTK_Multiple (which sequences children and needs one
 * direct recursive call; nesting is compile-bounded and shallow). */
void walk(const FlatTree& t, DTExecCtx& e, uint32_t n) {
    for (;;) {
        const FlatTree::Node& nd = t.nodes[n];
        switch (nd.kind) {
        case DTK_Random: {
            // parity: d — decision-tree contract: one uniform_01 per Random
            // node, then first cumulative key strictly greater (was:
            // branches.upper_bound(uniform_01()) on map<double>).
            uint32_t c = branch_select(t, nd, e.rng.uniform_01());
            assert(c != UINT32_MAX); // last key >= 1.0 > uniform_01()
            n = c;
            continue;
        }
        case DTK_Age:
            // age is that of human at start of time step (i.e. may be 0)
            n = branch_select(t, nd, e.ageYears);
            if (n == UINT32_MAX) throw TRACED_EXCEPTION("bad age-based decision tree switch", Error::PkPd);
            continue;
        case DTK_Multiple: {
            // parity: CMDTMultiple::exec quirk — treated ORs out of children
            // but screened does not leave the Multiple boundary: save/restore
            // around the children (nothing reads screened during descent, so
            // in-place accumulation with this restore is equivalent).
            bool sc = e.screened;
            for (uint32_t i = nd.aux; i != nd.end; ++i)
                walk(t, e, t.children[i]);
            e.screened = sc;
            return;
        }
        case DTK_CaseType:
            // Uses of this in complicated cases trigger an exception during
            // initialisation (compile).
            assert((e.pgState & SICK) && !(e.pgState & COMPLICATED));
            n = (e.pgState & SECOND_CASE) ? nd.c1 : nd.c0;
            continue;
        case DTK_InfectionOrigin: {
            auto o = e.human.withinHostModel->getInfectionOrigin();
            n = o == InfectionOrigin::Imported ? nd.c0 : o == InfectionOrigin::Introduced ? nd.c1 : nd.c2;
            continue;
        }
        case DTK_Diagnostic: {
            // parity: the diagnostic's own draws happen here, before descent.
            const bool pos = e.human.withinHostModel->diagnosticResult(e.rng, *t.diagnostics[nd.aux]);
            walk(t, e, pos ? nd.c0 : nd.c1);
            e.screened = true; // set after the child ran (was: exec result)
            return;
        }
        case DTK_Uncomplicated: {
            const Episode& latest = e.human.clinicalModel->latestReport;
            const SimTime memory = static_cast<SimTime>(nd.aux);
            const bool recent = ((e.pgState & SICK) && !(e.pgState & COMPLICATED)) || (e.pgState & MALARIA);
            n = (recent && latest.time + memory >= e.ctx.clock.nowOrTs0()) ? nd.c0 : nd.c1;
            continue;
        }
        case DTK_Severe:
            n = (e.pgState & COMPLICATED) ? nd.c0 : nd.c1;
            continue;
        case DTK_Cohort:
            // Rely on the health system memory to not count the same episode twice
            n = e.human.isInSubPop(ComponentId(nd.aux)) ? nd.c0 : nd.c1;
            continue;
        case DTK_NoTreatment:
            return;
        case DTK_TreatFailure:
            e.treated = true; // report treatment, no parasite effect
            return;
        case DTK_Report:
            for (uint32_t i = nd.aux; i != nd.end; ++i)
                recordEvent(measure("nCMDTReport"), e.human, 1, t.reports[i]);
            return;
        case DTK_TreatPKPD:
            for (uint32_t i = nd.aux; i != nd.end; ++i)
                e.human.withinHostModel->treatPkPd(t.pkpd[i].schedule, t.pkpd[i].dosage, e.ageYears, 0.0);
            e.treated = true;
            return;
        case DTK_TreatSimple: {
            bool bsTreatment = false;
            for (uint32_t i = nd.aux; i != nd.end; ++i)
                // last simple treatment's result is the outcome (parity:
                // CMDTTreatSimple::exec overwrote, did not OR)
                bsTreatment = e.human.withinHostModel->treatSimple(e.human, t.simple[i].first, t.simple[i].second);
            e.treated = bsTreatment;
            return;
        }
        case DTK_Deploy:
            //NOTE: it's not intuitively obvious what treated should mean for
            // intervention deployment; deploying reports nothing (CMDTOut(false)).
            t.deploys[nd.aux].deploy(e.ctx, e.human, Deploy::TREAT, VaccineLimits());
            return;
        }
    }
}

} // unnamed namespace

inline CMDTOut exec_decision_tree(const FlatTree& tree, DTExecCtx& e) {
    walk(tree, e, tree.root);
    CMDTOut out(e.treated);
    out.screened = e.screened;
    return out;
}

}
