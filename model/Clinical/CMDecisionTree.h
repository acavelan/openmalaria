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

#ifndef Hmod_ESDecisionTree
#define Hmod_ESDecisionTree

#include "core/RunContext.h"
#include "Host/Human.h"
#include "Host/WithinHost/Pathogenesis/State.h"
#include "Clinical/Episode.h"
#include "interventions/Interfaces.h"
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace scnXml { class DecisionTree; }

namespace OM {
namespace WithinHost { class Diagnostic; }
namespace Clinical {

/** All output data from the decision tree. */
struct CMDTOut {
    CMDTOut() : treated(false), screened(false) {}
    explicit CMDTOut(bool t) : treated(t), screened(false) {}
    bool treated;       // true iff some blood-stage treatment was administered
    bool screened;      // true iff some diagnostic was used
};

/** Node kinds of the compiled flat decision tree. See exec_decision_tree(). */
enum DTKind : uint8_t {
    // leaves (actions):
    DTK_NoTreatment,    // do nothing
    DTK_TreatFailure,   // report treatment without affecting parasites
    DTK_Report,         // record monitoring events for each output id
    DTK_TreatPKPD,      // deliver one or more treatments via the PK/PD model
    DTK_TreatSimple,    // deliver one or more simple treatments
    DTK_Deploy,         // deploy one or more interventions
    // branching nodes:
    DTK_Random,         // one uniform_01, upper-bound of cumulative-p keys
    DTK_Age,            // upper-bound of age-category keys (years)
    DTK_Multiple,       // execute all children in order
    DTK_CaseType,       // first/second line (uncomplicated only)
    DTK_InfectionOrigin,// imported/introduced/indigenous
    DTK_Diagnostic,     // run a diagnostic; positive/negative child
    DTK_Uncomplicated,  // recent uncomplicated report within memory
    DTK_Severe,         // pgState has COMPLICATED
    DTK_Cohort,         // human is in the component's sub-population
};

/** One keyed branch of a Random/Age node: select first key strictly > value. */
struct DTBranch { double key; uint32_t child; };

/** A decision tree compiled from XML into a flat node array plus side tables.
 * Built once per health-system load by compile_decision_tree() (xml layer);
 * owned by RunExtensions::decision_trees for the run's lifetime. Contains no
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
        uint32_t aux = 0, end = 0;  // side-table range [aux,end) or scalar
        uint32_t c0 = 0, c1 = 0, c2 = 0; // fixed children (node-array indices)
        uint8_t kind = 0;           // DTKind
    };
    struct PkPdTreat { size_t schedule; size_t dosage; };

    uint32_t root = 0;
    std::vector<Node> nodes;
    std::vector<DTBranch> branches;   // Random: cumulative p; Age: upper bounds
    std::vector<uint32_t> children;   // Multiple
    std::vector<int> reports;         // Report: monitoring output numbers
    std::vector<PkPdTreat> pkpd;      // TreatPKPD (delay_h is not applied here)
    std::vector<std::pair<SimTime,SimTime>> simple; // TreatSimple liver/blood
    /** THE clinical→interventions seam (phase 5, plan §5): compiled leaves
     * hold interventions by value (the only sanctioned direction; the
     * interventions module must not include Clinical/ — the last two
     * runtime edges live in the xml load layer: DecisionTreeComponent in
     * xml/HumanInterventionComponents.h and changeHS setHS in
     * xml/interventions.cpp; phase 7 owns moving those behind ctx.ext). */
    std::vector<interventions::HumanIntervention> deploys;
    std::vector<const WithinHost::Diagnostic*> diagnostics;
};

/** Per-execution state: one host, one walk, results accumulated in place.
 * Replaces the old CMHostData parameter-passing struct; the reference fields
 * must outlive the walk (they do: all call sites run it synchronously). */
struct DTExecCtx {
    DTExecCtx( RunContext& ctx_, Host::Human& human_, double ageYears_,
               Episode::State pgState_ ) :
        ctx( ctx_ ), human( human_ ), rng( human_.rng ),
        ageYears( ageYears_ ), pgState( pgState_ ) {}
    RunContext& ctx;
    Host::Human& human;
    util::LocalRng& rng;      // the host's stream; Random nodes draw from here
    double ageYears;          // human age at start of the time step
    Episode::State pgState;
    bool treated = false;     // accumulated outcome (see CMDTOut)
    bool screened = false;
};

/** Compile a decision tree from XML into the flat representation, resolving
 * drug/cohort/diagnostic names to ids, precomputing cumulative-p keys (in
 * map<double> numeric order, with tie-drop and identical-leaf-run folding)
 * and validating all constraints that the old node classes checked in their
 * constructors. Owned by ctx.ext->decision_trees; identical content compiled
 * twice is stored twice (the old decision-library dedup only saved memory).
 * @param isUC see the old CMDecisionTree::create docs: false rejects caseType. */
const FlatTree& compile_decision_tree( RunContext& ctx,
    const scnXml::DecisionTree& node, bool isUC );

/** Execute the decision tree: one free-function interpreter over the flat
 * node array, no virtual dispatch during treatment decisions.
 *
 * Reporting: use of diagnostics is reported (CMDTOut::screened); treatment
 * is not, but CMDTOut::treated says whether any took place. */
CMDTOut exec_decision_tree( const FlatTree& tree, DTExecCtx& e );

} }
#endif
