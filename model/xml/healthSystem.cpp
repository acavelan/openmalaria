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

// model/xml/ conversion layer — health system / clinical side of the host
// modules: non-malaria fevers (PathogenesisModel), neonatal mortality,
// diagnostic selection, within-host treatment options and the whole
// decision-tree / case-management construction from scnXml.
// Moved verbatim from Host/WithinHost/PathogenesisModel.cpp,
// Host/NeonatalMortality.cpp and Host/WithinHost/Treatments.cpp. State stays
// defined in the original files (extern re-declared below).

#include "Clock.h"
#include "RunContext.h"
#include "Parameters.h"
#include "Host/NeonatalMortality.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/WHInterface.h"
#include "Clinical/CM5DayCommon.h"
#include "Clinical/CMDecisionTree.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/CM5DayCommon.h"
#include "Clinical/Episode.h"
#include "Clinical/EventScheduler.h"
#include "PkPd/LSTMMedicate.h"
#include "interventions/Interfaces.h"
#include "interventions/InterventionManager.h"
#include "mon/Monitoring.h"
#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "schema/healthSystem.h"
#include "schema/scenario.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace std;

namespace OM {
namespace WithinHost {
using namespace OM::util;
namespace Pathogenesis {

// State defined in PathogenesisModel.cpp (external linkage):
extern double pg_indirRiskCoFactor, pg_severeMalThreshold, pg_comorbIntercept, pg_inv_critAgeComorb;
extern AgeGroupInterpolator pg_NMF_incidence;
extern bool opt_predetermined_episodes, opt_mueller_pres_model;

void PathogenesisModel::init(RunContext& ctx, const Parameters& parameters, const scnXml::Clinical& clinical,
                             bool nmfOnly) {
    if (util::ModelOptions::option(util::NON_MALARIA_FEVERS)) {
        if (!clinical.getNonMalariaFevers().present()) {
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc = clinical.getNonMalariaFevers().get();
        pg_NMF_incidence.set(ctx.clock, nmfDesc.getIncidence(), "incidence");
    }
    if (nmfOnly)
        return;

    pg_indirRiskCoFactor = 1 - exp(-parameters[Parameter::INDIRECT_RISK_COFACTOR]);
    pg_severeMalThreshold = parameters[Parameter::SEVERE_MALARIA_THRESHHOLD] + 1;
    pg_comorbIntercept = 1 - exp(-parameters[Parameter::COMORBIDITY_INTERCEPT]);
    pg_inv_critAgeComorb = 1 / parameters[Parameter::CRITICAL_AGE_FOR_COMORBIDITY];

    if (util::ModelOptions::option(util::PREDETERMINED_EPISODES)) {
        opt_predetermined_episodes = true;
        //no separate init:
        PathogenesisModel::initPyrogen(ctx, parameters);
    } else {
        if (util::ModelOptions::option(util::MUELLER_PRESENTATION_MODEL)) {
            opt_mueller_pres_model = true;
            PathogenesisModel::initMueller(ctx, parameters);
        } else {
            PathogenesisModel::initPyrogen(ctx, parameters);
        }
    }
}

}

TreatmentId Treatments::addTreatment(RunContext& ctx, const scnXml::TreatmentOption& desc) {
    Treatments treatment(ctx, desc);
    TreatmentId id(ctx.whTreatments->size());
    ctx.whTreatments->push_back(std::move(treatment));
    return id;
}

Treatments::Stages stageFromString(const std::string& str) {
    if (str == "liver")
        return Treatments::LIVER;
    if (str == "blood")
        return Treatments::BLOOD;
    if (str == "both")
        return Treatments::BOTH;
    throw util::format_error(std::string("stage must be liver, blood or both, not ").append(str));
}
Treatments::Treatments(RunContext& ctx, const scnXml::TreatmentOption& elt)
    : TriggeredDeployments(ctx, elt), timeLiver(sim::zero()), timeBlood(sim::zero()) {
    for (auto it = elt.getClearInfections().begin(), end = elt.getClearInfections().end(); it != end; ++it) {
        try {
            //NOTE: if changing XSD, this should not be called "timesteps" or have a default unit
            SimTime len = UnitParse::readShortDuration(ctx.clock, it->getTimesteps(), UnitParse::STEPS);
            if (len < -ctx.clock.one_TS() || len == sim::zero()) {
                throw util::format_error("timesteps must be ≥ 1 or have special value -1");
            }
            Stages stage = stageFromString(it->getStage());
            if (stage & LIVER) {
                if (timeLiver != sim::zero()) // existing treatment configuration
                    throw util::format_error("multiple specification of liver stage effect");
                timeLiver = len;
            }
            if (stage & BLOOD) {
                if (timeBlood != sim::zero()) // existing treatment configuration
                    throw util::format_error("multiple specification of blood stage effect");
                timeBlood = len;
            }
        } catch (const util::format_error& e) {
            throw util::xml_scenario_error(string(".../clearInfections/timesteps: ").append(e.message()));
        }
    }
}

}

namespace Host {
using namespace OM::util;

// State defined in NeonatalMortality.cpp (external linkage):
extern std::vector<double> prevByGestationalAge;
extern const WithinHost::Diagnostic* neonatalDiagnostic;

void NeonatalMortality::init(RunContext& ctx, const scnXml::Clinical& clinical) {
    SimTime fiveMonths = sim::fromDays(5 * 30);
    prevByGestationalAge.assign(ctx.clock.inSteps(fiveMonths), 0.0);

    if (clinical.getNeonatalMortality().present()) {
        neonatalDiagnostic = &WithinHost::diagnostics::get(ctx, clinical.getNeonatalMortality().get().getDiagnostic());
    } else {
        //NOTE: this is a compatibility option for older scenarios
        if (CommandLine::option(CommandLine::DEPRECATION_WARNINGS)) {
            cerr << "Deprecation warning: specification of the diagnostic "
                    "used by the Neonatal Mortality model is recommended "
                    "(model/clinical/neonatalMortality)"
                 << endl;
        }
        neonatalDiagnostic = &WithinHost::diagnostics::monitoringDiagnostic(ctx);
    }
}
}
}

// ———  clinical / health-system case management: decision-tree compile  ———

namespace OM {
namespace Clinical {

// NOTE: the decision-tree node classes formerly here (CMDTMultiple … CMDTAge,
// moved verbatim from Clinical/CMDecisionTree.cpp in phase 2) are replaced by
// the compile step below: XML → FlatTree (flat node array), executed by
// exec_decision_tree() in Clinical/CMDecisionTree.cpp with no virtual
// dispatch. Compile-stage canonicalisation: names→ids (diagnostics, cohorts,
// drug/schedule/dosage tables, deploy components), cumulative-p precompute,
// XSD-sequence-rank child order, map<double> key-order reproduction and
// consecutive-identical-leaf folding.
//
// Compile-time ORDER of side effects (component-id registration via
// getComponentId, and which validation error fires first) mirrors the
// reference binary exactly: the old node classes built their subtrees in
// constructor ARGUMENTS, which gcc 16 x86-64 evaluates right-to-left (verified
// empirically — probe: f(0),f(1),f(2) logs 2,1,0), and ran getComponentId in
// the constructor BODY (i.e. after all children). Every compile_* helper below
// reproduces that visit order explicitly; T1 (63/63) is the proof.

using WithinHost::diagnostics;
using namespace OM::util;
using std::vector;

// State defined in ClinicalModel.cpp (external linkage):
extern bool opt_event_scheduler, opt_imm_outcomes;
extern double oddsRatioThreshold;
// State defined in EventScheduler.cpp (external linkage):
extern bool opt_non_malaria_fevers;
// State defined in DecisionTree5Day.cpp (external linkage):
extern const FlatTree *treeUCOfficial, *treeUCSelfTreat;

namespace {

struct CompileCtx {
    RunContext& ctx;
    FlatTree& t;
    bool isUC;
    uint32_t add(DTKind kind) {
        t.nodes.push_back(FlatTree::Node());
        t.nodes.back().kind = kind;
        return static_cast<uint32_t>(t.nodes.size() - 1);
    }
};

uint32_t compile_node(CompileCtx& c, const scnXml::DecisionTree& node);
uint32_t compile_case_type(CompileCtx&, const scnXml::DTCaseType&);
uint32_t compile_infection_origin(CompileCtx&, const scnXml::DTInfectionOrigin&);
uint32_t compile_diagnostic(CompileCtx&, const scnXml::DTDiagnostic&);
uint32_t compile_uncomplicated(CompileCtx&, const scnXml::DTUncomplicated&);
uint32_t compile_severe(CompileCtx&, const scnXml::DTSevere&);
uint32_t compile_random(CompileCtx&, const scnXml::DTRandom&);
uint32_t compile_age(CompileCtx&, const scnXml::DTAge&);
uint32_t compile_cohort(CompileCtx&, const scnXml::DTCohort&);
uint32_t compile_report(CompileCtx&, const scnXml::DecisionTree::ReportSequence&);
uint32_t compile_treat_pkpd(CompileCtx&, const scnXml::DecisionTree::TreatPKPDSequence&);
uint32_t compile_treat_simple(CompileCtx&, const scnXml::DecisionTree::TreatSimpleSequence&);
uint32_t compile_deploy(CompileCtx&, const scnXml::DecisionTree::DeploySequence&);

/** Sort+normalise a keyed branch table to map<double,*> semantics, then fold
 * and append to the tree. Keys arrive in XSD order; the old map kept entries
 * in NUMERIC key order — reproduced by std::stable_sort with operator<:
 *  - ties: map::insert drops a key already present, keeping the FIRST
 *    inserted (XSD-earliest) child; stable_sort keeps equal keys adjacent in
 *    XSD order, so dropping all but the first of each tie group is identical;
 *  - NaN keys: impossible to sort (no strict weak ordering) and were undefined
 *    behaviour in the old map (broke its invariants silently). xs:double
 *    without a range facet parses "NaN", and the old sum check
 *    (NaN < 1.0 || NaN > 1.001) passed NaN through, so NaN keys WERE reachable
 *    as UB. Both callers now reject non-finite inputs with a clean
 *    xml_scenario_error before this point, so plain < sorting is exact.
 *    (±0 keys are unreachable: keys are running sums starting +0.0, and
 *    +0.0 + p never yields −0.0 for finite p in [0,1].)
 *
 * Leaf-run folding — legality proof per rule: adjacent (equal-key-free)
 * branches selecting the SAME child collapse to the run's last (largest) key.
 * upper_bound semantics: any value v either has the first key > v inside the
 * run (both versions select that child), or v >= every key of the run except
 * possibly reaching it (folded version skips straight to the run's last key =
 * same child), or v < all keys of the run (both select the run's first key
 * entry's child — the same child, since the run is child-identical). No draw
 * is consumed during selection: DTK_Random draws exactly ONE uniform_01
 * BEFORE selecting (contract §3d) and DTK_Age consumes no draws, so folding
 * cannot skip or duplicate a consumed draw. */
void add_branch_table(CompileCtx& c, uint32_t n, vector<DTBranch>&& b) {
    // parity: d — map<double> numeric key order == stable_sort by operator< with
    // XSD-first tie-drop (see file comment; plan §3d):
    std::stable_sort(b.begin(), b.end(), [](const DTBranch& x, const DTBranch& y) { return x.key < y.key; });
    // drop tied keys (keep first inserted, == map::insert behaviour):
    size_t w = 0;
    for (const DTBranch& br : b) {
        if (w && br.key == b[w - 1].key)
            continue;
        b[w++] = br;
    }
    b.resize(w);
    // fold adjacent identical children (keep run's last key):
    w = 0;
    for (const DTBranch& br : b) {
        if (w && b[w - 1].child == br.child) {
            b[w - 1] = br;
            continue;
        }
        b[w++] = br;
    }
    FlatTree::Node& nd = c.t.nodes[n];
    nd.aux = static_cast<uint32_t>(c.t.branches.size());
    c.t.branches.insert(c.t.branches.end(), b.begin(), b.begin() + w);
    nd.end = static_cast<uint32_t>(c.t.branches.size());
}

uint32_t compile_multiple(CompileCtx& c, const scnXml::DTMultiple& node) {
    const uint32_t n = c.add(DTK_Multiple);
    vector<uint32_t> kids;
    // parity: d — child chain order == XSD sequence of decision-tree element
    // choices (plan §3d): this fixed enumeration IS the schema-sequence rank
    // (the XSD codegen getters group by element; the old CMDTMultiple::create
    // walked them in this same rank order, not document order).
    for (const scnXml::DTCaseType& sn : node.getCaseType())
        kids.push_back(compile_case_type(c, sn));
    for (const scnXml::DTInfectionOrigin& sn : node.getInfectionOrigin())
        kids.push_back(compile_infection_origin(c, sn));
    for (const scnXml::DTDiagnostic& sn : node.getDiagnostic())
        kids.push_back(compile_diagnostic(c, sn));
    for (const scnXml::DTUncomplicated& sn : node.getUncomplicated())
        kids.push_back(compile_uncomplicated(c, sn));
    for (const scnXml::DTSevere& sn : node.getSevere())
        kids.push_back(compile_severe(c, sn));
    for (const scnXml::DTRandom& sn : node.getRandom())
        kids.push_back(compile_random(c, sn));
    for (const scnXml::DTAge& sn : node.getAge())
        kids.push_back(compile_age(c, sn));
    for (const scnXml::DTCohort& sn : node.getCohort())
        kids.push_back(compile_cohort(c, sn));
    if (node.getTreatPKPD().size())
        kids.push_back(compile_treat_pkpd(c, node.getTreatPKPD()));
    if (node.getTreatSimple().size())
        kids.push_back(compile_treat_simple(c, node.getTreatSimple()));
    if (node.getReport().size())
        kids.push_back(compile_report(c, node.getReport()));
    if (node.getDeploy().size())
        kids.push_back(compile_deploy(c, node.getDeploy()));
    FlatTree::Node& nd = c.t.nodes[n];
    nd.aux = static_cast<uint32_t>(c.t.children.size());
    c.t.children.insert(c.t.children.end(), kids.begin(), kids.end());
    nd.end = static_cast<uint32_t>(c.t.children.size());
    return n;
}

uint32_t compile_case_type(CompileCtx& c, const scnXml::DTCaseType& node) {
    if (!c.isUC) {
        throw util::xml_scenario_error("decision tree: caseType can only be used for uncomplicated cases");
    }
    const uint32_t n = c.add(DTK_CaseType);
    // gcc arg order: secondLine evaluated first (see file comment)
    const uint32_t second = compile_node(c, node.getSecondLine());
    const uint32_t first = compile_node(c, node.getFirstLine());
    c.t.nodes[n].c0 = first;
    c.t.nodes[n].c1 = second;
    return n;
}

uint32_t compile_infection_origin(CompileCtx& c, const scnXml::DTInfectionOrigin& node) {
    const uint32_t n = c.add(DTK_InfectionOrigin);
    // gcc arg order: indigenous, introduced, imported
    const uint32_t ind = compile_node(c, node.getIndigenous());
    const uint32_t intro = compile_node(c, node.getIntroduced());
    const uint32_t imp = compile_node(c, node.getImported());
    c.t.nodes[n].c0 = imp;
    c.t.nodes[n].c1 = intro;
    c.t.nodes[n].c2 = ind;
    return n;
}

uint32_t compile_diagnostic(CompileCtx& c, const scnXml::DTDiagnostic& node) {
    const uint32_t n = c.add(DTK_Diagnostic);
    // gcc arg order: negative, positive, then diagnostics::get (lookup)
    const uint32_t neg = compile_node(c, node.getNegative());
    const uint32_t pos = compile_node(c, node.getPositive());
    c.t.diagnostics.push_back(&diagnostics::get(c.ctx, node.getDiagnostic()));
    c.t.nodes[n].aux = static_cast<uint32_t>(c.t.diagnostics.size() - 1);
    c.t.nodes[n].c0 = pos;
    c.t.nodes[n].c1 = neg;
    return n;
}

uint32_t compile_uncomplicated(CompileCtx& c, const scnXml::DTUncomplicated& node) {
    const uint32_t n = c.add(DTK_Uncomplicated);
    // gcc arg order: negative, positive, then readShortDuration(memory)
    const uint32_t neg = compile_node(c, node.getNegative());
    const uint32_t pos = compile_node(c, node.getPositive());
    // (readShortDuration may throw format_error; like the old ctor argument,
    // it escapes uncaught)
    const SimTime memory = UnitParse::readShortDuration(c.ctx.clock, node.getMemory(), UnitParse::STEPS);
    if (memory > ClinicalModel::hsMemory())
        throw util::xml_scenario_error(
            "<uncomplicated> memory parameter must be less than or equal to the healthsystem memory (hsmemory parameter)");
    c.t.nodes[n].aux = static_cast<uint32_t>(memory); // two's-complement, C++20
    c.t.nodes[n].c0 = pos;
    c.t.nodes[n].c1 = neg;
    return n;
}

uint32_t compile_severe(CompileCtx& c, const scnXml::DTSevere& node) {
    const uint32_t n = c.add(DTK_Severe);
    const uint32_t neg = compile_node(c, node.getNegative());
    const uint32_t pos = compile_node(c, node.getPositive());
    c.t.nodes[n].c0 = pos;
    c.t.nodes[n].c1 = neg;
    return n;
}

uint32_t compile_cohort(CompileCtx& c, const scnXml::DTCohort& node) {
    const uint32_t n = c.add(DTK_Cohort);
    // gcc arg order: negative, positive; getComponentId ran in the ctor body
    const uint32_t neg = compile_node(c, node.getNegative());
    const uint32_t pos = compile_node(c, node.getPositive());
    c.t.nodes[n].aux = interventions::InterventionManager::getComponentId(c.ctx, node.getComponent()).id;
    c.t.nodes[n].c0 = pos;
    c.t.nodes[n].c1 = neg;
    return n;
}

// action leaves

uint32_t compile_report(CompileCtx& c, const scnXml::DecisionTree::ReportSequence& seq) {
    const uint32_t n = c.add(DTK_Report);
    FlatTree::Node& nd = c.t.nodes[n];
    nd.aux = static_cast<uint32_t>(c.t.reports.size());
    for (const scnXml::DTReport& reportElt : seq)
        c.t.reports.push_back(reportElt.getOutputNumber());
    nd.end = static_cast<uint32_t>(c.t.reports.size());
    assert(nd.aux != nd.end); // CMDTReport should not be used empty
    return n;
}

uint32_t compile_treat_pkpd(CompileCtx& c, const scnXml::DecisionTree::TreatPKPDSequence& seq) {
    const uint32_t n = c.add(DTK_TreatPKPD);
    FlatTree::Node& nd = c.t.nodes[n];
    nd.aux = static_cast<uint32_t>(c.t.pkpd.size());
    for (const scnXml::DTTreatPKPD& treatElt : seq) {
        // delay_h deliberately not stored: the old exec() ignored it, passing
        // 0.0 to treatPkPd (the delay-carrying variant is the intervention
        // component TreatPKPDComponent, unaffected here).
        c.t.pkpd.push_back(FlatTree::PkPdTreat{c.ctx.treatments.findSchedule(treatElt.getSchedule()),
                                               c.ctx.treatments.findDosages(treatElt.getDosage())});
    }
    nd.end = static_cast<uint32_t>(c.t.pkpd.size());
    assert(nd.aux != nd.end); // CMDTTreatPKPD should not be used empty
    return n;
}

uint32_t compile_treat_simple(CompileCtx& c, const scnXml::DecisionTree::TreatSimpleSequence& seq) {
    const uint32_t n = c.add(DTK_TreatSimple);
    FlatTree::Node& nd = c.t.nodes[n];
    nd.aux = static_cast<uint32_t>(c.t.simple.size());
    for (const scnXml::DTTreatSimple& treatElt : seq) {
        //NOTE: this code is currently identical to that in SimpleTreatComponent
        try {
            SimTime durL = UnitParse::readShortDuration(c.ctx.clock, treatElt.getDurationLiver(), UnitParse::NONE),
                    durB = UnitParse::readShortDuration(c.ctx.clock, treatElt.getDurationBlood(), UnitParse::NONE);
            SimTime neg1 = -c.ctx.clock.one_TS();
            if (durL < neg1 || durB < neg1) {
                throw util::xml_scenario_error("treatSimple: cannot have durationBlood or durationLiver less than -1");
            }
            c.t.simple.push_back(std::make_pair(durL, durB));
        } catch (const util::format_error& e) {
            throw util::xml_scenario_error(string("treatSimple: ").append(e.message()));
        }
    }
    nd.end = static_cast<uint32_t>(c.t.simple.size());
    return n;
}

uint32_t compile_deploy(CompileCtx& c, const scnXml::DecisionTree::DeploySequence& seq) {
    // HumanIntervention resolves component names to (deduplicated)
    // components and applies the standard stable_sort by Component::Type
    // (reproducibility contract b) — same construction as the old CMDTDeploy.
    const uint32_t n = c.add(DTK_Deploy);
    c.t.deploys.emplace_back(c.ctx, seq);
    c.t.nodes[n].aux = static_cast<uint32_t>(c.t.deploys.size() - 1);
    return n;
}

uint32_t compile_random(CompileCtx& c, const scnXml::DTRandom& node) {
    const uint32_t n = c.add(DTK_Random);
    vector<DTBranch> b;
    double cum_p = 0.0;
    for (const scnXml::Outcome& outcome : node.getOutcome()) {
        // parity: d — cumulative keys summed in XSD outcome order; selection
        // then follows map<double> numeric key order (reproduced by
        // add_branch_table; plan §3d); the last key must reach 1.
        const double p = outcome.getP();
        // xs:double parses "NaN", which the old sum check below passed and
        // the old map<double> stored (UB: breaks the map's ordering
        // invariants). Sorted tables cannot carry it, so reject it here —
        // behaviour changes only for inputs that were undefined before.
        if (p != p)
            throw util::xml_scenario_error("decision tree (random node): outcome probability is not a number");
        cum_p += p;
        b.push_back(DTBranch{cum_p, compile_node(c, outcome)});
    }
    // Test cum_p is approx. 1.0 in case the input tree is wrong. We require no
    // less than one to make sure generated random numbers are not greater than
    // the last option.
    if (cum_p < 1.0 || cum_p > 1.001) {
        throw util::xml_scenario_error("decision tree (random node): expected probability sum to be 1.0 but found " +
                                       to_string(cum_p));
    }
    add_branch_table(c, n, std::move(b));
    return n;
}

uint32_t compile_age(CompileCtx& c, const scnXml::DTAge& node) {
    const uint32_t n = c.add(DTK_Age);
    vector<DTBranch> b;
    double lastAge = numeric_limits<double>::quiet_NaN();
    uint32_t lastNode = UINT32_MAX; // == old null: never inserted for age[0]
    for (const scnXml::Age& age : node.getAge()) {
        if (lastAge != lastAge) {
            if (age.getLb() != 0.0)
                throw util::xml_scenario_error("decision tree age switch must have first lower bound equal 0");
        } else {
            if (age.getLb() <= lastAge) {
                throw util::xml_scenario_error("decision tree age switch must list age groups in increasing order");
            }
            b.push_back(DTBranch{age.getLb(), lastNode});
        }
        if (age.getLb() != age.getLb())
            throw util::xml_scenario_error("decision tree age switch: age bound is not a number");
        lastNode = compile_node(c, age);
        lastAge = age.getLb();
    }
    b.push_back(DTBranch{numeric_limits<double>::infinity(), lastNode});
    add_branch_table(c, n, std::move(b));
    return n;
}

uint32_t compile_node(CompileCtx& c, const scnXml::DecisionTree& node) {
    if (node.getMultiple().present())
        return compile_multiple(c, node.getMultiple().get());
    // branching nodes
    if (node.getCaseType().present())
        return compile_case_type(c, node.getCaseType().get());
    if (node.getInfectionOrigin().present())
        return compile_infection_origin(c, node.getInfectionOrigin().get());
    if (node.getDiagnostic().present())
        return compile_diagnostic(c, node.getDiagnostic().get());
    if (node.getUncomplicated().present())
        return compile_uncomplicated(c, node.getUncomplicated().get());
    if (node.getSevere().present())
        return compile_severe(c, node.getSevere().get());
    if (node.getRandom().present())
        return compile_random(c, node.getRandom().get());
    if (node.getAge().present())
        return compile_age(c, node.getAge().get());
    if (node.getCohort().present())
        return compile_cohort(c, node.getCohort().get());

    // action nodes
    if (node.getNoTreatment().present())
        return c.add(DTK_NoTreatment);
    if (node.getReport().size())
        return compile_report(c, node.getReport());
    if (node.getTreatFailure().present())
        return c.add(DTK_TreatFailure);
    if (node.getTreatPKPD().size())
        return compile_treat_pkpd(c, node.getTreatPKPD());
    if (node.getTreatSimple().size())
        return compile_treat_simple(c, node.getTreatSimple());
    if (node.getDeploy().size())
        return compile_deploy(c, node.getDeploy());
    throw xml_scenario_error("unterminated decision tree");
}

} // unnamed namespace

const FlatTree& compile_decision_tree(RunContext& ctx, const scnXml::DecisionTree& node, bool isUC) {
    std::unique_ptr<FlatTree> tree(new FlatTree());
    CompileCtx c{ctx, *tree, isUC};
    tree->root = compile_node(c, node);
    FlatTree& r = *tree;
    // owned for the run (frees at end of run, like the old decision_library).
    // The old dedup-by-equality is gone: identical content compiled twice
    // behaves identically, dedup only saved memory.
    ctx.decision_trees->push_back(std::move(tree));
    return r;
}

void ClinicalEventScheduler::init(RunContext& ctx, const Parameters& parameters, const scnXml::Clinical& clinical) {
    if (ctx.clock.one_TS() != sim::oneDay())
        throw util::xml_scenario_error("ClinicalEventScheduler is only designed for a 1-day time step.");

    opt_non_malaria_fevers = util::ModelOptions::option(util::NON_MALARIA_FEVERS);

    alpha = exp(-parameters[Parameter::CFR_NEG_LOG_ALPHA]);
    if (!(0.0 <= alpha && alpha <= 1.0)) {
        throw util::xml_scenario_error("Clinical outcomes: propDeathsFirstDay should be within range [0,1]");
    }
    neg_v = -parameters[Parameter::CFR_SCALE_FACTOR];

    if (util::ModelOptions::option(util::NON_MALARIA_FEVERS)) {
        if (!clinical.getNonMalariaFevers().present()) {
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc2 = clinical.getNonMalariaFevers().get();
        if (!nmfDesc2.getPrNeedTreatmentMF().present() || !nmfDesc2.getPrNeedTreatmentNMF().present()) {
            throw util::xml_scenario_error("prNeedTreatmentMF and "
                                           "prNeedTreatmentNMF elements required in model->clinical->"
                                           "NonMalariaFevers");
        }
        NMF_need_antibiotic.set(ctx.clock, nmfDesc2.getPrNeedTreatmentNMF().get(), "prNeedTreatmentNMF");
        MF_need_antibiotic.set(ctx.clock, nmfDesc2.getPrNeedTreatmentMF().get(), "prNeedTreatmentMF");
    }
}

void ClinicalEventScheduler::setParameters(const SimClock& clock, const scnXml::HSEventScheduler& esData) {
    const scnXml::ClinicalOutcomes& coData = esData.getClinicalOutcomes();

    maxUCSeekingMemory = sim::fromDays(coData.getMaxUCSeekingMemory());
    uncomplicatedCaseDuration = sim::fromDays(coData.getUncomplicatedCaseDuration());
    complicatedCaseDuration = sim::fromDays(coData.getComplicatedCaseDuration());
    extraDaysAtRisk = sim::fromDays(coData.getComplicatedRiskDuration()) - complicatedCaseDuration;
    if (uncomplicatedCaseDuration < sim::fromDays(1) || complicatedCaseDuration < sim::fromDays(1) ||
        maxUCSeekingMemory < sim::zero() ||
        extraDaysAtRisk + complicatedCaseDuration < sim::fromDays(1) // at risk at least 1 day
        || extraDaysAtRisk > sim::zero()                             // at risk longer than case duration
    ) {
        throw util::xml_scenario_error(
            "Clinical outcomes: constraints on case/risk/memory duration not met (see documentation)");
    }

    cumDailyPrImmUCTS.reserve(coData.getDailyPrImmUCTS().size());
    double cumP = 0.0;
    for (auto it = coData.getDailyPrImmUCTS().begin(); it != coData.getDailyPrImmUCTS().end(); ++it) {
        cumP += *it;
        cumDailyPrImmUCTS.push_back(cumP);
    }
    if (cumP < 0.999 || cumP > 1.001) {
        throw util::xml_scenario_error("Event scheduler: dailyPrImmUCTS seq must add up to 1");
    }
    cumDailyPrImmUCTS.back() = 1.0;

    ClinicalModel::caseFatalityRate.scale(alpha);

    if (util::ModelOptions::option(util::NON_MALARIA_FEVERS)) {
        if (!esData.getNonMalariaFevers().present()) {
            throw util::xml_scenario_error("NonMalariaFevers element of healthSystem->EventScheduler required");
        }
        const scnXml::HSESNMF& nmfDesc = esData.getNonMalariaFevers().get();

        double pT = nmfDesc.getPrTreatment();
        logOddsAbBase = log(pT / (1.0 - pT));
        logOddsAbNegTest = log(nmfDesc.getEffectNegativeTest());
        logOddsAbPosTest = log(nmfDesc.getEffectPositiveTest());
        logOddsAbNeed = log(nmfDesc.getEffectNeed());
        logOddsAbInformal = log(nmfDesc.getEffectInformal());
        oneMinusEfficacyAb = 1.0 - nmfDesc.getTreatmentEfficacy();
        severeNmfMortality.set(clock, nmfDesc.getCFR(), "CFR");
    }
}
void ESCaseManagement::setHealthSystem(RunContext& ctx, const scnXml::HSEventScheduler& esData) {
    escm_uncomplicated = &compile_decision_tree(ctx, esData.getUncomplicated(), true);
    escm_complicated = &compile_decision_tree(ctx, esData.getComplicated(), false);

    // Calling our parent class like this is messy. Changing this would require
    // moving change-of-health-system handling into ClinicalModel.
    ClinicalEventScheduler::setParameters(ctx.clock, esData);
}
// ———  static, utility functions  ———

double getHealthSystemACRByName(const scnXml::TreatmentDetails& td, const string& drug) {
    auto val = [](const auto& d) { return d.present() ? d.get().getValue() : 0; };
    double v;
    if (drug == "CQ")
        v = val(td.getCQ());
    else if (drug == "SP")
        v = val(td.getSP());
    else if (drug == "AQ")
        v = val(td.getAQ());
    else if (drug == "SPAQ")
        v = val(td.getSPAQ());
    else if (drug == "ACT")
        v = val(td.getACT());
    else if (drug == "QN")
        v = val(td.getQN());
    else if (drug == "selfTreatment")
        v = td.getSelfTreatment().getValue();
    else
        throw util::xml_scenario_error("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");

    if (!(v >= 0.0 && v <= 1.0)) {
        throw util::xml_scenario_error(string("healthSystem initialACR/compliance/nonCompliersEffective: ")
                                           .append(drug)
                                           .append(" value must be in range [0,1]"));
    }
    return v;
}

// note: "CQOptional" is the same type as for other drugs
const scnXml::TreatmentActions::CQOptional& getHealthSystemTreatmentOptional(const scnXml::TreatmentActions& treatActs,
                                                                             const string& drug) {
    if (drug == "CQ")
        return treatActs.getCQ();
    if (drug == "SP")
        return treatActs.getSP();
    if (drug == "AQ")
        return treatActs.getAQ();
    if (drug == "SPAQ")
        return treatActs.getSPAQ();
    if (drug == "ACT")
        return treatActs.getACT();
    if (drug == "QN")
        return treatActs.getQN();
    throw util::xml_scenario_error("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
}
WithinHost::TreatmentId getHealthSystemTreatmentByName(RunContext& ctx, const scnXml::TreatmentActions& treatActs,
                                                       const string& drug) {
    const scnXml::TreatmentActions::CQOptional& elt = getHealthSystemTreatmentOptional(treatActs, drug);
    if (!elt.present())
        throw util::xml_scenario_error(string("healthSystem.treatmentActions: description required for ").append(drug));
    return WithinHost::WHInterface::addTreatment(ctx, elt.get());
}

// ———  static, init  ———

// Access probabilities: ImmediateOutcomes and DecisionTree5Day computed these
// from identically-named getters, byte-identical code (merged; the checks and
// error string travel verbatim):
template <class HS> static void set_access_probabilities(const HS& hsDescription) {
    double accessUCOfficial1 = hsDescription.getPSeekOfficialCareUncomplicated1().getValue();
    double accessUCOfficial2 = hsDescription.getPSeekOfficialCareUncomplicated2().getValue();
    // Note: this asymmetry is historical, and probably matters little:
    CM5DayCommon::accessUCSelfTreat[CM5DayCommon::FirstLine] = hsDescription.getPSelfTreatUncomplicated().getValue();
    CM5DayCommon::accessUCSelfTreat[CM5DayCommon::SecondLine] = 0.0;
    CM5DayCommon::accessUCAny[CM5DayCommon::FirstLine] =
        accessUCOfficial1 + CM5DayCommon::accessUCSelfTreat[CM5DayCommon::FirstLine];
    CM5DayCommon::accessUCAny[CM5DayCommon::SecondLine] =
        accessUCOfficial1 + CM5DayCommon::accessUCSelfTreat[CM5DayCommon::SecondLine];
    CM5DayCommon::accessSevere = hsDescription.getPSeekOfficialCareSevere().getValue();

    if (accessUCOfficial1 < 0.0 || CM5DayCommon::accessUCSelfTreat[CM5DayCommon::FirstLine] < 0.0 ||
        CM5DayCommon::accessUCAny[CM5DayCommon::FirstLine] > 1.0 || accessUCOfficial2 < 0.0 ||
        CM5DayCommon::accessUCSelfTreat[CM5DayCommon::SecondLine] < 0.0 ||
        CM5DayCommon::accessUCAny[CM5DayCommon::SecondLine] > 1.0 || CM5DayCommon::accessSevere < 0.0 ||
        CM5DayCommon::accessSevere > 1.0) {
        throw util::xml_scenario_error("healthSystem: "
                                       "pSeekOfficialCareUncomplicated1 and pSelfTreatUncomplicated must be"
                                       " at least 0 and their sum must be at most 1, and"
                                       "pSeekOfficialCareUncomplicated2 and pSeekOfficialCareSevere must "
                                       "be in range [0,1]");
    }
}
// liverStageDrug handling was identical in both setHealthSystem paths:
static void set_liver_stage(const xsd::cxx::tree::optional<scnXml::LiverStageDrug>& elt) {
    if (ModelOptions::option(util::VIVAX_SIMPLE_MODEL)) {
        WithinHost::WHVivax::setHSParameters(elt.present() ? &elt.get() : 0);
    } else if (elt.present()) {
        throw util::xml_scenario_error("health-system's liverStageDrug element only supported by vivax");
    }
}

void ImmediateOutcomes::setHealthSystem(RunContext& ctx, const scnXml::HSImmediateOutcomes& hsDescription) {
    const string &firstLine = hsDescription.getDrugRegimen().getFirstLine(),
                 &secondLine = hsDescription.getDrugRegimen().getSecondLine(),
                 &inpatient = hsDescription.getDrugRegimen().getInpatient();

    // ———  calculate probability of getting treatment  ———

    set_access_probabilities(hsDescription);

    // ———  calculate probability of clearing parasites  ———

    const double complianceFirstLine = getHealthSystemACRByName(hsDescription.getCompliance(), firstLine);
    const double complianceSecondLine = getHealthSystemACRByName(hsDescription.getCompliance(), secondLine);

    const double cureRateFirstLine = getHealthSystemACRByName(hsDescription.getInitialACR(), firstLine);
    const double cureRateSecondLine = getHealthSystemACRByName(hsDescription.getInitialACR(), secondLine);

    const double nonCompliersEffectiveFirstLine =
        getHealthSystemACRByName(hsDescription.getNonCompliersEffective(), firstLine);
    const double nonCompliersEffectiveSecondLine =
        getHealthSystemACRByName(hsDescription.getNonCompliersEffective(), secondLine);

    const double complianceSelfTreatment = hsDescription.getCompliance().getSelfTreatment().getValue();
    const double cureRateSelfTreatment = hsDescription.getInitialACR().getSelfTreatment().getValue();
    if (!(complianceSelfTreatment >= 0.0 && complianceSelfTreatment <= 1.0 && cureRateSelfTreatment >= 0.0 &&
          cureRateSelfTreatment <= 1.0)) {
        throw util::xml_scenario_error(
            "healthSystem initialACR/compliance/nonCompliersEffective: pSelfTreatment must be in range [0,1]");
    }

    cureRateUCOfficial[CM5DayCommon::FirstLine] =
        complianceFirstLine * cureRateFirstLine + (1 - complianceFirstLine) * nonCompliersEffectiveFirstLine;
    cureRateUCOfficial[CM5DayCommon::SecondLine] =
        complianceSecondLine * cureRateSecondLine + (1 - complianceSecondLine) * nonCompliersEffectiveSecondLine;
    cureRateUCSelfTreat[CM5DayCommon::FirstLine] = complianceSelfTreatment * cureRateSelfTreatment +
                                                   (1 - complianceSelfTreatment) * nonCompliersEffectiveFirstLine;
    cureRateUCSelfTreat[CM5DayCommon::SecondLine] = 0.0; // as with access, this is a historical asymmetry

    CM5DayCommon::cureRateSevere = getHealthSystemACRByName(hsDescription.getInitialACR(), inpatient);

    treatmentUC[CM5DayCommon::FirstLine] =
        getHealthSystemTreatmentByName(ctx, hsDescription.getTreatmentActions(), firstLine);

    if (secondLine == firstLine)
        treatmentUC[CM5DayCommon::SecondLine] = treatmentUC[CM5DayCommon::FirstLine];
    else
        treatmentUC[CM5DayCommon::SecondLine] =
            getHealthSystemTreatmentByName(ctx, hsDescription.getTreatmentActions(), secondLine);

    if (inpatient == firstLine)
        CM5DayCommon::treatmentSevere = treatmentUC[CM5DayCommon::FirstLine];
    else if (inpatient == secondLine)
        CM5DayCommon::treatmentSevere = treatmentUC[CM5DayCommon::SecondLine];
    else
        CM5DayCommon::treatmentSevere =
            getHealthSystemTreatmentByName(ctx, hsDescription.getTreatmentActions(), inpatient);

    useDiagnosticUC = hsDescription.getUseDiagnosticUC();

    set_liver_stage(hsDescription.getLiverStageDrug());
}
void DecisionTree5Day::setHealthSystem(RunContext& ctx, const scnXml::HSDT5Day& hsDescription) {
    set_access_probabilities(hsDescription);

    treeUCOfficial = &compile_decision_tree(ctx, hsDescription.getTreeUCOfficial(), true);
    treeUCSelfTreat = &compile_decision_tree(ctx, hsDescription.getTreeUCSelfTreat(), true);

    CM5DayCommon::cureRateSevere = hsDescription.getCureRateSevere().getValue();
    CM5DayCommon::treatmentSevere = WHInterface::addTreatment(ctx, hsDescription.getTreatmentSevere());

    set_liver_stage(hsDescription.getLiverStageDrug());
}
void ClinicalModel::init(RunContext& ctx, const Parameters& parameters, const scnXml::Scenario& scenario) {
    const scnXml::Clinical& clinical = scenario.getModel().getClinical();
    try {
        indirectMortBugfix = util::ModelOptions::option(util::INDIRECT_MORTALITY_FIX);
        //NOTE: if changing XSD, this should not have a default unit:
        healthSystemMemory =
            UnitParse::readShortDuration(ctx.clock, clinical.getHealthSystemMemory(), UnitParse::STEPS);
        oddsRatioThreshold = exp(parameters[Parameter::LOG_ODDS_RATIO_CF_COMMUNITY]);
        InfantMortality::init(ctx, parameters);
    } catch (const util::format_error& e) {
        throw util::xml_scenario_error(string("model/clinical/healthSystemMemory: ").append(e.message()));
    }

    if (util::ModelOptions::option(util::CLINICAL_EVENT_SCHEDULER)) {
        opt_event_scheduler = true;
        ClinicalEventScheduler::init(ctx, parameters, clinical);
    } else {
        if (scenario.getHealthSystem().getImmediateOutcomes().present()) {
            opt_imm_outcomes = true;

            if (util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS)) {
                cerr
                    << "Deprecation warning: healthSystem: use of ImmediateOutcomes can be replaced by the more flexible DecisionTree5Day (optional)"
                    << endl;
            }
        }
        // else: decision tree 5 day

        CM5DayCommon::init();
    }
}

void ClinicalModel::setHS(RunContext& ctx, const scnXml::HealthSystem& healthSystem) {
    caseFatalityRate.set(ctx.clock, healthSystem.getCFR(), "CFR");
    pSequelaeInpatient.set(ctx.clock, healthSystem.getPSequelaeInpatient(), "pSequelaeInpatient");
    // Presence check + message for the three sections; bytes composed from the
    // section name exactly as the three former blocks spelled them out:
    auto need = [](const auto& opt, const char* name) -> decltype(auto) {
        if (!opt.present())
            throw util::xml_scenario_error(string("Expected ") + name +
                                           " section in healthSystem data (initial or intervention)");
        return (opt.get()); // parenthesised: stay a reference
    };
    if (opt_event_scheduler)
        ESCaseManagement::setHealthSystem(ctx, need(healthSystem.getEventScheduler(), "EventScheduler"));
    else if (opt_imm_outcomes)
        ImmediateOutcomes::setHealthSystem(ctx, need(healthSystem.getImmediateOutcomes(), "ImmediateOutcomes"));
    else
        DecisionTree5Day::setHealthSystem(ctx, need(healthSystem.getDecisionTree5Day(), "DecisionTree5Day"));
}
}
}
