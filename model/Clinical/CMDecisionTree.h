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

#include "Host/Human.h"
#include "Host/WithinHost/Pathogenesis/State.h"
#include "Clinical/Episode.h"
#include "interventions/Interfaces.h"

#include <map>
#include <optional>

namespace scnXml { class DecisionTree; }

namespace OM {
namespace WithinHost {
    class WHInterface;
    class Diagnostic;
}
namespace Clinical {
using WithinHost::WHInterface;

/** All data which needs to be passed to the decision tree evaluators. */
struct  CMHostData {
    /** Initialise, from a human h with age aY (years), and pathogenesis state
     * pS. Note that pS is only needed for uncomplicated cases
     * (see CMDecisionTree::create()). */
    CMHostData (Host::Human& h, double aY, Episode::State pS) :
        human(h), ageYears(aY), pgState(pS) {}
    Host::Human& human;
    double ageYears;    // only cached to save recalculating
    Episode::State pgState;
    inline WHInterface& withinHost(){ return *human.withinHostModel; }
};
/** All output data from the decision tree. */
struct CMDTOut {
    CMDTOut() : treated(false), screened(false) {}
    explicit CMDTOut(bool t) : treated(t), screened(false) {}
    bool treated;       // true iff some blood-stage treatment was administered
    bool screened;  // true iff some diagnostic was used
};


/**
 * A decision tree node: either a decision node (first/second line case, a
 * diagnostic with positive/negative outcome, a random decision) or an action.
 * 
 * Which fields are used depends on kind:
 *  - Multiple: children
 *  - CaseType: firstLine, secondLine
 *  - InfectionOrigin: imported, introduced, indigenous
 *  - Diagnostic: diagnostic, positive, negative
 *  - Uncomplicated: memory, positive, negative
 *  - Severe, Cohort: positive, negative (Cohort also component)
 *  - Random, Age: branches (keys: cumulative probabilities / age upper bounds)
 *  - Report: outIds
 *  - TreatPKPD: treatments
 *  - TreatSimple: timeLiver, timeBlood
 *  - Deploy: deployment
 *  - NoTreatment, TreatFailure: nothing
 *****************************************************************************/
class CMDecisionTree {
public:
    enum class Kind {
        Multiple, CaseType, InfectionOrigin, Diagnostic, Uncomplicated, Severe,
        Random, Age, Cohort, NoTreatment, Report, TreatFailure, TreatPKPD,
        TreatSimple, Deploy
    };
    
    explicit CMDecisionTree(Kind kind) : kind(kind) {}
    
    /** Create a user-configured decision from an XML node.
     * 
     * Memory management is handled internally (statically).
     * 
     * @param node XML element describing the tree
     * @param isUC If isUC is false and a "case type" decision is created, an
     *  xml_scenario_error exception is thrown. CMHostData::pgState is only
     *  used by the "case type" decision, so if isUC is false when creating the
     *  tree, pgState does not need to be set when executing the tree. */
    static const CMDecisionTree& create( const ::scnXml::DecisionTree& node, bool isUC,
        const interventions::InterventionRegistry& reg );
    
    /** Execute the decision tree.
     * 
     * Reporting: use of diagnostics is reported. Treatment is not, but the
     * output may be used to determine whether any treatment took place. */
    CMDTOut exec( CMHostData hostData ) const;

public:
    /** One treatment to be delivered via the PK/PD model. */
    struct TreatInfo {
        TreatInfo( const std::string& s, const std::string& d, double h );
        size_t schedule;        // index of the schedule
        size_t dosage;          // name of the dosage table
        double delay_h;         // delay in hours
    };
    
    Kind kind;
    // Branch targets (all owned by the decision library):
    std::vector<const CMDecisionTree*> children;    // Multiple
    const CMDecisionTree *firstLine = nullptr, *secondLine = nullptr;   // CaseType
    const CMDecisionTree *imported = nullptr, *introduced = nullptr,
                         *indigenous = nullptr;     // InfectionOrigin
    const CMDecisionTree *positive = nullptr, *negative = nullptr;  // several
    std::map<double, const CMDecisionTree*> branches;   // Random, Age
    const WithinHost::Diagnostic* diagnostic = nullptr; // Diagnostic
    SimTime memory = sim::zero();                 // Uncomplicated
    std::vector<int> outIds;                      // Report
    std::vector<TreatInfo> treatments;            // TreatPKPD
    std::vector<SimTime> timeLiver, timeBlood;    // TreatSimple
    std::optional<interventions::ComponentId> component;  // Cohort
    std::optional<interventions::HumanIntervention> deployment; // Deploy
};

} }
#endif
