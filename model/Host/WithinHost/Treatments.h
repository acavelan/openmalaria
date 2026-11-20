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

#ifndef Hmod_WithinHost_Treatments
#define Hmod_WithinHost_Treatments

#include "Clock.h"
#include "core/RunContext.h"
#include "interventions/Interfaces.h"

using namespace std;

namespace scnXml{
    class TreatmentOption;
}

namespace OM {
namespace WithinHost {

/**
 * Type used to select a treatment option.
 * 
 * Pass by value; it just hides an integer.
 * 
 * Note: this struct and the Treatments class offer a sub-set of the
 * functionality offered by CMDecisionTree, and thus is technically redundant.
 */
struct TreatmentId{
    inline bool operator==( const TreatmentId that ){ return id == that.id; }
    inline bool operator!=( const TreatmentId that ){ return id != that.id; }
    
    /// Default constructor: construct to an initial value. Don't pass this value to WHInterface::treatment()!
    TreatmentId() : id(numeric_limits<uint32_t>::max()) {}
    explicit TreatmentId(uint32_t id) : id(id) {}
    uint32_t id;
};

/**
 * Objects of this class describe effects of a treatment, after selecting a
 * compliance/adherence/... group.
 * 
 * The static part provides a way to find objects.
 * 
 * For use within WithinHost only.
 */
class Treatments : public interventions::TriggeredDeployments {
public:
    /// @brief Static methods
    //@{
    /** Configure a new treatment option, and return the code used to select
     * that option later. */
    static TreatmentId addTreatment( RunContext& ctx, const scnXml::TreatmentOption& desc );
    
    /** Return the corresponding treatment description. */
    static const Treatments& select( RunContext& ctx, TreatmentId treatId );
    //@}
    
    /// @brief Types used by the non-static methods
    //@{
    /// Stages effected
    enum Stages{
        NONE /*i.e. no effect*/,
        LIVER = 1,
        BLOOD = 2,
        BOTH = LIVER | BLOOD
    };
    //@}
    
    /// @brief Non-static methods
    //@}

    /// Liver/blood stage action: 0 = no action, <0 = retrospective action,
    /// n>0 = treatment for next n time steps.
    SimTime timeLiver, timeBlood;

private:
    // non-static:
    Treatments( RunContext& ctx, const scnXml::TreatmentOption& elt );
};

}
}
#endif
