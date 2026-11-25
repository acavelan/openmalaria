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

#ifndef OM_INTERVENTIONS_GVI
#define OM_INTERVENTIONS_GVI

#include "core/RunContext.h"
#include "Transmission/VectorIntervention.h"
#include "interventions/Interfaces.h"
#include "interventions/IvState.h"
#include "util/DecayFunction.h"
#include "util/sampler.h"
namespace scnXml
{
    class GVIDescription;
    class AnophelesParams11;
}

namespace OM { namespace interventions {
    using util::DecayFunction;
    using util::LocalRng;

/** Constant parameters for generic vector intervention model. */
class GVIComponent : public Transmission::HumanVectorInterventionComponent {
public:
    /** Initialise parameters
     * 
     * @param elt Component description from XML
     * @param species_name_map Map of species names to indices.
     */
    GVIComponent( RunContext& ctx, ComponentId, const scnXml::GVIDescription& elt,
               const map<string,size_t>& species_name_map );
    
    void deploy( RunContext& ctx, Host::Human& human, mon::Deploy::Method method, VaccineLimits )const;
    
    virtual Component::Type componentType() const;
    
    virtual void print_details( std::ostream& out )const;
    
    // parameters are read by the iv_* record functions (interventions/IvState)
    friend void iv_deploy( RunContext&, Host::Human&, const Transmission::HumanVectorInterventionComponent& );
    friend double iv_rel_attractiveness( RunContext&, const IvRecord&, size_t );
    friend double iv_preprandial_survival( RunContext&, const IvRecord&, size_t );
    friend double iv_postprandial_survival( RunContext&, const IvRecord&, size_t );
    friend double iv_rel_fecundity( RunContext&, const IvRecord&, size_t );
    friend void iv_checkpoint( RunContext&, std::vector<IvRecord>&, istream& );
    
private:
    /** Per mosquito-species parameters for generic vector intervention model. */
    class GVIAnopheles {
    public:
        GVIAnopheles() :
            proportionProtected( numeric_limits<double>::signaling_NaN() ),
            proportionUnprotected( numeric_limits<double>::signaling_NaN() ),
            deterrency( numeric_limits<double>::signaling_NaN() ),
            preprandialKilling( numeric_limits<double>::signaling_NaN() ),
            postprandialKilling( numeric_limits<double>::signaling_NaN() ),
            fecundityReduction( numeric_limits<double>::signaling_NaN() )
        {}
        void init(const scnXml::AnophelesParams11& elt, double proportionUse);
        
        /// Return x*proportionProtected + proportionUnprotected
        inline double byProtection(double x) const{
            return x*proportionProtected + proportionUnprotected;
        }
        
        double proportionProtected;
        double proportionUnprotected;
        double deterrency;
        double preprandialKilling;
        double postprandialKilling;
        double fecundityReduction;
        
        friend double iv_rel_attractiveness( RunContext&, const IvRecord&, size_t );
        friend double iv_preprandial_survival( RunContext&, const IvRecord&, size_t );
        friend double iv_postprandial_survival( RunContext&, const IvRecord&, size_t );
        friend double iv_rel_fecundity( RunContext&, const IvRecord&, size_t );
    };
    
    unique_ptr<DecayFunction> decay;
    vector<GVIAnopheles> species;  // vector specific params
    
    // This is sparse vector: only indexes corresponding to a GVI component are used
    // No memory management
};

} }

#endif
