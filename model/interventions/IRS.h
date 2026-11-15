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

#ifndef OM_INTERVENTIONS_IRS
#define OM_INTERVENTIONS_IRS

#include "util/DecayFunction.h"
#include "Transmission/PerHost.h"
#include "interventions/InterventionFactors.h"
namespace scnXml {
    class IRSDescription;
    class IRSDeterrency;
    class IRSKillingEffect;
    class AnophelesParams2;
}

namespace OM {
namespace interventions {
    using util::DecayFunction;
    using util::NormalSampler;
    using util::LognormalSampler;
    using util::LocalRng;

class IRSComponent : public Transmission::HumanVectorInterventionComponent {
public:
    IRSComponent( ComponentId id, const scnXml::IRSDescription& elt,
        const map<string,size_t>& species_name_map );
    
    virtual void deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits, SimTime /*now*/ )const;
    
    virtual Component::Type componentType() const;
    
    virtual void print_details( std::ostream& out )const;
    
public:
    /** Per mosquito-species parameters for extended IRS model. */
    class IRSAnopheles {
    public:
        IRSAnopheles() :
            proportionProtected( numeric_limits<double>::signaling_NaN() ),
            proportionUnprotected( numeric_limits<double>::signaling_NaN() )
        {}
        void init(const scnXml::AnophelesParams2& elt,
                 double proportionUse, double maxInsecticide);
        
        /// Get deterrency. See ComponentParams::effect for a more detailed description.
        inline double relativeAttractiveness( double insecticideContent )const{
            return _relativeAttractiveness.relativeAttractiveness( 0.0, insecticideContent );
        }
        /// Get survival effect on mosquitoes before feeding.
        /// See ComponentParams::effect for a more detailed description.
        inline double preprandialSurvivalFactor( double insecticideContent )const{
            return _preprandialKillingEffect.survivalFactor( 0.0, insecticideContent );
        }
        /// Get survival effect on mosquitoes after they've eaten.
        /// See ComponentParams::effect for a more detailed description.
        inline double postprandialSurvivalFactor( double insecticideContent )const{
            return _postprandialKillingEffect.survivalFactor( 0.0, insecticideContent );
        }
        /// Get fecundity effect on mosquitoes surviving feeding
        inline double fecundityEffect( double insecticideContent )const{
            return _fecundityEffect.survivalFactor( 0.0, insecticideContent );
        }
        
        /// Return x*proportionProtected + proportionUnprotected
        inline double byProtection(double x) const{
            return x*proportionProtected + proportionUnprotected;
        }
        
    public:
        double proportionProtected;
        double proportionUnprotected;
        factors::RelativeAttractiveness _relativeAttractiveness;
        factors::SurvivalFactor _preprandialKillingEffect;
        factors::SurvivalFactor _postprandialKillingEffect;
        factors::SurvivalFactor _fecundityEffect;
    };
    
    double sampleInitialInsecticide(LocalRng& rng) const;
    
    NormalSampler initialInsecticide;
    double maxInsecticide;              // maximum initial insecticide
    unique_ptr<DecayFunction> insecticideDecay;
    vector<IRSAnopheles> species; // vector specific params
};

} }

#endif
