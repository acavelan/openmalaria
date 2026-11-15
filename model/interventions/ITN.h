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

#ifndef OM_INTERVENTIONS_ITN
#define OM_INTERVENTIONS_ITN

#include "util/DecayFunction.h"
#include "Transmission/PerHost.h"
#include "util/sampler.h"
#include "interventions/InterventionFactors.h"
namespace scnXml {
    class ITNDescription;
    class AnophelesParams1;
}

namespace OM {
namespace interventions {
    using util::DecayFunction;
    using util::NormalSampler;
    using util::LognormalSampler;
    using util::LocalRng;


// —————  main, public classes  —————

class ITNComponent : public Transmission::HumanVectorInterventionComponent {
public:
    ITNComponent( ComponentId id, const scnXml::ITNDescription& elt,
               const map< string, size_t >& species_name_map );
    
    virtual void deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits, SimTime /*now*/ )const;
    
    virtual Component::Type componentType() const;
    
    virtual void print_details( std::ostream& out )const;
    
public:
    /** Per mosquito-species parameters for extended ITN model. */
    class ITNAnopheles {
    public:
        ITNAnopheles() :
            proportionProtected( numeric_limits<double>::signaling_NaN() ),
            proportionUnprotected( numeric_limits<double>::signaling_NaN() )
        {}
        void init(const scnXml::AnophelesParams1& elt,
                  double proportionUse, double maxInsecticide);
        
        /// Get deterrency. See ComponentParams::effect for a more detailed description.
        /// Range: ≥0 where 0=fullly deter, 1=no effect, >1 = attract
        inline double relativeAttractiveness( double holeIndex, double insecticideContent )const{
            return byProtection( relAttractiveness.relativeAttractiveness( holeIndex, insecticideContent ) );
        }
        /// Get killing effect on mosquitoes before feeding.
        /// See ComponentParams::effect for a more detailed description.
        inline double preprandialSurvivalFactor( double holeIndex, double insecticideContent )const{
            return byProtection( preprandialKillingEffect.survivalFactor( holeIndex, insecticideContent ) );
        }
        /// Get killing effect on mosquitoes after they've eaten.
        /// See ComponentParams::effect for a more detailed description.
        inline double postprandialSurvivalFactor( double holeIndex, double insecticideContent )const{
            return byProtection( postprandialKillingEffect.survivalFactor( holeIndex, insecticideContent ) );
        }
        /// Get mosquito fecundity multiplier
        inline double relFecundity( double holeIndex, double insecticideContent )const{
            return byProtection( relFecundityEffect.survivalFactor( holeIndex, insecticideContent ) );
        }
        
        /// Return x*proportionProtected + proportionUnprotected
        inline double byProtection(double x) const{
            return x*proportionProtected + proportionUnprotected;
        }
        
    public:
        double proportionProtected;
        double proportionUnprotected;
        factors::RelativeAttractiveness relAttractiveness;
        factors::SurvivalFactor preprandialKillingEffect;
        factors::SurvivalFactor postprandialKillingEffect;
        factors::SurvivalFactor relFecundityEffect;
    };

    NormalSampler initialInsecticide;
    unique_ptr<LognormalSampler> holeRate;	// holes per step
    unique_ptr<LognormalSampler> ripRate;	// rips per hole per step
    double maxInsecticide;		// maximum initial insecticide
    double ripFactor;			// factor expressing how significant rips are in comparison to holes
    unique_ptr<DecayFunction> insecticideDecay;
    unique_ptr<DecayFunction> attritionOfNets;
    vector<ITNAnopheles> species; // vector specific params
};

} }

#endif
