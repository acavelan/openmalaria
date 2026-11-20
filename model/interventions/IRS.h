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
#include "interventions/ComponentId.h"
#include "util/DecayFunction.h"
#include <cmath>
#include "R_nmath/qnorm.h"
#include "util/SpeciesIndexChecker.h"
#include "util/errors.h"
namespace scnXml
{
    class IRSDescription;
    class AnophelesParams2;
    class IRSDeterrency;
    class IRSKillingEffect;
}

namespace OM {
namespace interventions {
    using util::DecayFunction;
    using util::NormalSampler;
    using util::LocalRng;

// IRS component parameters (phase 7.2b: plain Component payload struct —
// cf. ITNParams). Members public: read by the iv_* record functions.
struct IRSParams {
public:
    IRSParams( RunContext& ctx, ComponentId id, const scnXml::IRSDescription& elt,
        const map<string,size_t>& species_name_map );
    
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
            return _relativeAttractiveness.relativeAttractiveness( insecticideContent );
        }
        /// Get survival effect on mosquitoes before feeding.
        /// See ComponentParams::effect for a more detailed description.
        inline double preprandialSurvivalFactor( double insecticideContent )const{
            return _preprandialKillingEffect.survivalFactor( insecticideContent );
        }
        /// Get survival effect on mosquitoes after they've eaten.
        /// See ComponentParams::effect for a more detailed description.
        inline double postprandialSurvivalFactor( double insecticideContent )const{
            return _postprandialKillingEffect.survivalFactor( insecticideContent );
        }
        /// Get fecundity effect on mosquitoes surviving feeding
        inline double fecundityEffect( double insecticideContent )const{
            return _fecundityEffect.survivalFactor( insecticideContent );
        }
        
        /// Return x*proportionProtected + proportionUnprotected
        inline double byProtection(double x) const{
            return x*proportionProtected + proportionUnprotected;
        }
        
    private:
        class RelativeAttractiveness {
        public:
            RelativeAttractiveness();
            
            /** Set parameters.
            * 
            * It is checked that input parameters lie in a range such that
            * the relative availability is always in the range (0,1] — that is,
            * the deterrent can never be perfect, but can have zero effect. */
            void init(const scnXml::IRSDeterrency& elt);
            
            /** Calculate effect. Positive is interpreted as having a positive effect
            * (thus decreasing availability or survival) and negative as having a
            * negative effect. Effect is not bounded, though it tends to
            * zero as holeIndex becomes large and insecticideContent tends to zero,
            * and parameters should be defined such that it is always in the
            * range [0,1]. */
            double relativeAttractiveness( double insecticideContent )const;
            
        protected:
            double lPF;      // log of insecticide factor
            double insecticideScaling;
        };
        class SurvivalFactor {
        public:
            SurvivalFactor();
            
            /** Set parameters.
            * 
            * It is checked that parameters lie in a suitible range, giving a
            * survival factor between 0 and 1. */
            void init(const scnXml::IRSKillingEffect& elt, bool postPrandial, double maxInsecticide);
            /// Initialise, such that factor returned is always 1 (i.e. no effect).
            void init1();
            
            /** Calculate additional survival factor imposed by IRS on pre-/post-
            * prandial killing. Should be bounded to [0,1] and tend to 1 as the
            * IRS ages. */
            double survivalFactor( double insecticideContent )const;
        private:
            double BF,PF;   // base and insecticide factors
            double insecticideScaling;
            double invBaseSurvival; // stored for performance only
        };
        double proportionProtected;
        double proportionUnprotected;
        RelativeAttractiveness _relativeAttractiveness;
        SurvivalFactor _preprandialKillingEffect;
        SurvivalFactor _postprandialKillingEffect;
        SurvivalFactor _fecundityEffect;
    };
    
    double sampleInitialInsecticide(LocalRng& rng) const;
    
    NormalSampler initialInsecticide;
    double maxInsecticide;              // maximum initial insecticide
    unique_ptr<DecayFunction> insecticideDecay;
    vector<IRSAnopheles> species; // vector specific params
    
    // This is sparse vector: only indexes corresponding to a IRS component are used
    // No memory management
};

// ———  implementation (was IRS.cpp)  ———



// (phase 7.2b: deploy/componentType/print_details virtuals are gone — the
// deploy_component switch calls iv_deploy + the mon report)

inline IRSParams::IRSAnopheles::RelativeAttractiveness::RelativeAttractiveness() :
    lPF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{}
inline IRSParams::IRSAnopheles::SurvivalFactor::SurvivalFactor() :
    BF( numeric_limits< double >::signaling_NaN() ),
    PF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() ),
    invBaseSurvival( numeric_limits< double >::signaling_NaN() )
{}
inline void IRSParams::IRSAnopheles::SurvivalFactor::init1(){
    BF = 0.0;
    PF = 0.0;
    insecticideScaling = 0.0;
    invBaseSurvival = 1.0;
}
inline double IRSParams::IRSAnopheles::RelativeAttractiveness::relativeAttractiveness(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lPF*insecticideComponent );
    assert( relAvail>=0.0 );
    return relAvail;
}
inline double IRSParams::IRSAnopheles::SurvivalFactor::survivalFactor(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double killingEffect = BF + PF*insecticideComponent;
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    assert( survivalFactor >= 0.0 );
    assert( survivalFactor <= 1.0 );
    return survivalFactor;
}

inline double IRSParams::sampleInitialInsecticide(LocalRng& rng) const{
    double sample = initialInsecticide.sample(rng);
    if( sample < 0.0 )
        sample = 0.0;       // avoid negative samples
    if( sample > maxInsecticide )
        sample = maxInsecticide;
    return sample;
}


// ———  per-human data  ———
// (phase 7.2a: the HumanIRS class is gone; its draws live in iv_deploy and
// its effects in iv_* — interventions/IvState.cpp)

} }

