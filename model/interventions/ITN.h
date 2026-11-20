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

#include "Context.h"
#include "interventions/ComponentId.h"
#include "util/DecayFunction.h"
#include "util/sampler.h"
#include "util/CommandLine.h"
#include <cmath>
#include "R_nmath/qnorm.h"
#include "util/random.h"
#include "util/SpeciesIndexChecker.h"
#include "util/errors.h"
namespace scnXml
{
    class ITNDescription;
    class AnophelesParams1;
    class ITNKillingEffect;
    class ITNEffectLogit;
    class ITNDeterrency;
    class TwoStageDeterrency;
}

namespace OM {
namespace interventions {
    using util::DecayFunction;
    using util::NormalSampler;
    using util::Sampler;
    using util::LocalRng;


// —————  utility classes (internal use only)  —————

namespace factors {
    class SurvivalFactor {
    public:
        /// Set parameters.
        /// 
        /// It is checked that parameters lie in a suitible range, giving a
        /// survival factor between 0 and 1.
        /// 
        /// @param elt  Element from XML
        /// @param maxInsecticide
        /// @param eltName  Element name (used in error messages)
        /// @param isDeterrent True when used to describe a deterrency effect, false when used to
        ///     describe a killing effect.
        void init(const scnXml::ITNKillingEffect& elt, double maxInsecticide,
                    const char* eltName, bool attractivenessConstraints);
        
        /// Variant to set paramters for the logit model.
        /// 
        /// @param elt  Element from XML
        /// @param maxHoleIndex Parameter from parent element in XML
        /// @param isDeterrent True when used to describe a deterrency effect, false when used to
        ///     describe a killing effect.
        void initLogit(const scnXml::ITNEffectLogit& elt,
                                                   double maxHoleIndex,
                                                   bool isDeterrent);
        
        /// Initialise the model to always return no effect (factor 1).
        void init1();
        
        /** Part of survival factor, used by new ITN deterrency model. */
        double rel_pAtt( double holeIndex, double insecticideContent )const;
        /** Calculate additional survival factor imposed by nets on pre-/post-
        * prandial killing. Should be bounded to [0,1] and tend to 1 as the
        * net ages. */
        double survivalFactor( double holeIndex, double insecticideContent )const;
        
    private:
        union {
            struct {
                double BF, HF, PF, IF;  // base, hole, insecticide and interaction factors
                double holeScaling, insecticideScaling;
                double invBaseSurvival; // stored for performance only; ≥1
            } a;
            struct {
                double BF, HF, PF, IF;
                double hMax;
                double invBaseEffect;   // 1 over base effect
            } b;
        };
        bool useLogitEqns;
    };

    class RelativeAttractiveness {
        enum DeterrencyModel {
            NO_MODEL, SINGLE_STAGE, TWO_STAGE, TWO_STAGE_LOGIT
        };
    public:
        RelativeAttractiveness() : model(0) {}
        
        /** Set parameters for single-stage deterrency model.
         * 
         * It is checked that input parameters lie in a range such that
         * the relative availability is always in the range (0,1] — that is,
         * the deterrent can never be perfect, but can have zero effect. */
        void initSingleStage(const scnXml::ITNDeterrency& elt, double maxInsecticide);
        
        /** Set parameters for two-stage deterrency model.
         * 
         * It is checked that input parameters lie in a range such that
         * the relative availability is always in the range (0,1] — that is,
         * the deterrent can never be perfect, but can have zero effect. */
        void initTwoStage(const scnXml::TwoStageDeterrency& elt,
                          double maxInsecticide,
                          double holeIndexMax, bool holeIndexMaxPresent = false);
        
        /** Calculate effect. Range of output is any value ≥ 0.
         * 
         * 0 implies a fully effective deterrent, 0.5 a 50% effective
         * deterrent, 1 has no effect, >1 attracts extra mosquitoes. */
        double relativeAttractiveness (double holeIndex, double insecticideContent) const;
        
    private:
        union {
            struct {
                double lHF, lPF, lIF;      // logs of hole, insecticide and interaction factors
                double holeScaling, insecticideScaling;
            } a;
            struct {
                double lPFEntering;      // log of insecticide factor
                double insecticideScalingEntering;
                double _pad;
                SurvivalFactor pAttacking;
            } b;
            struct {
                double entBaseFactor;
                double entInsecticideFactor;
                double pEnt0Inv;
                SurvivalFactor pAttacking;
            } c;
        };
        int model;
    };
}


// —————  main, public classes  —————

// ITN component parameters (phase 7.2b: plain Component payload struct; the
// HumanInterventionComponent/HumanVectorInterventionComponent class bases and
// the virtual deploy/componentType/print_details are gone — deployment is the
// deploy_component switch in model/xml/HumanInterventionComponents.h).
struct ITNParams {
public:
    ITNParams( Context& ctx, ComponentId id, const scnXml::ITNDescription& elt,
               const map< string, size_t >& species_name_map );
    
    // members are public: read by the iv_* record functions (interventions/IvState)
    
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
        
    private:
        double proportionProtected;
        double proportionUnprotected;
        factors::RelativeAttractiveness relAttractiveness;
        factors::SurvivalFactor preprandialKillingEffect;
        factors::SurvivalFactor postprandialKillingEffect;
        factors::SurvivalFactor relFecundityEffect;
    };

    NormalSampler initialInsecticide;
    unique_ptr<Sampler> holeRate;	// holes per step
    unique_ptr<Sampler> ripRate;	// rips per hole per step
    double maxInsecticide;		// maximum initial insecticide
    double ripFactor;			// factor expressing how significant rips are in comparison to holes
    unique_ptr<DecayFunction> insecticideDecay;
    unique_ptr<DecayFunction> attritionOfNets;
    vector<ITNAnopheles> species; // vector specific params
    
    // This is sparse vector: only indexes corresponding to ITN components are
    // used. No memory management.
};

// ———  implementation (was ITN.cpp)  ———




// —————  utility classes (internal use only)  —————


inline void factors::SurvivalFactor::init1() {
    useLogitEqns = false;
    a.invBaseSurvival = 1.0;
    a.BF = 0.0;
    a.HF = 0.0;
    a.PF = 0.0;
    a.IF = 0.0;
    a.insecticideScaling = 0.0;
    a.holeScaling = 0.0;
}

inline double factors::SurvivalFactor::rel_pAtt( double holeIndex, double insecticideContent )const {
    if (!useLogitEqns) {
        const double holeComponent = exp(-holeIndex * a.holeScaling);
        const double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        const double pAtt = a.BF
                + a.HF * holeComponent
                + a.PF * insecticideComponent
                + a.IF * holeComponent * insecticideComponent;
        assert( pAtt >= 0.0 );
        return pAtt / a.BF;
    } else {
        const double holeComponent = min(log(holeIndex + 1.0), b.hMax);
        const double insecticideComponent = log(insecticideContent + 1.0);
        const double x = exp(b.BF
                + b.HF * holeComponent
                + b.PF * insecticideComponent
                + b.IF * holeComponent * insecticideComponent);
        const double pAtt = x / (x + 1.0);
        
        // By construction, this is non-negative:
        return pAtt * b.invBaseEffect;
    }
}

inline double factors::SurvivalFactor::survivalFactor( double holeIndex, double insecticideContent )const {
    if (!useLogitEqns) {
        double holeComponent = exp(-holeIndex * a.holeScaling);
        double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        double killingEffect = a.BF
                + a.HF * holeComponent
                + a.PF * insecticideComponent
                + a.IF * holeComponent * insecticideComponent;
        
        double survivalFactor = (1.0 - killingEffect) * a.invBaseSurvival;
        assert( killingEffect <= 1.0 );
        // survivalFactor might be out of bounds due to precision error, see #49
        if (survivalFactor < 0.0)
            return 0.0;
        else if (survivalFactor > 1.0)
            return 1.0;
        return survivalFactor;
    } else {
        const double holeComponent = min(log(holeIndex + 1.0), b.hMax);
        const double insecticideComponent = log(insecticideContent + 1.0);
        const double x = exp(b.BF
                + b.HF * holeComponent
                + b.PF * insecticideComponent
                + b.IF * holeComponent * insecticideComponent);
        // K = x/(x+1), so (1-K) = 1/(x+1):
        const double surv = 1.0 / (x + 1.0);
        
        return surv * b.invBaseEffect;
    }
}


inline double factors::RelativeAttractiveness::relativeAttractiveness (
        double holeIndex, double insecticideContent) const
{
    if (model == SINGLE_STAGE) {
        double holeComponent = exp(-holeIndex * a.holeScaling);
        double insecticideComponent = 1.0 - exp(-insecticideContent * a.insecticideScaling);
        double relAvail = exp(a.lHF * holeComponent
                + a.lPF * insecticideComponent
                + a.lIF * holeComponent * insecticideComponent );
        assert( relAvail>=0.0 );
        return relAvail;
    }
    
    double factor = numeric_limits<double>::quiet_NaN();
    if (model == TWO_STAGE) {
        // This is essentially a combination of the relative attractiveness as used
        // by IRS and a killing factor.
        
        // Note that an alternative, simpler, model could have been used, but was
        // not for consistency with other models. Alternative (here we don't take
        // the logarithm of PF):
        // pEnt = 1 - PFEntering × insecticideComponent
        
        const double insecticideComponent = 1.0 - exp(-insecticideContent * b.insecticideScalingEntering);
        const double pEnt = exp(b.lPFEntering * insecticideComponent);
        // In this model, effect with 0 insectice, pEnt0 = exp(0) = 1, hence we don't need to
        // divide by a denominator like in the logit model:
        factor = pEnt;
    } else {
        assert (model == TWO_STAGE_LOGIT);
        const double p = log(insecticideContent + 1.0);
        // We directly take the exponential (this is exp(logit.pEnt0):
        const double q = exp(c.entBaseFactor + c.entInsecticideFactor * p);
        const double pEnt = q / (q + 1.0);
        factor = pEnt * c.pEnt0Inv;
    }
    assert( factor >= 0.0 );
    
    // Note: b.pAttacking and c.pAttacking overlap, so we can use either:
    const double rel_pAtt = b.pAttacking.rel_pAtt( holeIndex, insecticideContent );
    return factor * rel_pAtt;
}


// —————  main, public class  —————
// (phase 7.2b: deploy/componentType/print_details virtuals are gone — the
// deploy_component switch calls iv_deploy + the mon report)

} }

