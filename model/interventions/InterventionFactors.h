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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// Shared mosquito-response curves for ITNs and IRS (the IRS model uses the
// same equations with the hole terms zero and holeIndex fixed at 0).
// Definitions live in ITN.cpp.

#ifndef fhirdway_INTERVENTION_FACTORS_H
#define fhirdway_INTERVENTION_FACTORS_H

namespace scnXml {
    class ITNKillingEffect;
    class ITNEffectLogit;
    class ITNDeterrency;
    class TwoStageDeterrency;
}

namespace OM { namespace interventions {
namespace factors {
    class SurvivalFactor {
    public:
        /// Set parameters (net model; hole factors from XML).
        void init(const scnXml::ITNKillingEffect& elt, double maxInsecticide,
                    const char* eltName, bool attractivenessConstraints);
        
        /// Variant to set paramters for the logit model.
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
        
    public:
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
        /// ITN clamps the survival factor to [0,1]; IRS returns it unclamped.
        bool clamp = true;
    };

    class RelativeAttractiveness {
    public:
        enum DeterrencyModel {
            NO_MODEL, SINGLE_STAGE, TWO_STAGE, TWO_STAGE_LOGIT
        };
        RelativeAttractiveness() : model(0) {}
        
        /** Set parameters for single-stage deterrency model (net model). */
        void initSingleStage(const scnXml::ITNDeterrency& elt, double maxInsecticide);
        
        /** Set parameters for two-stage deterrency model (net model). */
        void initTwoStage(const scnXml::TwoStageDeterrency& elt,
                          double maxInsecticide,
                          double holeIndexMax, bool holeIndexMaxPresent = false);
        
        /** Calculate effect. Range of output is any value ≥ 0.
         * 
         * 0 implies a fully effective deterrent, 0.5 a 50% effective
         * deterrent, 1 has no effect, >1 attracts extra mosquitoes. */
        double relativeAttractiveness (double holeIndex, double insecticideContent) const;
        
    public:
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
} }
#endif
