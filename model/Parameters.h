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

#include <array>
#include <optional>
#include "util/ModelNameProvider.h"
#include "util/CommandLine.h"
#include "util/errors.h"

namespace scnXml {
class Parameters;
}
namespace xsd {
namespace cxx {
namespace tree {
template <typename T, bool> class optional; // bool spelled out in declarations below
}
}
}
namespace OpenMalaria {

/*
* Defines the names that clients of the Parameters class use in order to
* read the values of parameters.
*
* NOTE: enum order is the XML parameter order: IDs 1–28 map directly,
* deprecated ID 29 is skipped, and IDs 30–38 follow. Keep fromId() and the
* missing-value ID calculation in sync when adding or removing parameters.
*/
enum class Parameter {
    /// @b Infection incidence model parameters
    //@{
    NEG_LOG_ONE_MINUS_SINF,
    E_STAR,
    SIMM,
    X_STAR_P,
    GAMMA_P,
    //@}
    /// @b Immunity parameters, mostly on infections
    //@{
    SIGMA_I_SQ, ///< Host (not infection) parameter
    CUMULATIVE_Y_STAR,
    CUMULATIVE_H_STAR,
    NEG_LOG_ONE_MINUS_ALPHA_M,
    DECAY_M,
    //@}
    /// @b DescriptiveInfection specific
    //@{
    SIGMA0_SQ,
    X_NU_STAR,
    //@}
    /// @b Used in PathogenesisModel
    //@{
    Y_STAR_SQ,
    ALPHA,
    //@}
    DENSITY_BIAS_NON_GARKI,       ///< Used in Diagnostic
    BASELINE_AVAILABILITY_SHAPE,  ///< Used in InfectionIncidenceModel
    LOG_ODDS_RATIO_CF_COMMUNITY,  ///< Used in CaseManagementModel
    INDIRECT_RISK_COFACTOR,       ///< Used in PathogenesisModel
    NON_MALARIA_INFANT_MORTALITY, ///< Used in Summary
    DENSITY_BIAS_GARKI,           ///< Used in Diagnostic
    SEVERE_MALARIA_THRESHHOLD,    ///< Used in PathogenesisModel
    IMMUNITY_PENALTY,             ///< Used in WHFalciparum
    IMMUNE_EFFECTOR_DECAY,        ///< Used in WHFalciparum
    /// @b Used in PathogenesisModel
    //@{
    COMORBIDITY_INTERCEPT,
    Y_STAR_HALF_LIFE,
    Y_STAR_1,
    //@}
    ASEXUAL_IMMUNITY_DECAY, ///< Used in WHFalciparum
    /// @b Used in PathogenesisModel
    //@{
    Y_STAR_0,

    CRITICAL_AGE_FOR_COMORBIDITY,
    MUELLER_RATE_MULTIPLIER,
    MUELLER_DENSITY_EXPONENT,
    //@}
    /// EventScheduler: v in "Case Fatality Rate proposal"
    CFR_SCALE_FACTOR,

    /// @b Molineaux: sampling parameters (not pairwise mode only)
    MEAN_LOCAL_MAX_DENSITY,
    SD_LOCAL_MAX_DENSITY,
    MEAN_DIFF_POS_DAYS,
    SD_DIFF_POS_DAYS,

    /// EventScheduler: exp(-CFR_NEG_LOG_ALPHA) is the proportion of deaths occuring on the first day
    CFR_NEG_LOG_ALPHA,
    NUM_PARAMETERS,
};

/*
* This class' job is to make the values of model parameters available to clients when the simulation
* is running.
*
* From a user perspective, each parameter has these 3 parts:
*
* 1. The name of the parameter as written in XML.
*       OpenMalaria completely disregards this.
*
* 2. The numerical ID, "index" which identfies the parameter in XML.
*       OpenMalaria relies on this.
*
* 3. The value assigned to the parameter in XML.
*       OpenMalaria associates this with the numerical ID.
*
* From a developer perspective, each "parameter" can be thought of as having these 3 parts:
*
* 1. A parameter name (see Parameter enum class).
*       Clients of this class use a parameter name to refer to the parameter.
*
* 2. A numerical ID to identify the parameter.
*       This is used to map parameters in XML to parameters in OpenMalaria,
*       and for error reporting.
*
* 2. A (possible) floating point number representing the actual value assigned to the parameter.
*       If the user does not explicitly or implicitly set a value, an empty
*       value is assigned.  If a client of this class attempts to read the
*       parameter's value, but the value is empty, this is considered an
*       error - and this class is responsible for handling such errors.
*/
class Parameters {
  public:
    /*
    * Expects that the scenario XML either explicitly describes a collection of parameters and values
    * or describes the name of a model to use.
    */
    Parameters(const xsd::cxx::tree::optional<scnXml::Parameters, false>& xmlParameterValues, ModelNames modelName);

    /**
     * Get a parameter, using one of the Parameter codes.
     */
    double operator[](Parameter name) const {
        const auto& value = values.at(static_cast<size_t>(name));
        if (value)
            return *value;

        int id = static_cast<int>(name) + 1;
        if (id >= 29)
            ++id; // XML ID 29 is deprecated and absent from the enum.
        throw xml_scenario_error("Parameter with index " + to_string(id) + " required but not described.");
    }

  private:
    /*
    * Initializes parameters using explicit values specified in the input XML,
    * and does some validation on the specified parameters and values.
    * This method will be used if the user does not use the "base" model or any
    * other pre-set collection of parameter values.
    */
    void initializeParamsFromXML(const scnXml::Parameters& parameters);

    static std::optional<Parameter> fromId(int id) {
        if (id < 1 || id > 38 || id == 29)
            return std::nullopt;
        return static_cast<Parameter>(id - (id > 29 ? 2 : 1));
    }

    /*
    * Initializes some hardcoded values for some parameters, encapsulating the parameter values
    * that make up the base model.
    */
    void initializeParamsBaseModel() {
        static constexpr std::array<double, 29> defaults{
            0.050736,                // '-ln(1-Sinf)'
            0.03247,                 // Estar
            0.138161050830301,       // Simm
            1514.385853233699891,    // Xstar_p
            2.03692533424484,        // gamma_p
            10.173598698525799,      // sigma2i
            35158523.31132510304451, // CumulativeYstar
            97.334652723897705,      // CumulativeHstar
            2.33031045876193,        // '-ln(1-alpha_m)'
            2.53106547375805,        // decay_m
            0.655747311168152,       // sigma2_0
            0.916181104713054,       // Xstar_v
            6502.26335600001039,     // Ystar2
            142601.912520000012591,  // alpha
            0.177378570987455,       // Density bias (non Garki)
            1.0,                     // sigma2
            0.736202,                // log oddsr CF community
            0.018777338,             // Indirect risk cofactor
            49.539046599999999,      // Non-malaria infant mortality
            4.79610772546704,        // Density bias (Garki)
            784455.599999999976717,  // Severe Malaria Threshhold
            1,                       // Immunity Penalty
            0,                       // Immune effector decay
            0.0968,                  // comorbidity intercept
            0.275437402,             // Ystar half life
            0.596539864,             // Ystar1
            0,                       // Asexual immunity decay
            296.302437899999973,     // Ystar0
            0.117383,                // critical age for comorbidity (XML ID 30)
        };
        for (size_t i = 0; i < defaults.size(); ++i)
            values[i] = defaults[i];
    }

    std::array<std::optional<double>, static_cast<size_t>(Parameter::NUM_PARAMETERS)> values;
};

}
