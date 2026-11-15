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

/** \file XML parsing for Parameters; the header stays XML-free. */

#include "Parameters.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include <schema/scenario.h>
#include <iostream>
#include <set>

namespace OM {

Parameters::Parameters( const scnXml::Model& model, util::ModelNameProvider mnp )
{
    for (const std::pair<int, Parameter> pair : idCodeToNameMap)
    {
        const Parameter name = pair.second;
        nameToValueMap[name] = std::nullopt;
    }

    const util::ModelNames namedModelToUse = mnp.GetModelName();
    const bool useNamedModel = namedModelToUse != util::ModelNames::none;
    if (useNamedModel)
    {
        if (namedModelToUse == util::ModelNames::base)
        {
            initializeParamsBaseModel();
        }
        else
        {
            throw util::xml_scenario_error(
                "No pre-set parameter values are available for the specified model name.");
        }
    }

    // Get parameters specified expicitly in input XML, if any.
    // In the event where a named model was used to set parameter values earlier, it is intended
    // that any parameter values written explicitly in XML will override any values set by the
    // named model.
    const scnXml::Model::ParametersOptional& xmlParameterValues = model.getParameters();
    const bool useExplicitParamValues = xmlParameterValues.present();
    if (useExplicitParamValues)
    {
        initializeParamsFromXML(xmlParameterValues.get());
    }
}

/*
* Initializes parameters using explicit values specified in the input XML,
* and does some validation on the specified parameters and values.
* This method will be used if the user does not use the "base" model or any
* other pre-set collection of parameter values.
*/
void Parameters::initializeParamsFromXML(const scnXml::Parameters& parameters)
{
    const scnXml::Parameters::ParameterSequence& paramSeq = parameters.getParameter();

    // It's okay if a user overwrites the value of a parameter that was set before this
    // method was called.  E.g. if they used a named model, and then chose to manually
    // override the value of a parameter set by it.  It's not okay if a user sets a value
    // for the same parameter twice.  That represents a mistake in the input XML which
    // we need to handle here.  This has the potential to save users time debugging inputs.
    std::set<int> paramIdsSetByUser;

    for (auto param : paramSeq)
    {
        const int paramId = param.getNumber();
        const double paramValueFromXML = param.getValue();

        // C++20 offers "contains" which would be more expressive than "count" here.
        const bool paramIdIsValid = idCodeToNameMap.count(paramId);
        if (!paramIdIsValid)
        {
            if( util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS) ){
                std::cerr << "Deprecation warning: <parameter> index" << std::to_string(paramId) <<
                    " is no longer used" << std::endl;
            }

            // If we throw here, old scenarios containing the deprecated parameter
            // will need more manual work to migrate to new OpenMalaria versions.
            // This isn't unacceptable.  But it's simpler here to just skip over the
            // deprecated parameter.  In particular, many tests specify values for
            // deprecated some parameters.
            continue;
        }
        const Parameter nameOfParamToSet = idCodeToNameMap.at(paramId);

        const bool paramValueAlreadySet = !paramIdsSetByUser.insert(paramId).second;
        if (paramValueAlreadySet)
        {
            throw util::xml_scenario_error("Parameter with index " + std::to_string(paramId) +
                    " described twice in XML.");
        }

        nameToValueMap[nameOfParamToSet] = paramValueFromXML;
    }
}

}
