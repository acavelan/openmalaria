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

//parse the xml scenario file
//

#pragma once

#include "Clock.h"
#include "RunContext.h"
#include "util/errors.h"
#include <schema/scenario.h>
#include <fstream>
#include <string>
#include <memory>

namespace OpenMalaria {
static const int SCHEMA_VERSION = 50;

inline std::unique_ptr<scnXml::Scenario> loadScenario(std::string lXmlFile) {
    std::ifstream fileStream(lXmlFile.c_str(), ios::binary);
    if (!fileStream.good()) {
        std::string msg = "Error: unable to open " + lXmlFile;
        throw xml_scenario_error(msg);
    }
    std::unique_ptr<scnXml::Scenario> scenario = scnXml::parseScenario(fileStream);
    fileStream.close();
    const int scenarioVersion = scenario->getSchemaVersion();
    if (scenarioVersion < SCHEMA_VERSION) {
        std::cerr << "Warning: " << lXmlFile << " uses an old schema version (latest is " << SCHEMA_VERSION << ")."
                  << endl;
    } else if (scenarioVersion > SCHEMA_VERSION)
        throw xml_scenario_error("Error: new schema version unsupported");
    return scenario;
}
}

// ———  post-load checks (was xml/XMLChecker.h)  ———

namespace OpenMalaria {
class XMLChecker {
  public:
    XMLChecker() {}
    ~XMLChecker() {}

    /*
            * Should be called after the input XML has been loaded and before all subsequent
            * initialisation. I.e. well before the simulation starts running.
            *
            * Performs custom checks on the input XML which are not enforced in the schema itself.
            *
            * It would be better to handle as many such issues as possible in the schema itself
            * However, some required checks are not possible to do using the version of the XSD
            * spec supported by the XML validation library used by OpenMalaria.
            *
            * Throws iff a check fails.
            *
            * The purpose of these XML checks to:
            *  - identify certain problems in the input XML as early as possible, and
            *  - enable the user to obtain a more informative error message.
            */
    void PerformPostValidationChecks(const scnXml::Scenario& scenario) { CheckModelOptionsAndParams(scenario); }

  private:
    /*
            * Verifies that, if no model name is written in input XML, then both parameters and
            * model options and written explicitly.
            */
    void CheckModelOptionsAndParams(const scnXml::Scenario& scenario) {
        // For each relevant element, determine whether it's specified explicitly in XML.
        const bool modelName = scenario.getModel().getModelName().present();
        const bool parameters = scenario.getModel().getParameters().present();
        const bool modelOptions = scenario.getModel().getModelOptions().present();
        if (!modelName && (!parameters || !modelOptions)) {
            throw xml_scenario_error(
                "If a model name is not specified then both <ModelOptions> and <parameters> must be specified");
        }
        if (parameters) {
            // Same rule ModelOptions::compatCheck() applies later, checked before the
            // clock is initialised (a bogus interval would divide by zero there).
            const int interval = scenario.getModel().getParameters().get().getInterval();
            if (interval != 1 && interval != 5) {
                throw xml_scenario_error("Time step set to " + std::to_string(interval) +
                                         " days but only 1 and 5 days are supported.");
            }
        }
    }
};
}

// ———  shared latentP read (was xml/latentP.h; defined in xml/scenario.cpp)  ———

namespace OpenMalaria {

/// Default 15d, else model/parameters/latentP; wraps parse errors as
/// "model/parameters/latentP: <msg>" exactly as both callers did.
SimTime read_latentP(const SimClock& clock, const scnXml::Model& model);

} // namespace OpenMalaria
