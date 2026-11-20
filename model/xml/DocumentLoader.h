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

#ifndef Hmod_util_DocumentLoader
#define Hmod_util_DocumentLoader

#include "Clock.h"
#include "util/errors.h"
#include <schema/scenario.h>
#include <fstream>
#include <string>
#include <memory>

namespace OM
{ 
    namespace util
    {
        static const int SCHEMA_VERSION = 50;

        inline std::unique_ptr<scnXml::Scenario> loadScenario(std::string lXmlFile)
        {
            std::ifstream fileStream(lXmlFile.c_str(), ios::binary);
            if (!fileStream.good())
            {
                std::string msg = "Error: unable to open " + lXmlFile;
                throw util::xml_scenario_error(msg);
            }
            std::unique_ptr<scnXml::Scenario> scenario = scnXml::parseScenario (fileStream);
            fileStream.close ();
            const int scenarioVersion = scenario->getSchemaVersion();
            if (scenarioVersion < SCHEMA_VERSION) {
                std::cerr<<"Warning: "<<lXmlFile<<" uses an old schema version (latest is "<<SCHEMA_VERSION<<")."<<endl;
            }
            else if (scenarioVersion > SCHEMA_VERSION)
                throw util::xml_scenario_error ("Error: new schema version unsupported");
            return scenario;
        }
    }
}

#endif
