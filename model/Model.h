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

#ifndef Hmod_Model
#define Hmod_Model

#include "Global.h"
#include "util/random.h"
#include "Population.h"
#include "Host/NeonatalMortality.h"
#include "Transmission/TransmissionModel.h"
#include "Clinical/ClinicalParams.h"
#include "Clinical/InfantMortality.h"
#include "interventions/InterventionManager.h"

namespace scnXml { class Entomology; class AgeGroupValues; }

namespace OM {

/** Everything the simulation owns, in construction order. Members refer to
 * earlier members; never move or copy a Model. */
struct Model {
    /** Master RNG: seeds all per-human and vector model RNGs. */
    util::MasterRng master_rng;
    /** Owns species params and the age-availability curve. */
    std::unique_ptr<Transmission::TransmissionModel> transmission;
    /** Clinical/health-system parameters (read by the human update chain;
     * reloaded by ChangeHS interventions). */
    Clinical::ClinicalParams clinical_params;
    /** Infection-incidence parameters (scenario parameters + model options). */
    Host::InfectionIncidenceParams incidence_params;
    /** Number of new infections since the last continuous report. */
    int ctsNewInfections = 0;
    /** Human hosts. */
    Population humans;
    /** Neonatal mortality state. */
    Host::NeonatalMortality neonatal;
    /** Infant mortality reporting counters. */
    Clinical::InfantMortality infant_mortality;
    /** Intervention deployment state. */
    interventions::InterventionManager interventions;

    Model(size_t popSize, uint64_t iseed, const scnXml::Entomology& entomology,
          const Clinical::ClinicalParams& clinical_params,
          const Host::InfectionIncidenceParams& incidence_params,
          const scnXml::AgeGroupValues& availabilityToMosquitoes);
    ~Model(); // defined where the deployment types are complete
};

}
#endif
