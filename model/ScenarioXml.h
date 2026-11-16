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
#include "ModelParameters.h"
#include "util/DecayFunction.h"

namespace scnXml {
class Monitoring;
class Mosq;
class Entomology;
class VectorIntervention;
class VectorTrap;
class NonHumanHostsSpeciesIntervention;
class NonHumanHostsIntervention;
class NonHumanHosts2;
class Demography;
class Clinical;
class DecayFunction;
class AgeGroupValues;
class SampledValueCV;
class SampledValueLN;
class SampledValueN;
}
namespace OM {
class Parameters;
namespace Transmission { class VectorModel; class TransmissionModel; }
namespace Transmission::Anopheles { struct NhhIntervention; struct NhhDeploymentParameters; }
namespace mon::Continuous { struct State; }
namespace xml {
void initContinuous(mon::Continuous::State&, const scnXml::Monitoring&,
    const string& filename, bool isCheckpoint);
unique_ptr<Transmission::TransmissionModel> createTransmission(const scnXml::Entomology&,
    int populationSize, Transmission::HostTransmissionParameters&, util::MasterRng&);
void initVectorIntervention(Transmission::VectorModel&, const scnXml::VectorIntervention&, size_t instance);
void initVectorTrap(Transmission::VectorModel&, const scnXml::VectorTrap&, size_t instance);
void initNhhIntervention(Transmission::VectorModel&, const scnXml::NonHumanHostsIntervention&, size_t instance);
vector<Transmission::Anopheles::NhhDeploymentParameters> readNhhDeployment(
    const scnXml::NonHumanHosts2&, const map<string, size_t>& speciesIndex);
Transmission::Anopheles::NhhIntervention readNhhIntervention(
    const scnXml::NonHumanHostsSpeciesIntervention&, const scnXml::DecayFunction&);
vector<DemographyAgeGroup> readDemography(const scnXml::Demography&);
WithinHost::Diagnostic readNeonatalDiagnostic(const scnXml::Clinical&);
util::DecayFunction readDecay(const scnXml::DecayFunction&, const char* name);
Host::IncidenceConfiguration readIncidenceConfiguration(const Parameters&);
Transmission::HostSpeciesConfiguration readHostSpeciesConfiguration(const scnXml::Mosq&);
util::DistributionParameters readDistribution(const scnXml::SampledValueCV&);
util::LognormalSampler readSampler(const scnXml::SampledValueLN&);
util::NormalSampler readNormalSampler(const scnXml::SampledValueN&);
util::AgeGroupValues readAgeGroups(const scnXml::AgeGroupValues&);
void setAgeGroups(util::AgeGroupInterpolator&, const scnXml::AgeGroupValues&, const char* name);
}
}
