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

/** Deployment descriptions: plain data.
 *
 * What a deployment does is decided while reading the XML and never changes,
 * so this replaces the old TimedDeployment class hierarchy with a variant of
 * per-action data structs. Behaviour lives in the free functions declared at
 * the bottom (defined in InterventionManager.cpp). */

#ifndef OM_INTERVENTIONS_DEPLOYMENTS
#define OM_INTERVENTIONS_DEPLOYMENTS

#include "Global.h"
#include "interventions/Interfaces.h"
#include "Transmission/TransmissionModel.h"
#include "schema/interventions.h"

#include <map>
#include <memory>
#include <variant>

namespace OM {
    namespace Host {
        class Human;
    }
    namespace Transmission {
        class VectorSpeciesParams;
    }

namespace interventions {

/// Part shared by all deployments targeting humans.
struct HumanDeployment {
    HumanDeployment(const scnXml::DeploymentBase& deploy,
                    std::shared_ptr<const HumanIntervention> intervention,
                    ComponentId subPop, bool complement);

    bool copula = false;
    double coverage;    // proportion coverage within group meeting restrictions
    double coverageCorr = 0.0, coverageVar = 0.0;
    VaccineLimits vaccLimits;
    ComponentId subPop;  // ComponentId::wholePop() if deployment isn't restricted
    bool complement;
    std::shared_ptr<const HumanIntervention> intervention;
};

/** An intervention deployed at a fixed time: when, and what it does. */
struct TimedDeployment {
    struct Human {
        Human(const scnXml::MassDeployment& mass,
              std::shared_ptr<const HumanIntervention> intervention,
              ComponentId subPop, bool complement);
        HumanDeployment base;
        SimTime minAge, maxAge;
        SimTime minAvailability, maxAvailability;  // mosquito-availability percentiles
    };
    struct CumulativeHuman {  // brings the target group up to coverage
        CumulativeHuman(const scnXml::MassDeployment& mass,
                        std::shared_ptr<const HumanIntervention> intervention,
                        ComponentId subPop, bool complement, ComponentId cumCovId);
        Human human;
        ComponentId cumCovId;
    };
    struct ChangeHS { std::unique_ptr<scnXml::HealthSystem> newHS; };
    struct ChangeEIR { std::unique_ptr<scnXml::NonVector> newEIR; };
    struct UninfectVectors {};
    struct VectorPop { size_t instance; };
    struct Trap { size_t instance; double ratio; SimTime lifespan; };
    struct ReduceNhh { size_t instance; std::string intervName; };
    struct AddNhh {
        struct SpeciesParams {  // read from XML, keyed by mosquito species name
            double mosqRelativeAvailabilityHuman;
            double mosqProbBiting;
            double mosqProbFindingRestSite;
            double mosqProbResting;
            double hostFecundityFactor;
        };
        std::string intervName;
        SimTime lifespan = sim::never();
        std::map<std::string, SpeciesParams> params;
    };

    SimTime date = sim::never();
    std::variant<std::monostate /*list-end marker*/, ChangeHS, ChangeEIR,
                 UninfectVectors, Human, CumulativeHuman, VectorPop, Trap,
                 ReduceNhh, AddNhh> action;
};

/** Deployment once a human reaches the target age. */
struct ContinuousHumanDeployment {
    ContinuousHumanDeployment(SimTime begin, SimTime end,
                              const scnXml::ContinuousDeployment& elt,
                              std::shared_ptr<const HumanIntervention> intervention,
                              ComponentId subPop, bool complement);

    HumanDeployment base;
    SimTime begin, end;  // first time step active and first time step no-longer active
    SimTime deployAge = sim::never();
    SimTime minAvailability, maxAvailability;
};

/** Check a deployment date read from XML; warns if after the last survey.
 * Called once per deployment as it is read. */
SimTime checkDeployDate(SimTime date);

/** Do one timed deployment (now). */
void apply (TimedDeployment& dep, vector<Host::Human>& population,
            Transmission::TransmissionModel& transmission);

/** Apply non-human-host reduction effects from XML to the vector model
 * (at init; deployment itself only re-applies the decayed effect). */
void initReduceNhh (size_t instance, const std::string& intervName,
                    const scnXml::Description2::AnophelesSequence list,
                    const scnXml::DecayFunction& decay,
                    Transmission::TransmissionModel& transmission);

/** Read added-non-human-host params from XML (at init). */
TimedDeployment::AddNhh readAddNhh (const std::string& intervName, SimTime lifespan,
                                    const scnXml::Description3::AnophelesSequence list,
                                    Transmission::TransmissionModel& transmission);

/** Apply filters and potentially deploy.
 *
 * @returns false iff this deployment (and thus all later ones in the
 *  ordered list) happens in the future for this human. */
bool filterAndDeploy (const ContinuousHumanDeployment& dep, Host::Human& human,
                      const Transmission::VectorSpeciesParams& species_params);

void print_details (std::ostream& out, const TimedDeployment& dep);
void print_details (std::ostream& out, const ContinuousHumanDeployment& dep);

} }
#endif
