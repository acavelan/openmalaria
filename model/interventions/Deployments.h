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

// Note: this file is included by several sources (InterventionManager.cpp,
// xml/interventions.cpp and run.cpp). All definitions here must
// stay inline/template or header-only.

// Phase 7.1: the polymorphic TimedDeployment class tree is gone — one flat
// TimedDeploymentKind + variant payload, built at load, date-cursor driven.
// Processing order is unchanged: same single list, same stable_sort by date
// with insertion-order tie-break (contract b), and the deploy bodies are
// moved verbatim (contract a: the RNG draw stays the LAST filter test).

// The includes here are more for documentation than required.
#pragma once

#include "Clock.h"
#include "Context.h"
#include "util/errors.h"
#include "Population.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/PerHost.h"
#include "util/random.h"
#include <map>
#include <string>
#include <variant>
namespace scnXml {
class DeploymentBase;
class MassDeployment;
class ContinuousDeployment;
class HealthSystem;
class NonVector;
}

#include "util/SpeciesIndexChecker.h"
#include "Host/Human.h"

namespace OM {
namespace interventions {

using Host::Human;

/** Flat replacement of the TimedDeployment subclass tree: kind selects the
 * variant payload and the deploy/print switch case (no per-timestep clones). */
enum class TimedDeploymentKind : uint8_t {
    DUMMY,            // after-last sentinel; never deploys
    HUMAN,            // timed mass deployment to a human intervention
    CUMULATIVE_HUMAN, // bring coverage of cumCovId up to target
    UNINFECT_VECTORS,
    VECTOR, // vector-population intervention instance
    TRAP,
    ADD_NHH,  // add non-human hosts
    NHH_MODS, // modify non-human hosts
    CHANGE_HS,
    CHANGE_EIR
};

/// True for the kinds re-applied at checkpoint restore (they mutate model
/// state, unlike per-host deployments; was TimedDeployment::replayedOnRestore).
inline bool replayed_on_restore(TimedDeploymentKind kind) {
    return kind == TimedDeploymentKind::CHANGE_HS || kind == TimedDeploymentKind::CHANGE_EIR;
}

/** Old TimedDeployment ctor date checks; run at load, before payload reading
 * (same throw/warn order as the old base ctor). */
inline void check_timed_date(Context& ctx, SimTime deployDate) {
    if (deployDate < ctx.clock.start_date) {
        throw util::xml_scenario_error("timed intervention deployment: may not be negative");
    } else if (deployDate >= ctx.clock.end_date) {
        cerr << "Warning: timed intervention deployment at date " << deployDate;
        cerr << " happens after last survey" << endl;
    }
}

/** Deployment restrictions shared by timed and continuous human deployments
 * (old HumanDeploymentBase; its Context member is gone — the deploy entry
 * points already receive ctx). read() is defined in xml/interventions.cpp. */
struct HumanDeployParams {
    static HumanDeployParams read(Context& ctx, const scnXml::DeploymentBase& deploy,
                                  shared_ptr<const HumanIntervention> intervention, ComponentId subPop,
                                  bool complement);

    inline void deployToHuman(Context& ctx, Host::Human& human, mon::Deploy::Method method) const {
        intervention->deploy(ctx, human, method, vaccLimits);
    }

    bool copula = false;
    double coverage = 0.0; // proportion coverage within group meeting above restrictions
    double coverageCorr = 0.0, coverageVar = 0.0;
    VaccineLimits vaccLimits;
    ComponentId subPop = ComponentId::wholePop(); // wholePop() if not sub-population-restricted
    bool complement = false;
    shared_ptr<const HumanIntervention> intervention;
};

/// Old TimedHumanDeployment filter bounds (read from XML at load).
struct TimedHumanParams {
    static TimedHumanParams read(Context& ctx, SimTime date, const scnXml::MassDeployment& mass,
                                 shared_ptr<const HumanIntervention> intervention, ComponentId subPop, bool complement);

    HumanDeployParams d;
    SimTime minAge = sim::never(), maxAge = sim::never();
    SimTime minAvailability = sim::never(), maxAvailability = sim::never();
};

/// Old TimedCumulativeHumanDeployment: coverage target against cumCovId.
struct CumulativeHumanParams {
    TimedHumanParams h;
    ComponentId cumCovId = ComponentId::wholePop();
};

struct VectorDeployParams {
    size_t instance;
};
struct TrapDeployParams {
    size_t instance;
    double ratio;
    SimTime lifespan = sim::never();
};
struct NhhModsDeployParams {
    size_t instance;
    std::string intervName;
};
struct AddNhhDeployParams {
    struct NhhSpec { // old NhhParamsInterv
        double mosqRelativeAvailabilityHuman = 0.0;
        double mosqProbBiting = 0.0;
        double mosqProbFindingRestSite = 0.0;
        double mosqProbResting = 0.0;
        double hostFecundityFactor = 0.0;
    };
    std::string intervName;
    SimTime lifespan = sim::never();
    std::map<std::string, NhhSpec> bySpecies;
};
/// ChangeHS/ChangeEIR clones are consumed (deleted) on deploy. Raw pointers:
/// never-deployed clones leak at exit — exactly the old classes' behaviour
/// (their destructors never freed these).
struct ChangeHsDeployParams {
    scnXml::HealthSystem* newHS = nullptr;
};
struct ChangeEirDeployParams {
    scnXml::NonVector* newEIR = nullptr;
};

// Definitions in xml/interventions.cpp (they need the schema/clinical includes
// this header must not pull in):
void timed_deploy_change_hs(Context&, Transmission::TransmissionModel&, ChangeHsDeployParams&);
void timed_deploy_change_eir(Context&, Transmission::TransmissionModel&, ChangeEirDeployParams&);
void timed_deploy_add_nhh(Context&, std::vector<Host::Human>&, Transmission::TransmissionModel&,
                          AddNhhDeployParams&);
void timed_deploy_nhh_mods(Context&, std::vector<Host::Human>&, Transmission::TransmissionModel&,
                           const NhhModsDeployParams&);

// ———  timed human deployments  ———

/// Old TimedHumanDeployment::deploy body, verbatim (filters then bernoulli-last).
inline void timed_human_deploy(Context& ctx, std::vector<Host::Human>& population, const TimedHumanParams& th) {
    const HumanDeployParams& d = th.d;
    for (Human& human : population) {
        const SimTime age = human.age(ctx.clock.now());
        if (age >= th.minAge && age < th.maxAge) {
            double availability = 0.0;
            for (size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(ctx); ++i)
                availability += human.perHostTransmission.anoph_ento_avail()[i];

            double probability = d.coverage;

            if (d.copula && d.coverage > 0) {
                try {
                    // Beta distribution for intervention (Beta distributed)
                    const double beta_mean = d.coverage;   // Assuming this value is given
                    const double beta_var = d.coverageVar; // Given as well
                    const double alpha = ((1.0 - beta_mean) / beta_var - 1.0 / beta_mean) * (beta_mean * beta_mean);
                    const double beta = alpha * (1.0 / beta_mean - 1.0);

                    if (Transmission::PerHostAnophParams::numSpecies(ctx) > 1)
                        throw std::runtime_error("only supports one mosquito species");

                    const double ux = Transmission::PerHostAnophParams::get(ctx, 0).entoAvailability->cdf(
                        human.perHostTransmission.anoph_ento_raw()[0]);
                    if (ux < 1) {
                        // parity: a — per-human heterogeneity lognormal/beta draws (gsl_cdf pure quantiles, no stream cost) then bernoulli-last; per-candidate draw order/count frozen (plan §3a)
                        double xx = gsl_cdf_ugaussian_Pinv(ux); // Unit interval to Normal

                        // Apply the Gaussian copula transformation for correlation
                        const ComponentId cid = d.subPop;

                        /** human.rng().gauss(0.0, 1.0) is calculated once per component and per human. Re-deployment of
                         * the same component on the same human must use the same gaussian sample. Therefore we store
                         * this value in the human the first time it is calculated (copula_g column, 7.4). */
                        auto& cg = human.cols->copula_g[human.row];
                        const auto it = cg.find(cid);
                        double g = 0.0;
                        if (it != cg.end())
                            g = it->second;
                        else {
                            g = human.rng().gauss(0.0, 1.0);
                            cg[cid] = g;
                        }
                        const double yy = d.coverageCorr * xx + g * sqrt(1 - d.coverageCorr * d.coverageCorr);

                        // Transform back to unit interval
                        const double uy = gsl_cdf_ugaussian_P(yy);

                        if (alpha < 0 || beta < 0)
                            throw std::runtime_error("resulting alpha and beta parameters must be positive");

                        // Get the final probability from Beta distribution
                        probability = gsl_cdf_beta_Pinv(uy, alpha, beta);
                    } else // Gamma with CV=0
                        probability = human.rng().beta(alpha, beta);
                } catch (const std::exception& e) {
                    std::ostringstream oss;
                    oss << "[Gaussian copula]: computing correlated intervention probability using Gaussian copula and Beta distribution: "
                        << e.what()
                        << ". Possible reason: coverage might be too high and/or population size might be too small.";
                    throw std::runtime_error(oss.str());
                }
            }

            if (availability >=
                    Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, th.minAvailability) &&
                availability <=
                    Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, th.maxAvailability)) {
                if (d.subPop == ComponentId::wholePop() || (human.isInSubPop(d.subPop) != d.complement)) {
                    // parity: a — deploy short-circuit: RNG draw is the LAST test of the deploy condition; filter/age/cohort checks must precede it (plan §3a)
                    if (human.rng().bernoulli(probability)) {
                        d.deployToHuman(ctx, human, mon::Deploy::TIMED);
                    }
                }
            }
        }
    }
}

/// Old TimedCumulativeHumanDeployment::deploy body, verbatim.
inline void timed_cumulative_deploy(Context& ctx, std::vector<Host::Human>& population,
                                    const CumulativeHumanParams& cum) {
    const TimedHumanParams& th = cum.h;
    const HumanDeployParams& d = th.d;
    // Cumulative case: bring target group's coverage up to target coverage
    std::vector<Host::Human*> unprotected;
    size_t total = 0; // number of humans within age bound and optionally subPop
    for (Host::Human& human : population) {
        SimTime age = human.age(ctx.clock.now());
        if (age >= th.minAge && age < th.maxAge) {
            double availability = 0.0;
            for (size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(ctx); ++i)
                availability += human.perHostTransmission.anoph_ento_avail()[i];
            if (availability >=
                    Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, th.minAvailability) &&
                availability <=
                    Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, th.maxAvailability)) {
                if (d.subPop == ComponentId::wholePop() || (human.isInSubPop(d.subPop) != d.complement)) {
                    total += 1;
                    if (!human.isInSubPop(cum.cumCovId))
                        unprotected.push_back(&human);
                }
            }
        }
    }

    if (total == 0)
        return; // no humans to deploy to; avoid divide by zero
    double propProtected = static_cast<double>(total - unprotected.size()) / static_cast<double>(total);
    if (propProtected < d.coverage) {
        // Proportion propProtected are already covered, so need to
        // additionally cover the proportion (coverage - propProtected),
        // selected from the list unprotected.
        double additionalCoverage = (d.coverage - propProtected) / (1.0 - propProtected);
        cerr << "cum deployment: prop protected " << propProtected << "; additionalCoverage " << additionalCoverage
             << "; total " << total << endl;
        for (Human* human : unprotected) {
            // parity: a — deploy short-circuit: RNG draw is the LAST test of the deploy condition; filter/age/cohort checks must precede it (plan §3a)
            if (human->rng().uniform_01() < additionalCoverage) {
                d.deployToHuman(ctx, *human, mon::Deploy::TIMED);
            }
        }
    }
}

/** One timed deployment (old TimedDeployment + derivatives). List sorted by
 * date, ties keeping load order (contract b); consumed via nextTimed cursor.
 * Move-only (changeHS/changeEIR payloads are consumed on deploy). */
struct TimedDeployment {
    SimTime date = sim::never();
    TimedDeploymentKind kind = TimedDeploymentKind::DUMMY;
    std::variant<std::monostate, TimedHumanParams, CumulativeHumanParams, VectorDeployParams, TrapDeployParams,
                 AddNhhDeployParams, NhhModsDeployParams, ChangeHsDeployParams, ChangeEirDeployParams>
        p;

    TimedDeployment() = default;
    TimedDeployment(TimedDeployment&&) = default;
    TimedDeployment& operator=(TimedDeployment&&) = default;
    TimedDeployment(const TimedDeployment&) = delete;
    TimedDeployment& operator=(const TimedDeployment&) = delete;

    void deploy(Context& ctx, std::vector<Host::Human>& population, Transmission::TransmissionModel& transmission);
};

inline void TimedDeployment::deploy(Context& ctx, std::vector<Host::Human>& population,
                                    Transmission::TransmissionModel& transmission) {
    switch (kind) {

    case TimedDeploymentKind::DUMMY:
        break;
    case TimedDeploymentKind::HUMAN:
        timed_human_deploy(ctx, population, std::get<TimedHumanParams>(p));
        break;
    case TimedDeploymentKind::CUMULATIVE_HUMAN:
        timed_cumulative_deploy(ctx, population, std::get<CumulativeHumanParams>(p));
        break;
    case TimedDeploymentKind::UNINFECT_VECTORS:
        transmission.uninfectVectors();
        break;
    case TimedDeploymentKind::VECTOR:
        if (transmission.isVector())
            transmission.deployVectorPopInterv(std::get<VectorDeployParams>(p).instance);
        break;
    case TimedDeploymentKind::TRAP: {
        const TrapDeployParams& t = std::get<TrapDeployParams>(p);
        if (transmission.isVector()) {
            double number = population.size() * t.ratio;
            transmission.deployVectorTrap(t.instance, number, t.lifespan);
        }
        break;
    }
    case TimedDeploymentKind::ADD_NHH:
        timed_deploy_add_nhh(ctx, population, transmission, std::get<AddNhhDeployParams>(p));
        break;
    case TimedDeploymentKind::NHH_MODS:
        timed_deploy_nhh_mods(ctx, population, transmission, std::get<NhhModsDeployParams>(p));
        break;
    case TimedDeploymentKind::CHANGE_HS:
        timed_deploy_change_hs(ctx, transmission, std::get<ChangeHsDeployParams>(p));
        break;
    case TimedDeploymentKind::CHANGE_EIR:
        timed_deploy_change_eir(ctx, transmission, std::get<ChangeEirDeployParams>(p));
        break;
    default:
        throw SWITCH_DEFAULT_EXCEPTION;
    }
}

/// Builders (TimedDeployment declares move/copy ops, so isn't brace-constructible).
inline TimedDeployment timed_deployment(SimTime date, TimedDeploymentKind kind) {
    TimedDeployment td;
    td.date = date;
    td.kind = kind;
    return td;
}

template <class P> inline TimedDeployment timed_deployment(SimTime date, TimedDeploymentKind kind, P&& payload) {
    TimedDeployment td = timed_deployment(date, kind);
    td.p = std::forward<P>(payload);
    return td;
}

/// PRINT_INTERVENTIONS listing; texts byte-identical to the old print_details.
inline void print_timed_details(std::ostream& out, const TimedDeployment& td) {
    switch (td.kind) {
    case TimedDeploymentKind::DUMMY:
        out << "Dummy";
        break;
    case TimedDeploymentKind::UNINFECT_VECTORS:
        out << "UninfectVectors";
        break;
    case TimedDeploymentKind::VECTOR:
        out << "Vector";
        break;
    case TimedDeploymentKind::CHANGE_HS:
        out << "ChangeHS";
        break;
    case TimedDeploymentKind::CHANGE_EIR:
        out << "ChangeEIR";
        break;
    case TimedDeploymentKind::HUMAN:
    case TimedDeploymentKind::CUMULATIVE_HUMAN: {
        // old TimedCumulativeHumanDeployment inherited the plain Human text
        const TimedHumanParams& th = td.kind == TimedDeploymentKind::HUMAN ? std::get<TimedHumanParams>(td.p)
                                                                           : std::get<CumulativeHumanParams>(td.p).h;
        const HumanDeployParams& d = th.d;
        out << "Human:" << endl;
        out << "\tage: " << sim::inYears(th.minAge) << "y\tmax age: " << sim::inYears(th.maxAge) << "y\tsubpop: ";
        if (d.subPop == ComponentId::wholePop())
            out << "whole";
        else
            out << d.subPop.id;
        out << "\tcomplement: " << d.complement << "\tcoverage: " << d.coverage << "\t";
        d.intervention->print_details(out);
        break;
    }
    case TimedDeploymentKind::TRAP: {
        const TrapDeployParams& t = std::get<TrapDeployParams>(td.p);
        out << "VectorTrap:" << endl;
        out << "\tratio: " << t.ratio << "\tlifespan: " << t.lifespan;
        break;
    }
    case TimedDeploymentKind::ADD_NHH:
        out << "AddNonHumanHosts(" << std::get<AddNhhDeployParams>(td.p).intervName << "):" << endl;
        out << "\tlifespan: " << std::get<AddNhhDeployParams>(td.p).lifespan;
        break;
    case TimedDeploymentKind::NHH_MODS:
        out << "NonHumanHosts(" << std::get<NhhModsDeployParams>(td.p).intervName << ")";
        break;
    default:
        throw SWITCH_DEFAULT_EXCEPTION;
    }
}

// ———  ContinuousHumanDeployment  ———

/** Continuous deployment of an intervention (not polymorphic; kept as a value
 * struct, sorted by deployAge). read() defined in xml/interventions.cpp. */
struct ContinuousHumanDeployment {
    static ContinuousHumanDeployment read(Context& ctx, SimTime begin, SimTime end,
                                          const ::scnXml::ContinuousDeployment& elt,
                                          shared_ptr<const HumanIntervention> intervention, ComponentId subPop,
                                          bool complement);

    /** Apply filters and potentially deploy.
     * 
     * @returns false iff this deployment (and thus all later ones in the
     *  ordered list) happens in the future. */
    bool filterAndDeploy(Context& ctx, Host::Human& human) const {
        SimTime age = human.age(ctx.clock.now());
        if (deployAge > age) {
            // stop processing continuous deployments for this
            // human for now because remaining ones happen in the future
            return false;
        } else if (deployAge == age) {
            double availability = 0.0;
            for (size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(ctx); ++i)
                availability += human.perHostTransmission.anoph_ento_raw()[i];

            if (availability >= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, minAvailability) &&
                availability <= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(ctx, maxAvailability)) {
                auto now = ctx.clock.interv_date();
                if (begin <= now && now < end &&
                    (d.subPop == ComponentId::wholePop() || (human.isInSubPop(d.subPop) != d.complement)) &&
                    // parity: a — deploy short-circuit: RNG draw is the LAST test of the deploy condition; filter/age/cohort checks must precede it (plan §3a)
                    human.rng().uniform_01() < d.coverage) // RNG call should be last test
                {
                    d.deployToHuman(ctx, human, mon::Deploy::CTS);
                }
            }
        } //else: for some reason, a deployment age was missed; ignore it
        return true;
    }

    inline void print_details(std::ostream& out) const {
        out << "Human:" << endl;
        out << "\tbegin: " << begin << "\tend: ";
        if (end == sim::future())
            out << "(none)";
        else
            out << end;
        out << "\tdeployAge: " << sim::inYears(deployAge) << "y\tsubpop: ";
        if (d.subPop == ComponentId::wholePop())
            out << "whole";
        else
            out << d.subPop.id;
        out << "\tcomplement: " << d.complement << "\tcoverage: " << d.coverage << "\t";
        d.intervention->print_details(out);
    }

    HumanDeployParams d;
    SimTime begin = sim::never(), end = sim::never(); // first active/first inactive step
    SimTime deployAge = sim::never();
    SimTime minAvailability = sim::never(), maxAvailability = sim::never();
};

struct ByDeployTime {
    bool operator()(const TimedDeployment& a, const TimedDeployment& b) const { return a.date < b.date; }

    bool operator()(const ContinuousHumanDeployment& a, const ContinuousHumanDeployment& b) const {
        return a.deployAge < b.deployAge;
    }
};
// inline variable: Deployments.h is now included by more than one source
inline const ByDeployTime byDeployTime{};

}
}
