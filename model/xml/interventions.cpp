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

// model/xml/ conversion layer — interventions module (start): imported
// infection rates. Moved verbatim from Host/ImportedInfections.h (the
// inline ImportedInfections::init became an out-of-line definition).

#include "Clock.h"
#include "RunContext.h"
#include "Host/ImportedInfections.h"
#include "Population.h"
#include "interventions/InterventionManager.h"
#include "interventions/HumanComponents.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "interventions/GVI.h"
#include "interventions/IRS.h"
#include "interventions/ITN.h"
#include "interventions/Vaccine.h"
#include "xml/HumanInterventionComponents.h"
#include "xml/entomology.h"
#include "interventions/Deployments.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/PerHost.h"
#include "Clinical/ClinicalModel.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "R_nmath/qnorm.h"
#include "util/errors.h"
#include "schema/interventions.h"
#include "schema/healthSystem.h"

#include <algorithm>
#include <string>

using namespace std;

namespace OM {

// UnitParse::readDate with the STEPS fallback (STEPS default is only for
// backwards compatibility; former per-call comment, now stated once).
static SimTime readDate_steps(const SimClock& clock, const std::string& str) {
    return UnitParse::readDate(clock, str, UnitParse::STEPS);
}

namespace Host {

void ImportedInfections::init(RunContext& ctx, const scnXml::ImportedInfections& iiElt) {
    const scnXml::ImportedInfections::TimedType& tElt = iiElt.getTimed();
    try {
        //NOTE: if changing XSD, this should not have a default unit:
        period = UnitParse::readDuration(ctx.clock, tElt.getPeriod(), UnitParse::STEPS);
        if (period < sim::zero()) {
            throw util::format_error("cannot be negative");
        }
    } catch (const util::format_error& e) {
        throw util::xml_scenario_error(string("interventions/importedInfections/timed/period: ").append(e.message()));
    }
    rate.reserve(tElt.getRate().size());
    try {
        for (const auto& rateElt : tElt.getRate()) {
            SimTime date = readDate_steps(ctx.clock, rateElt.getTime());
            // convert to per-time-step, per-person
            double rateVal = rateElt.getValue() * ctx.clock.years_per_step * (1.0 / 1000.0);
            rate.push_back(Rate(date - ctx.clock.start_date, rateVal));
        }
        sort(rate.begin(), rate.end());
        if (rate.size() > 0) {
            if (period != sim::zero() && rate[0].time != sim::zero()) {
                throw util::xml_scenario_error("must specify rate at time zero when period is not zero");
            }
            // remove useless repeated entries from list
            double lastRate = rate[0].value;
            for (size_t i = 1; i < rate.size();) {
                if (rate[i].value == lastRate) {
                    rate.erase(rate.begin() + i);
                } else {
                    lastRate = rate[i].value;
                    ++i;
                }
            }
        }
    } catch (const util::format_error& e) {
        throw util::xml_scenario_error(string("interventions/importedInfections/timed/time: ").append(e.message()));
    }
}

}
}

namespace OM {
namespace interventions {

using namespace util; // moved code (ITN/IRS/GVI/Vaccine) uses unqualified util names

// Frozen error texts duplicated verbatim in the parse ladders below; single
// definitions here, every throw site emits the identical bytes (message
// multiset and first-error order unchanged — punch-list #16).
static const string min_avail_err = "timed intervention must have 0 <= minAvailability <= maxAvailability <= 100";
static const string nonvec_err = "vector model interventions can not be used with the non-vector model";

// 7.1: changeHS/changeEIR timed-deployment classes gone; the payload is the
// schema clone in the TimedDeployment variant, consumed (deleted) on deploy.
// Never-deployed clones leak at exit — exactly the old classes' behaviour
// (their destructors never freed these).
void timed_deploy_change_hs(RunContext& ctx, Transmission::TransmissionModel&, ChangeHsDeployParams& p) {
    Clinical::ClinicalModel::setHS(ctx, *p.newHS);
    delete p.newHS;
    p.newHS = 0;
}

void timed_deploy_change_eir(RunContext& ctx, Transmission::TransmissionModel& transmission, ChangeEirDeployParams& p) {
    transmission.changeEIRIntervention(*p.newEIR);
    delete p.newEIR;
    p.newEIR = 0;
}

static int vacc_nonneg(int v, const char* err) {
    if (v < 0)
        throw util::xml_scenario_error(err);
    return v;
}
void VaccineLimits::set(const scnXml::DeploymentBase& deploy) {
    if (deploy.getVaccMinPrevDoses().present())
        minPrevDoses = vacc_nonneg(deploy.getVaccMinPrevDoses().get(), "vaccMinPrevDoses: min value is 0");
    if (deploy.getVaccMaxCumDoses().present())
        maxCumDoses = vacc_nonneg(deploy.getVaccMaxCumDoses().get(), "vaccMaxCumDoses: min value is 0");
}

HumanDeployParams HumanDeployParams::read(RunContext&, const scnXml::DeploymentBase& deploy,
                                          shared_ptr<const HumanIntervention> intervention, ComponentId subPop,
                                          bool complement) {
    HumanDeployParams d;
    d.coverage = deploy.getCoverage();
    d.subPop = subPop;
    d.complement = complement;
    d.intervention = intervention;
    if (!(d.coverage >= 0.0 && d.coverage <= 1.0)) {
        throw util::xml_scenario_error("intervention deployment coverage must be in range [0,1]");
    }
    // If coverageCorr or coverageVar is present, then the user is asking for the copula
    if (deploy.getCoverageCorr().present() || deploy.getCoverageVar().present()) {
        if (!(deploy.getCoverageCorr().present() && deploy.getCoverageVar().present()))
            throw std::runtime_error("[Gaussian copula]: both coverageCorr and coverageVar must be present");

        if (deploy.getCoverageVar().get() < 0)
            throw std::runtime_error("[Gaussian copula]: negative variance is not allowed");

        // if variance is 0, then no copula and fall back to probability = coverage
        if (deploy.getCoverageVar().get() > 0) {
            d.coverageCorr = deploy.getCoverageCorr().get();
            d.coverageVar = deploy.getCoverageVar().get();
            d.copula = true;
        }
    }

    d.vaccLimits.set(deploy);
    return d;
}

/// max/minAvailability present-reads + range check, shared by the timed and
/// continuous human deployment readers (getMax before getMin, same message
/// and check position as the former duplicated blocks).
template <class Elt> static void read_avail_limits(const Elt& elt, SimTime& mn, SimTime& mx) {
    if (elt.getMaxAvailability().present())
        mx = elt.getMaxAvailability().get();
    if (elt.getMinAvailability().present())
        mn = elt.getMinAvailability().get();
    if (mn < sim::zero() || mx < mn)
        throw util::xml_scenario_error(min_avail_err);
}

TimedHumanParams TimedHumanParams::read(RunContext& ctx, SimTime date, const scnXml::MassDeployment& mass,
                                        shared_ptr<const HumanIntervention> intervention, ComponentId subPop,
                                        bool complement) {
    check_timed_date(ctx, date); // was the TimedDeployment base ctor: runs first
    TimedHumanParams th;
    th.d = HumanDeployParams::read(ctx, mass, intervention, subPop, complement);
    th.minAge = ctx.clock.from_years_n(mass.getMinAge());
    th.maxAge = sim::future();
    th.minAvailability = 0;
    th.maxAvailability = 100;

    if (mass.getMaxAge().present())
        th.maxAge = ctx.clock.from_years_n(mass.getMaxAge().get());

    if (th.minAge < sim::zero() || th.maxAge < th.minAge)
        throw util::xml_scenario_error(min_avail_err);

    read_avail_limits(mass, th.minAvailability, th.maxAvailability);
    return th;
}

ContinuousHumanDeployment ContinuousHumanDeployment::read(RunContext& ctx, SimTime begin, SimTime end,
                                                          const ::scnXml::ContinuousDeployment& elt,
                                                          shared_ptr<const HumanIntervention> intervention,
                                                          ComponentId subPop, bool complement) {
    ContinuousHumanDeployment c;
    c.d = HumanDeployParams::read(ctx, elt, intervention, subPop, complement);
    c.begin = begin;
    c.end = end;
    c.deployAge = ctx.clock.from_years_n(elt.getTargetAgeYrs());
    c.minAvailability = 0;
    c.maxAvailability = 100;
    if (begin < ctx.clock.start_date || end < c.begin) {
        throw util::xml_scenario_error("continuous intervention must have startDate <= begin <= end");
    }
    if (c.deployAge <= sim::zero()) {
        ostringstream msg;
        msg << "continuous intervention with target age " << elt.getTargetAgeYrs();
        msg << " years corresponds to time step " << ctx.clock.inSteps(c.deployAge);
        msg << "; must be at least 1.";
        throw util::xml_scenario_error(msg.str());
    }
    if (c.deployAge > ctx.clock.max_human_age) {
        ostringstream msg;
        msg << "continuous intervention must have target age no greater than ";
        msg << sim::inYears(ctx.clock.max_human_age);
        throw util::xml_scenario_error(msg.str());
    }
    read_avail_limits(elt, c.minAvailability, c.maxAvailability);
    return c;
}

/// One <anopheles> effect slot of nhh_mods_setup: present-check, initial-effect
/// range check [lo,hi] (both bounds compare-reject like the former per-type
/// checks; NaN passes like before), then decay set-up. Messages verbatim.
template <class M, class Opt>
static void nhh_slot(const SimClock& clock, M& m, const std::string& intervName, size_t instance, const Opt& opt,
                     const scnXml::DecayFunction& decay, double lo, double hi, const char* rangeErr, const char* tag) {
    if (!opt.present())
        return;
    const double v = opt.get().getInitial();
    if (v < lo || v > hi)
        throw util::xml_scenario_error(rangeErr);
    m[intervName][instance].set(clock, v, decay, tag);
}

/// Old TimedNonHumanHostsDeployment ctor body (decay set-up): runs at load,
/// once per timed deploy entry (idempotent, as before). check_timed_date is
/// the caller's (was the TimedDeployment base ctor, which ran first).
static void nhh_mods_setup(RunContext& ctx, size_t instance, string intervName,
                           const scnXml::Description2::AnophelesSequence list, const scnXml::DecayFunction& decay,
                           Transmission::TransmissionModel& transmission) {
    Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
    if (vectorModel) {
        util::SpeciesIndexChecker checker(intervName, vectorModel->speciesIndex);
        for (const scnXml::NonHumanHostsSpeciesIntervention& anoph : list) {
            const string& mosqName = anoph.getMosquito();
            size_t i = checker.getIndex(mosqName);

            Transmission::Anopheles::AnophelesModel* anophModel = vectorModel->species[i].get();

            // Resize every effect list first, as before (operator[] insertion of
            // intervName happens even for absent effects; loop order == former
            // textual order):
            for (auto* m : {&anophModel->reduceNhhAvailability, &anophModel->reduceP_B_I, &anophModel->reduceP_C_I,
                            &anophModel->reduceP_D_I, &anophModel->reduceFecundity}) {
                auto& slot = (*m)[intervName];
                if (slot.size() <= instance)
                    slot.resize(instance + 1);
            }

            nhh_slot(ctx.clock, anophModel->reduceNhhAvailability, intervName, instance,
                     anoph.getAvailabilityReduction(), decay, -numeric_limits<double>::infinity(), 1.0,
                     "availabilityReduction intervention: initial effect must be <= 1", "availabilityReduction");
            nhh_slot(ctx.clock, anophModel->reduceP_B_I, intervName, instance, anoph.getPreprandialKillingEffect(),
                     decay, 0.0, 1.0, "PreprandialKillingEffect intervention: initial effect must be between 0 and 1",
                     "reduceP_B_I");
            nhh_slot(ctx.clock, anophModel->reduceP_C_I, intervName, instance, anoph.getPostprandialKillingEffect(),
                     decay, 0.0, 1.0, "PostprandialKillingEffect intervention: initial effect must be between 0 and 1",
                     "reduceP_C_I");
            nhh_slot(ctx.clock, anophModel->reduceP_D_I, intervName, instance, anoph.getRestingKillingEffect(), decay,
                     0.0, 1.0, "RestingKillingEffect intervention: initial effect must be be between 0 and 1",
                     "reduceP_D_I");
            nhh_slot(ctx.clock, anophModel->reduceFecundity, intervName, instance, anoph.getFecundityReduction(), decay,
                     0.0, 1.0, "FecundityReduction intervention: initial effect must be be between 0 and 1",
                     "reduceFecundity");
        }
        checker.checkNoneMissed();
    }
}

/// Old TimedNonHumanHostsDeployment::deploy body, verbatim.
void timed_deploy_nhh_mods(RunContext& ctx, std::vector<Host::Human>&, Transmission::TransmissionModel& transmission,
                           const NhhModsDeployParams& params) {
    const size_t instance = params.instance;
    const string& intervName = params.intervName;
    Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
    if (vectorModel) {
        for (size_t i = 0; i < vectorModel->speciesIndex.size(); ++i) {
            Transmission::Anopheles::AnophelesModel* anophModel = vectorModel->species[i].get();

            if (anophModel->nhhInstances.count(intervName) == 0)
                throw util::xml_scenario_error("non human hosts type " + intervName +
                                               " not deployed during non human hosts intervention deployment");

            anophModel->reduceNhhAvailability[intervName][instance].deploy(vectorModel->m_rng, ctx.clock.now());
            anophModel->reduceP_B_I[intervName][instance].deploy(vectorModel->m_rng, ctx.clock.now());
            anophModel->reduceP_C_I[intervName][instance].deploy(vectorModel->m_rng, ctx.clock.now());
            anophModel->reduceP_D_I[intervName][instance].deploy(vectorModel->m_rng, ctx.clock.now());
            anophModel->reduceFecundity[intervName][instance].deploy(vectorModel->m_rng, ctx.clock.now());
        }
    }
}
/// Old TimedAddNonHumanHostsDeployment ctor body (per-species parameter read);
/// check_timed_date (was the TimedDeployment base ctor) runs first, as before.
static AddNhhDeployParams read_add_nhh(RunContext& ctx, SimTime date, const string& intervName, SimTime lifespan,
                                       const scnXml::Description3::AnophelesSequence list,
                                       Transmission::TransmissionModel& transmission) {
    check_timed_date(ctx, date);
    AddNhhDeployParams a;
    a.intervName = intervName;
    a.lifespan = lifespan;
    Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
    if (vectorModel) {
        util::SpeciesIndexChecker checker(intervName, vectorModel->speciesIndex);
        for (const scnXml::NonHumanHostsVectorSpecies& anoph : list) {
            const string& mosqName = anoph.getMosquito();
            checker.getIndex(mosqName);

            AddNhhDeployParams::NhhSpec nhh;
            nhh.mosqRelativeAvailabilityHuman = anoph.getMosqRelativeAvailabilityHuman().getValue();
            nhh.mosqProbBiting = anoph.getMosqProbBiting().getValue();
            nhh.mosqProbFindingRestSite = anoph.getMosqProbFindRestSite().getValue();
            nhh.mosqProbResting = anoph.getMosqProbResting().getValue();
            nhh.hostFecundityFactor = anoph.getHostFecundityFactor().getValue();
            a.bySpecies[mosqName] = nhh;
        }
        checker.checkNoneMissed();
    }
    return a;
}

/// Old TimedAddNonHumanHostsDeployment::deploy body, verbatim (nhhParams is
/// the old member of the same name).
void timed_deploy_add_nhh(RunContext& ctx, vector<Host::Human>& population,
                          Transmission::TransmissionModel& transmission, AddNhhDeployParams& nhhParams) {
    const string& intervName = nhhParams.intervName;
    const SimTime lifespan = nhhParams.lifespan;
    Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
    if (vectorModel) {
        double popSize = population.size();
        if (vectorModel->interventionMode != Transmission::SimulationMode::dynamicEIR) {
            throw util::xml_scenario_error(Transmission::vec_mode_err);
        }
        for (size_t i = 0; i < vectorModel->speciesIndex.size(); ++i) {
            Transmission::Anopheles::AnophelesModel* anophModel = vectorModel->species[i].get();

            if (anophModel->nhhInstances.count(intervName) != 0)
                throw util::xml_scenario_error("non human hosts type " + intervName +
                                               " already deployed during non human hosts deployment");

            const AddNhhDeployParams::NhhSpec& p = nhhParams.bySpecies[anophModel->mosq.name];
            const double adultAvail = Transmission::PerHostAnophParams::get(ctx, i).entoAvailability->mean();
            const double avail_i = Transmission::PerHostAnophParams::get(ctx, i).entoAvailabilityFactor * popSize *
                                   adultAvail * p.mosqRelativeAvailabilityHuman;
            Transmission::Anopheles::Nhh nhh;

            nhh.avail_i = avail_i;
            nhh.P_B_I = p.mosqProbBiting;
            nhh.P_C_I = p.mosqProbFindingRestSite;
            nhh.P_D_I = p.mosqProbResting;
            nhh.rel_fecundity = p.hostFecundityFactor;
            nhh.expiry = ctx.clock.now() + lifespan;

            // add the nhh to the active nhh instances
            anophModel->nhhInstances[intervName] = nhh;
        }
    }
}

// ———  ITN/IRS/GVI/Vaccine parameter readers  ———

// Shared property readers for the three vector interventions (their former
// duplicated blocks; present()/get() call order and messages travel verbatim).
template <class Usage> static double read_propUse(const Usage& usage, const char* rangeErr) {
    const double u = usage.present() ? usage.get().getValue() : 1.0;
    if (!(u >= 0.0 && u <= 1.0))
        throw util::xml_scenario_error(rangeErr);
    return u;
}
static void set_protected(double proportionUse, double propActive, double& prot, double& unprot) {
    if (propActive != 1.0 && propActive != 0.0)
        throw util::xml_scenario_error("GVI.description.propActive: must be either 0 or 1");
    assert(proportionUse >= 0.0 && proportionUse <= 1.0);
    assert(propActive >= 0.0 && propActive <= 1.0);
    prot = proportionUse * propActive;
    unprot = 1.0 - prot;
}
// Per-species parameter loop: resize, index-check each <anopheles>, init via
// the type's init(anoph, args...), then checkNoneMissed (former 3 copies).
template <class Seq, class V>
static void per_species(std::vector<V>& species, const char* intervName,
                        const std::map<std::string, size_t>& species_name_map, const Seq& ap, auto&&... initArgs) {
    species.resize(species_name_map.size());
    util::SpeciesIndexChecker checker(intervName, species_name_map);
    for (const auto& it : ap)
        species[checker.getIndex(it.getMosquito())].init(it, initArgs...);
    checker.checkNoneMissed();
}

void factors::SurvivalFactor::init(const scnXml::ITNKillingEffect& elt, double maxInsecticide, const char* eltName,
                                   bool isDeterrent) {
    useLogitEqns = false;
    a.BF = elt.getBaseFactor();
    a.HF = elt.getHoleFactor();
    a.PF = elt.getInsecticideFactor();
    a.IF = elt.getInteractionFactor();
    a.holeScaling = elt.getHoleScalingFactor();
    a.insecticideScaling = elt.getInsecticideScalingFactor();
    a.invBaseSurvival = 1.0 / (1.0 - a.BF);
    if (!(a.BF >= 0.0 && a.BF < 1.0)) {
        ostringstream msg;
        msg << eltName << ": expected baseFactor to be in range [0,1]";
        throw util::xml_scenario_error(msg.str());
    }
    if (!(a.holeScaling >= 0.0 && a.insecticideScaling >= 0.0)) {
        ostringstream msg;
        msg << eltName << ": expected scaling factors to be non-negative";
        throw util::xml_scenario_error(msg.str());
    }
    // see below
    double pmax = 1.0 - exp(-maxInsecticide * a.insecticideScaling);

    if (isDeterrent) {
        // Note: the following argument is a modification of the one below
        // (when !isDeterrent). The original may make more sense.
        /* We want K ≥ 0 where K is the killing factor:
    K=BF+HF×h+PF×p+IF×h×p, with h and p defined as:
    h=exp(-holeIndex×holeScalingFactor),
    p=1−exp(-insecticideContent×insecticideScalingFactor). 
    
    By their nature, holeIndex ≥ 0 and insecticideContent ≥ 0. We restrict:
        holeScalingFactor ≥ 0
        insecticideScalingFactor ≥ 0
    Which implies both h and p lie in the range [0,1]. We also know 0 ≤ BF ≤ 1.
    
    We need K ≥ 0 or:
        BF + HF×h + PF×p + IF×h×p ≥ 0   (1)
    
    Lets derive some limits on HF, PF and IF such that the above inequality (1)
    is satisfied.
    
    A net can theoretically be unholed (holeIndex=0 ⇒ h=1) and have no
    insecticide (thus have p=0). Substituting these values in (1) yields:
        BF + HF ≥ 0     (2)
    
    The maximum value for p depends on the maximum insecticide content; denote
    pmax = max(p). Note that holeIndex has no finite maximum; thus, although
    for any finite value of holeIndex, h > 0, there is no h₀ > 0 s.t. for all
    values of holeIndex h ≥ h₀. For the limiting case of a tattered but
    insecticide-saturated net our parameters are therefore p=pmax, h=0:
        BF + PF×pmax ≥ 0        (3)
    
    Consider a net saturated with insecticide (p=pmax) and without holes (h=1):
        BF + HF + (PF+IF)×pmax ≥ 0      (4)
    
    The opposite extreme (the limiting case of a decayed net with no remaining
    insecticide and a large number of holes) yields only BF ≥ 0.
    
    Some of the above examples of nets may be unlikely, but there is only one
    restriction in our model making any of these cases impossible: some
    insecticide must have been lost by the time any holes occur. We ignore this
    since its effect is likely small, and thus all of the above are required to
    keep the factor in the range [0,1]. Further, the inequalities (2) - (3) are
    sufficient to keep the factor within [0,1] since h and p are non-negative
    and act linearly in (1).
    
    From the definition of p, we always have pmax ≤ 1, so substituting pmax=1
    in (5) - (8) gives us bounds which imply our requirement, however if pmax
    is finite they are stricter than necessary. Since insecticideScalingFactor
    is constant, max(p) coincides with max(insecticideContent) which, since
    insecticide content only decays over time, coincides with the maximum
    initial insecticide content, Pmax. Since the initial insecticide content is
    sampled from a normal distribution in our model it should have no finite
    maximum, thus implying we cannot achieve more relaxed bounds than (5) - (8)
    when pmax=1 (unless the standard deviation of our normal distribution is 1).
    We would however like to impose less strict bounds than these, thus we
    impose a maximum value on the initial insecticide content, Pmax, such that
    the probability of sampling a value from our parameterise normal
    distribution greater than Pmax is 0.001. */
        if (!(a.BF + a.HF >= 0.0 && a.BF + a.PF * pmax >= 0.0 && a.BF + a.HF + (a.PF + a.IF) * pmax >= 0.0)) {
            ostringstream msg;
            msg << eltName << ": bounds not met:";
            if (!(a.BF + a.HF >= 0.0))
                msg << " baseFactor+holeFactor≥0";
            if (!(a.BF + a.PF * pmax >= 0.0))
                msg << " baseFactor+" << pmax << "×insecticideFactor≥0";
            if (!(a.PF + a.HF + (a.PF + a.IF) * pmax >= 0.0))
                msg << " baseFactor+holeFactor+" << pmax << "×(insecticideFactor+interactionFactor)≥0";
            throw util::xml_scenario_error(msg.str());
        }
    } else {
        /* We want the calculated survival factor (1−K)/(1−BF) to be in the range
    [0,1] where K is the killing factor: K=BF+HF×h+PF×p+IF×h×p, with h and p
    defined as: h=exp(-holeIndex×holeScalingFactor),
    p=1−exp(-insecticideContent×insecticideScalingFactor). 
    
    By their nature, holeIndex ≥ 0 and insecticideContent ≥ 0. We restrict:
        holeScalingFactor ≥ 0
        insecticideScalingFactor ≥ 0
    Which implies both h and p lie in the range [0,1]. We also know the base
    survival factor, 1−BF, is in the range [0,1].
    
    To make sure the survival factor is not negative we need (1−K)/(1−BF) ≥ 0.
    Since 1−BF > 0 we need 1−K ≥ 0, which, substituting K, gives us
        BF + HF×h + PF×p + IF×h×p ≤ 1	(1)
    We also want to make sure the survival factor is not greater than one (since
    nets shouldn't increase mosquito survival), (1−K)/(1−BF) ≤ 1 or equivalently
    1-K ≤ 1-BF or K ≥ BF, which, substituting K, yields
        HF×h + PF×p + IF×h×p ≥ 0		(2)
    
    Lets derive some limits on HF, PF and IF such that the above inequalities
    (1) and (2) are satisfied.
    
    A net can theoretically be unholed (holeIndex=0 ⇒ h=1) and have no
    insecticide (thus have p=0). Substituting these values in (1) and (2) yields:
        BF + HF ≤ 1	(3)
        HF ≥ 0		(4)
    
    The maximum value for p depends on the maximum insecticide content; denote
    pmax = max(p). Note that holeIndex has no finite maximum; thus, although
    for any finite value of holeIndex, h > 0, there is no h₀ > 0 s.t. for all
    values of holeIndex h ≥ h₀. For the limiting case of a tattered but
    insecticide-saturated net our parameters are therefore p=pmax, h=0:
        BF + PF×pmax ≤ 1	(5)
        PF×pmax ≥ 0			(6)
    (Assuming pmax > 0, (6) is equivalent to PF ≥ 0.)
    
    Consider a net saturated with insecticide (p=pmax) and without holes (h=1):
        BF + HF + (PF+IF)×pmax ≤ 1	(7)
        HF + (PF+IF)×pmax ≥ 0		(8)
    
    The opposite extreme (the limiting case of a decayed net with no remaining
    insecticide and a large number of holes) yields only BF ≤ 1 which we already
    know.
    
    Some of the above examples of nets may be unlikely, but there is only one
    restriction in our model making any of these cases impossible: some
    insecticide must have been lost by the time any holes occur. We ignore this
    since its effect is likely small, and thus all of the above are required to
    keep the survival factor in the range [0,1]. Further, these six inequalities
    (3) - (8) are sufficient to keep the survival factor within [0,1] since h
    and p are non-negative and act linearly in (1) and (2).
    
    From the definition of p, we always have pmax ≤ 1, so substituting pmax=1
    in (5) - (8) gives us bounds which imply our requirement, however if pmax
    is finite they are stricter than necessary. Since insecticideScalingFactor
    is constant, max(p) coincides with max(insecticideContent) which, since
    insecticide content only decays over time, coincides with the maximum
    initial insecticide content, Pmax. Since the initial insecticide content is
    sampled from a normal distribution in our model it should have no finite
    maximum, thus implying we cannot achieve more relaxed bounds than (5) - (8)
    when pmax=1 (unless the standard deviation of our normal distribution is 1).
    We would however like to impose less strict bounds than these, thus we
    impose a maximum value on the initial insecticide content, Pmax, such that
    the probability of sampling a value from our parameterise normal
    distribution greater than Pmax is 0.001. */
        if (!(a.BF + a.HF <= 1.0 && a.HF >= 0.0 && a.BF + a.PF * pmax <= 1.0 && a.PF >= 0.0 &&
              a.BF + a.HF + (a.PF + a.IF) * pmax <= 1.0 && a.HF + (a.PF + a.IF) * pmax >= 0.0)) {
            ostringstream msg;
            msg << eltName << ": bounds not met:";
            if (!(a.BF + a.HF <= 1.0))
                msg << " baseFactor+holeFactor≤1";
            if (!(a.HF >= 0.0))
                msg << " holeFactor≥0";
            if (!(a.BF + a.PF * pmax <= 1.0))
                msg << " baseFactor+" << pmax << "×insecticideFactor≤1";
            if (!(a.PF >= 0.0))
                msg << " insecticideFactor≥0"; // if this fails, we know pmax>0 (since it is in any case non-negative) — well, or an NaN
            if (!(a.PF + a.HF + (a.PF + a.IF) * pmax <= 1.0))
                msg << " baseFactor+holeFactor+" << pmax << "×(insecticideFactor+interactionFactor)≤1";
            if (!(a.HF + (a.PF + a.IF) * pmax >= 0.0))
                msg << " holeFactor+" << pmax << "×(insecticideFactor+interactionFactor)≥0";
            throw util::xml_scenario_error(msg.str());
        }
    }
}

void factors::SurvivalFactor::initLogit(const scnXml::ITNEffectLogit& elt, double holeIndexMax, bool isDeterrent) {
    useLogitEqns = true;
    b.BF = elt.getBaseFactor();
    b.HF = elt.getHoleFactor();
    b.PF = elt.getInsecticideFactor();
    b.IF = elt.getInteractionFactor();
    b.hMax = log(holeIndexMax + 1.0);
    const double x = exp(b.BF + b.HF * b.hMax);
    if (isDeterrent) {
        // pAtt0 =x/(x+1); this is 1/pAtt0:
        b.invBaseEffect = (x + 1.0) / x;
    } else {
        // K0 = x/(x+1), so 1/(1-K0) = x+1:
        b.invBaseEffect = x + 1.0;

        // We expect K >= K0 and deduce these "advisory limits":
        if (b.PF < 0.0 || b.PF + b.IF * b.hMax < 0.0) {
            cerr << "ITN.description.anophelesParams.*KillingEffectLogit: expected"
                 << "\n(insecticide) P >= 0, found P = " << b.PF
                 << "\n(interaction) P+I*log(holeIndexMax+1) >= 0, found " << b.PF + b.IF * b.hMax << endl;
        }
    }
}

void factors::RelativeAttractiveness::initSingleStage(const scnXml::ITNDeterrency& elt, double maxInsecticide) {
    model = SINGLE_STAGE;
    double HF = elt.getHoleFactor();
    double PF = elt.getInsecticideFactor();
    double IF = elt.getInteractionFactor();
    a.holeScaling = elt.getHoleScalingFactor();
    a.insecticideScaling = elt.getInsecticideScalingFactor();
    if (!(a.holeScaling >= 0.0 && a.insecticideScaling >= 0.0)) {
        throw util::xml_scenario_error(
            "ITN.description.anophelesParams.deterrency: expected scaling factors to be non-negative");
    }

    /* We need to ensure the relative availability is non-negative. However,
     * since it's an exponentiated value, it always will be.
     * 
     * If don't want nets to be able to increase transmission, the following
     * limits could also be applied. In general, however, there is no reason
     * nets couldn't make individuals more attractive to mosquitoes.
     * 
     * To ensure relative availability is at most one: relative availability is
     *  exp( log(HF)*h + log(PF)*p + log(IF)*h*p )
     * where HF, PF and IF are the hole, insecticide and interaction factors
     * respectively, with h and p defined as:
     *  h=exp(-holeIndex*holeScalingFactor),
     *  p=1−exp(-insecticideContent*insecticideScalingFactor).
     * We therefore need to ensure that:
     *  log(HF)*h + log(PF)*p + log(IF)*h*p ≤ 0
     * 
     * As with the argument below concerning limits of the killing effect
     * parameters, h and p will always be in the range [0,1] and p ≤ pmax.
     * We can then derive some bounds for HF and PF:
     *  log(HF) ≤ 0
     *  log(PF)×pmax ≤ 0
     *  log(HF) + (log(PF)+log(IF))×pmax = log(HF×(PF×IF)^pmax) ≤ 0
     * or equivalently (assuming pmax>0):
     *  HF ∈ (0,1]
     *  PF ∈ (0,1]
     *  HF×(PF×IF)^pmax ∈ (0,1]
     *
     * Weaker limits would not be sufficient, as with the argument for the
     * limits of killing effect arguments below. */
    // Print out a warning if nets may increase transmission.
    double pmax = 1.0 - exp(-maxInsecticide * a.insecticideScaling);
    if (!(HF > 0.0 && PF > 0.0 && IF > 0.0 && HF <= 1.0 && PF <= 1.0 && HF * pow(PF * IF, pmax) <= 1.0)) {
        cerr << "Note: since the following bounds are not met, the ITN could make humans more\n";
        cerr << "attractive to mosquitoes than they would be without a net.\n";
        cerr << "ITN.description.anophelesParams.deterrency: bounds not met:\n";
        if (!(HF > 0.0))
            cerr << "  holeFactor>0\n";
        if (!(PF > 0.0))
            cerr << "  insecticideFactor>0\n";
        if (!(IF > 0.0))
            cerr << "  interactionFactor>0\n";
        if (!(HF <= 1.0))
            cerr << "  holeFactor≤1\n";
        if (!(PF <= 1.0))
            cerr << "  insecticideFactor≤1\n";
        if (!(HF * pow(PF * IF, pmax) <= 1.0))
            cerr << "  holeFactor×(insecticideFactor×interactionFactor)^" << pmax << "≤1\n";
        cerr.flush();
    }
    a.lHF = log(HF);
    a.lPF = log(PF);
    a.lIF = log(IF);
}

void factors::RelativeAttractiveness::initTwoStage(const scnXml::TwoStageDeterrency& elt, double maxInsecticide,
                                                   double holeIndexMax, bool holeIndexMaxPresent) {
    b.lPFEntering = numeric_limits<double>::quiet_NaN();
    if (elt.getEntering().present()) {
        model = TWO_STAGE;
        const double PF = elt.getEntering().get().getInsecticideFactor();
        b.insecticideScalingEntering = elt.getEntering().get().getInsecticideScalingFactor();
        if (!(PF > 0.0)) {
            // we take the log of PF, so it must be positive
            throw util::xml_scenario_error("ITN.description.anophelesParams.twoStageDeterrency.entering: "
                                           "insecticideFactor must be positive since we take its logarithm.");
        }

        /* We need to ensure the relative availability is non-negative. However,
        * since it's an exponentiated value, it always will be.
        * 
        * If we don't want ITNs to be able to increase transmission, the following
        * limits could also be applied. In general, however, there is no reason
        * ITNs couldn't make individuals more attractive to mosquitoes.
        * 
        * To ensure relative availability is at most one: relative availability is
        *  exp( log(PF)*p ) = PF^p
        * where PF is the insecticide factor, with p∈[0,1] defined as:
        *  p=1−exp(-insecticideContent*insecticideScalingFactor).
        * We therefore just need PF ≤ 1. */
        // Print out a warning if ITNs may increase transmission.
        if (!(PF <= 1.0)) {
            cerr << "Note: since the following bounds are not met, the IRS could make humans more\n";
            cerr << "attractive to mosquitoes than they would be without IRS.\n";
            cerr << "IRS.description.anophelesParams.deterrency: bounds not met:\n";
            cerr << "  0<insecticideFactor≤1\n";
            cerr.flush();
        }
        b.lPFEntering = log(PF);
    } else {
        assert(elt.getEnteringLogit().present());
        model = TWO_STAGE_LOGIT;
        c.entBaseFactor = elt.getEnteringLogit().get().getBaseFactor();
        c.entInsecticideFactor = elt.getEnteringLogit().get().getInsecticideFactor();
        if (c.entInsecticideFactor > 0.0) {
            cerr << "ITN.description.anophelesParams.twoStageDeterrency.enteringLogit: \n"
                 << "insecticideFactor should be negative to for nets to reduce chance of entering." << endl;
        }
        // pre-calculate for efficieny:
        c.pEnt0Inv = (exp(c.entBaseFactor) + 1.0) / exp(c.entBaseFactor);
    }

    // Note: b.pAttacking and c.pAttacking overlap, so we can use either:
    if (elt.getAttacking().present()) {
        b.pAttacking.init(elt.getAttacking().get(), maxInsecticide,
                          "ITN.description.anophelesParams.twoStageDeterrency.attacking", true);
    } else {
        assert(elt.getAttackingLogit().present());
        if (!holeIndexMax) {
            throw util::xml_scenario_error(
                "ITN.description.anophelesParams: holeIndexMax required when using logit attacking deterrency");
        }
        b.pAttacking.initLogit(elt.getAttackingLogit().get(), holeIndexMax, true);
    }
}

ITNParams::ITNParams(RunContext& ctx, ComponentId /*id*/, const scnXml::ITNDescription& elt,
                     const map<string, size_t>& species_name_map)
    : ripFactor(numeric_limits<double>::signaling_NaN()) {
    initialInsecticide.setParams(elt.getInitialInsecticide());
    const double maxProp = 0.999; //NOTE: this could be exposed in XML, but probably doesn't need to be
    maxInsecticide = R::qnorm5(maxProp, initialInsecticide.getMu(), initialInsecticide.getSigma(), true, false);
    holeRate = util::createSampler(elt.getHoleRate()); // per year
    holeRate->scaleMean(ctx.clock.years_per_step);     // convert to per step
    ripRate = util::createSampler(elt.getRipRate());
    ripRate->scaleMean(ctx.clock.years_per_step);
    ripFactor = elt.getRipFactor().getValue();
    insecticideDecay =
        DecayFunction::makeObject(ctx.clock, elt.getInsecticideDecay(), "ITNDescription.insecticideDecay");
    attritionOfNets = DecayFunction::makeObject(ctx.clock, elt.getAttritionOfNets(), "ITNDescription.attritionOfNets");
    // assume usage modifier is 100% if none is specified
    const double propUse = read_propUse(elt.getUsage(), "ITN.description.proportionUse: must be within range [0,1]");
    per_species(species, "ITN", species_name_map, elt.getAnophelesParams(), propUse, maxInsecticide);
}

void ITNParams::ITNAnopheles::init(const scnXml::ITNDescription::AnophelesParamsType& elt, double proportionUse,
                                   double maxInsecticide) {
    bool holeIndexMaxPresent = false;

    double holeIndexMax = 0.0;
    if (elt.getHoleIndexMax().present()) {
        holeIndexMax = elt.getHoleIndexMax().get().getValue();
        holeIndexMaxPresent = true;
    }

    if (elt.getDeterrency().present()) {
        relAttractiveness.initSingleStage(elt.getDeterrency().get(), maxInsecticide);
    } else {
        assert(elt.getTwoStageDeterrency().present());
        relAttractiveness.initTwoStage(elt.getTwoStageDeterrency().get(), maxInsecticide, holeIndexMax,
                                       holeIndexMaxPresent);
    }
    // pre/postprandial killing: identical ladders (same holeIndexMax message),
    // differing only in target/accessor/name — present()/get() call order kept.
    auto init_kill = [&](factors::SurvivalFactor& target, const auto& eltOpt, const auto& logitOpt,
                         const char* eltName) {
        if (eltOpt.present()) {
            target.init(eltOpt.get(), maxInsecticide, eltName, false);
        } else {
            assert(logitOpt.present());
            if (!holeIndexMax) {
                throw util::xml_scenario_error(
                    "ITN.description.anophelesParams: holeIndexMax required when using logit killing effect");
            }
            target.initLogit(logitOpt.get(), holeIndexMax, false);
        }
    };
    init_kill(preprandialKillingEffect, elt.getPreprandialKillingEffect(), elt.getPreprandialKillingEffectLogit(),
              "ITN.description.anophelesParams.preprandialKillingFactor");
    init_kill(postprandialKillingEffect, elt.getPostprandialKillingEffect(), elt.getPostprandialKillingEffectLogit(),
              "ITN.description.anophelesParams.postprandialKillingFactor");
    if (elt.getFecundityReduction().present()) {
        relFecundityEffect.init(elt.getFecundityReduction().get(), maxInsecticide,
                                "ITN.description.anophelesParams.fecundityReduction", false);
    } else if (elt.getFecundityReductionLogit().present()) {
        if (!holeIndexMax) {
            throw util::xml_scenario_error(
                "ITN.description.anophelesParams: holeIndexMax required when using logit fecundity effect");
        }
        relFecundityEffect.initLogit(elt.getFecundityReductionLogit().get(), holeIndexMax, false);
    } else {
        relFecundityEffect.init1();
    }
    // Nets only affect people while they're using the net. NOTE: we may want
    // to revise this at some point (heterogeneity, seasonal usage patterns).
    set_protected(proportionUse, elt.getPropActive(), proportionProtected, proportionUnprotected);
}

IRSParams::IRSParams(RunContext& ctx, ComponentId /*id*/, const scnXml::IRSDescription& elt,
                     const map<string, size_t>& species_name_map) {
    initialInsecticide.setParams(elt.getInitialInsecticide());
    const double maxProp = 0.999; //NOTE: this could be exposed in XML, but probably doesn't need to be
    maxInsecticide = R::qnorm5(maxProp, initialInsecticide.getMu(), initialInsecticide.getSigma(), true, false);
    insecticideDecay = DecayFunction::makeObject(ctx.clock, elt.getInsecticideDecay(),
                                                 "interventions.human.IRS.description.insecticideDecay");
    // assume usage modifier is 100% if none is specified
    const double propUse = read_propUse(elt.getUsage(), "ITN.description.usage: must be within range [0,1]");
    per_species(species, "IRS intervention", species_name_map, elt.getAnophelesParams(), propUse, maxInsecticide);
}

void IRSParams::IRSAnopheles::init(const scnXml::IRSDescription::AnophelesParamsType& elt, double proportionUse,
                                   double maxInsecticide) {
    _relativeAttractiveness.init(elt.getDeterrency());
    _preprandialKillingEffect.init(elt.getPreprandialKillingEffect(), false, maxInsecticide);
    _postprandialKillingEffect.init(elt.getPostprandialKillingEffect(), true, maxInsecticide);
    if (elt.getFecundityReduction().present()) {
        _fecundityEffect.init(elt.getFecundityReduction().get(), false /*TODO: err msg*/, maxInsecticide);
    } else {
        _fecundityEffect.init1();
    }
    // Simpler version of ITN usage/action:
    set_protected(proportionUse, elt.getPropActive(), proportionProtected, proportionUnprotected);
}

void IRSParams::IRSAnopheles::RelativeAttractiveness::init(const scnXml::IRSDeterrency& elt) {
    double PF = elt.getInsecticideFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    if (!(PF > 0.0)) {
        throw util::xml_scenario_error(
            "IRS.description.anophelesParams.relativeAvailability: expected insecticideFactor to be positive.");
    }

    /* We need to ensure the relative availability is non-negative. However,
     * since it's an exponentiated value, it always will be.
     * 
     * If we don't want IRS to be able to increase transmission, the following
     * limits could also be applied. In general, however, there is no reason
     * IRS couldn't make individuals more attractive to mosquitoes.
     * 
     * To ensure relative availability is at most one: relative availability is
     *  exp( log(PF)*p ) = PF^p
     * where PF is the insecticide factor, with p∈[0,1] defined as:
     *  p=1−exp(-insecticideContent*insecticideScalingFactor).
     * We therefore just need PF ≤ 1. */
    //     if( !( PF <= 1.0 ) ) {
    // Potentially warn about this... but not necessary since making humans
    // more attractive isn't really an issue.
    //     }
    assert((std::isnan)(lPF)); // double init
    lPF = log(PF);
}

void IRSParams::IRSAnopheles::SurvivalFactor::init(const scnXml::IRSKillingEffect& elt, bool postPrandial,
                                                   double maxInsecticide) {
    BF = elt.getBaseFactor();
    PF = elt.getInsecticideFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    invBaseSurvival = 1.0 / (1.0 - BF);
    if (!(BF >= 0.0 && BF < 1.0)) {
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial ? "post" : "pre")
            << "killingFactor: expected baseFactor to be in range [0,1]";
        throw util::xml_scenario_error(msg.str());
    }
    if (!(insecticideScaling >= 0.0)) {
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial ? "post" : "pre")
            << "killingFactor: expected scaling factor to be non-negative";
        throw util::xml_scenario_error(msg.str());
    }

    /* We want the calculated survival factor (1−K)/(1−BF) to be in the range
    [0,1] where K is the killing factor: K=BF+PF×p, with p defined as:
    p=1−exp(-insecticideContent×insecticideScalingFactor).
    
    By its nature, insecticideContent ≥ 0. We restrict
        insecticideScalingFactor ≥ 0
    which implies that p lies in the range [0,1]. We also know the base
    survival factor, 1−BF, is in the range [0,1].
    
    To make sure the survival factor is not negative we need (1−K)/(1−BF) ≥ 0.
    Since 1−BF > 0 we need 1−K ≥ 0, which, substituting K, gives us
        BF + PF×p ≤ 1	(1)
    We also want to make sure the survival factor is not greater than one (since
    IRS shouldn't increase mosquito survival), (1−K)/(1−BF) ≤ 1 or equivalently
    1-K ≤ 1-BF or K ≥ BF, which, substituting K, yields
        PF×p ≥ 0		(2)
    Since p ≥ 0, ensuring PF is non-negative is enough for (2).
    
    The maximum value for p depends on the maximum insecticide content; denote
    pmax = max(p). In this extreme case (1) becomes:
        BF + PF×pmax ≤ 1	(3)
    
    As with the ITN model, we impose a maximum value on the initial insecticide
    content, Pmax, such that the probability of sampling a value from our
    parameterise normal distribution greater than Pmax is 0.001. */
    double pmax = 1.0 - exp(-maxInsecticide * insecticideScaling);
    if (!(PF >= 0.0 && BF + PF * pmax <= 1.0)) {
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial ? "post" : "pre")
            << "killingFactor: expected insecticideFactor≥0, baseFactor+" << pmax << "×insecticideFactor≤1";
        throw util::xml_scenario_error(msg.str());
    }
}

GVIParams::GVIParams(RunContext& ctx, ComponentId /*id*/, const scnXml::GVIDescription& elt,
                     const map<string, size_t>& species_name_map) {
    decay = DecayFunction::makeObject(ctx.clock, elt.getDecay(), "interventions.human.vector.decay");
    // assume usage modifier is 100% if none is specified
    const double propUse = read_propUse(elt.getUsage(), "GVI.description.usage: must be within range [0,1]");
    per_species(species, "GVI intervention", species_name_map, elt.getAnophelesParams(), propUse);
}

void GVIParams::GVIAnopheles::init(const scnXml::GVIDescription::AnophelesParamsType& elt, double proportionUse) {
    assert((std::isnan)(deterrency)); // double init
    deterrency = elt.getDeterrency().present() ? elt.getDeterrency().get().getValue() : 0.0;
    preprandialKilling =
        elt.getPreprandialKillingEffect().present() ? elt.getPreprandialKillingEffect().get().getValue() : 0.0;
    postprandialKilling =
        elt.getPostprandialKillingEffect().present() ? elt.getPostprandialKillingEffect().get().getValue() : 0.0;
    fecundityReduction = elt.getFecundityReduction().present() ? elt.getFecundityReduction().get().getValue() : 0.0;
    // Simpler version of ITN usage/action:
    set_protected(proportionUse, elt.getPropActive(), proportionProtected, proportionUnprotected);
}

VaccineParams::VaccineParams(RunContext& ctx, ComponentId component, const scnXml::VaccineDescription& vd,
                             Vaccine::Types type)
    : type(type), decayFunc(DecayFunction::makeObject(ctx.clock, vd.getDecay(), "decay")),
      efficacyB(vd.getEfficacyB().getValue()) {
    if (type == Vaccine::BSV && ModelOptions::option(util::VIVAX_SIMPLE_MODEL))
        throw util::unimplemented_exception("blood stage vaccines (BSV) cannot be used with vivax model");

    opt_vaccine_genotype = util::ModelOptions::option(util::VACCINE_GENOTYPE);

    if (ctx.vaccineReportComponent == ComponentId::wholePop() /* the initial value */) {
        // set to the first type described
        ctx.vaccineReportComponent = component;
    }

    const scnXml::VaccineDescription::InitialEfficacySequence ies = vd.getInitialEfficacy();
    const scnXml::VaccineDescription::PhenotypeSequence ps = vd.getPhenotype();

    //cout << ies.size() << " " << ps.size() << " " << ctx.genotypes.N() << endl;

    if (opt_vaccine_genotype == false && ies.size() == 0)
        throw util::xml_scenario_error("Vaccine: initialEfficacy is required but is missing");

    if (opt_vaccine_genotype == true && ies.size() != 0)
        throw util::xml_scenario_error(
            "Vaccine: initialEfficacy outside of phenotype definition with VACCINE_GENOTYPE option enabled is not allowed");

    /*    if(opt_vaccine_genotype == true && ies.size() == 0 && ps.size() != ctx.genotypes.N())
        throw util::xml_scenario_error( "Vaccine: Phenotype sequence must reference all genotypes" );
    */
    if (opt_vaccine_genotype == true) {
        perGenotypeInitialMeanEfficacy.resize(ctx.genotypes.N());

        for (size_t i = 0; i < ps.size(); ++i) {
            //cout << "Phenotype: " << ps[i].getName() << endl;

            const scnXml::VaccinePhenotype::RestrictionSequence restrictions = ps[i].getRestriction();
            const scnXml::VaccinePhenotype::InitialEfficacySequence vpies = ps[i].getInitialEfficacy();

            for (size_t j = 0; j < restrictions.size(); ++j) {
                //cout << "\tonLocus: " << restrictions[j].getOnLocus() << ", " << "toAllelle: " << restrictions[j].getToAllele() << endl;
                //cout << "\tmatching allelle code: " << ctx.genotypes.findAlleleCode(restrictions[j].getOnLocus(), restrictions[j].getToAllele()) << endl;

                uint32_t allele_code =
                    ctx.genotypes.findAlleleCode(restrictions[j].getOnLocus(), restrictions[j].getToAllele());

                if (allele_code == numeric_limits<uint32_t>::max()) {
                    throw util::xml_scenario_error("phenotype has "
                                                   "restriction on locus " +
                                                   string(restrictions[j].getOnLocus()) + " allele " +
                                                   string(restrictions[j].getToAllele()) +
                                                   " but this locus/allele "
                                                   "has not been defined in parasiteGenetics section");
                }

                if (allele_code >= ctx.genotypes.N())
                    throw util::xml_scenario_error(string("Vaccine phenotype: Invalid Locus ") +
                                                   restrictions[j].getOnLocus() + string(" and allele ") +
                                                   restrictions[j].getToAllele());

                if (perGenotypeInitialMeanEfficacy[allele_code].size() != 0)
                    throw util::xml_scenario_error(string(
                        "Vaccine phenotype: there are multiple restrictions locus " +
                        string(restrictions[j].getOnLocus() + string(" and allele ") + restrictions[j].getToAllele())));

                if (vpies.size() == 0)
                    throw util::xml_scenario_error(
                        string("Vaccine phenotype: initialEfficacy is required but is missing for phenotype: " +
                               string(restrictions[j].getOnLocus() + string(" ") + restrictions[j].getToAllele())));

                perGenotypeInitialMeanEfficacy[allele_code].resize(vpies.size());

                for (size_t k = 0; k < vpies.size(); ++k)
                    perGenotypeInitialMeanEfficacy[allele_code][k] = vpies[k].getValue();
            }
        }

        // Check that all genotypes have been adressed
        for (size_t i = 0; i < ctx.genotypes.N(); ++i) {
            if (perGenotypeInitialMeanEfficacy[i].size() == 0)
                throw util::xml_scenario_error(
                    "Vaccine phenotype: missing restriction for at least one locus and allele pair");
        }
    } else {
        perGenotypeInitialMeanEfficacy.resize(1);
        perGenotypeInitialMeanEfficacy[0].resize(ies.size());
        for (size_t i = 0; i < ies.size(); ++i)
            perGenotypeInitialMeanEfficacy[0][i] = ies[i].getValue();
    }

    /* initialMeanEfficacy.resize (ies.size());
    for(size_t i = 0; i < initialMeanEfficacy.size(); ++i)
        initialMeanEfficacy[i] = ies[i].getValue();*/
}

// ———  InterventionManager  ———

/// Register a per-type component parameter payload under id.
template <typename T> static void reg_at(std::vector<const T*>& v, ComponentId id, const T* p) {
    if (v.size() <= id.id)
        v.resize(id.id + 1, nullptr);
    assert(v[id.id] == nullptr);
    v[id.id] = p;
}

void InterventionManager::init(RunContext& ctx, const scnXml::Interventions& intervElt,
                               Transmission::TransmissionModel& transmission) {
    ctx.nextTimed = 0;

    // changeHS and changeEIR timed deployments: identical loops differing only
    // in payload kind and error prefix (former duplicated blocks, merged).
    auto timed_clone = [&ctx](const auto& seq, TimedDeploymentKind kind, const char* errPrefix, const auto& mk) {
        for (const auto& it : seq) {
            try {
                SimTime date = readDate_steps(ctx.clock, it.getTime());
                check_timed_date(ctx, date);
                ctx.timed->push_back(timed_deployment(date, kind, mk(it)));
            } catch (const util::format_error& e) {
                throw util::xml_scenario_error(string(errPrefix).append(e.message()));
            }
        }
    };
    if (intervElt.getChangeHS().present())
        timed_clone(intervElt.getChangeHS().get().getTimedDeployment(), TimedDeploymentKind::CHANGE_HS,
                    "interventions/changeHS/timedDeployment/time: ", [](const auto& it) {
                        return ChangeHsDeployParams{it._clone()};
                    });
    if (intervElt.getChangeEIR().present())
        timed_clone(intervElt.getChangeEIR().get().getTimedDeployment(), TimedDeploymentKind::CHANGE_EIR,
                    "interventions/changeEIR/timedDeployment/time: ", [](const auto& it) {
                        return ChangeEirDeployParams{it._clone()};
                    });
    // species_index_map is not available with the non-vector model or
    // non-dynamic mode, so setting it (lazily) also checks sim mode:
    map<string, size_t> species_index_map;

    {
        Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
        if (vectorModel != nullptr)
            species_index_map = vectorModel->speciesIndex;
    }

    if (intervElt.getHuman().present()) {
        const scnXml::HumanInterventions& human = intervElt.getHuman().get();

        // 1. Read components
        for (const scnXml::HumanInterventionComponent& component : human.getComponent()) {
            if (ctx.identifierMap.count(component.getId()) > 0) {
                ostringstream msg;
                msg << "The id attribute of interventions.human.component elements must be unique; found \""
                    << component.getId() << "\" twice.";
                throw util::xml_scenario_error(msg.str());
            }
            ComponentId id(ctx.humanComponents->size()); // i.e. index of next item
            ctx.identifierMap.insert(make_pair(component.getId(), id));

            SimTime expireAfter = sim::future();
            if (component.getSubPopRemoval().present()) {
                const scnXml::SubPopRemoval& removeOpts = component.getSubPopRemoval().get();
                if (removeOpts.getOnFirstBout()) {
                    ctx.removeAtIds[SubPopRemove::ON_FIRST_BOUT].push_back(id);
                }
                if (removeOpts.getOnFirstInfection()) {
                    ctx.removeAtIds[SubPopRemove::ON_FIRST_INFECTION].push_back(id);
                }
                if (removeOpts.getOnFirstTreatment()) {
                    ctx.removeAtIds[SubPopRemove::ON_FIRST_TREATMENT].push_back(id);
                }
                if (removeOpts.getAfterYears().present()) {
                    expireAfter = ctx.clock.from_years_n(removeOpts.getAfterYears().get());
                }
            }

            ComponentConfig comp{id, Component::RECRUIT_ONLY, expireAfter, RecruitOnlyParams{}};
            // Each branch: present() then get(), same call sequence as before; first child wins.
            if (auto& e = component.getScreen(); e.present())
                comp.set<ScreenParams>(Component::SCREEN, ctx, e.get());
            else if (auto& e = component.getTreatSimple(); e.present())
                comp.set<TreatSimpleParams>(Component::TREAT_SIMPLE, ctx, e.get());
            else if (auto& e = component.getTreatPKPD(); e.present())
                comp.set<TreatPKPDParams>(Component::TREAT_PKPD, ctx, e.get());
            else if (auto& e = component.getDecisionTree(); e.present())
                comp.set<DecisionTreeParams>(Component::CM_DT, ctx, e.get());
            else if (auto& e = component.getPEV(); e.present())
                comp.set<VaccineParams>(Component::PEV, ctx, id, e.get(), Vaccine::PEV);
            else if (auto& e = component.getBSV(); e.present())
                comp.set<VaccineParams>(Component::BSV, ctx, id, e.get(), Vaccine::BSV);
            else if (auto& e = component.getTBV(); e.present())
                comp.set<VaccineParams>(Component::TBV, ctx, id, e.get(), Vaccine::TBV);
            else if (auto& e = component.getITN(); e.present())
                comp.set<ITNParams>(Component::ITN, ctx, id, e.get(), species_index_map);
            else if (auto& e = component.getIRS(); e.present())
                comp.set<IRSParams>(Component::IRS, ctx, id, e.get(), species_index_map);
            else if (auto& e = component.getGVI(); e.present())
                comp.set<GVIParams>(Component::GVI, ctx, id, e.get(), species_index_map);
            else if (component.getRecruitmentOnly().present()) {
            } else if (component.getClearImmunity().present())
                comp.kind = Component::CLEAR_IMMUNITY; // no payload (cf. RecruitOnly)
            else
                throw util::xml_scenario_error("expected intervention.human.component element to have a "
                                               "child, didn't find it (perhaps I need updating)");
            ctx.humanComponents->push_back(std::move(comp));
        }

        // Fill per-type parameter registries. Payload pointers stay valid:
        // humanComponents is append-only and never shrinks (7.2b).
        for (auto& c : (*ctx.humanComponents)) {
            switch (c.kind) {
            case Component::ITN:
                reg_at(ctx.itnComponents, c.id, &std::get<ITNParams>(c.data));
                break;
            case Component::IRS:
                reg_at(ctx.irsComponents, c.id, &std::get<IRSParams>(c.data));
                break;
            case Component::GVI:
                reg_at(ctx.gviComponents, c.id, &std::get<GVIParams>(c.data));
                break;
            case Component::PEV:
            case Component::BSV:
            case Component::TBV:
                reg_at(ctx.vaccineComponents, c.id, &std::get<VaccineParams>(c.data));
                break;
            default:
                break;
            }
        }

        // 2. Read the list of deployments
        // Sub-population restriction; read identically in the continuous and
        // timed blocks below (getComponentId before getComplement, as before):
        auto sub_pop = [&ctx](const auto& dep) {
            ComponentId subPop = ComponentId::wholePop();
            bool complement = false;
            if (dep.getRestrictToSubPop().present()) {
                const string& subPopStr = dep.getRestrictToSubPop().get().getId();
                subPop = getComponentId(ctx, subPopStr);
                complement = dep.getRestrictToSubPop().get().getComplement();
            }
            return std::make_pair(subPop, complement);
        };
        for (const scnXml::Deployment& elt : human.getDeployment()) {
            // 2.a intervention components
            auto intervention =
                shared_ptr<HumanIntervention>(new HumanIntervention(ctx, elt.getComponent(), elt.getCondition()));

            // 2.b intervention deployments
            for (auto ctsIt = elt.getContinuous().begin(); ctsIt != elt.getContinuous().end(); ++ctsIt) {
                auto [subPop, complement] = sub_pop(*ctsIt);
                const scnXml::ContinuousList::DeploySequence& ctsSeq = ctsIt->getDeploy();
                for (auto it2 = ctsSeq.begin(), end2 = ctsSeq.end(); it2 != end2; ++it2) {
                    try {
                        SimTime begin = ctx.clock.start_date;
                        if (it2->getBegin().present()) {
                            begin = readDate_steps(ctx.clock, it2->getBegin().get());
                        }
                        SimTime end = sim::future();
                        if (it2->getEnd().present()) {
                            end = readDate_steps(ctx.clock, it2->getEnd().get());
                        }
                        ctx.continuous->push_back(
                            ContinuousHumanDeployment::read(ctx, begin, end, *it2, intervention, subPop, complement));
                    } catch (const util::format_error& e) {
                        throw util::xml_scenario_error(
                            string("interventions/human/deployment/continuous/deploy: ").append(e.message()));
                    }
                }
            }
            for (auto timedIt = elt.getTimed().begin(); timedIt != elt.getTimed().end(); ++timedIt) {
                auto [subPop, complement] = sub_pop(*timedIt);
                try {
                    multimap<SimTime, const scnXml::MassDeployment*> deployTimes;
                    for (auto it2 = timedIt->getDeploy().begin(), end2 = timedIt->getDeploy().end(); it2 != end2;
                         ++it2) {
                        SimTime date = readDate_steps(ctx.clock, it2->getTime());

                        if (it2->getRepeatStep().present() != it2->getRepeatEnd().present()) {
                            throw util::xml_scenario_error("deploy: use of repeatStep or repeatEnd without other");
                        }
                        if (it2->getRepeatStep().present()) {
                            SimTime step =
                                UnitParse::readDuration(ctx.clock, it2->getRepeatStep().get(), UnitParse::NONE);
                            if (step < ctx.clock.one_TS()) {
                                throw util::xml_scenario_error("deploy: repeatStep must be >= 1");
                            }
                            SimTime end = UnitParse::readDate(ctx.clock, it2->getRepeatEnd().get(), UnitParse::NONE);
                            while (date < end) {
                                deployTimes.insert(make_pair(date, &*it2));
                                date = date + step;
                            }
                        } else {
                            deployTimes.insert(make_pair(date, &*it2));
                        }
                    }
                    SimTime lastDate = sim::never();
                    for (auto deploy = deployTimes.begin(), dpEnd = deployTimes.end(); deploy != dpEnd; ++deploy) {
                        if (deploy->first == lastDate) {
                            ostringstream msg;
                            msg << "Timed deployment of components ";
                            bool first = true;
                            for (auto cp = elt.getComponent().begin(), cpEnd = elt.getComponent().end(); cp != cpEnd;
                                 ++cp) {
                                if (!first) {
                                    msg << ", ";
                                } else {
                                    first = false;
                                }
                                msg << cp->getId();
                            }
                            msg << " has multiple deplyoments at date " << lastDate << " (step "
                                << ctx.clock.inSteps(lastDate - ctx.clock.start_date) << ')';
                            throw util::xml_scenario_error(msg.str());
                        }
                        lastDate = deploy->first;
                    }

                    if (timedIt->getCumulativeCoverage().present()) {
                        const scnXml::CumulativeCoverage& cumCov = timedIt->getCumulativeCoverage().get();
                        ComponentId cumCovComponent = getComponentId(ctx, cumCov.getComponent());
                        for (auto deploy = deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy) {
                            ctx.timed->push_back(timed_deployment(
                                deploy->first, TimedDeploymentKind::CUMULATIVE_HUMAN,
                                CumulativeHumanParams{TimedHumanParams::read(ctx, deploy->first, *deploy->second,
                                                                             intervention, subPop, complement),
                                                      cumCovComponent}));
                        }
                    } else {
                        for (auto deploy = deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy) {
                            ctx.timed->push_back(
                                timed_deployment(deploy->first, TimedDeploymentKind::HUMAN,
                                                 TimedHumanParams::read(ctx, deploy->first, *deploy->second,
                                                                        intervention, subPop, complement)));
                        }
                    }
                } catch (const util::format_error& e) {
                    throw util::xml_scenario_error(
                        string("interventions/human/deployment/timed/deploy/time: ").append(e.message()));
                }
            }
        }
    }
    if (intervElt.getImportedInfections().present()) {
        const scnXml::ImportedInfections& ii = intervElt.getImportedInfections().get();
        ctx.importedInfections->init(ctx, ii);
    }
    // Must come after vaccines are initialised:
    if (intervElt.getInsertR_0Case().present()) {
        throw util::xml_scenario_error("R_0 code is disabled to reduce "
                                       "maintenance. If you need it, either use an old OpenMalaria "
                                       "version (pre 32) or request it be reenabled.");
        // code disabled (former #if 0 timed-deployment block deleted; search for R_0 to find more blocks like this)
    }
    if (intervElt.getUninfectVectors().present()) {
        const scnXml::UninfectVectors& elt = intervElt.getUninfectVectors().get();
        // timed deployments:
        for (const auto& it : elt.getTimedDeployment()) {
            SimTime date = readDate_steps(ctx.clock, it.getTime());
            check_timed_date(ctx, date);
            ctx.timed->push_back(timed_deployment(date, TimedDeploymentKind::UNINFECT_VECTORS));
        }
    }
    if (intervElt.getVectorPop().present()) {
        Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
        if (vectorModel) {
            typedef scnXml::VectorPop::InterventionSequence SeqT;
            const SeqT& seq = intervElt.getVectorPop().get().getIntervention();
            size_t instance = 0;
            for (const scnXml::VectorIntervention& elt : seq) {
                if (elt.getTimed().present()) {
                    initVectorInterv(*vectorModel, elt.getDescription().getAnopheles(), instance, elt.getName());

                    const scnXml::TimedBaseList::DeploySequence& seq = elt.getTimed().get().getDeploy();
                    for (auto it = seq.begin(); it != seq.end(); ++it) {
                        SimTime date = readDate_steps(ctx.clock, it->getTime());
                        check_timed_date(ctx, date);
                        ctx.timed->push_back(
                            timed_deployment(date, TimedDeploymentKind::VECTOR, VectorDeployParams{instance}));
                    }
                    instance++;
                }
            }
        } else
            throw util::xml_scenario_error(nonvec_err);
    }
    if (intervElt.getAddNonHumanHosts().present()) {
        Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
        if (vectorModel) {
            typedef scnXml::AddNonHumanHosts::NonHumanHostsSequence SeqT;
            const SeqT& seq = intervElt.getAddNonHumanHosts().get().getNonHumanHosts();
            for (const scnXml::NonHumanHosts2& elt : seq) {
                if (elt.getTimed().present()) {
                    for (const scnXml::Deploy2& deploy : elt.getTimed().get().getDeploy()) {
                        SimTime date =
                            readDate_steps(ctx.clock, deploy.getTime());
                        SimTime lifespan = UnitParse::readDuration(ctx.clock, deploy.getLifespan(), UnitParse::NONE);
                        ctx.timed->push_back(
                            timed_deployment(date, TimedDeploymentKind::ADD_NHH,
                                             read_add_nhh(ctx, date, elt.getName(), lifespan,
                                                          elt.getDescription().getAnopheles(), transmission)));
                    }
                }
            }
        } else
            throw util::xml_scenario_error(nonvec_err);
    }
    if (intervElt.getNonHumanHostsModifications().present()) {
        Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
        if (vectorModel) {
            typedef scnXml::NonHumanHostsModifications::InterventionSequence SeqT;
            const SeqT& seq = intervElt.getNonHumanHostsModifications().get().getIntervention();
            size_t instance = 0;
            for (const scnXml::NonHumanHostsIntervention& elt : seq) {
                if (elt.getTimed().present()) {
                    const scnXml::DecayFunction& decay = elt.getDecay();
                    const scnXml::TimedBaseList::DeploySequence& seq = elt.getTimed().get().getDeploy();
                    for (auto it = seq.begin(); it != seq.end(); ++it) {
                        SimTime date = readDate_steps(ctx.clock, it->getTime());
                        check_timed_date(ctx, date);
                        nhh_mods_setup(ctx, instance, elt.getNonHumanHostsName(), elt.getDescription().getAnopheles(),
                                       decay, transmission);
                        ctx.timed->push_back(
                            timed_deployment(date, TimedDeploymentKind::NHH_MODS,
                                             NhhModsDeployParams{instance, elt.getNonHumanHostsName()}));
                    }
                    instance++;
                }
            }
        } else
            throw util::xml_scenario_error(nonvec_err);
    }
    if (intervElt.getVectorTrap().present()) {
        Transmission::TransmissionModel* vectorModel = transmission.isVector() ? &transmission : nullptr;
        if (vectorModel) {
            size_t instance = 0;
            for (const scnXml::VectorTrap& trap : intervElt.getVectorTrap().get().getIntervention()) {
                initVectorTrap(*vectorModel, trap.getDescription(), instance, trap.getName());
                if (trap.getTimed().present()) {
                    for (const scnXml::Deploy1& deploy : trap.getTimed().get().getDeploy()) {
                        SimTime date = UnitParse::readDate(ctx.clock, deploy.getTime(), UnitParse::STEPS);
                        double ratio = deploy.getRatioToHumans();
                        SimTime lifespan = UnitParse::readDuration(ctx.clock, deploy.getLifespan(), UnitParse::NONE);
                        check_timed_date(ctx, date);
                        ctx.timed->push_back(timed_deployment(date, TimedDeploymentKind::TRAP,
                                                              TrapDeployParams{instance, ratio, lifespan}));
                    }
                }
                instance += 1;
            }
        } else
            throw util::xml_scenario_error(nonvec_err);
    }

    // lists must be sorted, increasing
    // For reproducability, we need to use stable_sort, not sort.
    // parity: b — stable_sort order = reproducibility contract; comparator ties resolve by insertion order
    stable_sort(ctx.continuous->begin(), ctx.continuous->end(), byDeployTime);
    stable_sort(ctx.timed->begin(), ctx.timed->end(), byDeployTime);

    // make sure the list ends with something always in the future, so we don't
    // have to check ctx.nextTimed is within range:
    ctx.timed->push_back(timed_deployment(sim::future(), TimedDeploymentKind::DUMMY));

    if (util::CommandLine::option(util::CommandLine::PRINT_INTERVENTIONS)) {
        cout << "Continuous deployments:" << endl << "begin\tend\tage\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for (auto it = ctx.continuous->begin(); it != ctx.continuous->end(); ++it) {
            it->print_details(std::cout);
            cout << endl;
        }
        cout << "Timed deployments:" << endl << "time\tmin age\tmax age\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for (auto& deploy : (*ctx.timed)) {
            print_timed_details(std::cout, deploy);
            cout << endl;
        }
        cout << "Human components:" << endl;
        for (auto& component : (*ctx.humanComponents)) {
            print_component(cout, component);
            cout << endl;
        }
    }
}

}
}
