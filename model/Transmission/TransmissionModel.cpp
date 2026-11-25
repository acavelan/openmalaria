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

#include "RunContext.h"
#include "Transmission/TransmissionModel.h"
#include "Population.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Genotypes.h"
#include "interventions/IvState.h" // history-ordered iv folds (E5 hoisted rows)
#include "mon/Continuous.h"
#include "util/vectors.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"

#include <fstream>
#include <map>
#include <cmath>
#include <set>
#include <numeric>

namespace OpenMalaria {

void TransmissionModel::changeEIRIntervention(const EirSeries& eir) {
    // parity: xml/entomology.cpp changeEIR variants — preserve validation and
    // mutation order. Each sum is a left fold in day order; no RNG draws.
    if (isVector() && simulationMode != forcedEIR && simulationMode != transientEIRknown)
        throw xml_scenario_error("changeEIR intervention can only be used in forced EIR mode");
    simulationMode = transientEIRknown;
    if (isVector()) {
        if (speciesIndex.size() != 1)
            throw xml_scenario_error("changeEIR intervention can only be used with a single anopheles species");
    } else if (nSpore != eir.eipDuration) {
        throw xml_scenario_error("change-of-EIR intervention cannot change EIP duration");
    }

    const auto& daily = eir.daily;
    if (daily.empty()) throw xml_scenario_error("no EIRDaily values given");
    const double minEIR = min_EIR_mult * (std::accumulate(daily.begin(), daily.end(), 0.0) / daily.size());
    if (isVector() && daily.size() < static_cast<size_t>(oneYear()))
        throw xml_scenario_error("entomology.anopheles.seasonality.dailyValues insufficient daily data for a year");

    auto& values = isVector() ? species[0]->interventionEIR : interventionEIR;
    vector<int> nDays(ctx.clock.inSteps(fromDays(daily.size() - 1)) + 1, 0);
    values.assign(nDays.size(), 0.0);
    size_t required_days = static_cast<size_t>((ctx.clock.end_date - ctx.clock.start_date) + 1);
    if (daily.size() < required_days) {
        if (!isVector())
            cerr << "Days: " << daily.size() << "\nIntervals: " << nDays.size() << "\nRequired: " << required_days
                 << endl;
        throw xml_scenario_error("Insufficient intervention phase EIR values provided");
    }
    for (SimTime day = zero(), end = fromDays(daily.size()); day < end; day = day + oneDay()) {
        double value = std::max(daily[day], minEIR);
        size_t step = ctx.clock.inSteps(day);
        ++nDays[step];
        values[step] += value;
    }
    for (size_t i = 0; i < values.size(); ++i)
        values[i] *= ctx.clock.one_TS() / static_cast<double>(nDays[i]);
    if (!isVector()) annualEIR = numeric_limits<double>::quiet_NaN();
}

void TransmissionModel::ctsCbN_v0(ostream& stream) {
    for (auto& s : species)
        stream << '\t' << s->getLastN_v0();
}
void TransmissionModel::ctsCbVecStat(VecStat vs, ostream& stream) {
    for (auto& s : species)
        stream << '\t' << s->getLastVecStat(vs);
}
namespace {
/// Mean over the population of a per-host term: left fold in host order, the
/// same FP sequence the three separate callbacks accumulated (cts bytes).
template <typename Term> double meanOver(const vector<Human>& population, Term term) {
    double total = 0.0;
    for (const Human& human : population)
        total += term(human);
    return total / population.size();
}

/** E5: per-host view of the transmission row containers, bound once per host.
 * Rows stay hoisted; the history-ordered folds are shared with PerHost via
 * IvState.h (parity: c). Bind after any iv-deploying call (e.g. PerHost::update). */
struct IvRows {
    explicit IvRows(const PerHost& host)
        : ctx(*host.ctx), iv(host.cols->iv[host.row]), anoph(host.anoph()) {}
    double entoAvailabilityHetVecItv(size_t species) const {
        return iv_rel_attractiveness(ctx, iv, anoph[species].ento_avail, species);
    }
    double probMosqBiting(size_t species) const {
        return iv_preprandial_survival(ctx, iv, anoph[species].prob_biting, species);
    }
    double probMosqResting(size_t species) const {
        return iv_postprandial_survival(ctx, iv, anoph[species].prob_resting, species);
    }
    double relMosqFecundity(size_t species) const { return iv_rel_fecundity(ctx, iv, 1.0, species); }

    RunContext& ctx;
    const std::vector<IvRecord>& iv;
    std::span<const HostSpecies> anoph;
};
}

void TransmissionModel::ctsCbAlpha(const vector<Human>& population, ostream& stream) {
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << meanOver(population, [&](const Human& human) {
            return human.perHostTransmission.entoAvailabilityFull(i, inYears(human.age(ctx.clock.now())));
        });
}
void TransmissionModel::ctsCbP_B(const vector<Human>& population, ostream& stream) {
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t'
               << meanOver(population, [&](const Human& human) { return human.perHostTransmission.probMosqBiting(i); });
}
void TransmissionModel::ctsCbP_CD(const vector<Human>& population, ostream& stream) {
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << meanOver(population, [&](const Human& human) {
            return human.perHostTransmission.probMosqResting(i);
        });
}
void TransmissionModel::ctsCbResAvailability(ostream& stream) {
    for (auto& s : species)
        stream << '\t' << s->getResAvailability();
}
void TransmissionModel::ctsCbResRequirements(ostream& stream) {
    for (auto& s : species)
        stream << '\t' << s->getResRequirements();
}

const string& reverseLookup(const map<string, size_t>& m, size_t i) {
    auto it = std::find_if(m.begin(), m.end(), [i](const auto& p) { return p.second == i; });
    if (it == m.end())
        throw TRACED_EXCEPTION_DEFAULT("reverseLookup: key not found"); // shouldn't ever happen
    return it->first;
}

TransmissionModel::TransmissionModel(RunContext& ctx, vector<double> initEIR, int interventionMode,
                                     vector<std::unique_ptr<AnophelesModel>> speciesList,
                                     vector<std::unique_ptr<EmergenceRateEstimator>> emergenceRateEstimatorsList,
                                     map<string, size_t> speciesIndexList)
    : TransmissionModel(ctx, std::move(initEIR), interventionMode, Kind::Vector, std::move(speciesList),
                        std::move(emergenceRateEstimatorsList), std::move(speciesIndexList)) {
    // set VACCINE_GENOTYPE option
    opt_vaccine_genotype = ModelOptions::option(VACCINE_GENOTYPE);

    // Calculate total annual EIR
    annualEIR = vectors::sum(initialisationEIR);

    // -----  Continuous reporting  -----
    // One column per species, in the order specified in XML, so outputting is
    // just iterating the list. Unfortunately, we need to reverse-lookup the name.
    auto ctsHeads = [this](const char* prefix) {
        ostringstream s;
        for (size_t i = 0; i < species.size(); ++i)
            s << prefix << reverseLookup(speciesIndex, i) << ")";
        return s.str();
    };

    Continuous::registerCallback(ctx, "N_v0", ctsHeads("\tN_v0("), std::bind(&TransmissionModel::ctsCbN_v0, this, _1));
    Continuous::registerCallback(ctx, "P_A", ctsHeads("\tP_A("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, PA, _1));
    Continuous::registerCallback(ctx, "P_Amu", ctsHeads("\tP_Amu("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, PAmu, _1));
    Continuous::registerCallback(ctx, "P_A1", ctsHeads("\tP_A1("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, PA1, _1));
    Continuous::registerCallback(ctx, "P_Ah", ctsHeads("\tP_Ah("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, PAh, _1));
    Continuous::registerCallback(ctx, "P_df", ctsHeads("\tP_df("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, PDF, _1));
    Continuous::registerCallback(ctx, "P_dif", ctsHeads("\tP_dif("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, PDIF, _1));
    Continuous::registerCallback(ctx, "N_v", ctsHeads("\tN_v("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, NV, _1));
    Continuous::registerCallback(ctx, "O_v", ctsHeads("\tO_v("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, OV, _1));
    Continuous::registerCallback(ctx, "S_v", ctsHeads("\tS_v("),
                                 std::bind(&TransmissionModel::ctsCbVecStat, this, SV, _1));

    // availability to mosquitoes relative to other humans, excluding age factor
    Continuous::registerCallback(ctx, "alpha", ctsHeads("\talpha_i("),
                                 std::bind(&TransmissionModel::ctsCbAlpha, this, _1, _2));
    Continuous::registerCallback(ctx, "P_B", ctsHeads("\tP_B("), std::bind(&TransmissionModel::ctsCbP_B, this, _1, _2));
    Continuous::registerCallback(ctx, "P_C*P_D", ctsHeads("\tP_C*P_D("),
                                 std::bind(&TransmissionModel::ctsCbP_CD, this, _1, _2));

    Continuous::registerCallback(ctx, "resource availability", ctsHeads("\tres avail("),
                                 std::bind(&TransmissionModel::ctsCbResAvailability, this, _1));
    Continuous::registerCallback(ctx, "resource requirements", ctsHeads("\tres req("),
                                 std::bind(&TransmissionModel::ctsCbResRequirements, this, _1));
}

void TransmissionModel::init2_vector(const vector<Human>& population) {
    double sumRelativeAvailability = 0.0;
    for (const Human& human : population) {
        sumRelativeAvailability +=
            human.perHostTransmission.relativeAvailabilityAge(inYears(human.age(ctx.clock.now())));
    }
    int popSize = population.size();
    // value should be unimportant when no humans are available, though inf/nan is not acceptable
    if (popSize > 0) {
        meanAvail = sumRelativeAvailability / popSize; // mean-rel-avail
    }

    for (size_t i = 0; i < speciesIndex.size(); ++i) {
        double sum_avail = 0.0;
        double sigma_f = 0.0;
        double sigma_df = 0.0;
        double sigma_dff = 0.0;

        for (const Human& human : population) {
            const PerHost& host = human.perHostTransmission;
            double prod = host.entoAvailabilityFull(i, inYears(human.age(ctx.clock.now())));
            sum_avail += prod;
            prod *= host.probMosqBiting(i);
            sigma_f += prod;
            sigma_df += prod * host.probMosqResting(i);
            sigma_dff += prod * host.probMosqResting(i) * host.relMosqFecundity(i);
        }

        species[i]->init2(population.size(), meanAvail, sum_avail, sigma_f, sigma_df, sigma_dff);
    }
    simulationMode = forcedEIR; // now we should be ready to start
}

void TransmissionModel::scaleEIR_vector(double factor) {
    for (auto& s : species)
        s->scaleEIR(factor);
    vectors::scale(initialisationEIR, factor);
    annualEIR = vectors::sum(initialisationEIR);
}

SimTime TransmissionModel::initIterate_vector() {
    if (interventionMode != dynamicEIR) {
        // allow forcing equilibrium mode like with non-vector model
        return zero(); // no initialization to do
    }

    // Fitting is done
    if (initIterations < 0) {
        simulationMode = dynamicEIR;
        return zero();
    }

    if (++initIterations > 30) {
        throw TRACED_EXCEPTION("Transmission warmup exceeded 30 iterations!", Error::VectorWarmup);
    }

    bool needIterate = false;
    for (size_t i = 0; i < speciesIndex.size(); ++i) {
        needIterate = emergenceRateEstimators[i]->estimate(*species[i], laggedKappa, meanAvail);
        species[i]->initIterate();
        if (needIterate) break;
    }

    if (needIterate) {
        // stabilization + 5 years data-collection time:
        return oneYear() + fromYearsI(5);
    } else {
        // One year stabilisation, then we're finished:
        initIterations = -1;
        return oneYear();
    }
}

void TransmissionModel::calculateEIR_vector(Human& human, double ageYears, vector<double>& EIR_i,
                                            vector<double>& EIR_l) const {

    PerHost& host = human.perHostTransmission;
    host.update(human);
    const IvRows iv{host}; // E5: rows bound once per host (see IvRows)
    if (simulationMode == transientEIRknown) {
        EIR_i.assign(1, 0.0);
        EIR_l.assign(1, 0.0);

        for (size_t i = 0; i < speciesIndex.size(); ++i) {
            /* Calculates EIR per individual (hence N_i == 1).
             *
             * See comment in AnophelesModel::advancePeriod for method. */
            double eir =
                species[i]->getInterventionEIR() * host.relativeAvailabilityAge(ageYears) * iv.entoAvailabilityHetVecItv(i);
            recordStat(measure("innoculationsPerAgeGroup"), human, eir);
            recordStat(measure("innoculationsPerVector"), human, eir, i);
            // cout << host.availBite(i) << endl;
            // cout << "species: " << eir << " " << host.entoAvailabilityHetVecItv(i) << " " << host.entoAvailabilityHetVecItv(i) * (1000*0.00112608) << endl;
            EIR_l[0] += eir;
        }
    } else if (simulationMode == forcedEIR) {
        // double eir = initialisationEIR[ctx.clock.moduloYearSteps(ctx.clock.ts0())] * host.relativeAvailabilityHetAge(ageYears);

        // cout << "global: " << eir << " " <<  endl;
        // recordStat(measure("innoculationsPerAgeGroup"), human, eir);
        // EIR.assign(1, eir);

        EIR_i.assign(1, 0.0);
        EIR_l.assign(1, 0.0);

        for (size_t i = 0; i < speciesIndex.size(); ++i) {
            /* Calculates EIR per individual (hence N_i == 1).
             *
             * See comment in AnophelesModel::advancePeriod for method. */
            double eir =
                species[i]->getInitPartialEIR() * host.relativeAvailabilityAge(ageYears) * iv.entoAvailabilityHetVecItv(i);
            recordStat(measure("innoculationsPerAgeGroup"), human, eir);
            recordStat(measure("innoculationsPerVector"), human, eir, i);
            // cout << host.availBite(i) << endl;
            // cout << "species: " << eir << " " << host.entoAvailabilityHetVecItv(i) << " " << host.entoAvailabilityHetVecItv(i) * (1000*0.00112608) << endl;
            EIR_l[0] += eir;
        }
    } else {
        EIR_i.assign(ctx.genotypes.N(), 0.0);
        EIR_l.assign(ctx.genotypes.N(), 0.0);
        const double ageFactor = host.relativeAvailabilityAge(ageYears);
        for (size_t i = 0; i < speciesIndex.size(); ++i) {
            const vector<double>& partialEIR_i = species[i]->partialEIR_i;
            const vector<double>& partialEIR_l = species[i]->partialEIR_l;

            assert(EIR_i.size() == partialEIR_i.size());
            assert(EIR_l.size() == partialEIR_l.size());

            // if ((std::isnan)(vectors::sum(partialEIR_i))) { cerr << "partialEIR is not a number; " << i << endl; }

            /* Calculates EIR per individual (hence N_i == 1).
             *
             * See comment in AnophelesModel::advancePeriod for method. */
            double entoFactor = ageFactor * (iv.entoAvailabilityHetVecItv(i) * iv.probMosqBiting(i));
            for (size_t g = 0; g < EIR_i.size(); ++g) {
                double eir_i = partialEIR_i[g] * entoFactor;
                double eir_l = partialEIR_l[g] * entoFactor;
                recordStat(measure("innoculationsPerAgeGroup"), human, eir_i + eir_l, 0, g);
                recordStat(measure("innoculationsPerVector"), human, eir_i + eir_l, i, g);
                EIR_i[g] += eir_i;
                EIR_l[g] += eir_l;
            }
        }
    }
}

// Every Global::interval days:
// Single-pass host-major accumulation (old upstream loop shape): each accumulator
// receives its addends in ascending host index order, exactly as the pass-2 ordered
// reductions of the two-pass slab version did ⇒ bit-exact.
void TransmissionModel::vectorUpdate_vector(const vector<Human>& population) {
    const size_t nGenotypes = ctx.genotypes.N();
    const size_t nSpec = speciesIndex.size();
    std::vector<double> probTransmission_i, probTransmission_l;
    std::vector<double> sum_avail(nSpec);
    std::vector<double> sigma_df(nSpec);
    std::vector<double> sigma_dff(nSpec);
    std::vector<std::vector<double>> sigma_dif_i(nSpec, std::vector<double>(nGenotypes)),
        sigma_dif_l(nSpec, std::vector<double>(nGenotypes));

    for (const Human& human : population) {
        const PerHost& host = human.perHostTransmission;
        WHInterface& whm = *human.withinHostModel;

        // p7b B4 hoist: vaccine_factor is pure over the host's vac column (frozen
        // while the population is scanned) and returns 1.0 when it is empty — the
        // per-species × per-genotype calls collapse to one check per host.
        const bool anyVac = human.has_vaccines();
        probTransmission_i.assign(nGenotypes, 0.0);
        probTransmission_l.assign(nGenotypes, 0.0);
        whm.probTransmissionToMosquito(probTransmission_i, probTransmission_l);

        // NOTE: calculate availability relative to age at end of time step;
        // not my preference but consistent with TransmissionModel::getEIR().
        // TODO: even stranger since probTransmission comes from the previous time step
        // PR-480 hoist: the age factor is species-independent (was re-evaluated per
        // species inside entoAvailabilityFull). Operand order kept (PerHost.h:167:
        // entoAvailabilityHetVecItv(s) * relativeAvailabilityAge(ageYears)) ⇒ FP-exact.
        const double ageFactor = host.relativeAvailabilityAge(inYears(human.age(ctx.clock.ts1())));
        // E5 hoist: bind this host's row containers once per host instead of
        // per species x getter. Each out-of-line getter re-did cols->col[row]
        // (bounds-checked vector<vector<double>> index) plus its own pass over
        // cols->iv[row]; per host that's 5 row lookups x nSpec. The iv product
        // runs below keep list order and multiply order, so each value is the
        // same FP expression as the getter returned (parity: c).
        const IvRows iv{host};
        for (size_t s = 0; s < nSpec; ++s) {
            const double avail = iv.entoAvailabilityHetVecItv(s) * ageFactor;
            sum_avail[s] += avail;
            const double df = avail * iv.probMosqBiting(s) * iv.probMosqResting(s);
            sigma_df[s] += df;
            for (size_t g = 0; g < nGenotypes; ++g) {
                const double tbv = anyVac ? human.vaccine_factor(Vaccine::TBV, opt_vaccine_genotype ? g : 0) : 1.0;
                sigma_dif_i[s][g] += df * probTransmission_i[g] * tbv;
                sigma_dif_l[s][g] += df * probTransmission_l[g] * tbv;
            }
            sigma_dff[s] += df * iv.relMosqFecundity(s);
        }
    }

    for (size_t s = 0; s < nSpec; ++s) {
        species[s]->advancePeriod(sum_avail[s], sigma_df[s], sigma_dif_i[s], sigma_dif_l[s], sigma_dff[s],
                                  simulationMode == dynamicEIR);
    }
}

extern const string vec_mode_err = "vector interventions can only be used in "
                                   "dynamic transmission mode (mode=\"dynamic\")";

void TransmissionModel::deployVectorPopInterv(size_t instance) {
    if (interventionMode != dynamicEIR) throw xml_scenario_error(vec_mode_err);
    for (auto& s : species)
        s->deployVectorPopInterv(m_rng, instance);
}
void TransmissionModel::deployVectorTrap(size_t instance, double number, SimTime lifespan) {
    if (interventionMode != dynamicEIR) throw xml_scenario_error(vec_mode_err);
    for (size_t i = 0; i < speciesIndex.size(); ++i) {
        species[i]->deployVectorTrap(m_rng, i, instance, number, lifespan);
    }
}

void TransmissionModel::uninfectVectors_vector() {
    for (auto& s : species)
        s->uninfectVectors();
}

void TransmissionModel::summarize_vector(RunContext& ctx) {

    for (size_t i = 0; i < speciesIndex.size(); ++i) {
        species[i]->summarize(ctx, i);
    }
}

} // namespace OpenMalaria
