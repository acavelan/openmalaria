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
#include "mon/Continuous.h"
#include "util/vectors.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"

#include <fstream>
#include <map>
#include <cmath>
#include <set>

namespace OM
{
namespace Transmission
{
using namespace OM::util;
using Anopheles::AnophelesModel;

void TransmissionModel::ctsCbN_v0(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastN_v0();
}
void TransmissionModel::ctsCbP_A(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PA);
}
void TransmissionModel::ctsCbP_Amu(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PAmu);
}
void TransmissionModel::ctsCbP_A1(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PA1);
}
void TransmissionModel::ctsCbP_Ah(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PAh);
}
void TransmissionModel::ctsCbP_df(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PDF);
}
void TransmissionModel::ctsCbP_dif(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PDIF);
}
void TransmissionModel::ctsCbN_v(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::NV);
}
void TransmissionModel::ctsCbO_v(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::OV);
}
void TransmissionModel::ctsCbS_v(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::SV);
}
void TransmissionModel::ctsCbAlpha(const vector<Host::Human> &population, ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
    {
        double total = 0.0;
        for (const Host::Human &human : population)
        {
            total += human.perHostTransmission.entoAvailabilityFull(i, sim::inYears(human.age(ctx.clock.now())));
        }
        stream << '\t' << total / population.size();
    }
}
void TransmissionModel::ctsCbP_B(const vector<Host::Human> &population, ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
    {
        double total = 0.0;
        for (const Host::Human &human : population)
        {
            total += human.perHostTransmission.probMosqBiting(i);
        }
        stream << '\t' << total / population.size();
    }
}
void TransmissionModel::ctsCbP_CD(const vector<Host::Human> &population, ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
    {
        double total = 0.0;
        for (const Host::Human &human : population)
        {
            total += human.perHostTransmission.probMosqResting(i);
        }
        stream << '\t' << total / population.size();
    }
}
void TransmissionModel::ctsCbResAvailability(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getResAvailability();
}
void TransmissionModel::ctsCbResRequirements(ostream &stream)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getResRequirements();
}

const string &reverseLookup(const map<string, size_t> &m, size_t i)
{
    for (auto it = m.begin(); it != m.end(); ++it)
    {
        if (it->second == i) return it->first;
    }
    throw TRACED_EXCEPTION_DEFAULT("reverseLookup: key not found"); // shouldn't ever happen
}

TransmissionModel::TransmissionModel(RunContext& ctx, vector<double> initEIR, int interventionMode, vector<std::unique_ptr<Anopheles::AnophelesModel>> speciesList,
                         vector<std::unique_ptr<Anopheles::EmergenceRateEstimator>> emergenceRateEstimatorsList, map<string, size_t> speciesIndexList)
    : TransmissionModel(ctx, std::move(initEIR), interventionMode, Kind::Vector,
        std::move(speciesList), std::move(emergenceRateEstimatorsList), std::move(speciesIndexList))
{
    // set VACCINE_GENOTYPE option
    opt_vaccine_genotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    // Calculate total annual EIR
    annualEIR = vectors::sum(initialisationEIR);

    // -----  Continuous reporting  -----
    ostringstream ctsNv0, ctsPA, ctsPAmu, ctsPA1, ctsPAh, ctsPdf, ctsPdif, ctsNv, ctsOv, ctsSv, ctsAlpha, ctsPB, ctsPCD, ctsIRSEffects,
        ctsRA, ctsRR;

    size_t numSpecies = species.size();

    // Output in order of species so that (1) we can just iterate through this
    // list when outputting and (2) output is in order specified in XML.
    for (size_t i = 0; i < numSpecies; ++i)
    {
        // Unfortunately, we need to reverse-lookup the name.
        const string &name = reverseLookup(speciesIndex, i);
        ctsNv0 << "\tN_v0(" << name << ")";
        ctsPA << "\tP_A(" << name << ")";
        ctsPAmu << "\tP_Amu(" << name << ")";
        ctsPA1 << "\tP_A1(" << name << ")";
        ctsPAh << "\tP_Ah(" << name << ")";
        ctsPdf << "\tP_df(" << name << ")";
        ctsPdif << "\tP_dif(" << name << ")";
        ctsNv << "\tN_v(" << name << ")";
        ctsOv << "\tO_v(" << name << ")";
        ctsSv << "\tS_v(" << name << ")";
        ctsAlpha << "\talpha_i(" << name << ")";
        ctsPB << "\tP_B(" << name << ")";
        ctsPCD << "\tP_C*P_D(" << name << ")";
        ctsIRSEffects << "\tIRS rel attr (" << name << ")"
                      << "\tIRS preprand surv factor (" << name << ")"
                      << "\tIRS postprand surv factor (" << name << ")";
        ctsRA << "\tres avail(" << name << ")";
        ctsRR << "\tres req(" << name << ")";
    }

    mon::Continuous::registerCallback(ctx, "N_v0", ctsNv0.str(), std::bind( &TransmissionModel::ctsCbN_v0, this, _1));
    mon::Continuous::registerCallback(ctx, "P_A", ctsPA.str(), std::bind( &TransmissionModel::ctsCbP_A, this, _1));
    mon::Continuous::registerCallback(ctx, "P_Amu", ctsPAmu.str(), std::bind( &TransmissionModel::ctsCbP_Amu, this, _1));
    mon::Continuous::registerCallback(ctx, "P_A1", ctsPA1.str(), std::bind( &TransmissionModel::ctsCbP_A1, this, _1));
    mon::Continuous::registerCallback(ctx, "P_Ah", ctsPAh.str(), std::bind( &TransmissionModel::ctsCbP_Ah, this, _1));
    mon::Continuous::registerCallback(ctx, "P_df", ctsPdf.str(), std::bind( &TransmissionModel::ctsCbP_df, this, _1));
    mon::Continuous::registerCallback(ctx, "P_dif", ctsPdif.str(), std::bind( &TransmissionModel::ctsCbP_dif, this, _1));
    mon::Continuous::registerCallback(ctx, "N_v", ctsNv.str(), std::bind( &TransmissionModel::ctsCbN_v, this, _1));
    mon::Continuous::registerCallback(ctx, "O_v", ctsOv.str(), std::bind( &TransmissionModel::ctsCbO_v, this, _1));
    mon::Continuous::registerCallback(ctx, "S_v", ctsSv.str(), std::bind( &TransmissionModel::ctsCbS_v, this, _1));

    // availability to mosquitoes relative to other humans, excluding age factor
    mon::Continuous::registerCallback(ctx, "alpha", ctsAlpha.str(), std::bind( &TransmissionModel::ctsCbAlpha, this, _1, _2));
    mon::Continuous::registerCallback(ctx, "P_B", ctsPB.str(), std::bind( &TransmissionModel::ctsCbP_B, this, _1, _2));
    mon::Continuous::registerCallback(ctx, "P_C*P_D", ctsPCD.str(), std::bind( &TransmissionModel::ctsCbP_CD, this, _1, _2));

    mon::Continuous::registerCallback(ctx, "resource availability", ctsRA.str(), std::bind( &TransmissionModel::ctsCbResAvailability, this, _1));
    mon::Continuous::registerCallback(ctx, "resource requirements", ctsRR.str(), std::bind( &TransmissionModel::ctsCbResRequirements, this, _1));
}

void TransmissionModel::init2_vector(const vector<Host::Human> &population)
{
    double sumRelativeAvailability = 0.0;
    for (const Host::Human &human : population)
    {
        sumRelativeAvailability += human.perHostTransmission.relativeAvailabilityAge(sim::inYears(human.age(ctx.clock.now())));
    }
    int popSize = population.size();
    // value should be unimportant when no humans are available, though inf/nan is not acceptable
    if (popSize > 0)
    {
        meanAvail = sumRelativeAvailability / popSize; // mean-rel-avail
    }

    for (size_t i = 0; i < speciesIndex.size(); ++i)
    {
        double sum_avail = 0.0;
        double sigma_f = 0.0;
        double sigma_df = 0.0;
        double sigma_dff = 0.0;

        for (const Host::Human &human : population)
        {
            const OM::Transmission::PerHost &host = human.perHostTransmission;
            double prod = host.entoAvailabilityFull(i, sim::inYears(human.age(ctx.clock.now())));
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

void TransmissionModel::scaleEIR_vector(double factor)
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        species[i]->scaleEIR(factor);
    vectors::scale(initialisationEIR, factor);
    annualEIR = vectors::sum(initialisationEIR);
}

SimTime TransmissionModel::initIterate_vector()
{
    if (interventionMode != dynamicEIR)
    {
        // allow forcing equilibrium mode like with non-vector model
        return sim::zero(); // no initialization to do
    }

    // Fitting is done
    if (initIterations < 0)
    {
        simulationMode = dynamicEIR;
        return sim::zero();
    }

    if (++initIterations > 30) { throw TRACED_EXCEPTION("Transmission warmup exceeded 30 iterations!", util::Error::VectorWarmup); }

    bool needIterate = false;
    for (size_t i = 0; i < speciesIndex.size(); ++i)
    {
        needIterate = emergenceRateEstimators[i]->estimate(*species[i], laggedKappa, meanAvail);
        species[i]->initIterate();
        if (needIterate) break;
    }

    if (needIterate)
    {
        // stabilization + 5 years data-collection time:
        return sim::oneYear() + sim::fromYearsI(5);
    }
    else
    {
        // One year stabilisation, then we're finished:
        initIterations = -1;
        return sim::oneYear();
    }
}

void TransmissionModel::calculateEIR_vector(Host::Human &human, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l) const
{

    PerHost &host = human.perHostTransmission;
    host.update(human);
    if(simulationMode == transientEIRknown)
    {
        EIR_i.assign(1, 0.0);
        EIR_l.assign(1, 0.0);

        assert(simulationMode == forcedEIR);
        for (size_t i = 0; i < speciesIndex.size(); ++i)
        {
            /* Calculates EIR per individual (hence N_i == 1).
             *
             * See comment in AnophelesModel::advancePeriod for method. */
            double eir = species[i]->getInterventionEIR() * host.relativeAvailabilityAge(ageYears) * host.entoAvailabilityHetVecItv(i);
            mon::recordStat(mon::measure("innoculationsPerAgeGroup"), human, eir);
            mon::recordStat(mon::measure("innoculationsPerVector"), human, eir, i);
            // cout << host.availBite(i) << endl;
            // cout << "species: " << eir << " " << host.entoAvailabilityHetVecItv(i) << " " << host.entoAvailabilityHetVecItv(i) * (1000*0.00112608) << endl;
            EIR_l[0] += eir;
        }
    }
    else if (simulationMode == forcedEIR)
    {
        // double eir = initialisationEIR[ctx.clock.moduloYearSteps(ctx.clock.ts0())] * host.relativeAvailabilityHetAge(ageYears);

        // cout << "global: " << eir << " " <<  endl;
        // mon::recordStat(mon::measure("innoculationsPerAgeGroup"), human, eir);
        // EIR.assign(1, eir);

        EIR_i.assign(1, 0.0);
        EIR_l.assign(1, 0.0);

        assert(simulationMode == forcedEIR);
        for (size_t i = 0; i < speciesIndex.size(); ++i)
        {
            /* Calculates EIR per individual (hence N_i == 1).
             *
             * See comment in AnophelesModel::advancePeriod for method. */
            double eir = species[i]->getInitPartialEIR() * host.relativeAvailabilityAge(ageYears) * host.entoAvailabilityHetVecItv(i);
            mon::recordStat(mon::measure("innoculationsPerAgeGroup"), human, eir);
            mon::recordStat(mon::measure("innoculationsPerVector"), human, eir, i);
            // cout << host.availBite(i) << endl;
            // cout << "species: " << eir << " " << host.entoAvailabilityHetVecItv(i) << " " << host.entoAvailabilityHetVecItv(i) * (1000*0.00112608) << endl;
            EIR_l[0] += eir;
        }
    }
    else
    {
        assert(simulationMode == dynamicEIR);
        EIR_i.assign(ctx.genotypes.N(), 0.0);
        EIR_l.assign(ctx.genotypes.N(), 0.0);
        const double ageFactor = host.relativeAvailabilityAge(ageYears);
        for (size_t i = 0; i < speciesIndex.size(); ++i)
        {
            const vector<double> &partialEIR_i = species[i]->partialEIR_i;
            const vector<double> &partialEIR_l = species[i]->partialEIR_l;

            assert(EIR_i.size() == partialEIR_i.size());
            assert(EIR_l.size() == partialEIR_l.size());
            
            // if ((std::isnan)(vectors::sum(partialEIR_i))) { cerr << "partialEIR is not a number; " << i << endl; }

            /* Calculates EIR per individual (hence N_i == 1).
             *
             * See comment in AnophelesModel::advancePeriod for method. */
            double entoFactor = ageFactor * host.availBite(i);
            for (size_t g = 0; g < EIR_i.size(); ++g)
            {
                double eir_i = partialEIR_i[g] * entoFactor;
                double eir_l = partialEIR_l[g] * entoFactor;
                mon::recordStat(mon::measure("innoculationsPerAgeGroup"), human, eir_i + eir_l, 0, g);
                mon::recordStat(mon::measure("innoculationsPerVector"), human, eir_i + eir_l, i, g);
                EIR_i[g] += eir_i;
                EIR_l[g] += eir_l;
            }
        }
    }
}

// Every Global::interval days:
// 7c-C6 (kept per user ruling: uniform stage+reduce shape across the transmission
// sweeps — same two-pass pattern as calculateEIR_vector/updateKappa — even though
// the win, −0.56 %, sits in the ±1 % unknown band). Pass 1 writes pure elementwise
// per-host quantities into reused slabs; pass 2 reduces in host order, SCALAR
// (operand order = old host-major accumulation ⇒ bit-exact; SIMD reductions are
// phase 8). Slabs are function-local statics, capacity reused never shrunk
// (precedent: B1a m_y_lag).
void TransmissionModel::vectorUpdate_vector(const vector<Host::Human> &population)
{
    const size_t nGenotypes = ctx.genotypes.N();
    const size_t nSpec = speciesIndex.size();
    const size_t nHost = population.size();
    static std::vector<double> s_avail, s_df, s_fec, s_pti, s_ptl, s_tbv;
    s_avail.resize(nHost * nSpec);
    s_df.resize(nHost * nSpec);
    s_fec.resize(nHost * nSpec);
    s_pti.resize(nHost * nGenotypes);
    s_ptl.resize(nHost * nGenotypes);
    s_tbv.resize(nHost * nGenotypes);
    std::vector<double> probTransmission_i, probTransmission_l;
    std::vector<double> sum_avail(nSpec);
    std::vector<double> sigma_df(nSpec);
    std::vector<double> sigma_dff(nSpec);
    std::vector<std::vector<double>> sigma_dif_i(nSpec, std::vector<double>(nGenotypes)), sigma_dif_l(nSpec, std::vector<double>(nGenotypes));

    size_t h = 0;
    for (const Host::Human &human : population)
    {
        const OM::Transmission::PerHost &host = human.perHostTransmission;
        WithinHost::WHInterface &whm = *human.withinHostModel;

        // p7b B4 hoist: vaccine_factor is pure over the host's vac column (frozen
        // while the population is scanned) and returns 1.0 when it is empty — the
        // per-species × per-genotype calls collapse to one check per host.
        const bool anyVac = human.has_vaccines();
        probTransmission_i.assign(nGenotypes, 0.0);
        probTransmission_l.assign(nGenotypes, 0.0);
        whm.probTransmissionToMosquito(probTransmission_i, probTransmission_l);
        std::copy(probTransmission_i.begin(), probTransmission_i.end(), &s_pti[h * nGenotypes]);
        std::copy(probTransmission_l.begin(), probTransmission_l.end(), &s_ptl[h * nGenotypes]);
        for (size_t g = 0; g < nGenotypes; ++g)
            s_tbv[h * nGenotypes + g] = anyVac ? human.vaccine_factor(interventions::Vaccine::TBV, opt_vaccine_genotype? g : 0) : 1.0;

        for (size_t s = 0; s < nSpec; ++s)
        {
            // NOTE: calculate availability relative to age at end of time step;
            // not my preference but consistent with TransmissionModel::getEIR().
            // TODO: even stranger since probTransmission comes from the previous time step
            const double avail = host.entoAvailabilityFull(s, sim::inYears(human.age(ctx.clock.ts1())));
            const double df = avail * host.probMosqBiting(s) * host.probMosqResting(s);
            s_avail[h * nSpec + s] = avail;
            s_df[h * nSpec + s] = df;
            s_fec[h * nSpec + s] = df * host.relMosqFecundity(s);
        }
        ++h;
    }

    // pass 2: ordered reductions (same addend sequence as the old host-major loop)
    for (size_t s = 0; s < nSpec; ++s)
    {
        for (h = 0; h < nHost; ++h)
        {
            sum_avail[s] += s_avail[h * nSpec + s];
            sigma_df[s] += s_df[h * nSpec + s];
            sigma_dff[s] += s_fec[h * nSpec + s];
        }
        for (size_t g = 0; g < nGenotypes; ++g)
        {
            for (h = 0; h < nHost; ++h)
            {
                const double tbv = s_tbv[h * nGenotypes + g];
                sigma_dif_i[s][g] += s_df[h * nSpec + s] * s_pti[h * nGenotypes + g] * tbv;
                sigma_dif_l[s][g] += s_df[h * nSpec + s] * s_ptl[h * nGenotypes + g] * tbv;
            }
        }
    }

    for (size_t s = 0; s < nSpec; ++s)
    {
        species[s]->advancePeriod(sum_avail[s], sigma_df[s], sigma_dif_i[s], sigma_dif_l[s], sigma_dff[s], simulationMode == dynamicEIR);
    }
}

const string vec_mode_err = "vector interventions can only be used in "
                            "dynamic transmission mode (mode=\"dynamic\")";

void TransmissionModel::deployVectorPopInterv(size_t instance)
{
    if (interventionMode != dynamicEIR) { throw xml_scenario_error(vec_mode_err); }
    for (auto it = species.begin(); it != species.end(); ++it)
    {
        (*it)->deployVectorPopInterv(m_rng, instance);
    }
}
void TransmissionModel::deployVectorTrap(size_t instance, double number, SimTime lifespan)
{
    if (interventionMode != dynamicEIR) { throw xml_scenario_error(vec_mode_err); }
    for (size_t i = 0; i < speciesIndex.size(); ++i)
    {
        species[i]->deployVectorTrap(m_rng, i, instance, number, lifespan);
    }
}


void TransmissionModel::uninfectVectors_vector()
{
    for (size_t i = 0; i < speciesIndex.size(); ++i)
        species[i]->uninfectVectors();
}

void TransmissionModel::summarize_vector(RunContext& ctx)
{

    for (size_t i = 0; i < speciesIndex.size(); ++i)
    {
        species[i]->summarize(ctx, i);
    }
}

} // namespace Transmission
} // namespace OM
