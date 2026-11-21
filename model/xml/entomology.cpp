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

// model/xml/ conversion layer — entomology: vector/non-vector transmission
// model construction and EIR/change-EIR intervention reading from scnXml.
// Moved verbatim from Transmission/transmission.h (the inline factory
// functions; readMode stays), TransmissionModel.h (averageEIR),
// NonVectorModel.h (constructor and changeEIRIntervention bodies),
// PerHost.cpp, VectorModel.cpp (per-species intervention loops became the
// free initVectorInterv/initVectorTrap functions, called via
// model/xml/entomology.h) and Anopheles/AnophelesModel.cpp.

#include "Clock.h"
#include "RunContext.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/transmission.h"
#include "Transmission/PerHost.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/Anopheles/EmergenceRateEstimator.h"
#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/SpeciesIndexChecker.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "util/sampler.h"
#include "util/vectors.h"
#include "schema/scenario.h"
#include "schema/interventions.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

using namespace std;

namespace OpenMalaria {

// Arithmetic mean of a whole daily-EIR sequence (both overloads read it from
// the .XML file the same way); noneMsg is the pedantic empty-vector check text.
template <class Seq> static double averageEIR_seq(const Seq& daily, const char* noneMsg) {
    double valaverageEIR = 0.0;
    size_t i = 0;
    for (; i < daily.size(); ++i)
        valaverageEIR += (double)daily[i];
    if (i == 0)
        throw xml_scenario_error(noneMsg);
    return valaverageEIR / i;
}

inline double averageEIR(const scnXml::NonVector& nonVectorData) {
    return averageEIR_seq(nonVectorData.getEIRDaily(), "no EIRDaily values given");
}

inline double averageEIR(const scnXml::DailyValues& vectorData) {
    return averageEIR_seq(vectorData.getValue(), "no daily values given");
}

TransmissionModel* createNonVectorModel(RunContext& ctx, const scnXml::Entomology& entoData) {
    const scnXml::NonVector& nonVectorData = entoData.getNonVector().get();

    int interventionMode = readMode(entoData.getMode());
    int eipDuration = nonVectorData.getEipDuration();

    vector<double> initialisationEIR;
    initialisationEIR.assign(ctx.clock.steps_per_year, 0.0);

    return new TransmissionModel(ctx, initialisationEIR, interventionMode, nonVectorData, eipDuration);
}

bool anophelesCompare(const scnXml::AnophelesParams& a1, const scnXml::AnophelesParams& a2) {
    return (a1.getSeasonality().getAnnualEIR().get() > a2.getSeasonality().getAnnualEIR().get());
}

AnophelesModel* createAnophelesModel(RunContext& ctx, size_t i, const scnXml::AnophelesParams& anoph,
                                     vector<double>& initialisationEIR, int populationSize, int interventionMode) {
    AnophelesModel* anophModel;

    if (ModelOptions::option(VECTOR_LIFE_CYCLE_MODEL)) {
        throw xml_scenario_error("VECTOR_LIFE_CYCLE_MODEL not yet "
                                 "implemented. Use VECTOR_SIMPLE_MPD_MODEL instead.");
        // TODO
        //  * Note: this model is older than SimpleMPD and more complicated.
        //  * Difficulties are in parameterisation and estimation of resources.
        // if (!lcOpt.present())
        //     throw xml_scenario_error(
        //         "VECTOR_LIFE_CYCLE_MODEL: requires <lifeCycle> element with "
        //         "model parameters for each anopheles species");
        // emergence = unique_ptr<EmergenceModel>( new LCEmergence() );
        // emergence->initLifeCycle( lcOpt.get() );
    } else if (ModelOptions::option(VECTOR_SIMPLE_MPD_MODEL)) {
        if (!anoph.getSimpleMPD().present())
            throw xml_scenario_error("VECTOR_SIMPLE_MPD_MODEL: requires <simpleMPD> element with "
                                     "model parameters for each anopheles species");

        const scnXml::SimpleMPD& smpd = anoph.getSimpleMPD().get();

        SimTime developmentDuration = fromDays(smpd.getDevelopmentDuration().getValue());
        if (!(developmentDuration > zero()))
            throw xml_scenario_error("entomology.vector.simpleMPD.developmentDuration: "
                                     "must be positive");
        double probPreadultSurvival = smpd.getDevelopmentSurvival().getValue();
        if (!(0.0 <= probPreadultSurvival && probPreadultSurvival <= 1.0))
            throw xml_scenario_error("entomology.vector.simpleMPD.developmentSurvival: "
                                     "must be a probability (in range [0,1]");
        double fEggsLaidByOviposit = smpd.getFemaleEggsLaidByOviposit().getValue();
        if (!(fEggsLaidByOviposit > 0.0))
            throw xml_scenario_error("entomology.vector.simpleMPD.femaleEggsLaidByOviposit: "
                                     "must be positive");

        anophModel = new AnophelesModel(ctx);
        anophModel->initSimpleMPD(developmentDuration, probPreadultSurvival, fEggsLaidByOviposit);
    } else
        anophModel = new AnophelesModel(ctx);

    const scnXml::Seasonality& seasonality = anoph.getSeasonality();

    if (seasonality.getInput() != "EIR")
        throw xml_scenario_error("entomology.anopheles.seasonality.input: must be EIR (for now)");

    vector<double> initEIR365;
    vector<double> FSCoeffic;
    double EIRRotateAngle = 0;
    double targetEIR = seasonality.getAnnualEIR().get();
    if (seasonality.getFourierSeries().present()) {
        const scnXml::FourierSeries& seasFC = seasonality.getFourierSeries().get();
        const scnXml::FourierSeries::CoefficSequence& fsCoeffic = seasFC.getCoeffic();

        FSCoeffic.reserve(2 * fsCoeffic.size() + 1);

        FSCoeffic.push_back(0.0); // value doesn't matter; EIR will be scaled
        for (auto it = fsCoeffic.begin(); it != fsCoeffic.end(); ++it) {
            FSCoeffic.push_back(it->getA());
            FSCoeffic.push_back(it->getB());
        }

        EIRRotateAngle = 0;
        if (seasFC.getEIRRotateAngle().present())
            EIRRotateAngle = seasFC.getEIRRotateAngle().get();

        // EIR for this species, with index 0 refering to value over first interval
        initEIR365.resize(oneYear(), 0.0);

        // Now we rescale to get an EIR of targetEIR.
        // Calculate current sum as is usually done.
        vectors::expIDFT(initEIR365, FSCoeffic, EIRRotateAngle);

        // And scale (also acts as a unit conversion):
        FSCoeffic[0] += log(targetEIR / vectors::sum(initEIR365));

        // Calculate forced EIR for pre-intervention phase from FSCoeffic:
        vectors::expIDFT(initEIR365, FSCoeffic, EIRRotateAngle);
    } else if (seasonality.getMonthlyValues().present()) {
        const scnXml::MonthlyValues& seasM = seasonality.getMonthlyValues().get();
        if (seasM.getSmoothing() == "fourier") {
            const size_t N_m = 12;
            const scnXml::MonthlyValues::ValueSequence seq = seasM.getValue();
            assert(seq.size() == N_m); // enforced by schema
            vector<double> months(N_m);
            double maxEIR = months[0];
            for (size_t i = 0; i < N_m; ++i) {
                months[i] = seq[i];
                if (months[i] > maxEIR)
                    maxEIR = months[i];
            }
            // arbitrary minimum we allow (cannot have zeros since we take the logarithm)
            double min = maxEIR / 100.0;
            for (size_t i = 0; i < N_m; ++i) {
                if (months[i] < min)
                    months[i] += min;
            }

            FSCoeffic.assign(5, 0.0);
            // TODO: determine whether to use Fourier Series Coefficient method or
            // Discrete Fourier Transform. Former is designed for integrals (which
            // roughly what we have?), the latter for discrete values. DFT doesn't
            // work well when number of intervals changes?
            vectors::logFourierCoefficients(months, FSCoeffic);

            // Assume the values are for the start of the month
            EIRRotateAngle = 0;
            if (seasM.getEIRRotateAngle().present())
                EIRRotateAngle = seasM.getEIRRotateAngle().get();

            // EIR for this species, with index 0 refering to value over first interval
            initEIR365.resize(oneYear());

            // Now we rescale to get an EIR of targetEIR.
            // Calculate current sum as is usually done.
            vectors::expIDFT(initEIR365, FSCoeffic, EIRRotateAngle);
            // And scale (also acts as a unit conversion):
            FSCoeffic[0] += log(targetEIR / vectors::sum(initEIR365));

            // Calculate forced EIR for pre-intervention phase from FSCoeffic:
            vectors::expIDFT(initEIR365, FSCoeffic, EIRRotateAngle);
        } else if (seasM.getSmoothing() == "none") {
            if (interventionMode != forcedEIR) {
                throw xml_scenario_error(
                    "entomology.anopheles.seasonality.monthlyValues.smoothing: smoothing mode \"none\" is not allowed with monthly EIR values");
            }
        } else {
            if (interventionMode != forcedEIR) {
                throw xml_scenario_error(
                    "entomology.anopheles.seasonality.monthlyValues.smoothing: unknown smoothing mode");
            }
        }
    } else {
        assert(seasonality.getDailyValues().present()); // XML loading code should enforce this

        if (interventionMode != forcedEIR) {
            throw xml_scenario_error(
                "entomology.anopheles.seasonality.dailyValues: daily values are only allowed with forced EIR");
        }

        const scnXml::DailyValues& seasM = seasonality.getDailyValues().get();
        const scnXml::DailyValues::ValueSequence& daily = seasM.getValue();

        // The minimum EIR allowed in the array. The product of the average EIR and a constant.
        double minEIR = 0.01 * averageEIR(seasM);

        if (daily.size() < static_cast<size_t>(oneYear()))
            throw xml_scenario_error("entomology.anopheles.seasonality.dailyValues insufficient daily data for a year");

        // EIR for this species, with index 0 refering to value over first interval
        initEIR365.resize(oneYear(), 0.0);
        vector<int> nDays(oneYear(), 0.0);

        for (SimTime d = zero(), endDay = daily.size(); d < endDay; d = d + oneDay()) {
            double EIRdaily = std::max(static_cast<double>(daily[d]), minEIR);

            // Index 0 of initialisationEIR refers to the EIR affecting the
            // first day(s) of the year. Correspondingly, the first 1 or 5 values
            // of EIRDaily affect this (1- or 5-day) time-step.
            size_t i = mod_nn(d, oneYear());

            nDays[i] += 1;
            initEIR365[i] += EIRdaily;
        }

        // Calculate total annual EIR
        // divide by number of records assigned to each interval (usually one per day)
        for (SimTime d = zero(); d < oneYear(); d = d + oneDay())
            initEIR365[d] *= 1.0 / static_cast<double>(nDays[d]);
    }

    if (!seasonality.getAnnualEIR().present())
        throw xml_scenario_error("entomology.anopheles.seasonality.annualEIR is required at the moment");

    // -----  Set model variables  -----
    const scnXml::Mosq& mosq = anoph.getMosq();

    MosquitoParams mosqParams;
    mosqParams.name = anoph.getMosquito();
    mosqParams.laidEggsSameDayProportion = mosq.getMosqLaidEggsSameDayProportion().getValue();
    mosqParams.survivalFeedingCycleProbability = mosq.getMosqSurvivalFeedingCycleProbability().getValue();
    mosqParams.humanBloodIndex = mosq.getMosqHumanBloodIndex().getValue();
    mosqParams.probBiting = mosq.getMosqProbBiting().getMean();
    mosqParams.probFindRestSite = mosq.getMosqProbFindRestSite().getMean();
    mosqParams.probResting = mosq.getMosqProbResting().getMean();
    mosqParams.probOvipositing = mosq.getMosqProbOvipositing().getValue();
    mosqParams.seekingDuration = mosq.getMosqSeekingDuration().getValue();
    mosqParams.probMosqSurvivalOvipositing = mosq.getMosqProbOvipositing().getValue();
    mosqParams.restDuration = fromDays(mosq.getMosqRestDuration().getValue());
    mosqParams.EIPDuration = fromDays(mosq.getExtrinsicIncubationPeriod().getValue());
    mosqParams.minInfectedThreshold = mosq.getMinInfectedThreshold();

    if (oneDay() > mosqParams.restDuration || mosqParams.restDuration * 2 >= mosqParams.EIPDuration) {
        // TODO: limit was EIPDuration >= mosqRestDuration >= 1
        // but in usage of ftauArray this wasn't enough. Check why.
        throw xml_scenario_error("Code expects EIPDuration > 2*mosqRestDuration >= 2");
    }

    const scnXml::AnophelesParams::NonHumanHostsSequence& xmlSeqNNHs = anoph.getNonHumanHosts();

    vector<NhhParams> nhhs;
    for (const scnXml::NonHumanHosts& xmlNNH : xmlSeqNNHs) {
        NhhParams nhh;
        nhh.mosqRelativeEntoAvailability = xmlNNH.getMosqRelativeEntoAvailability().getValue();
        nhh.mosqProbBiting = xmlNNH.getMosqProbBiting().getValue();
        nhh.mosqProbFindRestSite = xmlNNH.getMosqProbFindRestSite().getValue();
        nhh.mosqProbResting = xmlNNH.getMosqProbResting().getValue();
        nhh.hostFecundityFactor =
            xmlNNH.getHostFecundityFactor().present() ? xmlNNH.getHostFecundityFactor().get().getValue() : 1.0;
        nhh.name = xmlNNH.getName();
        nhhs.push_back(nhh);
    }

    // Initial estimate of the proportion of mosquitoes which are infectious, s: S_v(t) = s*N_v(t). Used as a starting value and then fit.
    double propInfectious = 0.021;
    // Initial guess of the proportion of mosquitoes which are infected, o: O_v(t) = o*N_v(t). Only used as a starting value.
    double propInfected = 0.078;

    anophModel->initialise(mosqParams);
    anophModel->initAvailability(i, nhhs, populationSize);
    anophModel->initEIR(initEIR365, FSCoeffic, EIRRotateAngle, propInfectious, propInfected);

    for (SimTime i = zero(); i < oneYear(); i = i + oneDay())
        initialisationEIR[mod_nn(ctx.clock.inSteps(i), ctx.clock.steps_per_year)] += initEIR365[i];

    return anophModel;
}

TransmissionModel* createVectorModel(RunContext& ctx, const scnXml::Entomology& entoData, int populationSize) {
    const scnXml::Vector& vectorData = entoData.getVector().get();
    scnXml::Vector::AnophelesSequence anophelesList = vectorData.getAnopheles();

    int interventionMode = readMode(entoData.getMode());

    vector<double> initialisationEIR;
    initialisationEIR.assign(ctx.clock.steps_per_year, 0.0);

    vector<std::unique_ptr<AnophelesModel>> species;
    vector<std::unique_ptr<EmergenceRateEstimator>> emergenceRateEstimators;
    map<string, size_t> speciesIndex;

    size_t numSpecies = anophelesList.size();

    if (numSpecies < 1)
        throw xml_scenario_error("Can't use Vector model without data for at least one anopheles species!");

    // Sort Anopheles by Decreasing EIR
    sort(anophelesList.begin(), anophelesList.end(), anophelesCompare);

    for (size_t i = 0; i < numSpecies; ++i) {
        auto anoph = anophelesList[i];

        PerHostAnophParams::init(ctx, anoph.getMosq());

        AnophelesModel* anophModel =
            createAnophelesModel(ctx, i, anoph, initialisationEIR, populationSize, interventionMode);
        std::unique_ptr<EmergenceRateEstimator> estimator = std::make_unique<EmergenceRateEstimator>(
            ModelOptions::option(USE_EXACT_NV0_SOLVER) ? EmergenceRateEstimator::Kind::Solver
                                                       : EmergenceRateEstimator::Kind::AdaptiveFitter,
            *anophModel, populationSize);

        species.push_back(std::unique_ptr<AnophelesModel>(anophModel));
        emergenceRateEstimators.push_back(std::move(estimator));
        speciesIndex[anophModel->mosq.name] = i;
    }

    // Note: in forcedEIR mode species/emergenceRateEstimators are not needed
    // anymore (we have initialisationEIR), but interventions still use them.

    return new TransmissionModel(ctx, initialisationEIR, interventionMode, std::move(species),
                                 std::move(emergenceRateEstimators), speciesIndex);
}

TransmissionModel* createTransmissionModel(RunContext& ctx, const scnXml::Entomology& entoData, int populationSize) {
    // Entomology contains either a list of at least one anopheles or a list of at
    // least one EIRDaily.
    TransmissionModel* model;
    if (entoData.getVector().present())
        model = createVectorModel(ctx, entoData, populationSize);
    else {
        if (!entoData.getNonVector().present()) // should be a validation error, but anyway...
            throw xml_scenario_error("Neither vector nor non-vector data present in the XML!");
        if (ModelOptions::option(VECTOR_LIFE_CYCLE_MODEL) || ModelOptions::option(VECTOR_SIMPLE_MPD_MODEL))
            throw xml_scenario_error(
                "VECTOR_*_MODEL is only compatible with the vector model (and non-vector data is present).");
        model = createNonVectorModel(ctx, entoData);
    }

    if (entoData.getScaledAnnualEIR().present()) {
        model->scaleEIR(entoData.getScaledAnnualEIR().get() / model->annualEIR);
        assert(vectors::approxEqual(model->annualEIR, entoData.getScaledAnnualEIR().get()));
    }

    if (CommandLine::option(CommandLine::PRINT_ANNUAL_EIR)) {
        // Note: after internal scaling (which doesn't imply exit)
        // but before external scaling.
        cout << "Total annual EIR: " << model->annualEIR << endl;
    }

    return model;
}

TransmissionModel::TransmissionModel(RunContext& ctx, vector<double> initEIR, int interventionMode,
                                     const scnXml::NonVector& nonVectorData, int eipDuration)
    : TransmissionModel(ctx, std::move(initEIR), interventionMode, Kind::NonVector, {}, {}, {}) {
    nSpore = fromDays(eipDuration);
    laggedKappa.resize(ctx.clock.inSteps(nSpore) + 1, 0.0);

    vector<int> nDays(ctx.clock.steps_per_year, 0);
    // The minimum EIR allowed in the array. The product of the average EIR and a constant.
    double minEIR = min_EIR_mult * averageEIR(nonVectorData);

    const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
    if (daily.size() < static_cast<size_t>(oneYear()))
        throw xml_scenario_error("insufficient EIRDaily data for a year");

    for (SimTime mpcday = zero(), endDay = fromDays(daily.size()); mpcday < endDay; mpcday = mpcday + oneDay()) {
        double EIRdaily = std::max(static_cast<double>(daily[mpcday]), minEIR);

        // Index 0 of initialisationEIR refers to the EIR affecting the
        // first day(s) of the year. Correspondingly, the first 1 or 5 values
        // of EIRDaily affect this (1- or 5-day) time-step.
        size_t i = mod_nn(ctx.clock.inSteps(mpcday), ctx.clock.steps_per_year);

        nDays[i] += 1;
        initialisationEIR[i] += EIRdaily;
    }

    // Calculate total annual EIR
    // divide by number of records assigned to each interval (usually one per day)
    for (size_t indTS = 0; indTS < ctx.clock.steps_per_year; indTS += 1) {
        initialisationEIR[indTS] *= ctx.clock.one_TS() / static_cast<double>(nDays[indTS]);
        annualEIR += initialisationEIR[indTS];
    }

    initialKappa.assign(ctx.clock.inSteps(fromYearsI(5)), 0.0); // 5y warmup data (was nYearsWarmupData)
}

// Shared loop tail of both changeEIRIntervention variants: fills
// interventionEIR with the per-interval mean of max(daily, minEIR). The
// days_cerr flag reproduces the extra diagnostic of the non-vector variant.
static void fill_intervention_eir(RunContext& ctx, vector<double>& interventionEIR,
                                  const scnXml::NonVector::EIRDailySequence& daily, double minEIR, bool days_cerr) {
    vector<int> nDays(ctx.clock.inSteps(fromDays(daily.size() - 1)) + 1, 0);
    interventionEIR.assign(nDays.size(), 0.0);
    size_t required_days = static_cast<size_t>((ctx.clock.end_date - ctx.clock.start_date) + 1);
    if (daily.size() < required_days) {
        if (days_cerr)
            cerr << "Days: " << daily.size() << "\nIntervals: " << nDays.size() << "\nRequired: " << required_days
                 << endl;
        throw xml_scenario_error("Insufficient intervention phase EIR values provided");
    }
    // The minimum EIR allowed in the array. The product of the average EIR and a constant.
    for (SimTime mpcday = zero(), endDay = fromDays(daily.size()); mpcday < endDay; mpcday = mpcday + oneDay()) {
        double EIRdaily = std::max(static_cast<double>(daily[mpcday]), minEIR);

        // istep is the time period to which the day is assigned.
        size_t istep = ctx.clock.inSteps(mpcday);
        nDays[istep]++;
        interventionEIR[istep] += EIRdaily;
    }
    // divide by number of records assigned to each interval (usually one per day)
    for (size_t i = 0; i < interventionEIR.size(); ++i) {
        interventionEIR[i] *= ctx.clock.one_TS() / static_cast<double>(nDays[i]);
    }
}

void TransmissionModel::changeEIRIntervention_nonVector(const scnXml::NonVector& nonVectorData) {
    // Note: requires ctx.clock.interv >= zero(), but this can only be
    // called in intervention period anyway.
    simulationMode = transientEIRknown;

    if (nSpore != fromDays(nonVectorData.getEipDuration()))
        throw xml_scenario_error("change-of-EIR intervention cannot change EIP duration");

    // EIRDaily is schema-required (minOccurs=1), so computing the average
    // before the insufficient-days check cannot change which error is seen.
    fill_intervention_eir(ctx, interventionEIR, nonVectorData.getEIRDaily(), min_EIR_mult * averageEIR(nonVectorData),
                          true);

    // I've no idea what this should be, so until someone asks it can be NaN.
    // It was -9.99 and later 0.0. It could of course be recalculated from interventionEIR.
    annualEIR = numeric_limits<double>::quiet_NaN();
}

PerHostAnophParams::PerHostAnophParams(const scnXml::Mosq& mosq) {
    if (mosq.getAvailability().getDistr() != "const" &&
        (ModelOptions::option(NEGATIVE_BINOMIAL_MASS_ACTION) || ModelOptions::option(LOGNORMAL_MASS_ACTION) ||
         ModelOptions::option(COMORB_TRANS_HET)))
        throw xml_scenario_error(
            "PerHostAnophParams::PerHostAnophParams(): ModelOptions NEGATIVE_BINOMIAL_MASS_ACTION, LOGNORMAL_MASS_ACTION, COMORB_TRANS_HET are not compatible with ento-availability heterogeneity distributions other than \"const\"");

    entoAvailability = createSampler(1.0, mosq.getAvailability());
    probMosqBiting.setParams(mosq.getMosqProbBiting());
    probMosqFindRestSite.setParams(mosq.getMosqProbFindRestSite());
    probMosqSurvivalResting.setParams(mosq.getMosqProbResting());
}

void PerHost::init(RunContext& ctx, const scnXml::AgeGroupValues& availabilityToMosquitoes) {
    relAvailAge.set(ctx.clock, availabilityToMosquitoes, "availabilityToMosquitoes");
}

void initVectorInterv(TransmissionModel& vm, const scnXml::Description::AnophelesSequence& list, size_t instance,
                      const string& name) {
    SpeciesIndexChecker checker(name, vm.speciesIndex);
    for (const scnXml::VectorSpeciesIntervention& anoph : list) {
        const string& mosq = anoph.getMosquito();
        vm.species[checker.getIndex(mosq)]->initVectorInterv(anoph, instance);
    }
    checker.checkNoneMissed();
}
void initVectorTrap(TransmissionModel& vm, const scnXml::VectorTrap::DescriptionSequence list, size_t instance,
                    const scnXml::VectorTrap::NameOptional name_opt) {
    stringstream name_ss;
    if (name_opt.present())
        name_ss << name_opt.get();
    else
        name_ss << "vector trap intervention " << (instance + 1);
    string name = name_ss.str();
    SpeciesIndexChecker checker(name, vm.speciesIndex);
    for (const scnXml::Description1& anoph : list) {
        const string& mosq = anoph.getMosquito();
        vm.species[checker.getIndex(mosq)]->initVectorTrap(anoph, instance);
    }
    checker.checkNoneMissed();
}

void TransmissionModel::changeEIRIntervention_vector(const scnXml::NonVector& nonVectorData) {
    if ((simulationMode != forcedEIR) && (simulationMode != transientEIRknown))
        throw xml_scenario_error("changeEIR intervention can only be used in forced EIR mode");

    simulationMode = transientEIRknown;

    if (speciesIndex.size() != 1)
        throw xml_scenario_error("changeEIR intervention can only be used with a single anopheles species");

    species[0]->changeEIRIntervention(nonVectorData);
}

// Shared range-check + set of the four initVectorInterv decay slots (the
// resizes stay at the call sites, before the present() checks, as before;
// the bound pair reproduces each original check exactly -- NaN initial
// fails none of them either way).
template <class V, class E2>
static void set_interv_decay(SimClock& clock, V& slot, size_t instance, const E2& elt2, double lo, double hi,
                             const char* rangeMsg, const char* tag) {
    double initial = elt2.getInitial();
    if (initial < lo || initial > hi)
        throw xml_scenario_error(rangeMsg);
    slot[instance].set(clock, initial, elt2.getDecay(), tag);
}

void AnophelesModel::initVectorInterv(const scnXml::VectorSpeciesIntervention& elt, size_t instance) {
    if (emergenceReduction.size() <= instance)
        emergenceReduction.resize(instance + 1);
    if (elt.getEmergenceReduction().present())
        set_interv_decay(ctx.clock, emergenceReduction, instance, elt.getEmergenceReduction().get(),
                         -numeric_limits<double>::infinity(), 1.0,
                         "emergenceReduction intervention: initial effect must be ≤ 1", "emergenceReduction");
    if (probAdditionalDeathSugarFeedingIntervs.size() <= instance)
        probAdditionalDeathSugarFeedingIntervs.resize(instance + 1);
    if (seekingDeathRateIntervs.size() <= instance)
        seekingDeathRateIntervs.resize(instance + 1);
    if (probDeathOvipositingIntervs.size() <= instance)
        probDeathOvipositingIntervs.resize(instance + 1);

    if (elt.getProbAdditionalDeathSugarFeeding().present())
        set_interv_decay(ctx.clock, probAdditionalDeathSugarFeedingIntervs, instance,
                         elt.getProbAdditionalDeathSugarFeeding().get(), -1.0, numeric_limits<double>::infinity(),
                         "probAdditionalDeathSugarFeeding intervention: initial effect must be ≥ -1",
                         "seekingDeathRateIncrease");
    if (elt.getSeekingDeathRateIncrease().present())
        set_interv_decay(ctx.clock, seekingDeathRateIntervs, instance, elt.getSeekingDeathRateIncrease().get(), -1.0,
                         numeric_limits<double>::infinity(),
                         "seekingDeathRateIncrease intervention: initial effect must be ≥ -1",
                         "seekingDeathRateIncrease");
    if (elt.getProbDeathOvipositing().present())
        set_interv_decay(ctx.clock, probDeathOvipositingIntervs, instance, elt.getProbDeathOvipositing().get(), 0.0,
                         1.0, "probDeathOvipositing intrevention: initial effect must be in range [0,1]",
                         "probDeathOvipositing");
}
void AnophelesModel::initVectorTrap(const scnXml::Description1& desc, size_t instance) {
    assert(trapParams.size() == instance); // if triggered, this is a code error not XML
    TrapParams params;
    params.relAvail = desc.getRelativeAvailability().getValue();
    params.availDecay = DecayFunction::makeObject(ctx.clock, desc.getDecayOfAvailability(), "decayOfAvailability");
    trapParams.push_back(std::move(params));
}

void AnophelesModel::changeEIRIntervention(const scnXml::NonVector& nonVectorData) {
    // assert(nonVectorData.getEIRDaily().present()); // XML loading code should enforce this
    const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();

    // The minimum EIR allowed in the array. The product of the average EIR and a constant.
    double minEIR = 0.01 * averageEIR(nonVectorData);

    if (daily.size() < static_cast<size_t>(oneYear()))
        throw xml_scenario_error("entomology.anopheles.seasonality.dailyValues insufficient daily data for a year");

    fill_intervention_eir(ctx, interventionEIR, daily, minEIR, false);
}

}
