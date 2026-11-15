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

#include "schema/entomology.h"
#include "Transmission/NonVectorModel.h"

namespace OM
{
namespace Transmission
{

NonVectorModel::NonVectorModel(vector<double> initEIR, int interventionMode, const scnXml::NonVector &nonVectorData, int eipDuration)
    : TransmissionModel(std::move(initEIR), interventionMode, 1 /*this model doesn't support multiple genotypes*/)
    , nSpore(sim::fromDays(eipDuration))
{
    laggedKappa.resize(sim::inSteps(nSpore) + 1, 0.0);

    vector<int> nDays(sim::stepsPerYear(), 0);
    // The minimum EIR allowed in the array. The product of the average EIR and a constant.
    double minEIR = min_EIR_mult * averageEIR(nonVectorData);

    const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily();
    if (daily.size() < static_cast<size_t>(sim::oneYear()))
        throw util::xml_scenario_error("insufficient EIRDaily data for a year");

    for (SimTime mpcday = sim::zero(), endDay = sim::fromDays(daily.size()); mpcday < endDay; mpcday = mpcday + sim::oneDay())
    {
        double EIRdaily = std::max(static_cast<double>(daily[mpcday]), minEIR);

        // Index 0 of initialisationEIR refers to the EIR affecting the
        // first day(s) of the year. Correspondingly, the first 1 or 5 values
        // of EIRDaily affect this (1- or 5-day) time-step.
        size_t i = mod_nn(sim::inSteps(mpcday), sim::stepsPerYear());

        nDays[i] += 1;
        initialisationEIR[i] += EIRdaily;
    }

    // Calculate total annual EIR
    // divide by number of records assigned to each interval (usually one per day)
    for (size_t indTS = 0; indTS < sim::stepsPerYear(); indTS += 1)
    {
        initialisationEIR[indTS] *= sim::oneTS() / static_cast<double>(nDays[indTS]);
        annualEIR += initialisationEIR[indTS];
    }

    initialKappa.assign(sim::inSteps(sim::fromYearsI(nYearsWarmupData)), 0.0);
}

void NonVectorModel::changeEIRIntervention(const scnXml::NonVector &nonVectorData)
{
    // Note: requires sim::intervTime() >= sim::zero(), but this can only be
    // called in intervention period anyway.
    simulationMode = transientEIRknown;

    if (nSpore != sim::fromDays(nonVectorData.getEipDuration()))
        throw util::xml_scenario_error("change-of-EIR intervention cannot change EIP duration");

    const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily();
    vector<int> nDays(sim::inSteps(sim::fromDays(daily.size() - 1)) + 1, 0);
    interventionEIR.assign(nDays.size(), 0.0);
    size_t required_days = static_cast<size_t>((sim::endDate() - sim::startDate()) + 1);
    if (daily.size() < required_days)
    {
        cerr << "Days: " << daily.size() << "\nIntervals: " << nDays.size() << "\nRequired: " << required_days << endl;
        throw util::xml_scenario_error("Insufficient intervention phase EIR values provided");
    }
    // The minimum EIR allowed in the array. The product of the average EIR and a constant.
    double minEIR = min_EIR_mult * averageEIR(nonVectorData);
    for (SimTime mpcday = sim::zero(), endDay = sim::fromDays(daily.size()); mpcday < endDay; mpcday = mpcday + sim::oneDay())
    {
        double EIRdaily = std::max(static_cast<double>(daily[mpcday]), minEIR);

        // istep is the time period to which the day is assigned.
        size_t istep = sim::inSteps(mpcday);
        nDays[istep]++;
        interventionEIR[istep] += EIRdaily;
    }
    // divide by number of records assigned to each interval (usually one per day)
    for (size_t i = 0; i < interventionEIR.size(); ++i)
    {
        interventionEIR[i] *= sim::oneTS() / static_cast<double>(nDays[i]);
    }

    // I've no idea what this should be, so until someone asks it can be NaN.
    // It was -9.99 and later 0.0. It could of course be recalculated from interventionEIR.
    annualEIR = numeric_limits<double>::quiet_NaN();
}

} // namespace Transmission
} // namespace OM
