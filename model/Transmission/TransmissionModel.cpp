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

#include "Transmission/TransmissionModel.h"
#include "schema/interventions.h"
#include "schema/util.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Genotypes.h"

namespace OM
{
namespace Transmission
{

double averageEIR(const scnXml::NonVector &nonVectorData)
{
    double valaverageEIR = 0.0;
    size_t i = 0;
    for (const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily(); i < daily.size(); ++i)
        valaverageEIR += (double)daily[i];
    if (i == 0) throw util::xml_scenario_error("no EIRDaily values given"); // pedantic check
    return valaverageEIR / i;
}

double averageEIR(const scnXml::DailyValues &vectorData)
{
    double valaverageEIR = 0.0;
    size_t i = 0;
    for (const scnXml::DailyValues::ValueSequence &daily = vectorData.getValue(); i < daily.size(); ++i)
        valaverageEIR += (double)daily[i];
    if (i == 0) throw util::xml_scenario_error("no daily values given"); // pedantic check
    return valaverageEIR / i;
}

void TransmissionModel::setAgeAvailability(const scnXml::AgeGroupValues& availabilityToMosquitoes)
{
    ageAvailability.set(availabilityToMosquitoes, "availabilityToMosquitoes");
    adultAge = sim::fromYearsD(ageAvailability.firstGlobalMaximum());
}
double TransmissionModel::getEIR(Host::Human &human, SimTime age, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l,
                                  const SimClock& clock)
{
    startHumanUpdate(human, clock);
    /* For the NonVector model, the EIR should just be multiplied by the
     * availability. For the Vector model, the availability is also required
     * for internal calculations, but again the EIR should be multiplied by the
     * availability. */
    calculateEIR(human, ageYears, EIR_i, EIR_l, clock);

    double sum_EIR_i = util::vectors::sum(EIR_i);
    double sum_EIR_l = util::vectors::sum(EIR_l);

    double allEIR = sum_EIR_i + sum_EIR_l;
    if (age >= adultAge)
    {
        tsAdultEntoInocs_i += sum_EIR_i;
        tsAdultEntoInocs_l += sum_EIR_l;
        tsAdultEntoInocs += allEIR;
        tsNumAdults += 1;
    }
    return allEIR;
}

double TransmissionModel::updateKappa(const vector<Host::Human> &population, const SimClock& clock)
{
    // We calculate kappa for output and the non-vector model.
    double sumWt_kappa = 0.0;
    double sumWeight = 0.0;
    numTransmittingHumans = 0;

    for (const Host::Human &human : population)
    {
        // NOTE: calculate availability relative to age at end of time step;
        // not my preference but consistent with TransmissionModel::getEIR().
        const double avail = human.perHostTransmission.relativeAvailabilityHetAge(sim::inYears(human.age(clock.ts1)));
        sumWeight += avail;

        vector<double> probTransGenotype_i(WithinHost::Genotypes::N());
        vector<double> probTransGenotype_l(WithinHost::Genotypes::N());
        const double pTransmit = human.withinHostModel->probTransmissionToMosquito(probTransGenotype_i, probTransGenotype_l);

        double riskTrans = 0.0;

        // Only to be consistent with old simulation runs when set to false
        // Setting this option to true will only affect reporting
        if (opt_vaccine_genotype == false)
            riskTrans = avail * pTransmit * human.vaccine.getFactor(interventions::Vaccine::TBV, clock.ts1);
        else
        {
            for (size_t g = 0; g < WithinHost::Genotypes::N(); ++g)
                riskTrans += (probTransGenotype_i[g] + probTransGenotype_l[g]) * human.vaccine.getFactor(interventions::Vaccine::TBV, clock.ts1, g);
            riskTrans *= avail;
        }

        if (riskTrans > 0.0) ++numTransmittingHumans;

        sumWt_kappa += riskTrans;
    }

    size_t lKMod = sim::moduloSteps(clock.ts1, laggedKappa.size()); // now
    if (population.size() == 0)
    {                             // this is valid
        laggedKappa[lKMod] = 0.0; // no humans: no infectiousness
    }
    else
    {
        if (!(sumWeight > DBL_MIN * 10.0))
        { // if approx. eq. 0, negative or an NaN
            ostringstream msg;
            msg << "sumWeight is invalid: " << sumWeight << ", " << sumWt_kappa << ", " << population.size();
            throw TRACED_EXCEPTION(msg.str(), util::Error::SumWeight);
        }
        laggedKappa[lKMod] = sumWt_kappa / sumWeight;
    }

    size_t tmod = sim::moduloYearSteps(clock.now);

    // Calculate time-weighted average of kappa
    _sumAnnualKappa += laggedKappa[lKMod] * initialisationEIR[tmod];
    if (tmod == sim::stepsPerYear() - 1)
    {
        _annualAverageKappa = _sumAnnualKappa / annualEIR; // inf or NaN when annualEIR is 0
        _sumAnnualKappa = 0.0;
    }

    return laggedKappa[lKMod]; // kappa now
}
}
}
