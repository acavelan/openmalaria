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

#ifndef Hmod_NonVectorModel
#define Hmod_NonVectorModel

#include "core/RunContext.h"
#include "TransmissionModel.h"
#include "util/vectors.h"
#include "Host/Human.h"
#include "mon/Monitoring.h"

namespace scnXml
{
class Entomology;
class NonVector;
}

namespace OM
{
namespace Transmission
{

//! Base transmission model, as used in Phase A
class NonVectorModel : public TransmissionModel
{
public:
    NonVectorModel(RunContext& ctx, vector<double> initEIR, int interventionMode, const scnXml::Entomology &entoData, const scnXml::NonVector &nonVectorData, int eipDuration);

    virtual ~NonVectorModel() {}

    virtual void init2(const vector<Host::Human> &population) { simulationMode = forcedEIR; } // no set-up needed; just indicate we're ready to roll

    virtual void scaleEIR(double factor)
    {
        util::vectors::scale(initialisationEIR, factor);
        annualEIR = util::vectors::sum(initialisationEIR);
    }

    virtual SimTime initIterate()
    {
        simulationMode = interventionMode;
        if (simulationMode != dynamicEIR) { return sim::zero(); }

        // initialKappa is used in calculateEIR
        assert(initialKappa.size() >= ctx.clock.steps_per_year);
        assert(mod_nn(initialKappa.size(), ctx.clock.steps_per_year) == 0);
        for (size_t i = ctx.clock.steps_per_year; i < initialKappa.size(); ++i)
        {
            initialKappa[mod_nn(i, ctx.clock.steps_per_year)] += initialKappa[i];
        }
        double factor = static_cast<double>(ctx.clock.steps_per_year) / static_cast<double>(initialKappa.size());
        initialKappa.resize(ctx.clock.steps_per_year);
        for (size_t i = 0; i < initialKappa.size(); ++i)
        {
            initialKappa[i] *= factor;
            // error check:
            if (!(initialKappa[i] > 0.0)) // if not positive
                throw TRACED_EXCEPTION("initialKappa is invalid", util::Error::InitialKappa);
        }

        return sim::zero(); // nothing to do
    }

    virtual void changeEIRIntervention(const scnXml::NonVector &nonVectorData);

    virtual void uninfectVectors()
    {
        if (simulationMode != dynamicEIR) cerr << "Warning: uninfectVectors is not efficacious with forced EIR" << endl;
        // reset history of human infectivity, which scales dynamic EIR:
        laggedKappa.assign(laggedKappa.size(), 0.0);
    }

    virtual double updateKappa(const vector<Host::Human> &population)
    {
        double currentKappa = TransmissionModel::updateKappa(population);
        if (simulationMode == forcedEIR) { initialKappa[ctx.clock.moduloSteps(ctx.clock.ts1(), initialKappa.size())] = currentKappa; }
        return currentKappa;
    }

    virtual void calculateEIR(Host::Human &human, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l) const
    {
        EIR_i.resize(1);
        EIR_l.resize(1); // no support for per-genotype tracking in this model (possible, but we're lazy)
        // where the full model, with estimates of human mosquito transmission is in use, use this:
        if (simulationMode == forcedEIR) { EIR_l[0] = initialisationEIR[ctx.clock.moduloYearSteps(ctx.clock.ts0())]; }
        else if (simulationMode == transientEIRknown)
        {
            // where the EIR for the intervention phase is known, obtain this from
            // the interventionEIR array
            EIR_l[0] = interventionEIR[ctx.clock.inSteps(ctx.clock.interv)];
        }
        else if (simulationMode == dynamicEIR)
        {
            EIR_l[0] = initialisationEIR[ctx.clock.moduloYearSteps(ctx.clock.ts0())];
            if (ctx.clock.interv >= sim::zero())
            {
                // we modulate the initialization based on the human infectiousness time steps ago in the
                // simulation relative to infectiousness at the same time-of-year, pre-intervention.
                // nspore gives the sporozoite development delay.
                size_t t = ctx.clock.inSteps(ctx.clock.ts1() - nSpore);
                EIR_l[0] *= laggedKappa[mod_nn(t, laggedKappa.size())] / initialKappa[mod_nn(t, ctx.clock.steps_per_year)];
            }
        }
        else
        {
            throw util::xml_scenario_error("Invalid simulation mode");
        }
#ifndef NDEBUG
        if (!(std::isfinite)(EIR_l[0]))
        {
            size_t t = ctx.clock.inSteps(ctx.clock.ts1() - nSpore);
            ostringstream msg;
            msg << "Error: non-vect eir is: " << EIR_l[0] << "\nlaggedKappa:\t" << laggedKappa[mod_nn(t, laggedKappa.size())]
                << "\ninitialKappa:\t" << initialKappa[mod_nn(t, ctx.clock.steps_per_year)] << endl;
            throw TRACED_EXCEPTION(msg.str(), util::Error::InitialKappa);
        }
#endif
        EIR_l[0] *= human.perHostTransmission.relativeAvailabilityHetAge(ageYears);

        mon::recordStat(mon::measure("innoculationsPerAgeGroup"), human, EIR_l[0]);
        mon::recordStat(mon::measure("innoculationsPerVector"), human, EIR_l[0]);
    }

private:
    virtual void checkpoint(istream &stream)
    {
        TransmissionModel::checkpoint(stream);
        nSpore &stream;
        interventionEIR &stream;
        initialKappa &stream;
    }

    virtual void checkpoint(ostream &stream)
    {
        TransmissionModel::checkpoint(stream);
        nSpore &stream;
        interventionEIR &stream;
        initialKappa &stream;
    }

    //! multiplier used to calculate a positive EIR value where the measured value is zero
    /*
      0.01 was old pv(30) Now a constant. min_EIR_mult multiplies the average EIR to obtain a value used for the EIR during periods when
      it is too low to be measureable. The value of 0.01 was old pv(30) Now a constant. 0.01 was old pv(30) Now a constant.
    */
    const double min_EIR_mult = 0.01;

    /** @brief Variables set by constructor.
     *
     * There shouldn't be any need to checkpoint these, at least before
     * interventions take effect. */
    //@{
    //! Variance of Infection Rate according to fielddata
    const double totalInfectionrateVariance = 1.0;

    const int nYearsWarmupData = 5;

    //! The duration of sporogony in time steps
    // doesn't need checkpointing
    SimTime nSpore = sim::never();
    //@}

    /** EIR per time interval during the intervention period. Value at index
     * ctx.clock.interv.inSteps() used each time-step.
     *
     * Units: inoculations per adult per time step */
    vector<double> interventionEIR;

    /** When simulationMode == dynamicEIR, this is the annual cycle of kappa
     * from the warmup phase and has length 1 year (in time steps).
     *
     * When simulationMode == equilibriumMode, this may be multiple years long and
     * is used to collect values of kappa (human infectiousness).
     *
     * In either case, ctx.clock.ts0().moduloSteps(initialKappa.size()) is the index
     * for the current infectiousness during updates. */
    vector<double> initialKappa;
};
} // namespace Transmission
} // namespace OM
#endif
