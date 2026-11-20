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

#ifndef Hmod_TransmissionModel
#define Hmod_TransmissionModel

#include <string.h>
#include <fstream>

#include "Global.h"
#include "core/RunContext.h"
#include "util/errors.h"

namespace scnXml { class NonVector; class Entomology; class Vector; }

#include "Transmission/PerHost.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/Anopheles/EmergenceRateEstimator.h"
#include "Population.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Genotypes.h"
#include "mon/Continuous.h"
#include "mon/Monitoring.h"
#include "util/StreamValidator.h"
#include "util/CommandLine.h"
#include "util/vectors.h"
#include "util/ModelOptions.h"

#include <cmath>
#include <cfloat>
#include <gsl/gsl_vector.h>

using namespace std::placeholders;

namespace OM
{
class Summary;
class Population;
namespace Host
{
class Human;
}
namespace Transmission
{



/** Variable describing current simulation mode. */
enum SimulationMode
{
    /** Initial mode. Indicates that initialisation still needs to happen (i.e.
     * it is an error if this mode is still set when getEIR is called). */
    preInit,
    /** Equilibrium mode (i.e. just apply a fixed EIR)
     *
     * This is used for the warm-up period and if we want to separate direct
     * effect of an intervention from indirect effects via transmission
     * intensity. The seasonal pattern and intensity of the EIR do not change
     * over years.
     *
     * In this mode vector calculations aren't run. */
    forcedEIR,
    /** Transient EIR known
     *
     * This is used to simulate an intervention that changes EIR, and where we
     * have measurements of the EIR over time during the intervention period.
     */
    transientEIRknown,
    /** EIR changes
     *
     * The simulation is driven by the EIR which changes dynamically during the
     * intervention phase as a function of the characteristics of the
     * interventions.
     *
     * Dependending on whether the Vector or NonVector model is in use, this
     * EIR may be calculated from a mosquito emergence rate or be an input EIR
     * scaled by the relative infectiousness of the humans. */
    dynamicEIR,
};

/// Abstract base class, defines behaviour of transmission models
class TransmissionModel
{
public:
    /// Which behaviour this (single, final) model class implements.
    enum class Kind { Vector, NonVector };

    /// Vector-model constructor (was VectorModel::VectorModel).
    TransmissionModel(RunContext& ctx, vector<double> initEIR, int interventionMode,
        vector<std::unique_ptr<Anopheles::AnophelesModel>> speciesList,
        vector<std::unique_ptr<Anopheles::EmergenceRateEstimator>> emergenceRateEstimatorsList,
        map<string, size_t> speciesIndexList, int populationSize);

    /// Non-vector constructor (was NonVectorModel; defined in xml/entomology.cpp).
    TransmissionModel(RunContext& ctx, vector<double> initEIR, int interventionMode,
        const scnXml::Entomology &entoData, const scnXml::NonVector &nonVectorData, int eipDuration);

protected:
    /// The run this model belongs to (clock accessed as ctx.clock).
    RunContext& ctx;

    /// Reads all entomological parameters from the input datafile.
    /// @param entoData input configuration for model
    /// @param nGenotypes number of genotypes transmission model is using
    TransmissionModel(RunContext& ctx, vector<double> initialisationEIR, int interventionMode, size_t nGenotypes, Kind kind,
        vector<std::unique_ptr<Anopheles::AnophelesModel>> speciesList,
        vector<std::unique_ptr<Anopheles::EmergenceRateEstimator>> emergenceRateEstimatorsList,
        map<string, size_t> speciesIndexList)
        : ctx(ctx)
        , kind(kind)
        , initialisationEIR(std::move(initialisationEIR))
        , simulationMode(forcedEIR)
        , interventionMode(interventionMode)
        , laggedKappa(ctx.clock.steps_per_year, 0.0)
        , // if using non-vector model, it will resize this
        annualEIR(0.0)
        , _annualAverageKappa(numeric_limits<double>::signaling_NaN())
        , _sumAnnualKappa(0.0)
        , tsAdultEIR(0.0)
        , surveyInputEIR(0.0)
        , surveyInputEIR_i(0.0)
        , surveyInputEIR_l(0.0)
        , surveySimulatedEIR(0.0)
        , surveySimulatedEIR_i(0.0)
        , surveySimulatedEIR_l(0.0)
        , adultAge(PerHost::adultAge(ctx.clock))
        , numTransmittingHumans(0)
        // m_rng(0,0) draws nothing from master; the vector kind re-seeds below
        , m_rng(0, 0)
        , initIterations(0)
        , species(std::move(speciesList))
        , emergenceRateEstimators(std::move(emergenceRateEstimatorsList))
        , speciesIndex(std::move(speciesIndexList))
    {
        // parity: order — vector kind only: seeding LocalRng from master_RNG
        // consumes master draws (Xoshiro seeds via 4 gen_u64); the non-vector
        // kind must not touch master here (old NonVectorModel never had an m_rng).
        if (kind == Kind::Vector) m_rng = LocalRng(ctx.master_RNG);
        // Set VACCINE_GENOTYPE option
        opt_vaccine_genotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

        mon::Continuous::registerCallback(ctx, "input EIR", "\tinput EIR", std::bind(&TransmissionModel::ctsCbInputEIR, this, _1));
        mon::Continuous::registerCallback(ctx, "simulated EIR", "\tsimulated EIR", std::bind(&TransmissionModel::ctsCbSimulatedEIR, this, _1));
        mon::Continuous::registerCallback(ctx, "human infectiousness", "\thuman infectiousness", std::bind(&TransmissionModel::ctsCbKappa, this, _1));
        mon::Continuous::registerCallback(ctx, "num transmitting humans", "\tnum transmitting humans", std::bind(&TransmissionModel::ctsCbNumTransmittingHumans, this, _1));
    }

public:
    bool isVector() const { return kind == Kind::Vector; }

    /// Get the map of species names to indicies (vector kind).
    const map<string, size_t> &getSpeciesIndexMap() { return speciesIndex; }

    /// ~Deallocate memory for TransmissionModel parameters and clean up
    ~TransmissionModel() {}

    /** Extra initialisation when not loading from a checkpoint, requiring
     * information from the human population structure. */
    void init2(const vector<Host::Human> &population)
    {
        if (kind == Kind::Vector) { init2_vector(population); }
        else { simulationMode = forcedEIR; } // no set-up needed; just indicate we're ready to roll
    }

    /// Checkpointing
    template <class S>
    void operator&(S &stream)
    {
        checkpoint(stream);
    }
    //@}

    /** Set some summary items.
     *
     * Overriding functions should call this base version too. */
    void summarize(RunContext& ctx)
    {
        mon::recordStat(ctx, mon::measure("nTransmit"), laggedKappa[ctx.clock.moduloSteps(ctx.clock.now(), laggedKappa.size())]);
        mon::recordStat(ctx, mon::measure("annAvgK"), _annualAverageKappa);

        if (!mon::isReported(ctx)) return; // cannot use counters below when not reporting

        double duration = ctx.clock.inSteps(ctx.clock.now() - lastSurveyTime);
        if (duration > 0.0)
        {
            mon::recordStat(ctx, mon::measure("inputEIR"), surveyInputEIR / duration);
            mon::recordStat(ctx, mon::measure("simulatedEIR"), surveySimulatedEIR / duration);
            mon::recordStat(ctx, mon::measure("simulatedEIR_Introduced"), surveySimulatedEIR_i / duration);
            mon::recordStat(ctx, mon::measure("simulatedEIR_Indigenous"), surveySimulatedEIR_l / duration);
        }

        surveyInputEIR = 0.0;
        surveySimulatedEIR = 0.0;
        surveySimulatedEIR_i = 0.0;
        surveySimulatedEIR_l = 0.0;
        lastSurveyTime = ctx.clock.now();

        if (kind == Kind::Vector) summarize_vector(ctx);
    }

    /** Scale the EIR used by the model.
     *
     * EIR is scaled in memory (so will affect this simulation).
     * XML data is not touched. */
    void scaleEIR(double factor)
    {
        if (kind == Kind::Vector) { scaleEIR_vector(factor); return; }
        util::vectors::scale(initialisationEIR, factor);
        annualEIR = util::vectors::sum(initialisationEIR);
    }

    SimTime initIterate()
    {
        if (kind == Kind::Vector) return initIterate_vector();

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

    /** Needs to be called each step of the simulation before Human::update().
     *
     * when the vector model is used this updates mosquito populations. */
    void vectorUpdate(const vector<Host::Human> &population)
    {
        if (kind == Kind::Vector) vectorUpdate_vector(population);
    }

    void changeEIRIntervention(const scnXml::NonVector &nonVectorData)
    {
        if (kind == Kind::Vector) changeEIRIntervention_vector(nonVectorData);
        else changeEIRIntervention_nonVector(nonVectorData);
    }

    /** Does per-time-step updates and returns the EIR (inoculation rate per host
     * per time step). Should be called exactly once per time-step (at least,
     * during the intervention period when ITNs may be in use).
     *
     * Non-vector:
     * During the pre-intervention phase, the EIR is forced, using values from
     * the XML file. During the main simulation phase, it may be calculated or
     * obtained from data in the XML file.
     *
     * Vector:
     * During vector initialisation phase, EIR is forced based on the EIR given
     * in the XML file as a Fourier Series. After endVectorInitPeriod() is called
     * the simulation switches to using dynamic EIR. advanceStep _must_ be
     * called before this function in order to return the correct value.
     *
     * @param human A reference to the human who's EIR is being calculated.
     *    The human's "per host transmission" potentially needs updating.
     * @param age Age of the human in time units
     * @param ageYears Age of the human in years
     * @param EIR Out-vector of EIR per parasite genotype. The length is also set
     *    by the called function. Where genotype tracking is not supported (e.g.
     *    the non-vector model), the length is set to one.
     * @returns the sum of EIR across genotypes
     */
    double getEIR(Host::Human &human, SimTime age, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l)
    {
        /* For the NonVector model, the EIR should just be multiplied by the
         * availability. For the Vector model, the availability is also required
         * for internal calculations, but again the EIR should be multiplied by the
         * availability. */
        calculateEIR(human, ageYears, EIR_i, EIR_l);

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

    /** Remove all current infections to mosquitoes, such that without re-
     * infection, humans will then be exposed to zero EIR. */
    void uninfectVectors()
    {
        if (kind == Kind::Vector) { uninfectVectors_vector(); return; }
        if (simulationMode != dynamicEIR) cerr << "Warning: uninfectVectors is not efficacious with forced EIR" << endl;
        // reset history of human infectivity, which scales dynamic EIR:
        laggedKappa.assign(laggedKappa.size(), 0.0);
    }

    /** Needs to be called each time-step after Human::update() to update summary
     * statististics related to transmission. Also returns kappa (the average
     * human infectiousness weighted by availability to mosquitoes). */
    double updateKappa(const vector<Host::Human> &population)
    {
        // We calculate kappa for output and the non-vector model.
        double sumWt_kappa = 0.0;
        double sumWeight = 0.0;
        numTransmittingHumans = 0;

        // Scratch hoisted out of the per-human loop (one allocation for the
        // sweep; assign() re-zeroes each row — same bytes as the old fresh
        // vector(N), same pattern as vectorUpdate).
        vector<double> probTransGenotype_i, probTransGenotype_l;

        for (const Host::Human &human : population)
        {
            // NOTE: calculate availability relative to age at end of time step;
            // not my preference but consistent with TransmissionModel::getEIR().
            const double avail = human.perHostTransmission.relativeAvailabilityHetAge(sim::inYears(human.age(ctx.clock.ts1())));
            sumWeight += avail;

            probTransGenotype_i.assign(ctx.genotypes.N(), 0.0);
            probTransGenotype_l.assign(ctx.genotypes.N(), 0.0);
            const double pTransmit = human.withinHostModel->probTransmissionToMosquito(probTransGenotype_i, probTransGenotype_l);

            double riskTrans = 0.0;

            // Only to be consistent with old simulation runs when set to false
            // Setting this option to true will only affect reporting
            // p7b B4 hoist: pure + empty⇒1.0 (see vectorUpdate_vector).
            const bool anyVac = human.has_vaccines();
            if (opt_vaccine_genotype == false)
                riskTrans = avail * pTransmit * (anyVac ? human.vaccine_factor(interventions::Vaccine::TBV) : 1.0);
            else
            {
                for (size_t g = 0; g < ctx.genotypes.N(); ++g)
                    riskTrans += (probTransGenotype_i[g] + probTransGenotype_l[g]) * (anyVac ? human.vaccine_factor(interventions::Vaccine::TBV, g) : 1.0);
                riskTrans *= avail;
            }

            if (riskTrans > 0.0) ++numTransmittingHumans;

            sumWt_kappa += riskTrans;
        }

        size_t lKMod = ctx.clock.moduloSteps(ctx.clock.ts1(), laggedKappa.size()); // now
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

        size_t tmod = ctx.clock.moduloYearSteps(ctx.clock.ts0());

        // Calculate time-weighted average of kappa
        _sumAnnualKappa += laggedKappa[lKMod] * initialisationEIR[tmod];
        if (tmod == ctx.clock.steps_per_year - 1)
        {
            _annualAverageKappa = _sumAnnualKappa / annualEIR; // inf or NaN when annualEIR is 0
            _sumAnnualKappa = 0.0;
        }

        if (kind == Kind::NonVector && simulationMode == forcedEIR)
        { initialKappa[ctx.clock.moduloSteps(ctx.clock.ts1(), initialKappa.size())] = laggedKappa[lKMod]; }

        return laggedKappa[lKMod]; // kappa now
    }

    void surveyEIR()
    {
        tsAdultEIR = tsAdultEntoInocs / tsNumAdults;
        tsAdultEntoInocs = 0.0;

        tsAdultEIR_i = tsAdultEntoInocs_i / tsNumAdults;
        tsAdultEntoInocs_i = 0.0;

        tsAdultEIR_l = tsAdultEntoInocs_l / tsNumAdults;
        tsAdultEntoInocs_l = 0.0;

        tsNumAdults = 0;

        surveyInputEIR += initialisationEIR[ctx.clock.moduloYearSteps(ctx.clock.ts0())];

        surveySimulatedEIR += tsAdultEIR;
        surveySimulatedEIR_i += tsAdultEIR_i;
        surveySimulatedEIR_l += tsAdultEIR_l;
    }

protected:
    /** Calculates the EIR individuals are exposed to.
     *
     * Call once per time-step: updates ITNs in vector model.
     *
     * @param host Transmission data for the human to calculate EIR for.
     * @param ageGroupData Age group of this host for availablility data.
     * @param EIR Out-vector. Set to the age- and heterogeneity-specific EIR an
     *    individual human is exposed to, per parasite genotype, in units of
     *    inoculations per day. Length set by callee. 
     *    _i for imported infections and _l for local infections */
    void calculateEIR(Host::Human &human, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l) const
    {
        if (kind == Kind::Vector) { calculateEIR_vector(human, ageYears, EIR_i, EIR_l); return; }

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

    void checkpoint(istream &stream)
    {
        checkpoint_common(stream);
        if (kind == Kind::Vector) checkpoint_vector(stream);
        else { nSpore &stream; interventionEIR &stream; initialKappa &stream; }
    }

    void checkpoint(ostream &stream)
    {
        checkpoint_common(stream);
        if (kind == Kind::Vector) checkpoint_vector(stream);
        else { nSpore &stream; interventionEIR &stream; initialKappa &stream; }
    }

    template<class S> void checkpoint_common(S &stream)
    {
        simulationMode &stream;
        interventionMode &stream;
        initialisationEIR &stream;
        laggedKappa &stream;
        annualEIR &stream;
        _annualAverageKappa &stream;
        _sumAnnualKappa &stream;
        tsAdultEIR &stream;
        surveyInputEIR &stream;
        surveySimulatedEIR &stream;
        surveySimulatedEIR_i &stream;
        surveySimulatedEIR_l &stream;
        lastSurveyTime &stream;
        adultAge &stream;
        numTransmittingHumans &stream;
    }

private:
    // The times here should be for the last updated index of arrays:
    void ctsCbInputEIR(ostream &stream)
    {
        int prevStep = (ctx.clock.now() - ctx.clock.one_TS()) / ctx.clock.one_TS();
        // Note: prevStep may be negative, hence util::mod not mod_nn:
        stream << '\t' << initialisationEIR[util::mod(prevStep, ctx.clock.steps_per_year)];
    }
    void ctsCbSimulatedEIR(ostream &stream) { stream << '\t' << tsAdultEIR; }
    void ctsCbKappa(ostream &stream)
    {
        // The latest time-step's kappa:
        stream << '\t' << laggedKappa[ctx.clock.moduloSteps(ctx.clock.now(), laggedKappa.size())];
    }
    void ctsCbNumTransmittingHumans(ostream &stream) { stream << '\t' << numTransmittingHumans; }

public:
    /** Entomological inoculation rate for adults during the
     * pre-intervention phase.
     *
     * Length: time-steps per year
     *
     * Index sim::now_mod_steps_per_year() corresponds to the EIR
     * acting on the current time-step: i.e. total inoculations since the
     * previous time-step.
     * Since time-step 0 is not calculated, initialisationEIR[0] is actually the
     * last value used (to calculate the state at the start of the second year).
     *
     * Units: infectious bites per adult per time step
     *
     * Not checkpointed; doesn't need to be except when a changeEIR intervention
     * occurs. */
    vector<double> initialisationEIR;

    /** The type of EIR calculation. Checkpointed. */
    int simulationMode;
    /** New simulation mode during intervention period. Not checkpointed. */
    int interventionMode;

    /** The probability of infection of a mosquito at each bite.
     * It is calculated as the average infectiousness per human.
     *
     * The value in index ctx.clock.ts1().moduloSteps(initialKappa.size()) is the
     * kappa from this time step (i.e. the infectiousness of humans at the end of
     * this step). Length depends on entomological incubation period from
     * non-vector model.
     *
     * Checkpointed. */
    vector<double> laggedKappa;

    /** Total annual infectious bites per adult.
     *
     * Checkpointed. */
    double annualEIR;

private:
    /*! annAvgKappa is the overall proportion of mosquitoes that get infected
     * allowing for the different densities in different seasons (approximating
     * relative mosquito density with the EIR).
     *
     * Checkpointed. */
    double _annualAverageKappa;

    /*! Used to calculate annAvgKappa.
     *
     * Checkpointed. */
    double _sumAnnualKappa;

    /// Adult-only EIR over the last update
    double tsAdultEIR, tsAdultEIR_i, tsAdultEIR_l;

    /** Per-time-step input EIR summed over inter-survey period.
     * Units: infectious bites/adult/inter-survey period. */
    double surveyInputEIR, surveyInputEIR_i, surveyInputEIR_l;

    /** Per-time-step simulated EIR summed over inter-survey period.
     * Units: infectious bites/adult/inter-survey period. */
    double surveySimulatedEIR, surveySimulatedEIR_i, surveySimulatedEIR_l;

    /** Time of last survey. */
    SimTime lastSurveyTime = sim::never();

    /// age at which an individual is considered an adult
    SimTime adultAge = sim::never();

    /// For "num transmitting humans" cts output.
    int numTransmittingHumans;

    // Reporting data. Doesn't need checkpointing due to reset every time-step.
    // accumulator for time step EIR of adults
    double tsAdultEntoInocs = 0.0, tsAdultEntoInocs_i = 0.0, tsAdultEntoInocs_l = 0.0;
    int tsNumAdults = 0; // accumulator for time step adults requesting EIR

    bool opt_vaccine_genotype = false;

    // -----  vector kind (was VectorModel; unassigned for the NonVector kind)  -----
    void init2_vector(const vector<Host::Human> &population);
    void scaleEIR_vector(double factor);
    SimTime initIterate_vector();
    void vectorUpdate_vector(const vector<Host::Human> &population);
    void calculateEIR_vector(Host::Human &human, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l) const;
    void uninfectVectors_vector();
    void summarize_vector(RunContext& ctx);
    void changeEIRIntervention_vector(const scnXml::NonVector &nonVectorData);
    void checkpoint_vector(istream &stream);
    void checkpoint_vector(ostream &stream);

    void ctsCbN_v0(ostream &stream);
    void ctsCbP_A(ostream &stream);
    void ctsCbP_Amu(ostream &stream);
    void ctsCbP_A1(ostream &stream);
    void ctsCbP_Ah(ostream &stream);
    void ctsCbP_df(ostream &stream);
    void ctsCbP_dif(ostream &stream);
    void ctsCbN_v(ostream &stream);
    void ctsCbO_v(ostream &stream);
    void ctsCbS_v(ostream &stream);
    void ctsCbAlpha(const vector<Host::Human> &population, ostream &stream);
    void ctsCbP_B(const vector<Host::Human> &population, ostream &stream);
    void ctsCbP_CD(const vector<Host::Human> &population, ostream &stream);
    void ctsNetInsecticideContent(const vector<Host::Human> &population, ostream &stream);
    void ctsIRSInsecticideContent(const vector<Host::Human> &population, ostream &stream);
    void ctsIRSEffects(const vector<Host::Human> &population, ostream &stream);
    void ctsCbResAvailability(ostream &stream);
    void ctsCbResRequirements(ostream &stream);

public:
    /// RNG used by the transmission model (vector kind)
    LocalRng m_rng;

    /// Number of iterations performed during initialization (vector kind).
    ///
    /// Special cases: 0 during initial one-human-lifespan warmup, <0 after
    /// initialisation.
    int initIterations;

    double meanAvail = 1.0;

    /** Per anopheles species data (vector kind).
     *
     * Array will be recreated by constructor, but some members of AnophelesModel
     * need to be checkpointed. */
    vector<std::unique_ptr<Anopheles::AnophelesModel>> species;
    vector<std::unique_ptr<Anopheles::EmergenceRateEstimator>> emergenceRateEstimators;

    /** A map of anopheles species/variant name to an index in species
     * (vector kind; doesn't need checkpointing). */
    map<string, size_t> speciesIndex;

    /// Deploy a vector population-level intervention (vector kind only).
    void deployVectorPopInterv(size_t instance);
    /// Deploy some traps (vector kind only).
    void deployVectorTrap(size_t instance, double number, SimTime lifespan);

public:
    // -----  non-vector kind (was NonVectorModel; unused for the Vector kind)  -----
    void changeEIRIntervention_nonVector(const scnXml::NonVector &nonVectorData);

    /// multiplier used to calculate a positive EIR value where the measured value is zero
    const double min_EIR_mult = 0.01;

    //! Variance of Infection Rate according to fielddata
    const double totalInfectionrateVariance = 1.0;

    const int nYearsWarmupData = 5;

    //! The duration of sporogony in time steps (doesn't need checkpointing)
    SimTime nSpore = sim::never();

    /** EIR per time interval during the intervention period (non-vector kind).
     * Units: inoculations per adult per time step */
    vector<double> interventionEIR;

    /** (non-vector kind) Annual cycle of kappa from the warmup phase.
     * Checkpointed. */
    vector<double> initialKappa;

    /// Kind tag: which behaviour this instance implements (was the class split)
    const Kind kind;
};

} // namespace Transmission
} // namespace OM
#endif
