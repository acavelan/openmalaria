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

// Wire section: transmission model (kind state; vector species or
// non-vector EIR data). Config (species list, MPD params...) is rebuilt by
// init; only what the run mutates goes on the wire.

#include "sections.h"

#include "Transmission/TransmissionModel.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/Anopheles/EmergenceRateEstimator.h"
#include "util/SimpleDecayingValue.h"

namespace OM {

using checkpoint::Archive;

namespace util {
void SimpleDecayingValue::checkpoint (checkpoint::Archive& a) {
    a.value (deploy_t);
    // het = per-deployment RNG sample of the decay curve. Deployments are
    // replayed by interventions at restore; sampled clones aren't on the
    // wire yet (ponytail: extend with a hetFactor payload when a species
    // intervention needs round-tripping, see HANDOFF limits).
    bool hasSample = het != nullptr;
    a.boolean (hasSample);
    if (a.reading() && hasSample && het == nullptr)
        throw checkpoint::checkpoint_error (
            "checkpoint: sampled decay deployment not yet supported");
}
} // namespace util

namespace Transmission {
namespace Anopheles {

namespace {
/// Count must match the config-built size (species/interventions lists).
void ckSdvVector (Archive& a, std::vector<util::SimpleDecayingValue>& v) {
    uint64_t n = v.size();
    a.value (n);
    if (a.reading() && n != v.size())
        throw checkpoint::checkpoint_error (
            "species intervention count mismatch");
    for (auto& sdv : v) sdv.checkpoint (a);
}
} // namespace

void AnophelesModel::checkpoint (Archive& a) {
    // --- per-species dynamics state (order mirrors the old v2 wire) ---
    a.value (mosq.seekingDeathRate);
    a.value (mosq.seekingDuration);
    a.value (mosq.probMosqSurvivalOvipositing);
    ckSdvVector (a, seekingDeathRateIntervs);
    ckSdvVector (a, probDeathOvipositingIntervs);
    a.vector (partialEIR_i);
    a.vector (partialEIR_l);
    a.value (EIRRotateAngle);
    a.vector (FSCoeffic);
    a.vector (forcedS_v);
    a.value (initNvFromSv);
    a.value (initOvFromSv);
    ckSdvVector (a, emergenceReduction);
    a.vector (mosqEmergeRate);
    a.vector (quinquennialS_v);
    a.value (initNv0FromSv);
    a.value (initSvFromEIR);
    a.value (mosq.restDuration);
    a.value (mosq.EIPDuration);
    a.value (N_v_length);
    a.vector (P_A);
    a.vector (P_df);
    a.vector (P_dif_i);
    a.vector (P_dif_l);
    a.vector (P_dff);
    a.vector (N_v);
    a.vector (O_v_i);
    a.vector (O_v_l);
    a.vector (S_v_i);
    a.vector (S_v_l);
    a.vector (P_Amu);
    a.vector (P_A1);
    a.vector (P_Ah);
    a.vector (fArray);
    a.vector (ftauArray);
    a.vector (uninfected_v);
    a.value (timeStep_N_v0);
    bool mpd = simpleMPD != nullptr;
    a.boolean (mpd);
    if (a.reading() && mpd && !simpleMPD)
        throw checkpoint::checkpoint_error (
            "checkpoint has MPD state; scenario doesn't enable MPD");
    if (mpd) {
        a.vector (simpleMPD->quinquennialOvipositing);
        a.value (simpleMPD->developmentDuration);
        a.value (simpleMPD->probPreadultSurvival);
        a.value (simpleMPD->fEggsLaidByOviposit);
        a.vector (simpleMPD->invLarvalResources);
        a.vector (simpleMPD->nOvipositingDelayed);
    }
}

void EmergenceRateEstimator::checkpoint (Archive& a) {
    a.value (scaleFactor);
    a.value (shiftAngle);
    a.boolean (rotated);
    a.boolean (scaled);
}

} // namespace Anopheles

void TransmissionModel::checkpoint (Archive& a) {
    a.rng (m_rng);
    a.value (initIterations);
    a.value (meanAvail);
    // --- shared state (both kinds) ---
    a.value (simulationMode);
    a.value (interventionMode);
    a.vector (initialisationEIR);
    a.vector (laggedKappa);
    a.value (annualEIR);
    a.value (_annualAverageKappa);
    a.value (_sumAnnualKappa);
    a.value (tsAdultEIR);
    a.value (tsAdultEIR_i);
    a.value (tsAdultEIR_l);
    a.value (surveyInputEIR);
    a.value (surveyInputEIR_i);
    a.value (surveyInputEIR_l);
    a.value (surveySimulatedEIR);
    a.value (surveySimulatedEIR_i);
    a.value (surveySimulatedEIR_l);
    a.value (lastSurveyTime);
    a.value (adultAge);
    a.value (numTransmittingHumans);
    if (kind == Kind::Vector) {
        uint64_t n = species.size();
        a.value (n);
        if (a.reading() && n != species.size())
            throw checkpoint::checkpoint_error ("species count mismatch");
        for (auto& s : species) s->checkpoint (a);
        if (a.reading() && n != emergenceRateEstimators.size())
            throw checkpoint::checkpoint_error ("estimator count mismatch");
        for (auto& e : emergenceRateEstimators) e->checkpoint (a);
    } else {
        a.value (nSpore);
        a.vector (interventionEIR);
        a.vector (initialKappa);
    }
}

} // namespace Transmission

namespace checkpoint {

void ckTransmission (Archive& a, RunContext&,
        Transmission::TransmissionModel& transmission) {
    transmission.checkpoint (a);
}

} // namespace checkpoint
} // namespace OM
