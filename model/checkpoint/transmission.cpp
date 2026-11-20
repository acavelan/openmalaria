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
void checkpoint (checkpoint::Archive& a, SimpleDecayingValue& sdv) {
    a.value (sdv.deploy_t);
    // het = per-deployment RNG sample of the decay curve. Deployments are
    // replayed by interventions at restore; sampled clones aren't on the
    // wire yet (ponytail: extend with a hetFactor payload when a species
    // intervention needs round-tripping, see HANDOFF limits).
    bool hasSample = sdv.het != nullptr;
    a.boolean (hasSample);
    if (a.reading() && hasSample && sdv.het == nullptr)
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
    for (auto& sdv : v) util::checkpoint (a, sdv);
}
} // namespace

void checkpoint (Archive& a, AnophelesModel& m) {
    // --- per-species dynamics state (order mirrors the old v2 wire) ---
    a.value (m.mosq.seekingDeathRate);
    a.value (m.mosq.seekingDuration);
    a.value (m.mosq.probMosqSurvivalOvipositing);
    ckSdvVector (a, m.seekingDeathRateIntervs);
    ckSdvVector (a, m.probDeathOvipositingIntervs);
    a.vector (m.partialEIR_i);
    a.vector (m.partialEIR_l);
    a.value (m.EIRRotateAngle);
    a.vector (m.FSCoeffic);
    a.vector (m.forcedS_v);
    a.value (m.initNvFromSv);
    a.value (m.initOvFromSv);
    ckSdvVector (a, m.emergenceReduction);
    a.vector (m.mosqEmergeRate);
    a.vector (m.quinquennialS_v);
    a.value (m.initNv0FromSv);
    a.value (m.initSvFromEIR);
    a.value (m.mosq.restDuration);
    a.value (m.mosq.EIPDuration);
    a.value (m.N_v_length);
    a.vector (m.P_A);
    a.vector (m.P_df);
    a.vector (m.P_dif_i);
    a.vector (m.P_dif_l);
    a.vector (m.P_dff);
    a.vector (m.N_v);
    a.vector (m.O_v_i);
    a.vector (m.O_v_l);
    a.vector (m.S_v_i);
    a.vector (m.S_v_l);
    a.vector (m.P_Amu);
    a.vector (m.P_A1);
    a.vector (m.P_Ah);
    a.vector (m.fArray);
    a.vector (m.ftauArray);
    a.vector (m.uninfected_v);
    a.value (m.timeStep_N_v0);
    bool mpd = m.simpleMPD != nullptr;
    a.boolean (mpd);
    if (a.reading() && mpd && !m.simpleMPD)
        throw checkpoint::checkpoint_error (
            "checkpoint has MPD state; scenario doesn't enable MPD");
    if (mpd) {
        a.vector (m.simpleMPD->quinquennialOvipositing);
        a.value (m.simpleMPD->developmentDuration);
        a.value (m.simpleMPD->probPreadultSurvival);
        a.value (m.simpleMPD->fEggsLaidByOviposit);
        a.vector (m.simpleMPD->invLarvalResources);
        a.vector (m.simpleMPD->nOvipositingDelayed);
    }
}

void checkpoint (Archive& a, EmergenceRateEstimator& e) {
    a.value (e.scaleFactor);
    a.value (e.shiftAngle);
    a.boolean (e.rotated);
    a.boolean (e.scaled);
}

} // namespace Anopheles

void checkpoint (Archive& a, TransmissionModel& tm) {
    a.rng (tm.m_rng);
    a.value (tm.initIterations);
    a.value (tm.meanAvail);
    // --- shared state (both kinds) ---
    a.value (tm.simulationMode);
    a.value (tm.interventionMode);
    a.vector (tm.initialisationEIR);
    a.vector (tm.laggedKappa);
    a.value (tm.annualEIR);
    a.value (tm._annualAverageKappa);
    a.value (tm._sumAnnualKappa);
    a.value (tm.tsAdultEIR);
    a.value (tm.tsAdultEIR_i);
    a.value (tm.tsAdultEIR_l);
    a.value (tm.surveyInputEIR);
    a.value (tm.surveyInputEIR_i);
    a.value (tm.surveyInputEIR_l);
    a.value (tm.surveySimulatedEIR);
    a.value (tm.surveySimulatedEIR_i);
    a.value (tm.surveySimulatedEIR_l);
    a.value (tm.lastSurveyTime);
    a.value (tm.adultAge);
    a.value (tm.numTransmittingHumans);
    if (tm.kind == TransmissionModel::Kind::Vector) {
        uint64_t n = tm.species.size();
        a.value (n);
        if (a.reading() && n != tm.species.size())
            throw checkpoint::checkpoint_error ("tm.species count mismatch");
        for (auto& s : tm.species) Anopheles::checkpoint (a, *s);
        if (a.reading() && n != tm.emergenceRateEstimators.size())
            throw checkpoint::checkpoint_error ("estimator count mismatch");
        for (auto& e : tm.emergenceRateEstimators) Anopheles::checkpoint (a, *e);
    } else {
        a.value (tm.nSpore);
        a.vector (tm.interventionEIR);
        a.vector (tm.initialKappa);
    }
}

} // namespace Transmission

namespace checkpoint {

void ckTransmission (Archive& a, RunContext&,
        Transmission::TransmissionModel& transmission) {
    Transmission::checkpoint (a, transmission);
}

} // namespace checkpoint
} // namespace OM
