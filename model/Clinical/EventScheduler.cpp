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
#include "Clinical/EventScheduler.h"
#include "util/random.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"

#include <limits>

namespace OM {
namespace Clinical {
using namespace OM::util;

SimTime ClinicalEventScheduler::maxUCSeekingMemory(sim::never());
SimTime ClinicalEventScheduler::uncomplicatedCaseDuration(sim::never());
SimTime ClinicalEventScheduler::complicatedCaseDuration(sim::never());
SimTime ClinicalEventScheduler::extraDaysAtRisk(sim::never());
vector<double> ClinicalEventScheduler::cumDailyPrImmUCTS;
double ClinicalEventScheduler::neg_v;
double ClinicalEventScheduler::alpha;

double ClinicalEventScheduler::logOddsAbBase = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::logOddsAbNegTest = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::logOddsAbPosTest = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::logOddsAbNeed = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::logOddsAbInformal = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::oneMinusEfficacyAb = std::numeric_limits<double>::signaling_NaN();
AgeGroupInterpolator ClinicalEventScheduler::severeNmfMortality;

bool opt_non_malaria_fevers = false;

AgeGroupInterpolator ClinicalEventScheduler::NMF_need_antibiotic;
AgeGroupInterpolator ClinicalEventScheduler::MF_need_antibiotic;

// -----  static init  -----

// -----  construction, destruction and checkpointing  -----

void esClinicalUpdate(ClinicalModel& cm, RunContext& ctx, Human& human, double ageYears,
                      WithinHost::Pathogenesis::StatePair& pg) {
    WHInterface& withinHostModel = *human.withinHostModel;
    Episode::State newState = static_cast<Episode::State>(pg.state);
    util::streamValidate((newState << 16) & cm.esPgState);

    if (ctx.clock.ts0() == cm.es_time_recovery()) {
        if (cm.esPgState & Episode::DIRECT_DEATH) {
            // Human dies this time step (last day of risk of death)
            cm.doomed = DOOMED_COMPLICATED;

            cm.latestReport.update(ctx, human, cm.esPgState);
        } else {
            if (cm.esPgState & Episode::COMPLICATED) {
                const double pSequelae = ClinicalModel::pSequelaeInpatient.eval(ageYears);
                mon::recordStat(mon::measure("expectedSequelae"), human, pSequelae);
                // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
                if (human.rng().uniform_01() < pSequelae) {
                    cm.esPgState = Episode::State(cm.esPgState | Episode::SEQUELAE);
                } else {
                    cm.esPgState = Episode::State(cm.esPgState | Episode::RECOVERY);
                }
            } else {
                cm.esPgState = Episode::State(cm.esPgState | Episode::RECOVERY);
            }
            // report bout, at conclusion of episode:
            cm.latestReport.update(ctx, human, cm.esPgState);

            // Individual recovers (and is immediately susceptible to new cases)
            cm.esPgState = Episode::NONE; // recovery (reset to healthy state)

            // And returns to transmission (if was removed)
            human.perHostTransmission.outside_trans() = false;
        }
    }

    // The bug fixed here is that indirectMortality was not set if no event occurred
    bool indirectMortality = false;
    if (ClinicalModel::indirectMortBugfix)
        indirectMortality = pg.indirectMortality;

    if (newState & Episode::SICK) {
        // we have some new case: is it severe/complicated?
        if (newState & Episode::COMPLICATED) {
            if (cm.esPgState & Episode::COMPLICATED) {
                // previously severe: no events happen for course of medication
            } else {
                // previously healthy or UC: progress to severe
                cm.esPgState = Episode::State(cm.esPgState | newState | Episode::RUN_CM_TREE);
                indirectMortality = pg.indirectMortality;
                cm.es_case_start() = ctx.clock.ts0();
            }
        } else {
            // uncomplicated case (UC/UC2/NMF): is it new?
            if (cm.esPgState & Episode::SICK) {
                // previously UC; nothing to do
            } else {
                if (cm.es_case_start() < ctx.clock.ts0()) {
                    // new UC case
                    cm.esPgState = Episode::State(cm.esPgState | newState | Episode::RUN_CM_TREE);
                    indirectMortality = pg.indirectMortality;

                    // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
                    double uVariate = human.rng().uniform_01();
                    size_t i = 0; // units: days
                    for (; i < ClinicalEventScheduler::cumDailyPrImmUCTS.size(); ++i) {
                        if (uVariate < ClinicalEventScheduler::cumDailyPrImmUCTS[i]) {
                            goto gotDelay;
                        }
                    }
                    assert(false); // should have uVariate < 1 = ClinicalEventScheduler::cumDailyPrImmUCTS[len-1]
                gotDelay:
                    // set start time: current time plus length of delay (days)
                    cm.es_case_start() = ctx.clock.ts0() + sim::fromDays(i);
                }
            }
        }

        if (indirectMortality && cm.doomed == NOT_DOOMED)
            cm.doomed = -ctx.clock.one_TS(); // start indirect mortality countdown
    }

    if (cm.es_case_start() == ctx.clock.ts0() && (cm.esPgState & Episode::RUN_CM_TREE)) {
        // OK, we're about to run the CM tree
        cm.esPgState = Episode::State(cm.esPgState & ~Episode::RUN_CM_TREE);

        // If last treatment prescribed was in recent memory, consider second line.
        if (cm.es_time_treatment() + ClinicalModel::healthSystemMemory > ctx.clock.ts0()) {
            cm.esPgState = Episode::State(cm.esPgState | Episode::SECOND_CASE);
        }

        CMDTOut auxOut = ESCaseManagement::execute(ctx, human, ageYears, cm.esPgState);

        if (auxOut.treated) { // I.E. some treatment was given
            cm.es_time_treatment() = ctx.clock.ts0();
            if (cm.esPgState & Episode::COMPLICATED) {
                mon::recordEvent(mon::measure("nTreatments3"), human);
            } else {
                if (cm.esPgState & Episode::SECOND_CASE) {
                    mon::recordEvent(mon::measure("nTreatments2"), human);
                } else {
                    mon::recordEvent(mon::measure("nTreatments1"), human);
                }
            }
        }
        if (auxOut.screened) {
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
        }

        if (true /*FIXME auxOut.hospitalisation != CMAuxOutput::NONE*/) { // in hospital
            cm.esPgState = Episode::State(cm.esPgState | Episode::EVENT_IN_HOSPITAL);

            /*FIXME if( auxOut.hospitalisation == CMAuxOutput::DELAYED )
		++cm.es_case_start();*/
        }

        // Case fatality rate (first day of illness)
        // P(death) is some fixed input scaled by age-specific CFR.
        if ((cm.esPgState & Episode::COMPLICATED) && !(cm.esPgState & Episode::DIRECT_DEATH)) {
            double pDeath = ClinicalModel::caseFatalityRate.eval(ageYears);
            // NOTE: hospitalisation is hardwired (was: const bool inHospital = true,
            // FIXME auxOut.hospitalisation == CMAuxOutput::IMMEDIATE), so the
            // community-CFR (!inHospital) branch and the inHospital guard around
            // expectedHospitalDeaths were removed as dead (audit sweep).
            mon::recordStat(mon::measure("expectedDirectDeaths"), human, pDeath);
            mon::recordStat(mon::measure("expectedHospitalDeaths"), human, pDeath);
            // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
            if (human.rng().uniform_01() < pDeath) {
                cm.esPgState = Episode::State(cm.esPgState | Episode::DIRECT_DEATH | Episode::EVENT_FIRST_DAY);
                // Human is killed at end of time at risk
                //cm.es_time_recovery() += ClinicalEventScheduler::extraDaysAtRisk;	(no point updating; will be set later: ATORWD)
            }
            cm.es_prev_density() = withinHostModel.getTotalDensity();
        }

        if (opt_non_malaria_fevers) {
            if ((cm.esPgState & Episode::SICK) && !(cm.esPgState & Episode::COMPLICATED)) {
                // Have a NMF or UC malaria case

                /** Given a non-malaria fever, return the probability of it requiring
                 * treatment. */
                bool isMalarial = cm.esPgState & Episode::MALARIA;
                double pNeedTreat = isMalarial ? ClinicalEventScheduler::MF_need_antibiotic.eval(ageYears)
                                               : ClinicalEventScheduler::NMF_need_antibiotic.eval(ageYears);
                // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
                bool needTreat = human.rng().bernoulli(pNeedTreat);

                // Calculate chance of antibiotic administration:
                /*FIXME
                if( auxOut.AB_provider != CMAuxOutput::NO_AB ){
                    // Treatment seeking; may be at a health facility or from
                    // the informal sector.
                    double logOddsAb = ClinicalEventScheduler::logOddsAbBase - ClinicalEventScheduler::logOddsAbNeed * pNeedTreat;
                    if( auxOut.AB_provider == CMAuxOutput::FACILITY ){
                        if( auxOut.diagnostic == CMAuxOutput::NEGATIVE ){
                            logOddsAb += ClinicalEventScheduler::logOddsAbNegTest;
                        }else if( auxOut.diagnostic == CMAuxOutput::POSITIVE ){
                            logOddsAb += ClinicalEventScheduler::logOddsAbPosTest;
                        }
                        if( needTreat ){
                            logOddsAb += ClinicalEventScheduler::logOddsAbNeed;
                        }
                    }else{
                        assert( auxOut.AB_provider == CMAuxOutput::INFORMAL );
                        logOddsAb += ClinicalEventScheduler::logOddsAbInformal;
                    }
                    
                    double oddsTreatment = exp( logOddsAb );
                    pTreatment = oddsTreatment / (1.0 + oddsTreatment);
                }*/

                // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
                // Draw kept; the old `uniform_01() < pTreatment` branch was
                // provably never taken (pTreatment was hardwired 0.0, see the
                // FIXME above), so its empty body and treatmentEffectMult
                // (always 1.0) are removed.
                human.rng().uniform_01();

                // In a severe NMF case (only when not malarial), there is a
                // chance of death:
                if (needTreat) {
                    double pDeath = ClinicalEventScheduler::severeNmfMortality.eval(ageYears);
                    if (human.rng().uniform_01() < pDeath) {
                        cm.esPgState = Episode::State(cm.esPgState | Episode::DIRECT_DEATH);
                    }
                }
            }
        }
    } else {
        // No new event (haven't changed state this time step).

        // Case fatality rate (subsequent days)
        // Complicated case & at risk of death (note: ClinicalEventScheduler::extraDaysAtRisk <= 0)
        if ((cm.esPgState & Episode::COMPLICATED) && !(cm.esPgState & Episode::DIRECT_DEATH) &&
            (ctx.clock.ts0() < cm.es_time_recovery() + ClinicalEventScheduler::extraDaysAtRisk)) {
            // In complicated episodes, S(t), the probability of survival on
            // subsequent days t, is described by log(S(t)) = -v(Y(t)/Y(t-1)),
            // for parasite density Y(t). v_neg below is -v.
            // TODO: this model should be revisited, and if possible placed
            // within the WHFalciparum class (or a subclass).
            if (withinHostModel.getTotalDensity() > 0.0) { // avoid division by zero
                double parasiteReductionEffect = withinHostModel.getTotalDensity() / cm.es_prev_density();
                double pDeath = 1.0 - exp(ClinicalEventScheduler::neg_v * parasiteReductionEffect);
                // community fatality rate when not in hospital
                if (!(cm.esPgState & Episode::EVENT_IN_HOSPITAL))
                    pDeath = ClinicalModel::getCommunityCFR(pDeath);
                mon::recordStat(mon::measure("expectedDirectDeaths"), human, pDeath);
                if (cm.esPgState & Episode::EVENT_IN_HOSPITAL)
                    mon::recordStat(mon::measure("expectedHospitalDeaths"), human, pDeath);
                // parity: order — per-host Xoshiro stream: this call’s position in the draw sequence fixes the host RNG state
                if (human.rng().uniform_01() < pDeath) {
                    cm.esPgState = Episode::State(cm.esPgState | Episode::DIRECT_DEATH);
                    // Human is killed at end of time at risk
                    cm.es_time_recovery() = cm.es_time_recovery() +
                                            ClinicalEventScheduler::extraDaysAtRisk; // may be re-set later (see ATORWD)
                }
            }
            cm.es_prev_density() = withinHostModel.getTotalDensity();
        }
    }

    // Start of case. Not necessarily start of sickness due to treatment-seeking
    // delays and travel time.
    if (cm.es_case_start() == ctx.clock.ts0()) {
        // Patients in hospital are removed from the transmission cycle.
        // This should have an effect from the start of the next time step.
        // NOTE: This is not very accurate, but considered of little importance.
        if (cm.esPgState & Episode::EVENT_IN_HOSPITAL)
            human.perHostTransmission.outside_trans() = true;

        if (cm.esPgState & Episode::COMPLICATED) {
            // ClinicalEventScheduler::complicatedCaseDuration should to some respects be associated
            // with medication duration, however ongoing medications after
            // exiting hospital are OK and medications terminating before the
            // end of hospitalisation shouldn't matter too much if the person
            // can't recieve new infections due to zero transmission in hospital.
            cm.es_time_recovery() = ctx.clock.ts0() + ClinicalEventScheduler::complicatedCaseDuration;
            // Time should be adjusted to end of at-risk period when patient dies:
            if (cm.esPgState & Episode::DIRECT_DEATH) // death may already have been determined
                cm.es_time_recovery() =
                    cm.es_time_recovery() + ClinicalEventScheduler::extraDaysAtRisk; // ATORWD (search keyword)
        } else {
            cm.es_time_recovery() = ctx.clock.ts0() + ClinicalEventScheduler::uncomplicatedCaseDuration;
        }
    }

    // Remove on first models...
    if (cm.es_time_treatment() == ctx.clock.ts0())
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_TREATMENT);

    if (cm.esPgState & Episode::SICK)
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_BOUT);
}

}
}
