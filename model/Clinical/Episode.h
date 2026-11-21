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

#pragma once

#include "Clock.h"
#include "Host/WithinHost/PathogenesisModel.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include <ostream>

namespace OpenMalaria {
class Human;

/** Summary of clinical events during a caseManagementMemory period, in one individual.
 *
 * Terminology:
 * An "event" is an instantaeneous alteration of state.
 * A "bout" is a single fever or other sickness (falling sick to recovery).
 * An "episode" is a clinical-view of sickness caused by a malaria infection.
 * There's no reason an "episode" can't span multiple infections and multiple
 * bouts of sickness and recovery (the most severe is reported). */
class Episode {
  public:
    /**
     * State description flags. Note that values are carefully engineered to
     * work with State values.
     */
    using State = OpenMalaria::State;

    explicit Episode(RunContext& ctx)
        : time(never()), surveyPeriod(NOT_USED), ageGroup(0), cohortSet(0), state(NONE), rctx(&ctx) {
        healthSystemMemoryFix = ModelOptions::option(HEALTH_SYSTEM_MEMORY_FIX);
    };

    ~Episode();

    /// Report anything pending, as on destruction
    void flush();

    /** Report an episode, its severity, and any outcomes it entails.
     *
     * @param human The human whose info is being reported
     * @param newState The severity (diagnosis) and outcome.
     */
    void update(RunContext& ctx, const Human& human, Episode::State newState);

    /// Time of event, potentially never
    SimTime time = never();
    /// Survey during which the event occured
    size_t surveyPeriod;
    /// Age group of the individual when the episode's first bout occurred
    size_t ageGroup;
    /// Cohort membership
    uint32_t cohortSet;
    /// Descriptor of state, containing reporting info. Not all information will
    /// be reported (e.g. indirect deaths are reported independantly).
    Episode::State state;

    InfectionOrigin infectionType = InfectionOrigin::Indigenous;

  private:
    /** Report a clinical episode.
     *
     * From _state, an episode is reported based on severity (SICK,
     * MALARIA or COMPLICATED), and any outcomes are reported: RECOVERY (in
     * hospital, i.e. with EVENT_IN_HOSPITAL, only), SEQUELAE and DIRECT_DEATH
     * (both in and out of hospital). */
    void report();

    bool healthSystemMemoryFix = false;
    /// The run this episode reports to (set at construction).
    RunContext* rctx;
};

}
