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

/** Clock module (plan §2: one clock, no time-state in free helpers).
 *
 * Everything time-related lives here:
 *  - SimTime: a number of days (integer), global typedef.
 *  - sim: stateless constants/conversions (never/future/inYears/...) —
 *    plain namespace, no clock dependency (formerly the  static-only
 *    class of model/sim.h).
 *  - SimClock: the simulation clock (time-step interval, dates, t0/t1,
 *    update markers). Owned by RunContext; init/initInterval defined in
 *    run.cpp.
 **/
#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>

#include "util/mod.h"

// Old sim.h included its header chain before SimTime was needed;
// keep the whole-tree forward declarations the old sim.h provided:
namespace scnXml {
class Scenario;
}

/// SimTime is a number of days (integer).
typedef int SimTime;

namespace OpenMalaria {

class RunContext; // SimClock::init takes a reference (defined in RunContext.h)

class ModelNameProvider;

inline int floorToInt(double x) { return static_cast<int>(std::floor(x)); }

/** Stateless simulation-time constants and conversions. */
/// Number of days in a year; defined as 365 (leap years are not simulated).
constexpr int DAYS_IN_YEAR = 365;

/// Convert to years
inline double inYears(SimTime d) { return d * (1.0 / DAYS_IN_YEAR); }

/// Duration zero and the time at the start of the simulation.
inline SimTime zero() { return SimTime(0); }

inline SimTime origin() { return SimTime(0); }

/** Special value representing a time point always in the past, such that
     * never() + x < zero() and x - never() will not to overflow for all valid
     * simulation times x (including any value now() may take as well as
     * never() and future()). */
inline SimTime never() { return SimTime(-0x3FFFFFFF); }

/** Special value representing a time point always in the future, such that
     * now() < future() and now() + future() does not overflow. */
inline SimTime future() { return SimTime(0x3FFFFFFF); }

/// One day.
inline SimTime oneDay() { return SimTime(1); }

/// One year. See DAYS_IN_YEAR.
inline SimTime oneYear() { return SimTime(DAYS_IN_YEAR); }

/// Duration in days. Should be fast (currently no conversion required).
inline SimTime fromDays(int days) { return SimTime(days); }

/// Convert from a whole number of years.
inline SimTime fromYearsI(int years) { return SimTime(DAYS_IN_YEAR * years); }

class SimClock {
  public:
    /// Days per time step (1 or 5; set by initInterval).
    int interval = 0;
    /// Number of time steps in one year.
    size_t steps_per_year = 0;
    /// A cached value: one year divided by one time step.
    double years_per_step = 0.0;

    /// Maximum possible age of a human.
    SimTime max_human_age = never();
    /// The starting date of the simulation.
    SimTime start_date = never();
    /// The ending date of the simulation.
    SimTime end_date = never();

    SimTime t0 = never();
    SimTime t1 = never();
    /// Time relative to the start of the intervention period.
    SimTime interv = never();

#ifndef NDEBUG
    bool in_update = false; // only true during human/population/transmission update
#endif

    ///@brief Conversions to other types/units
    //@{
    /// Convert to time steps (rounding down)
    int inSteps(SimTime d) const { return d / interval; }

    /// Return this time in time steps modulo some positive integer.
    int moduloSteps(SimTime d, int denominator) const { return mod_nn(d / interval, denominator); }

    /// Return this time in time steps modulo some positive integer.
    int moduloYearSteps(SimTime d) const { return mod_nn(d / interval, steps_per_year); }

    /// One time step (currently either one or five days).
    SimTime one_TS() const { return SimTime(interval); }

    /// Convert from time steps.
    SimTime fromTS(int ts) const { return SimTime(int(one_TS()) * ts); }

    /** Convert from years to nearest time step. */
    SimTime from_years_n(double years) const { return roundToTSFromDays(DAYS_IN_YEAR * years); }

    /** Convert from years, rounding down to the next time step. */
    SimTime from_years_d(double years) const { return fromTS(floorToInt(steps_per_year * years)); }

    /** Round to the nearest time-step, where input is in days. */
    SimTime roundToTSFromDays(double days) const { return fromTS(floorToInt(days / interval + 0.5)); }
    //@}

    ///@brief Access simulation time variables
    //@{
    /** Time at the beginning of a time step update.
     *
     * This is what is mostly used during an update. It is never negative and
     * increases throughout the simulation. */
    SimTime ts0() const {
        assert(in_update); // should only be used during updates
        return t0;
    }
    /** Time at the end of a time step update.
     *
     * During an update, ts0() + one_TS() = ts1(). Neither this nor ts0 should
     * be used outside of updates. */
    SimTime ts1() const {
        assert(in_update); // should only be used during updates
        return t1;
    }
    /**
     * Time steps are mid-day to mid-day, and this is the time at mid-day (i.e.
     * this equals ts1 from the last step and ts0 from the next one).
     *
     * This is for monitoring and intervention deployment which happens between
     * updates. Cannot be used during human or vector update.
     */
    SimTime now() const {
        assert(!in_update); // only for use outside of step updates
        return t0;          // which is equal to t1 outside of updates, but that's a detail
    }
    /** During updates, this is ts0; between, this is now. */
    SimTime nowOrTs0() const { return t0; }
    /** During updates, this is ts1; between, this is now. */
    SimTime nowOrTs1() const { return t1; }
    /** During updates, this is ts0; between, it is now - 1. */
    SimTime latestTs0() const { return t1 - one_TS(); }
    //@}

    ///@brief Access intervention-time variables
    //@{
    /// The current date. Only valid during the intervention phase (the warmup
    /// duration is not known in advance; before that this is very negative).
    /// Intervention deployment times are relative to this date.
    SimTime interv_date() const { return start_date + interv; }
    //@}

    /// Initial set-up: called by run_simulation via Scenario::load.
    void init(RunContext& ctx, const scnXml::Scenario& scenario, ModelNameProvider mnp);

    /** Sets up time step duration: either a named model's interval or the
     * scenario's explicit interval. Must be 1 or 5 days (checked by
     * XMLChecker::PerformPostValidationChecks before seeding). */
    void initInterval(const scnXml::Scenario& scenario, ModelNameProvider mnp);

    /// Start of update: called by the update loop.
    void start_update() {
        t1 = t1 + one_TS();
#ifndef NDEBUG
        in_update = true;
#endif
    }
    /// End of update: called by the update loop.
    void end_update() {
#ifndef NDEBUG
        in_update = false;
#endif
        t0 = t1;
        interv = interv + one_TS();
    }
};

}
