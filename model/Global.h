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

/** The Global module is only for items wanted nearly everywhere.
 * It is not for miscelaneous stuff and should only include headers wanted nearly everywhere.
 */
#ifndef Hmod_Global
#define Hmod_Global

#ifdef _MSC_VER
#pragma warning(disable: 4290)  // disable some warnings on MSVC
// Avoid min/max as macros on windows (breaks numeric_limits<T>::max())
#define NOMINMAX
#endif

// unit32_t and similar, mostly for compatibility with MSVC:
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <memory>

// Checkpointing and time-step operations are used _everywhere_:
#include "util/checkpoint.h"
#include "util/mod.h"

#define PI 3.14159265358979323846

/******************************************************************************
 * SimTime is a number of days (integer). (Formerly model/sim.h.)
 *
 * Only stateless time helpers live here. Everything that depends on the time
 * step interval or on the simulation clock lives in OM::SimClock (owned by
 * OM::RunContext, see core/RunContext.h): the now/ts0/ts1/nowOrTs0/nowOrTs1/
 * latestTs0 accessors, start/end date, intervention time, oneTS, fromTS,
 * inSteps, moduloSteps, moduloYearSteps, fromYearsN, fromYearsD,
 * roundToTSFromDays, stepsPerYear, yearsPerStep and the init functions.
 *****************************************************************************/
typedef int SimTime;

// Forward declaration the old sim.h provided to the whole tree (global
// ::scnXml — the generated schema namespace).
namespace scnXml {
    class Scenario;
}

namespace OM {
	
using namespace util::checkpoint;
using util::mod;
using util::mod_nn;

namespace util {
    class ModelNameProvider;
}

inline int floorToInt( double x ){
	return static_cast<int>(std::floor(x));
}

/** Stateless simulation-time constants and conversions. */
class sim {
public:
    /// Number of days in a year; defined as 365 (leap years are not simulated).
    static const int DAYS_IN_YEAR = 365;

    /// Convert to years
    static inline double inYears(SimTime d) { return d * (1.0 / DAYS_IN_YEAR); }

    /// Duration zero and the time at the start of the simulation.
    static inline SimTime zero(){ return SimTime(0); }

    static inline SimTime origin(){ return SimTime(0); }

    /** Special value representing a time point always in the past, such that
     * never() + x < zero() and x - never() will not to overflow for all valid
     * simulation times x (including any value now() may take as well as
     * never() and future()). */
    static inline SimTime never(){ return SimTime(-0x3FFFFFFF); }

    /** Special value representing a time point always in the future, such that
     * now() < future() and now() + future() does not overflow. */
    static inline SimTime future(){ return SimTime(0x3FFFFFFF); }

    /// One day.
    static inline SimTime oneDay(){ return SimTime(1); }

    /// One year. See DAYS_IN_YEAR.
    static inline SimTime oneYear(){ return SimTime(DAYS_IN_YEAR); }

    /// Duration in days. Should be fast (currently no conversion required).
    static inline SimTime fromDays(int days){ return SimTime(days); }

    /// Convert from a whole number of years.
    static inline SimTime fromYearsI(int years){
        return SimTime(DAYS_IN_YEAR * years);
    }
};

}
#endif
