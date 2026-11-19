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

#ifndef Hmod_core_RunContext
#define Hmod_core_RunContext

#include "Global.h"
#include "util/mod.h"
#include "util/random.h"
#include "Host/WithinHost/Genotypes.h"
#include "PkPd/Drug/LSTMDrugType.h"
#include "PkPd/LSTMMedicate.h"
#include <memory>
#include <string>
#include <vector>

namespace scnXml {
class Scenario;
}

namespace OM {

struct RunExtensions;   // core/RunExtensions.h: registries needing complete RunContext

namespace util {
class ModelNameProvider;
}

/** The simulation clock and scenario time constants: everything the old
 * sim class stored in globals. Constructed fresh with each RunContext, so
 * independent runs never share clock state. Values and construction order
 * match the old sim:: statics exactly. */
class SimClock {
public:
    /// Days per time step (1 or 5; set by initInterval).
    int interval = 0;
    /// Number of time steps in one year.
    size_t steps_per_year = 0;
    /// A cached value: one year divided by one time step.
    double years_per_step = 0.0;

    /// Maximum possible age of a human.
    SimTime max_human_age = sim::never();
    /// The starting date of the simulation.
    SimTime start_date = sim::never();
    /// The ending date of the simulation.
    SimTime end_date = sim::never();

    SimTime t0 = sim::never();
    SimTime t1 = sim::never();
    /// Time relative to the start of the intervention period.
    SimTime interv = sim::never();

#ifndef NDEBUG
    bool in_update = false;     // only true during human/population/transmission update
#endif

    ///@brief Conversions to other types/units
    //@{
    /// Convert to time steps (rounding down)
    int inSteps(SimTime d) const { return d / interval; }

    /// Return this time in time steps modulo some positive integer.
    int moduloSteps(SimTime d, int denominator) const {
        return util::mod_nn(d / interval, denominator);
    }

    /// Return this time in time steps modulo some positive integer.
    int moduloYearSteps(SimTime d) const {
        return util::mod_nn(d / interval, steps_per_year);
    }

    /// One time step (currently either one or five days).
    SimTime one_TS() const { return SimTime(interval); }

    /// Convert from time steps.
    SimTime fromTS(int ts) const { return SimTime(int(one_TS()) * ts); }

    /** Convert from years to nearest time step. */
    SimTime from_years_n(double years) const {
        return roundToTSFromDays(sim::DAYS_IN_YEAR * years);
    }

    /** Convert from years, rounding down to the next time step. */
    SimTime from_years_d(double years) const {
        return fromTS( floorToInt(steps_per_year * years) );
    }

    /** Round to the nearest time-step, where input is in days. */
    SimTime roundToTSFromDays(double days) const {
        return fromTS( floorToInt( days / interval + 0.5 ));
    }
    //@}

    ///@brief Access simulation time variables
    //@{
    /** Time at the beginning of a time step update.
     *
     * This is what is mostly used during an update. It is never negative and
     * increases throughout the simulation. */
    SimTime ts0() const {
        assert(in_update);      // should only be used during updates
        return t0;
    }
    /** Time at the end of a time step update.
     *
     * During an update, ts0() + one_TS() = ts1(). Neither this nor ts0 should
     * be used outside of updates. */
    SimTime ts1() const {
        assert(in_update);      // should only be used during updates
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
        assert(!in_update);     // only for use outside of step updates
        return t0;   // which is equal to t1 outside of updates, but that's a detail
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
    void init( RunContext& ctx, const scnXml::Scenario& scenario, util::ModelNameProvider mnp );

    /** Sets up time step duration: either a named model's interval or the
     * scenario's explicit interval. Must be 1 or 5 days (checked by
     * XMLChecker::PerformPostValidationChecks before seeding). */
    void initInterval( const scnXml::Scenario& scenario, util::ModelNameProvider mnp );

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

/** Everything one simulation run owns: the clock now; further run state
 * (parameters, interventions, monitoring, RNGs, ...) moves in here through
 * the rest of phase 3. Owned by the caller of run_simulation(); never a
 * global. */
class RunContext {
public:
    ~RunContext();      // out-of-line: deletes ext (incomplete type here)
    
    SimClock clock;
    /** Master RNG. (0, 0) matches the old global's construction — the
     * percentile-stream seed depends on it; run_simulation reseeds from
     * --seed/Iseed after validation. */
    util::MasterRng master_RNG{0, 0};
    /// Genotype registry (built by Genotypes::init during scenario load).
    WithinHost::Genotypes genotypes;
    /// Drug-type registry (PkPd; built by LSTMDrugType::init).
    PkPd::DrugRegistry drugs;
    /// PK/PD treatment-course registry (LSTMTreatments::init).
    PkPd::TreatmentRegistry treatments;
    /// Registries whose types need a complete RunContext; created by
    /// ensure_extensions() (model/core/RunExtensions.h) at run start,
    /// owned and deleted by ~RunContext (raw pointer: gcc 16 instantiates
    /// unique_ptr's deleter on incomplete types at class scope).
    RunExtensions* ext = nullptr;
    /// Intervention states saved for checkpoints (checkpoint read/write).
};

}
#endif
