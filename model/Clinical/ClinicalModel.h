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

#ifndef Hmod_ClinicalModel
#define Hmod_ClinicalModel

#include "Host/Human.h"
#include "Episode.h"
#include "Clinical/ClinicalParams.h"
#include <memory>

namespace scnXml{
    class Scenario;
    class HealthSystem;
}
namespace OM { namespace interventions { struct InterventionRegistry; } }
namespace OM { namespace Host { struct NeonatalMortality; } }
namespace OM { namespace Clinical {
    struct InfantMortality;
    using Host::Human;

/** The clinical model models the effects of sickness dependant on malarial
 * parasite densities and administers anti-malaria treatments via the drug
 * model (or in a simpler case, directly clearing infections).
 * 
 * So far, sickness types include uncomplicated and severe malaria cases and
 * non-malaria sickness.
 * 
 * Patient outcomes include full recovery, recovery with sequelae and death.
 * 
 * Reporting includes patient outcome and potentially drug usage and use of
 * RDTs (Rapid Diagnostic Tests) for costing purposes. */\
class ClinicalModel
{
public:
    /** Return a new ClinicalModel of the kind selected by params.
     *
     * @param tSF treatment seeking factor, passed to CaseManagementModel */
    static unique_ptr<ClinicalModel> createClinicalModel (double tSF, const ClinicalParams& params);
    
    /// Destructor
    virtual ~ClinicalModel ();
    
    /** Returns true if the human has been killed by some means.
     * 
     * Also kills the human if he/she reaches the simulation age limit.
     */
    bool isDead( SimTime age );
    
    /** Run main part of the model: determine the sickness status and any
     * treatment for the human.
     * 
     * @param ageYears Age of human.
     * @param newBorn True if human age is one time step old */
    void update (Human& human, double ageYears, bool newBorn, Host::NeonatalMortality& neonatal,
                 const ClinicalParams& params);
    
    /** For infants, updates the infantIntervalsAtRisk and potentially
     * infantDeaths arrays. */
    void updateInfantDeaths( SimTime age, InfantMortality& infant_mortality );
    
    /** Special option to allow reports not to be delivered for existing cases
     * (within health-system-memory and not new cases). */
    virtual bool isExistingCase(const ClinicalParams& params) =0;
    
    /// Force all pending summaries to be reported. Should only be called when
    /// class is about to be destroyed anyway to avoid affecting output.
    inline void flushReports (){
        latestReport.flush();
    }

    inline const Episode &getLatestReport () const{
        return latestReport;
    }
    
    inline const Episode::State &getLatestState () const{
        return latestState;
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }

    /** Time of the last treatment (sim::never() if never treated). */
    SimTime m_tLastTreatment = sim::never();
    
public:
    /// Constructor.
    ClinicalModel ();
    
    /** Update for clinical model - new pathogenesis status, treatment, etc.
     *
     * @param withinHostModel WithinHostModel of human.
     * @param hostTransmission per-host transmission data of human.
     * @param ageYears Age of human.
     * @param ageGroup Survey age group of human. */
    virtual void doClinicalUpdate (Human& human, double ageYears, WithinHost::Pathogenesis::StatePair &pg,
                                   const ClinicalParams& params) =0;
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
    /** Last episode; report to survey pending a new episode or human's death. */
    Episode latestReport;
    Episode::State latestState;

    /** @brief Positive values of _doomed variable (codes). */
    enum {
        DOOMED_EXPIRED = -35,   // codes less than or equal to this mean "dead now"
        DOOMED_NEXT_TS = -30,   // will expire on next time step
        DOOMED_START_TIMER = -1,        // set on start of doomed timer
        NOT_DOOMED = 0, // all codes greater than this mean "already dead"; codes less than this mean a count-down to death has started
        DOOMED_TOO_OLD = 1,		///< died because reached age limit
        DOOMED_COMPLICATED = 4,	///< died from severe malaria or malaria with a coinfection
        DOOMED_NEONATAL = 6,	///< died due to mother's malaria infection
        DOOMED_INDIRECT = 7,	///< died indirectly from malaria (after a delay)
    };
    
    /** Can indicate that the individual is dead or about to die.
     *
     * If doomed < 0, the individual is doomed to die.
     * 
     * If doomed > 0, the individual is dead, and will be removed from the
     * population at the beginning of the next time step. NOTE:
     * updateInfantDeaths() counts deaths after the fact, thus cannot remove immediately. */
    int doomed;
};


} }
#endif
