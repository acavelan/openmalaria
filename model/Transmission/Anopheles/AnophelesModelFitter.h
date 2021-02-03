/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_AnophelesModelFitter
#define Hmod_AnophelesModelFitter

#include "Global.h"
#include "Transmission/Anopheles/PerHostAnoph.h"
#include "Transmission/PerHost.h"
#include "util/SimpleDecayingValue.h"
#include "util/vectors.h"
#include "util/vecDay.h"

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

#include "rotate.h"

#include <vector>
#include <limits>

#include "CalcNvO.h"

namespace OM {
namespace Transmission {
namespace Anopheles {
    using std::numeric_limits;
    using util::vector2D;
    using namespace OM::util;

class AnophelesModelFitter
{
public:
    AnophelesModelFitter(AnophelesModel &m) : scaleFactor(1.0), rotated(false), scaled(false)
    {
        // usually around 20 days; no real analysis for effect of changing EIPDuration or mosqRestDuration
        shiftAngle = m.EIRRotateAngle - (m.mosq.EIPDuration.inDays() + 10) / 365. * 2. *M_PI; 
    }

    
    // bool fit(AnophelesModel &m, vector<double> &laggedKappa, vector<double> &initialisationEIR)
    // {
    //     vecDay<double> avgAnnualS_v(SimTime::oneYear(), 0.0);
    //     for (SimTime i = SimTime::fromYearsI(4); i < SimTime::fromYearsI(5); i += SimTime::oneDay())
    //     {
    //         avgAnnualS_v[mod_nn(i, SimTime::oneYear())] = m.quinquennialS_v[i];
    //     }

    //     double factor = vectors::sum(m.forcedS_v) / vectors::sum(avgAnnualS_v);

    //     // cout << "check: " << vectors::sum(forcedS_v) << " " << vectors::sum(avgAnnualS_v) << endl;
    //     // cout << "Pre-calced Sv, dynamic Sv:\t"<<sumAnnualForcedS_v<<'\t'<<vectors::sum(annualS_v)<<endl;
    //     if (!(factor > 1e-6 && factor < 1e6))
    //     {
    //         if (factor > 1e6 && vectors::sum(m.quinquennialS_v) < 1e-3)
    //         {
    //             throw util::base_exception("Simulated S_v is approx 0 (i.e.\
    //  mosquitoes are not infectious, before interventions). Simulator cannot handle this; perhaps\
    //  increase EIR or change the entomology model.",
    //                                        util::Error::VectorFitting);
    //         }
    //         if (vectors::sum(m.forcedS_v) == 0.0)
    //         {
    //             return false; // no EIR desired: nothing to do
    //         }
    //         cerr << "Input S_v for this vector:\t" << vectors::sum(m.forcedS_v) << endl;
    //         cerr << "Simulated S_v:\t\t\t" << vectors::sum(m.quinquennialS_v) / 5.0 << endl;
    //         throw TRACED_EXCEPTION("factor out of bounds", util::Error::VectorFitting);
    //     }

    //     const double LIMIT = 0.1;

    //     if (fabs(factor - 1.0) > LIMIT)
    //     {
    //         scaled = false;
    //         double factorDiff = (scaleFactor * factor - scaleFactor) * 1.0;
    //         scaleFactor += factorDiff;
    //     }
    //     else
    //         scaled = true;

    //     double rAngle = findAngle(m.EIRRotateAngle, m.FSCoeffic, avgAnnualS_v);
    //     shiftAngle += rAngle;
    //     rotated = true;

    //     // cout << "EIRRotateAngle: " << EIRRotateAngle << " rAngle = " << rAngle << ", angle = " << shiftAngle << " scalefactor: " <<
    //     // scaleFactor << " , factor: " << factor << endl;

    //     // Compute forced_sv from the Fourrier Coeffs
    //     // shiftAngle rotate the vector to correct the offset between simulated and input EIR
    //     // shiftAngle is the offset between the
    //     vectors::expIDFT(m.mosqEmergeRate, m.FSCoeffic, -shiftAngle);
    //     // Scale the vector according to initNv0FromSv to get the mosqEmergerate
    //     // scaleFactor scales the vector to correct the ratio between simulated and input EIR
    //     vectors::scale(m.mosqEmergeRate, scaleFactor * m.initNv0FromSv);

    //     // initNvFromSv *= scaleFactor;     //(not currently used)

    //     // What factor exactly these should be scaled by isn't obvious; in any case
    //     // they should reach stable values quickly.
    //     m.scale(factor);

    //     m.initIterate();

    //     return !(scaled && rotated);
    // }

    bool fit(AnophelesModel &m, vector<double> &laggedKappa, vector<double> &initialisationEIR)
    {
        double *mosqEmergerateVector = m.mosqEmergeRate.internal().data();
        int daysInYear = m.mosqEmergeRate.internal().size(); // 365?
        int mosqRestDuration = m.mosq.restDuration.inDays();
        int EIPDuration = m.mosq.EIPDuration.inDays();
        int nHostTypesInit = 1; // + types of non human hosts
        int nMalHostTypesInit = 1; // ??
        double popSizeInit = m.populationSize; // ??
        double hostAvailabilityRateInit = PerHostAnophParams::get(m.species).entoAvailability.mean(); //meanAvail; // ?
        //double mosqSeekingDuration = mosqSeekingDuration;
        // double mosqProbBiting = mosq.getMosqProbBiting().getMean(); //PerHostAnophParams::get(species).ProbMosqbBiting;
        // double mosqProbFindRestSite = mosq.getMosqProbFindRestSite().getMean(); //PerHostAnophParams::get(species).ProbMosqbFindRestSite;
        // double mosqProbResting = mosq.getMosqProbResting().getMean(); //PerHostAnophParams::get(species).ProbMosqbResting;
        // double mosqProbOvipositing = mosq.getMosqProbOvipositing().getValue(); //PerHostAnophParams::get(species).ProbMosqbOvipositing;

        vector<double> Kvi(daysInYear, laggedKappa[0]);
        double *FHumanInfectivityInitVector = Kvi.data();
        double *FEIRInitVector = initialisationEIR.data();// m.ini m.getEIR().internal().data();

        // Copy of mosqEmergeRate
        std::vector<double> & Nv0guessRef = m.mosqEmergeRate.internal();
        std::vector<double> NvOguess(Nv0guessRef.begin(), Nv0guessRef.end());
        double *FMosqEmergeRateInitEstimateVector = NvOguess.data();

        cout << "daysInYear: " << daysInYear << endl;
        cout << "mosqRestDuration: " << mosqRestDuration << endl;
        cout << "EIPDuration: " << EIPDuration << endl;
        cout << "nHostTypesInit: " << nHostTypesInit << endl;
        cout << "nMalHostTypesInit: " << nMalHostTypesInit << endl;
        cout << "popSizeInit: " << popSizeInit << endl;
        cout << "hostAvailabilityRateInit: " << hostAvailabilityRateInit << endl;
        cout << "mosqSeekingDuration: " << m.mosq.seekingDuration << endl;
        cout << "mosqProbBiting: " << m.mosq.probBiting << endl;
        cout << "mosqProbFindRestSite: " << m.mosq.probFindRestSite << endl;
        cout << "mosqProbResting: " << m.mosq.probResting << endl;
        cout << "mosqProbOvipositing: " << m.mosq.probOvipositing << endl;

        cout << "FHumanInfectivityInitVector: " << Kvi.size() << endl;
        cout << "FEIRInitVector: " << initialisationEIR.size() << endl;
        cout << "FMosqEmergeRateInitEstimateVector: " << NvOguess.size() << endl;
        cout << "initNv0FromSv: " << m.initNv0FromSv << endl;

        vector<double> EIR;
        EIR.resize(365);

        for(unsigned int i=0; i<initialisationEIR.size(); i++)
        {
            EIR[i*5] = initialisationEIR[i];
            EIR[i*5+1] = initialisationEIR[i];
            EIR[i*5+2] = initialisationEIR[i];
            EIR[i*5+3] = initialisationEIR[i];
            EIR[i*5+4] = initialisationEIR[i];
        }

        // for(const auto &i : initialisationEIR)
        //     cout << i << endl;

        CalcInitMosqEmergeRate(mosqEmergerateVector, &daysInYear,
            &mosqRestDuration, &EIPDuration, &nHostTypesInit,
            &nMalHostTypesInit, &popSizeInit, 
            &hostAvailabilityRateInit, &m.mosq.seekingDeathRate,
            &m.mosq.seekingDuration, &m.mosq.probBiting,
            &m.mosq.probFindRestSite, &m.mosq.probResting,
            &m.mosq.probOvipositing, FHumanInfectivityInitVector,
            EIR.data(), FMosqEmergeRateInitEstimateVector);//, m.initNv0FromSv);

        return false;
    }

private:
    double scaleFactor, shiftAngle;
    bool rotated, scaled;
};

}
}
}
#endif
