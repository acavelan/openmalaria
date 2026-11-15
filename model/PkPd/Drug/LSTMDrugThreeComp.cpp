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

#include "PkPd/Drug/LSTMDrug.h"
#include "PkPd/Drug/LSTMDrugType.h"
#include "Host/WithinHost/Infection/CommonInfection.h"
#include "util/errors.h"
#include "util/StreamValidator.h"

#include <gsl/gsl_integration.h>
#include <limits>

using namespace std;

namespace OM { namespace PkPd { namespace lstmdrug { namespace threeComp {

/// Parameters for func_fC
struct Params_fC {
    double cA, cB, cC, cABC;    // concentration parameters
    double na, nb, ng, nka;     // decay parameters
    double n;       // slope: unitless
    double V;       // max killing rate: unitless
    double Kn;      // IC50^n: (mg/kg) ^ n
};
/** Function for calculating concentration and then killing function at time t
 * 
 * @param t The variable being integrated over (in this case, time since start
 *      of day or last dose, units days)
 * @param pp Pointer to a Params_fC struct
 * @return killing rate (unitless)
 */
static double func_fC( double t, void* pp ){
    const Params_fC& p = *static_cast<const Params_fC*>( pp );
    
    // exponential decay of drug concentration:
    const double concA = p.cA * exp(p.na * t);
    const double concB = p.cB * exp(p.nb * t);
    const double concC = p.cC * exp(p.ng * t);
    const double concABC = p.cABC * exp(p.nka * t);
    const double conc = concA + concB + concC - concABC;      // mg/l
    
    const double cn = pow(conc, p.n);        // (mg/l) ^ n
    const double fC = p.V * cn / (cn + p.Kn);       // unitless
    return fC;
}
const size_t GSL_INTG_MAX_ITER = 1000;     // 10 seems enough, but no harm in using a higher value
static gsl_integration_workspace *gsl_intgr_wksp = gsl_integration_workspace_alloc (GSL_INTG_MAX_ITER);
//NOTE: we "should" free, but mem-leaks at end of program aren't really important
// gsl_integration_workspace_free (gsl_intgr_wksp);

static double calculateFactor(const Params_fC& p, double duration){
    gsl_function F;
    F.function = &func_fC;
    // gsl_function doesn't accept const; we re-apply const later
    F.params = static_cast<void*>(const_cast<Params_fC*>(&p));
    
    // NOTE: tolerances are arbitrary, but seem to be sufficient
    const double abs_eps = 1e-4, rel_eps = 1e-4;
    // NOTE: 1 through 6 are different algorithms of increasing complexity
    const int qag_rule = 1;     // alg 1 seems to be good enough
    double intfC, err_eps;
    
    int r = gsl_integration_qag (&F, 0.0, duration, abs_eps, rel_eps,
                                 GSL_INTG_MAX_ITER, qag_rule, gsl_intgr_wksp, &intfC, &err_eps);
    if( r != 0 ){
        throw TRACED_EXCEPTION( "calculateFactor: error from gsl_integration_qag",util::Error::GSL );
    }
    if( err_eps > 5e-2 ){
        // This could be a warning, except that warnings tend to be ignored.
        ostringstream msg;
        msg << "calculateFactor: error epsilon is large: "<<err_eps<<" (integral is "<<intfC<<")";
        throw TRACED_EXCEPTION( msg.str(), util::Error::GSL );
    }
    
    return 1.0 / exp( intfC );  // drug factor
}

// Update cached constants. bm: body mass (kg)
static void updateCached(const LSTMDrug& d, double body_mass){
    if( d.last_bm == body_mass ) return;
    
    double k = d.elim_sample * pow(body_mass, d.typeData.neg_m_exponent());
    
    const double a0 = k*d.k21*d.k31;
    const double a1 = k*d.k31 + d.k21*d.k31 + d.k21*d.k13 + k*d.k21 + d.k31*d.k12;
    const double a2 = k + d.k12 + d.k13 + d.k21 + d.k31;
    
    const double third = 1.0 / 3.0;
    const double p = a1 - third * a2 * a2;
    const double q = (2.0 / 27.0) * a2 * a2 * a2 - third * a1 * a2 + a0;
    
    const double at = third * a2;
    const double rt = -third * p;
    const double r2 = 2.0 * sqrt(rt);
    
    const double phi = acos(-q / (rt * r2)) * third;
    const double pi23 = (2.0 / 3.0) * PI;
    
    // negative alpha, beta, gamma:
    d.na = r2 * cos(phi) - at;
    d.nb = r2 * cos(phi + pi23) - at;
    d.ng = r2 * cos(phi + 2.0 * pi23) -at;
    
    // A, B, C from Monolix 1.3.3 (p44):
    const double kaOv = -d.nka / (d.vol_dist * body_mass);
    d.A = kaOv * (d.k21 + d.na) * (d.k31 + d.na) / ((d.na - d.nka) * (d.nb - d.na) * (d.ng - d.na));
    d.B = kaOv * (d.k21 + d.nb) * (d.k31 + d.nb) / ((d.nb - d.nka) * (d.na - d.nb) * (d.ng - d.nb));
    d.C = kaOv * (d.k21 + d.ng) * (d.k31 + d.ng) / ((d.ng - d.nka) * (d.nb - d.ng) * (d.na - d.ng));
    
    d.last_bm = body_mass;
}
double getConcentration(const LSTMDrug& d, size_t index) {
    if( index == d.typeData.getIndex() ) return d.conc();
    else return 0.0;
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double drugFactor(const LSTMDrug& d, LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) {
    if( d.conc() == 0.0 && d.doses.size() == 0 ) return 1.0; // nothing to do
    updateCached(d, body_mass);
    
    Params_fC p;
    p.cA = d.concA;     p.cB = d.concB;   p.cC = d.concC;   p.cABC = d.concABC;
    p.na = d.na;        p.nb = d.nb;      p.ng = d.ng;      p.nka = d.nka;
    const LSTMDrugPD& pd = d.typeData.getPD(inf->genotype());
    p.n = pd.slope();   p.V = pd.max_killing_rate();
    p.Kn = pd.IC50_pow_slope(rng, d.typeData.getIndex(), inf);
    
    double time = 0.0;  // time since start of day
    double totalFactor = 1.0;   // survival factor for whole day
    
    typedef pair<double,double> TimeConc;
    for( const TimeConc& time_conc : d.doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                double duration = time_conc.first - time;
                totalFactor *= calculateFactor(p, duration);
                p.cA *= exp(p.na * duration);
                p.cB *= exp(p.nb * duration);
                p.cC *= exp(p.ng * duration);
                p.cABC *= exp(p.nka * duration);
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add dose:
            p.cA += d.A * time_conc.second;
            p.cB += d.B * time_conc.second;
            p.cC += d.C * time_conc.second;
            p.cABC += (d.A + d.B + d.C) * time_conc.second;
        }else /*i.e. tomorrow or later*/{
            // ignore
        }
    }
    if( time < 1.0 ){
        totalFactor *= calculateFactor(p, 1.0 - time);
    }
    
    return totalFactor;
}

void update(LSTMDrug& d, double body_mass ) {
    if( d.conc() == 0.0 && d.doses.size() == 0 ) return;        // nothing to do
    updateCached(d, body_mass);
    
    // exponential decay of existing quantities:
    //TODO(performance): is it faster to pre-calculate these and either store extra
    // parameters or adapt uses of alpha, beta, gamma, etc. below?
    d.concA *= exp(d.na);
    d.concB *= exp(d.nb);
    d.concC *= exp(d.ng);
    d.concABC *= exp(d.nka);
    
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    for( TimeConc& time_conc : d.doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            // add dose:
            const double qty = time_conc.second;
            d.concA += d.A * qty * exp(d.na * (1.0 - time_conc.first));
            d.concB += d.B * qty * exp(d.nb * (1.0 - time_conc.first));
            d.concC += d.C * qty * exp(d.ng * (1.0 - time_conc.first));
            d.concABC += (d.A + d.B + d.C) * qty * exp(d.nka * (1.0 - time_conc.first));
            doses_taken += 1;
        }else /*i.e. tomorrow or later*/{
            time_conc.first -= 1.0;
        }
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    d.doses.erase(d.doses.begin(), d.doses.begin() + doses_taken);
    
    util::streamValidate( d.conc() );
    if( d.conc() < d.typeData.getNegligibleConcentration() ){
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        d.concA = 0.0;
        d.concB = 0.0;
        d.concC = 0.0;
        d.concABC = 0.0;
    }
}

void checkpoint(LSTMDrug& d, ostream& stream){
    d.concA & stream;
    d.concB & stream;
    d.concC & stream;
    d.concABC & stream;
    d.elim_sample & stream;
    d.k12 & stream;
    d.k21 & stream;
    d.k13 & stream;
    d.k31 & stream;
    d.nka & stream;
    d.last_bm & stream;
    d.na & stream;
    d.nb & stream;
    d.ng & stream;
    d.A & stream;
    d.B & stream;
    d.C & stream;
}
void read(LSTMDrug& d, istream& stream){
    d.concA & stream;
    d.concB & stream;
    d.concC & stream;
    d.concABC & stream;
    d.elim_sample & stream;
    d.k12 & stream;
    d.k21 & stream;
    d.k13 & stream;
    d.k31 & stream;
    d.nka & stream;
    d.last_bm & stream;
    d.na & stream;
    d.nb & stream;
    d.ng & stream;
    d.A & stream;
    d.B & stream;
    d.C & stream;
}

} } } }
