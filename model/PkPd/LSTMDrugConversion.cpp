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

// LSTMDrug Conversion kind (was LSTMDrugConversion; hierarchy collapsed phase 4b):
// two 1-compartment models with conversion (parent typeData, metaboliteData).
#include "PkPd/LSTMDrug.h"
#include "util/checkpoint.h"
#include "PkPd/LSTMDrugType.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <gsl/gsl_integration.h>
#include <limits>

using namespace std;

namespace OM { namespace PkPd {

double LSTMDrug::getParentConcentration() const {
    return qtyP / (vol_dist * last_bm);
}
double LSTMDrug::getMetaboliteConcentration() const {
    return qtyM / (vol_dist_metabolite * last_bm);
}

/// Parameters for func_convFactor and LSTMDrug::calculateFactor_conv
struct Params_convFactor {
    // Quantities of parent in gut, parent in circulation and metabolite in
    // circulation, in mg (A, B, C in paper):
    double qtyG, qtyP, qtyM;
    // nka = -x (negative of absorption rate)
    // nkM = -k
    // nl = -(y + z)
    double nka, nkM, nl;
    // terms involving x, y, z and the molecular weight ratio:
    double f, g, h, i, j;
    
    double invVdP, invVdM;     // 1.0 / (Vd * body mass) for parent and metabolite; units: 1/l
    double nP, nM;       // slope: unitless
    double VP, VM;       // max killing rate: unitless
    double KnP, KnM;      // IC50^n: (mg/kg) ^ n
};

double calculateParentQuantity( struct Params_convFactor p, double expAbsorb, double expPLoss ) {
    return p.f * p.qtyG * expAbsorb
        + (p.qtyP - p.f * p.qtyG) * expPLoss;
}

double calculateParentDrugFactor( struct Params_convFactor p, double expAbsorb, double expPLoss ) {
    const double qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
    const double cP = qtyP * p.invVdP;                  // concentrations; mg/l*/
    const double cnP = pow(cP, p.nP);                   // (mg/l) ^ n
    const double fCP = p.VP * cnP / (cnP + p.KnP);      // unitless
    return fCP;
}

double calculateMetaboliteQuantity(struct Params_convFactor p, double expAbsorb, double expPLoss, double t) {
    return p.g * p.qtyG * expAbsorb
        + (p.h * p.qtyG - p.i * p.qtyP) * expPLoss
        + (p.j * p.qtyG + p.i * p.qtyP + p.qtyM) * exp(p.nkM * t);
}

double calculateMetaboliteDrugFactor( struct Params_convFactor p, double expAbsorb, double expPLoss, double t ) {
    const double qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, t);
    const double cM = qtyM * p.invVdM;              // concentrations; mg/l
    const double cnM = pow(cM, p.nM);               // (mg/l) ^ n
    const double fCM = p.VM * cnM / (cnM + p.KnM);  // unitless
    return fCM;
}

// size_t intg_steps = 0;
/** Function for calculating concentration and then killing function at time t
 * 
 * @param t The variable being integrated over (in this case, time since start
 *      of day or last dose, units days)
 * @param pp Pointer to a Params_convFactor struct
 * @return killing rate (unitless)
 */
double func_convFactor( double t, void* pp ){
//     intg_steps += 1;
    const Params_convFactor& p = *static_cast<const Params_convFactor*>( pp );
    const double expAbsorb = exp(p.nka * t), expPLoss = exp(p.nl * t);
    const double fCP = calculateParentDrugFactor( p, expAbsorb, expPLoss );
    const double fCM = calculateMetaboliteDrugFactor( p, expAbsorb, expPLoss, t );
    // use the most effective killing factor (from area under the drug kill curve), which is the one with the bigger number
    return max(fCP,fCM);
}

const size_t GSL_INTG_CONV_MAX_ITER = 1000;     // 10 seems enough, but no harm in using a higher value
gsl_integration_workspace *gsl_intgr_conv_wksp = gsl_integration_workspace_alloc (GSL_INTG_CONV_MAX_ITER);
//NOTE: we "should" free, but mem-leaks at end of program aren't really important
// gsl_integration_workspace_free (gsl_intgr_conv_wksp);
double LSTMDrug::calculateFactor_conv(const Params_convFactor& p, double duration) const{
    gsl_function F;
    F.function = &func_convFactor;
    // gsl_function doesn't accept const; we re-apply const later
    F.params = static_cast<void*>(const_cast<Params_convFactor*>(&p));
    
    // We use exp(-result), so small absolute differences can matter (but also
    // using smaller abs_eps is cheap). We likely don't need high rel precision.
    const double abs_eps = 1e-5, rel_eps = 1e-2;
    // NOTE: 1 through 6 are different algorithms of increasing complexity
    const int qag_rule = 1;     // alg 1 seems to be good enough
    double intfC, err_eps;      // intfC will carry our result; err_eps is a measure of accuracy of the result
    
//     intg_steps = 0;
    int r = gsl_integration_qag (&F, 0.0, duration, abs_eps, rel_eps,
                                 GSL_INTG_CONV_MAX_ITER, qag_rule, gsl_intgr_conv_wksp, &intfC, &err_eps);
    if( r != 0 ){
        throw TRACED_EXCEPTION( "calculateFactor: error from gsl_integration_qag",util::Error::GSL );
    }
    // Testing err_eps is redundant with GSL's built-in tests
//     cout << "integration steps: " << intg_steps << endl;
//     cout << "duration: " << duration << ", AUC: " << intfC << endl;
    return exp( -intfC );  // drug factor
}

void LSTMDrug::setConversionParameters(Params_convFactor& p,
        double body_mass) const
{
    p.qtyG = qtyG; p.qtyP = qtyP; p.qtyM = qtyM;
    
    // decay "constants" (dependent on body mass):
    const double nkP = nkP_sample * pow(body_mass, typeData->neg_m_exp);      // -y
    const double nconv = nconv_sample * pow(body_mass, typeData->neg_m_exp);  // -z
    p.nka = nka;        // -x
    p.nkM = nkM_sample * pow(body_mass, metaboliteData->neg_m_exp);  // -k
    p.nl = nkP + nconv;    // -(y + z)
    
    p.f = nka / (p.nl - nka);                           // -x / (x-y-z) = x / (y+z-x)
    // Let u = typeData->mwr = M_C / M_B:
    const double rz = typeData->mwr * nconv; // -u*z
    p.g = rz * nka / ((nka - p.nl) * (nka - p.nkM));    // u*z*x / ((y+z-x) * (k-x))
    p.h = rz * nka / ((nka - p.nl) * (p.nkM - p.nl));   // u*z*x / ((y+z-x) * (y+z-k))
    p.i = rz / (p.nl - p.nkM);                          // -u*z / (k - y - z) = u*z / (y+z-k)
    p.j = rz * nka / ((p.nkM - p.nl) * (p.nkM - nka));  // u*z*x / ((y+z-k) * (x-k))
    
    p.invVdP = 1.0 / (vol_dist * body_mass); p.invVdM = 1.0 / (vol_dist_metabolite * body_mass);
}

void LSTMDrug::setKillingParameters(LocalRng& rng, Params_convFactor& p,
        WithinHost::Infection *inf) const
{
    uint32_t genotype = inf->genotype();
    const LSTMDrugPD& pdP = typeData->getPD(genotype), &pdM =
            metaboliteData->getPD(genotype);
    p.nP = pdP.slope();
    p.VP = pdP.max_killing_rate();
    p.nM = pdM.slope();
    p.VM = pdM.max_killing_rate();
    
    // Use custom code here because we need to handle covariance
    auto pIndex = typeData->getIndex();
    if (inf->Kn().count(pIndex) != 0) {
        // Read cached values: IC50 ^ n
        p.KnP = inf->Kn().at(pIndex);
        p.KnM = inf->Kn().at(metaboliteData->getIndex());
    } else {
        // First usage for this infection / treatment: sample, optionally with correlation.
        auto zscore = NormalSample::generate(rng);
        p.KnP = pdP.IC50_pow_slope(zscore);
        inf->Kn()[pIndex] = p.KnP;
        
        auto metab_zscore = typeData->IC50_correlated_sample(zscore, rng);
        p.KnM = pdM.IC50_pow_slope(metab_zscore);
        inf->Kn()[metaboliteData->getIndex()] = p.KnM;
    }
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrug::drugFactor_conversion(LocalRng& rng, WithinHost::Infection *inf, double body_mass) const {
    if( qtyG == 0.0 && qtyP == 0.0 && qtyM == 0.0 && doses.size() == 0 ){
        return 1.0; // nothing to do
    }
    
    Params_convFactor p;
    setConversionParameters(p, body_mass);
    setKillingParameters(rng, p, inf);
    
    double time = 0.0;  // time since start of day
    double totalFactor = 1.0;   // survival factor for whole day
    
    typedef pair<double,double> TimeConc;
    for( const TimeConc& time_conc : doses ){
        // we iterate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                double duration = time_conc.first - time;
                totalFactor *= calculateFactor_conv(p, duration);
                const double expAbsorb = exp(nka * duration), expPLoss = exp(p.nl * duration);
                p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
                p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
                p.qtyG *= expAbsorb;
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add to quantity of drug in gut:
            p.qtyG += time_conc.second;   // units: mg
        }else /*i.e. tomorrow or later*/{
            // ignore
        }
    }
    if( time < 1.0 ){
        totalFactor *= calculateFactor_conv(p, 1.0 - time);
    }
    
    return totalFactor;
}

void LSTMDrug::updateConcentration_conversion( double body_mass ){
    if( qtyG == 0.0 && qtyP == 0.0 && qtyM == 0.0 && doses.size() == 0 ){
        return; // nothing to do
    }
    last_bm = body_mass;
    
    Params_convFactor p;
    setConversionParameters(p, body_mass);
    
    double time = 0.0, duration;
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    for( TimeConc& time_conc : doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( (duration = time_conc.first - time) > 0.0 ){
                const double expAbsorb = exp(nka * duration), expPLoss = exp(p.nl * duration);
                p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
                p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
                p.qtyG *= expAbsorb;
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add to quantity of drug in gut:
            p.qtyG += time_conc.second;   // units: mg
            doses_taken += 1;
        }else /*i.e. tomorrow or later*/{
            time_conc.first -= 1.0;
        }
    }
    if( time < 1.0 ){
        duration = 1.0 - time;
        const double expAbsorb = exp(nka * duration), expPLoss = exp(p.nl * duration);
        p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
        p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
        p.qtyG *= expAbsorb;
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    doses.erase(doses.begin(), doses.begin() + doses_taken);
    
    qtyM = p.qtyM;
    qtyP = p.qtyP;
    qtyG = p.qtyG;
    
    util::streamValidate( qtyM );
    if( qtyP < typeData->getNegligibleConcentration() * body_mass * vol_dist &&
            qtyM < metaboliteData->getNegligibleConcentration() * body_mass * vol_dist_metabolite )
    {
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        qtyG = qtyP = qtyM = 0.0;
    }
}

template<class S> void LSTMDrug::checkpoint_conversion(S& stream){
    qtyG & stream;
    qtyP & stream;
    qtyM & stream;
    nka & stream;
    nkP_sample & stream;
    nconv_sample & stream;
    nkM_sample & stream;
}
template void LSTMDrug::checkpoint_conversion(std::ostream&);
template void LSTMDrug::checkpoint_conversion(std::istream&);

}
}
