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
#include "util/vectors.h"

#include <gsl/gsl_integration.h>
#include <limits>

using namespace std;

namespace OM { namespace PkPd { namespace lstmdrug { namespace conversion {

/// Parameters for func_convFactor and calculateFactor
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

static double calculateParentQuantity( struct Params_convFactor p, double expAbsorb, double expPLoss ) {
    return p.f * p.qtyG * expAbsorb
        + (p.qtyP - p.f * p.qtyG) * expPLoss;
}

static double calculateParentDrugFactor( struct Params_convFactor p, double expAbsorb, double expPLoss ) {
    const double qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
    const double cP = qtyP * p.invVdP;                  // concentrations; mg/l*/
    const double cnP = pow(cP, p.nP);                   // (mg/l) ^ n
    const double fCP = p.VP * cnP / (cnP + p.KnP);      // unitless
    return fCP;
}

static double calculateMetaboliteQuantity(struct Params_convFactor p, double expAbsorb, double expPLoss, double t) {
    return p.g * p.qtyG * expAbsorb
        + (p.h * p.qtyG - p.i * p.qtyP) * expPLoss
        + (p.j * p.qtyG + p.i * p.qtyP + p.qtyM) * exp(p.nkM * t);
}

static double calculateMetaboliteDrugFactor( struct Params_convFactor p, double expAbsorb, double expPLoss, double t ) {
    const double qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, t);
    const double cM = qtyM * p.invVdM;              // concentrations; mg/l
    const double cnM = pow(cM, p.nM);               // (mg/l) ^ n
    const double fCM = p.VM * cnM / (cnM + p.KnM);  // unitless
    return fCM;
}

/** Function for calculating concentration and then killing function at time t
 * 
 * @param t The variable being integrated over (in this case, time since start
 *      of day or last dose, units days)
 * @param pp Pointer to a Params_convFactor struct
 * @return killing rate (unitless)
 */
static double func_convFactor( double t, void* pp ){
    const Params_convFactor& p = *static_cast<const Params_convFactor*>( pp );
    const double expAbsorb = exp(p.nka * t), expPLoss = exp(p.nl * t);
    const double fCP = calculateParentDrugFactor( p, expAbsorb, expPLoss );
    const double fCM = calculateMetaboliteDrugFactor( p, expAbsorb, expPLoss, t );
    // use the most effective killing factor (from area under the drug kill curve), which is the one with the bigger number
    return max(fCP,fCM);
}

const size_t GSL_INTG_CONV_MAX_ITER = 1000;     // 10 seems enough, but no harm in using a higher value
static gsl_integration_workspace *gsl_intgr_conv_wksp = gsl_integration_workspace_alloc (GSL_INTG_CONV_MAX_ITER);
//NOTE: we "should" free, but mem-leaks at end of program aren't really important
// gsl_integration_workspace_free (gsl_intgr_conv_wksp);

static double calculateFactor(const Params_convFactor& p, double duration){
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
    
    int r = gsl_integration_qag (&F, 0.0, duration, abs_eps, rel_eps,
                                 GSL_INTG_CONV_MAX_ITER, qag_rule, gsl_intgr_conv_wksp, &intfC, &err_eps);
    if( r != 0 ){
        throw TRACED_EXCEPTION( "calculateFactor: error from gsl_integration_qag",util::Error::GSL );
    }
    return exp( -intfC );  // drug factor
}

static void setConversionParameters(const LSTMDrug& d, Params_convFactor& p, double body_mass)
{
    p.qtyG = d.qtyG; p.qtyP = d.qtyP; p.qtyM = d.qtyM;
    
    // decay "constants" (dependent on body mass):
    const double nkP = d.nkP_sample * pow(body_mass, d.typeData.neg_m_exponent());      // -y
    const double nconv = d.nconv_sample * pow(body_mass, d.typeData.neg_m_exponent());  // -z
    p.nka = d.nka;        // -x
    p.nkM = d.nkM_sample * pow(body_mass, d.metaboliteType->neg_m_exponent());  // -k
    p.nl = nkP + nconv;    // -(y + z)
    
    p.f = d.nka / (p.nl - d.nka);                           // -x / (x-y-z) = x / (y+z-x)
    // Let u = parentType.molecular_weight_ratio() = M_C / M_B:
    const double rz = d.typeData.molecular_weight_ratio() * nconv; // -u*z
    p.g = rz * d.nka / ((d.nka - p.nl) * (d.nka - p.nkM));    // u*z*x / ((y+z-x) * (k-x))
    p.h = rz * d.nka / ((d.nka - p.nl) * (p.nkM - p.nl));   // u*z*x / ((y+z-x) * (y+z-k))
    p.i = rz / (p.nl - p.nkM);                          // -u*z / (k - y - z) = u*z / (y+z-k)
    p.j = rz * d.nka / ((p.nkM - p.nl) * (p.nkM - d.nka));  // u*z*x / ((y+z-k) * (x-k))
    
    p.invVdP = 1.0 / (d.vol_dist * body_mass); p.invVdM = 1.0 / (d.vol_dist_metabolite * body_mass);
}

static void setKillingParameters(const LSTMDrug& d, LocalRng& rng, Params_convFactor& p,
        WithinHost::CommonInfection *inf)
{
    uint32_t genotype = inf->genotype();
    const LSTMDrugPD& pdP = d.typeData.getPD(genotype), &pdM =
            d.metaboliteType->getPD(genotype);
    p.nP = pdP.slope();
    p.VP = pdP.max_killing_rate();
    p.nM = pdM.slope();
    p.VM = pdM.max_killing_rate();
    
    // Use custom code here because we need to handle covariance
    auto pIndex = d.typeData.getIndex();
    if (inf->Kn.count(pIndex) != 0) {
        // Read cached values: IC50 ^ n
        p.KnP = inf->Kn.at(pIndex);
        p.KnM = inf->Kn.at(d.metaboliteType->getIndex());
    } else {
        // First usage for this infection / treatment: sample, optionally with correlation.
        auto zscore = NormalSample::generate(rng);
        p.KnP = pdP.IC50_pow_slope(zscore);
        inf->Kn[pIndex] = p.KnP;
        
        auto metab_zscore = d.typeData.IC50_correlated_sample(zscore, rng);
        p.KnM = pdM.IC50_pow_slope(metab_zscore);
        inf->Kn[d.metaboliteType->getIndex()] = p.KnM;
    }
}

double getConcentration(const LSTMDrug& d, size_t index) {
    if( index == d.typeData.getIndex() ){
        return d.qtyP / (d.vol_dist * d.last_bm);
    }else if( index == d.metaboliteType->getIndex() ){
        return d.qtyM / (d.vol_dist_metabolite * d.last_bm);
    }else return 0.0;
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double drugFactor(const LSTMDrug& d, LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) {
    if( d.qtyG == 0.0 && d.qtyP == 0.0 && d.qtyM == 0.0 && d.doses.size() == 0 ){
        return 1.0; // nothing to do
    }
    
    Params_convFactor p;
    setConversionParameters(d, p, body_mass);
    setKillingParameters(d, rng, p, inf);
    
    double time = 0.0;  // time since start of day
    double totalFactor = 1.0;   // survival factor for whole day
    
    typedef pair<double,double> TimeConc;
    for( const TimeConc& time_conc : d.doses ){
        // we iterate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                double duration = time_conc.first - time;
                totalFactor *= calculateFactor(p, duration);
                const double expAbsorb = exp(p.nka * duration), expPLoss = exp(p.nl * duration);
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
        totalFactor *= calculateFactor(p, 1.0 - time);
    }
    
    return totalFactor;
}

void update(LSTMDrug& d, double body_mass ){
    if( d.qtyG == 0.0 && d.qtyP == 0.0 && d.qtyM == 0.0 && d.doses.size() == 0 ){
        return; // nothing to do
    }
    d.last_bm = body_mass;
    
    Params_convFactor p;
    setConversionParameters(d, p, body_mass);
    
    double time = 0.0, duration;
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    for( TimeConc& time_conc : d.doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( (duration = time_conc.first - time) > 0.0 ){
                const double expAbsorb = exp(p.nka * duration), expPLoss = exp(p.nl * duration);
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
        const double expAbsorb = exp(p.nka * duration), expPLoss = exp(p.nl * duration);
        p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
        p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
        p.qtyG *= expAbsorb;
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    d.doses.erase(d.doses.begin(), d.doses.begin() + doses_taken);
    
    d.qtyM = p.qtyM;
    d.qtyP = p.qtyP;
    d.qtyG = p.qtyG;
    
    util::streamValidate( d.qtyM );
    if( d.qtyP < d.typeData.getNegligibleConcentration() * body_mass * d.vol_dist &&
            d.qtyM < d.metaboliteType->getNegligibleConcentration() * body_mass * d.vol_dist_metabolite )
    {
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        d.qtyG = d.qtyP = d.qtyM = 0.0;
    }
}

void checkpoint(LSTMDrug& d, ostream& stream){
    d.qtyG & stream;
    d.qtyP & stream;
    d.qtyM & stream;
    d.nka & stream;
    d.nkP_sample & stream;
    d.nconv_sample & stream;
    d.nkM_sample & stream;
}
void read(LSTMDrug& d, istream& stream){
    d.qtyG & stream;
    d.qtyP & stream;
    d.qtyM & stream;
    d.nka & stream;
    d.nkP_sample & stream;
    d.nconv_sample & stream;
    d.nkM_sample & stream;
}

} } } }
