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

#include "PkPd/LSTMDrug.h"
#include "PkPd/LSTMDrugType.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <limits>
#include <gsl/gsl_integration.h>

using namespace std;

namespace OM {
namespace PkPd {

LSTMDrug::LSTMDrug (Kind kind, const LSTMDrugType& type,
                    const LSTMDrugType* metabolite, LocalRng& rng) :
    kind(kind), typeData(&type), metaboliteData(metabolite)
{
    switch (kind) {
    case Kind::OneComp:
        vol_dist = type.sample_Vd(rng);
        concentration = 0.0;
        neg_elim_sample = -type.sample_elim_rate(rng);
        break;
    case Kind::ThreeComp:
        vol_dist = type.sample_Vd(rng);
        concA = 0.0; concB = 0.0; concC = 0.0; concABC = 0.0;
        elim_sample = type.sample_elim_rate(rng);
        k12 = type.sample_k12(rng);
        k21 = type.sample_k21(rng);
        k13 = type.sample_k13(rng);
        k31 = type.sample_k31(rng);
        nka = -type.sample_ka(rng);
        last_bm = numeric_limits<double>::quiet_NaN();
        na = numeric_limits<double>::quiet_NaN();
        nb = numeric_limits<double>::quiet_NaN();
        ng = numeric_limits<double>::quiet_NaN();
        A = numeric_limits<double>::quiet_NaN();
        B = numeric_limits<double>::quiet_NaN();
        C = numeric_limits<double>::quiet_NaN();
        break;
    case Kind::Conversion:
        vol_dist = type.sample_Vd(rng);
        qtyG = 0.0; qtyP = 0.0; qtyM = 0.0;
        nka = -type.sample_ka(rng);
        nkP_sample = -type.sample_elim_rate(rng);
        nconv_sample = -type.sample_conv_rate(rng);
        nkM_sample = -metaboliteData->sample_elim_rate(rng);
        vol_dist_metabolite = metaboliteData->sample_Vd(rng);
        break;
    }
}

size_t LSTMDrug::getIndex() const {
    return typeData->index;   // Conversion: parent drug's index
}

double LSTMDrug::getConcentration(size_t index) const {
    switch (kind) {
    case Kind::OneComp:
        return index == typeData->index ? concentration : 0.0;
    case Kind::ThreeComp:
        return index == typeData->index ? conc() : 0.0;
    case Kind::Conversion:
        if( index == typeData->index ) return getParentConcentration();
        else if( index == metaboliteData->index ) return getMetaboliteConcentration();
        else return 0.0;
    }
    return 0.0;
}

double LSTMDrug::calculateDrugFactor(LocalRng& rng, WithinHost::Infection *inf, double body_mass) const {
    switch (kind) {
    case Kind::OneComp: return drugFactor_oneComp (rng, inf, body_mass);
    case Kind::ThreeComp: return drugFactor_threeComp (rng, inf, body_mass);
    default: return drugFactor_conversion (rng, inf, body_mass);
    }
}

void LSTMDrug::updateConcentration (double body_mass) {
    switch (kind) {
    case Kind::OneComp: updateConcentration_oneComp (body_mass); break;
    case Kind::ThreeComp: updateConcentration_threeComp (body_mass); break;
    case Kind::Conversion: updateConcentration_conversion (body_mass); break;
    }
}

bool comp(pair<double,double> lhs, pair<double,double> rhs){
    return lhs.first < rhs.first;
}

// Create our list of doses. Optimise for the case where we only have 1 or less
// per day, but be able to handle more.
// If two doses coincide, they can be combined but doing so is not essential.

void LSTMDrug::medicate(double time, double qty){
    // Insert in the right position to maintain sorting:
    auto elt = make_pair (time, qty);
    auto pos = lower_bound(doses.begin(), doses.end(), elt, comp);
    doses.insert(pos, std::move(elt));
    assert(is_sorted(doses.begin(), doses.end(), comp));
}

// ———  OneComp kind (was LSTMDrugOneComp.cpp)  ———


// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrug::drugFactor_oneComp(LocalRng& rng, WithinHost::Infection *inf, double body_mass) const {
    if( concentration == 0.0 && doses.size() == 0 ) return 1.0; // nothing to do
    
    /* Survival factor of the parasite (this multiplies the parasite density).
    Calculated below for each time interval. */
    double totalFactor = 1.0;
    
    // Make a copy of concetration and use that over today. Don't adjust concentration because this
    // function may be called multiple times (or not at all) in a day.
    double concentration_today = concentration; // mg / l
    double neg_elim_rate = neg_elim_sample * pow(body_mass, typeData->neg_m_exp);
    
    const LSTMDrugPD& drugPD = typeData->getPD(inf->genotype());
    const double Kn = drugPD.IC50_pow_slope(rng, typeData->index, inf);
    
    double time = 0.0;
    typedef pair<double,double> TimeConc;
    for( TimeConc time_conc : doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                totalFactor *= drugPD.calcFactor( Kn, neg_elim_rate, &concentration_today, time_conc.first - time );
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add dose (instantaneous absorption):
            concentration_today += time_conc.second / (vol_dist * body_mass);
        }else/*i.e. tomorrow or later*/{
            break;
        }
    }
    if( time < 1.0 ){
        totalFactor *= drugPD.calcFactor( Kn, neg_elim_rate, &concentration_today, 1.0 - time );
    }
    
    return totalFactor; // Drug effect per day per drug per parasite
}

void LSTMDrug::updateConcentration_oneComp( double body_mass ){
    if( concentration == 0.0 && doses.size() == 0 ) return;     // nothing to do
    
    // exponential decay of drug concentration (portion without new doses):
    //TODO: is it faster to pre-calculate this and either store an extra
    // parameter or adapt uses of neg_elim_rate, etc. below?
    double neg_elim_rate = neg_elim_sample * pow(body_mass, typeData->neg_m_exp);
    concentration *= exp(neg_elim_rate);
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    for( TimeConc& time_conc : doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            // calculate decayed dose and add:
            concentration += time_conc.second / (vol_dist * body_mass) * exp(neg_elim_rate * (1.0 - time_conc.first));
            doses_taken += 1;
        }else /*i.e. tomorrow or later*/{
            time_conc.first -= 1.0;
        }
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    doses.erase(doses.begin(), doses.begin() + doses_taken);
    
    util::streamValidate( concentration );
    if( concentration < typeData->getNegligibleConcentration() ){
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        concentration = 0.0;
    }
}

// ———  ThreeComp kind (was LSTMDrugThreeComp.cpp)  ———


void LSTMDrug::updateCached(double body_mass) const{
    if( last_bm == body_mass ) return;
    
    double k = elim_sample * pow(body_mass, typeData->neg_m_exp);
    
    const double a0 = k*k21*k31;
    const double a1 = k*k31 + k21*k31 + k21*k13 + k*k21 + k31*k12;
    const double a2 = k + k12 + k13 + k21 + k31;
    
    const double third = 1.0 / 3.0;
    const double p = a1 - third * a2 * a2;
    const double q = (2.0 / 27.0) * a2 * a2 * a2 - third * a1 * a2 + a0;
    
    const double at = third * a2;
    const double rt = -third * p;
    const double r2 = 2.0 * sqrt(rt);
    
    const double phi = acos(-q / (rt * r2)) * third;
    const double pi23 = (2.0 / 3.0) * 3.14159265358979323846;
    
    // negative alpha, beta, gamma:
    na = r2 * cos(phi) - at;
    nb = r2 * cos(phi + pi23) - at;
    ng = r2 * cos(phi + 2.0 * pi23) -at;
    
    // A, B, C from Monolix 1.3.3 (p44):
    const double kaOv = -nka / (vol_dist * body_mass);
    A = kaOv * (k21 + na) * (k31 + na) / ((na - nka) * (nb - na) * (ng - na));
    B = kaOv * (k21 + nb) * (k31 + nb) / ((nb - nka) * (na - nb) * (ng - nb));
    C = kaOv * (k21 + ng) * (k31 + ng) / ((ng - nka) * (nb - ng) * (na - ng));
    
    last_bm = body_mass;
}

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
double func_fC( double t, void* pp ){
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
gsl_integration_workspace *gsl_intgr_wksp = gsl_integration_workspace_alloc (GSL_INTG_MAX_ITER);
//NOTE: we "should" free, but mem-leaks at end of program aren't really important
// gsl_integration_workspace_free (gsl_intgr_wksp);

double LSTMDrug::calculateFactor_fC(const Params_fC& p, double duration) const{
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

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrug::drugFactor_threeComp(LocalRng& rng, WithinHost::Infection *inf, double body_mass) const {
    if( conc() == 0.0 && doses.size() == 0 ) return 1.0; // nothing to do
    updateCached(body_mass);
    
    Params_fC p;
    p.cA = concA;       p.cB = concB;   p.cC = concC;   p.cABC = concABC;
    p.na = na;  p.nb = nb;      p.ng = ng;      p.nka = nka;
    const LSTMDrugPD& pd = typeData->getPD(inf->genotype());
    p.n = pd.slope();   p.V = pd.max_killing_rate();
    p.Kn = pd.IC50_pow_slope(rng, typeData->index, inf);
    
    double time = 0.0;  // time since start of day
    double totalFactor = 1.0;   // survival factor for whole day
    
    typedef pair<double,double> TimeConc;
    for( const TimeConc& time_conc : doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                double duration = time_conc.first - time;
                totalFactor *= calculateFactor_fC(p, duration);
                p.cA *= exp(p.na * duration);
                p.cB *= exp(p.nb * duration);
                p.cC *= exp(p.ng * duration);
                p.cABC *= exp(p.nka * duration);
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add dose:
            p.cA += A * time_conc.second;
            p.cB += B * time_conc.second;
            p.cC += C * time_conc.second;
            p.cABC += (A + B + C) * time_conc.second;
        }else /*i.e. tomorrow or later*/{
            // ignore
        }
    }
    if( time < 1.0 ){
        totalFactor *= calculateFactor_fC(p, 1.0 - time);
    }
    
    return totalFactor;
}

void LSTMDrug::updateConcentration_threeComp (double body_mass) {
    if( conc() == 0.0 && doses.size() == 0 ) return;     // nothing to do
    updateCached(body_mass);
    
    // exponential decay of existing quantities:
    //TODO(performance): is it faster to pre-calculate these and either store extra
    // parameters or adapt uses of alpha, beta, gamma, etc. below?
    concA *= exp(na);
    concB *= exp(nb);
    concC *= exp(ng);
    concABC *= exp(nka);
    
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    for( TimeConc& time_conc : doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            // add dose:
            const double qty = time_conc.second;
            concA += A * qty * exp(na * (1.0 - time_conc.first));
            concB += B * qty * exp(nb * (1.0 - time_conc.first));
            concC += C * qty * exp(ng * (1.0 - time_conc.first));
            concABC += (A + B + C) * qty * exp(nka * (1.0 - time_conc.first));
            doses_taken += 1;
        }else /*i.e. tomorrow or later*/{
            time_conc.first -= 1.0;
        }
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    doses.erase(doses.begin(), doses.begin() + doses_taken);
    
    util::streamValidate( conc() );
    if( conc() < typeData->getNegligibleConcentration() ){
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        concA = 0.0;
        concB = 0.0;
        concC = 0.0;
        concABC = 0.0;
    }
}

// ———  Conversion kind (parent + metabolite) (was LSTMDrugConversion.cpp)  ———


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
    auto pIndex = typeData->index;
    if (inf->Kn().count(pIndex) != 0) {
        // Read cached values: IC50 ^ n
        p.KnP = inf->Kn().at(pIndex);
        p.KnM = inf->Kn().at(metaboliteData->index);
    } else {
        // First usage for this infection / treatment: sample, optionally with correlation.
        auto zscore = NormalSample::generate(rng);
        p.KnP = pdP.IC50_pow_slope(zscore);
        inf->Kn()[pIndex] = p.KnP;
        
        auto metab_zscore = typeData->IC50_correlated_sample(zscore, rng);
        p.KnM = pdM.IC50_pow_slope(metab_zscore);
        inf->Kn()[metaboliteData->index] = p.KnM;
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

}
}
