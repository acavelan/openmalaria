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

// model/xml/ conversion layer — scenario core: within-host model selection,
// host weight/heterogeneity, parasite genetics and diagnostics from scnXml.
// Moved verbatim from Host/WithinHost/WHInterface.cpp, WHFalciparum.cpp,
// CommonWithinHost.cpp, Genotypes.cpp and Diagnostic.cpp. State written by
// these functions stays defined in the original files (extern re-declared
// below). Declarations for the moved entry points stay in the original
// headers (scnXml forward-declared).

#include "Global.h"
#include "Parameters.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "Host/WithinHost/DescriptiveWithinHost.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Infection/DummyInfection.h"
#include "Host/WithinHost/Infection/EmpiricalInfection.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "Host/WithinHost/Infection/MolineauxInfection.h"
#include "Host/WithinHost/Infection/PennyInfection.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "Host/WithinHost/WHFalciparum.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "mon/Monitoring.h"
#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "schema/scenario.h"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace OM {
namespace WithinHost {
using namespace OM::util;

// State defined in WHInterface.cpp (external linkage):
extern bool reportInfectionsByGenotype;
extern bool opt_vivax_simple, opt_dummy_whm, opt_empirical_whm,
        opt_molineaux_whm, opt_penny_whm, opt_common_whm;

// State defined in WHFalciparum.cpp (statics made external for this layer):
extern double sigma_i, immPenalty_22, asexImmRemain, immEffectorRemain,
        invCumulativeYstar, invCumulativeHstar, alpha_m, decayM;

// State defined in CommonWithinHost.cpp (external linkage):
extern double hetMassMultStdDev, minHetMassMult;
extern util::AgeGroupInterpolator massByAge;
extern bool reportInfectedOrPatentInfected;

// State defined in Diagnostic.cpp (external linkage):
extern map<string, Diagnostic> diagnostic_set;

void WHInterface::init( const OM::Parameters& parameters, const scnXml::Scenario& scenario ) {
    reportInfectionsByGenotype = mon::isUsed(mon::measure("nInfectByGenotype")) ||
        mon::isUsed(mon::measure("nPatentByGenotype")) ||
        mon::isUsed(mon::measure("logDensByGenotype"));
    
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        opt_vivax_simple = true;
        WHVivax::init( parameters, scenario.getModel() );
    }else{
        WHFalciparum::init( parameters, scenario.getModel() );
        
        if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL)) {
            opt_dummy_whm = true;
            DummyInfection::init ();
        } else if (util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
            opt_empirical_whm = true;
            EmpiricalInfection::init();    // 1-day time step check
        } else if (util::ModelOptions::option (util::MOLINEAUX_WITHIN_HOST_MODEL)) {
            opt_molineaux_whm = true;
            MolineauxInfection::init( parameters );
        } else if (util::ModelOptions::option (util::PENNY_WITHIN_HOST_MODEL)) {
            opt_penny_whm = true;
            PennyInfection::init();
        } else {
            DescriptiveInfection::init( parameters );      // 5-day time step check
        }
        opt_common_whm = opt_dummy_whm || opt_empirical_whm
                || opt_molineaux_whm || opt_penny_whm;
        
        if( opt_common_whm ){
            CommonWithinHost::init( scenario );
        }else{
            DescriptiveWithinHostModel::initDescriptive();
        }
    }
}

void WHFalciparum::init( const OM::Parameters& parameters, const scnXml::Model& model ) {
    sigma_i=sqrt(parameters[Parameter::SIGMA_I_SQ]);
    immPenalty_22=1-exp(parameters[Parameter::IMMUNITY_PENALTY]);
    immEffectorRemain=exp(-parameters[Parameter::IMMUNE_EFFECTOR_DECAY]);
    asexImmRemain=exp(-parameters[Parameter::ASEXUAL_IMMUNITY_DECAY]);
    
    // calculate inverses here, so we can use multiplication later (faster):
    invCumulativeYstar = 1.0 / parameters[Parameter::CUMULATIVE_Y_STAR];
    invCumulativeHstar = 1.0 / parameters[Parameter::CUMULATIVE_H_STAR];
    alpha_m = 1.0 - exp(-parameters[Parameter::NEG_LOG_ONE_MINUS_ALPHA_M]);
    decayM = parameters[Parameter::DECAY_M];
    
    y_lag_len = sim::inSteps(sim::fromDays(20)) + 1;
    
    //NOTE: should also call cleanup() on the PathogenesisModel, but it only frees memory which the OS does anyway
    Pathogenesis::PathogenesisModel::init( parameters, model.getClinical(), false );
    
    // Default value to use if no <parameters> element is present in XML.
    SimTime latentP = UnitParse::readShortDuration("15d", UnitParse::DAYS);
    try
    {
        if (model.getParameters().present())
        {
            latentP = UnitParse::readShortDuration(
                model.getParameters().get().getLatentp(), UnitParse::NONE);
        }
    }
    catch( const util::format_error& e )
    {
            throw util::xml_scenario_error( string("model/parameters/latentP: ").append(e.message()) );
    }

    Infection::init( latentP );
}

void CommonWithinHost::init( const scnXml::Scenario& scenario ){
    const scnXml::Human human = scenario.getModel().getHuman();
    if( !human.getWeight().present() ){
        // Technically this is needed by the PK/PD and Molineaux models
        throw util::xml_scenario_error( "model->human->weight element required by certain models" );
    }
    massByAge.set( human.getWeight().get(), "weight" );
    hetMassMultStdDev = human.getWeight().get().getMultStdDev();
    // hetWeightMult must be large enough that birth weight is at least 0.5 kg:
    minHetMassMult = 0.5 / massByAge.eval( 0.0 );
    
    reportInfectedOrPatentInfected = mon::isUsed(mon::measure("totalInfs")) ||
        mon::isUsed(mon::measure("totalPatentInf"));
    
    PkPd::LSTMModel::init( scenario );
}

// ———  parasite genetics  ———

namespace GT /*for genotype impl details*/{
// State defined in Genotypes.cpp (external linkage):
extern map<double,uint32_t> cum_initial_freqs;
extern map<string, map<string, uint32_t> > alleleCodes;
extern uint32_t nextAlleleCode;
extern vector<Genotypes::Genotype> genotypes;

// NOTE: duplicated from Genotypes.cpp (identical definition; ODR-safe) —
// needed by Genotypes::init below.
enum SampleMode{
    SAMPLE_FIRST,      // always choose first genotype (essentially the off switch)
    SAMPLE_INITIAL,    // sample from initial probabilities
    SAMPLE_TRACKING    // sample from tracked success at genotype level (no recombination)
};
extern SampleMode current_mode, interv_mode;

// Represent a set of loci: all possible combinations of alleles
// This is actually just machinery to calculate the list of all genotypes
struct LocusSet{
    LocusSet( const scnXml::ParasiteLocus& elt_l ){
        alleles.reserve( elt_l.getAllele().size() );
        double cum_p = 0.0;
        for( const scnXml::ParasiteAllele& elt_a : elt_l.getAllele() ){
            uint32_t alleleCode = (nextAlleleCode++);
            alleleCodes[elt_l.getName()][elt_a.getName()] = alleleCode;
            cum_p += elt_a.getInitialFrequency();
            alleles.push_back( Genotypes::Genotype( alleleCode, elt_a.getInitialFrequency(), elt_a.getFitness(), elt_a.getHrp2_deletion() ) );
        }
        if( cum_p < 0.999 || cum_p > 1.001 )
            throw util::xml_scenario_error("expected sum of initial probabilities of alleles to be 1, but for the " 
                + to_string(alleles.size()) + " alleles under locus " + string(elt_l.getName()) + " this is " + to_string(cum_p));

        alleles[0].init_freq += 1.0 - cum_p;     // account for any small errors by adjusting the first frequency
    }
    void include( LocusSet& that ){
        vector<Genotypes::Genotype> newAlleles;
        newAlleles.reserve( alleles.size() * that.alleles.size() );
        for( size_t i = 0; i < alleles.size(); i += 1 ){
            for( size_t j = 0; j < that.alleles.size(); j += 1 ){
                // note that we generate new codes each time; we just waste the old ones
                newAlleles.push_back( alleles[i].cross(that.alleles[j]) );
            }
        }
        alleles.swap(newAlleles);
    }
    
    vector<Genotypes::Genotype> alleles;
};

}

// utility function: does a vector contain an element?
template<class T>
bool contains(const vector<T>& vec, const T& x){
    for( auto it = vec.begin(), end = vec.end(); it != end; ++it ){
        if( *it == x ) return true;
    }
    return false;
}

void Genotypes::init( const scnXml::Scenario& scenario ){
    if( scenario.getParasiteGenetics().present() ){
        const scnXml::ParasiteGenetics& genetics =
            scenario.getParasiteGenetics().get();
        
        GT::current_mode = GT::SAMPLE_INITIAL;      // turn on sampling
        if( genetics.getSamplingMode() == "initial" ){
            GT::interv_mode = GT::SAMPLE_INITIAL;
        }else if( genetics.getSamplingMode() == "tracking" ){
            GT::interv_mode = GT::SAMPLE_TRACKING;
            if( !scenario.getEntomology().getVector().present() ){
                throw util::xml_scenario_error( "incompatibility; either use "
                    "entomology/vector (not nonVector) or set "
                    "parasiteGenetics/samplingMode to \"initial\" (not \"tracking\")" );
            }
            if( scenario.getEntomology().getMode() != "dynamic" ){
                throw util::xml_scenario_error( "incompatibility; either set "
                    "entomology/mode to \"dynamic\" (not \"forced\") or set "
                    "parasiteGenetics/samplingMode to \"initial\" (not \"tracking\")" );
            }
        }else{
            throw util::xml_scenario_error( "parasiteGenetics/samplingMode: expected \"initial\" or \"tracking\"" );
        }
        
        // Build the list of all allele combinations by iterating over loci (plural of locus):
        GT::LocusSet loci( genetics.getLocus()[0] );
        for( size_t i = 1; i < genetics.getLocus().size(); i += 1 ){
            GT::LocusSet newLocus(genetics.getLocus()[i]);
            loci.include( newLocus );
        }
        GT::genotypes.swap( loci.alleles );
        N_genotypes = GT::genotypes.size();
        
        double cum_p = 0.0;
        for( size_t i = 0; i < GT::genotypes.size(); ++i ){
            cum_p += GT::genotypes[i].init_freq;
            uint32_t genotype_id = static_cast<uint32_t>(i);
            assert( genotype_id == i );
            GT::cum_initial_freqs.insert( make_pair(cum_p, genotype_id) );
        }
        
        // Test cum_p is approx. 1.0 in case the input tree is wrong.
        if (cum_p < 0.999 || cum_p > 1.001)
            throw util::xml_scenario_error ("decision tree (random node): expected probability sum to be 1.0 but found " + to_string(cum_p));

        // last cum_p might be slightless less than 1 due to arithmetic errors; add a failsafe:
        GT::cum_initial_freqs[1.0] = GT::genotypes.size() - 1;
    }else{
        initSingle();
        GT::cum_initial_freqs[1.0] = GT::genotypes.size() - 1;
    }
    
    if( util::CommandLine::option( util::CommandLine::PRINT_GENOTYPES ) ){
        // reorganise GT::alleleCodes so that we can look up codes, not names
        vector<pair<string,string> > allele_codes( GT::cum_initial_freqs.size() );
        for( auto i = GT::alleleCodes.begin(), iend = GT::alleleCodes.end(); i != iend; ++i ) {
            const string locus = i->first;
            for( auto j = i->second.begin(),
                jend = i->second.end(); j != jend; ++j )
            {
                uint32_t code = j->second;
                const string allele = j->first;
                assert( code < allele_codes.size() );
                allele_codes[code] = make_pair( locus, allele );
            }
        }
        
        // determine our columns
        map<string,uint32_t> longest;      // longest name in column; key is locus
        for( auto i = GT::genotypes.begin(),
            iend = GT::genotypes.end(); i != iend; ++i )
        {
            if( longest.size() == 0 ){
                for( auto j = i->alleles.begin();
                    j != i->alleles.end(); ++j ){
                    const string& locus = allele_codes[*j].first;

                    longest[locus] = locus.length();  // locus name is included in column
                }
            }else assert( longest.size() == i->alleles.size() );
            
            for( auto j = i->alleles.begin();
                j != i->alleles.end(); ++j ){
                auto it = longest.find(allele_codes[*j].first);
                assert( it != longest.end() );
                uint32_t len_allele = allele_codes[*j].second.length();
                if( len_allele > it->second ) it->second = len_allele;
            }
        }
        
        // find original loci order
        vector<string> loci;
        loci.reserve( longest.size() );
        for( size_t i = 0; loci.size() < longest.size(); ++i ){
            assert( i < allele_codes.size() );
            const string& locus = allele_codes[i].first;
            if( !contains(loci, locus) ){
                loci.push_back( locus );
            }
        }
        
//         string fm = "|%1$-12d|%|14t|%2$-12d|";
        cout << endl;
        stringstream fmt;
        fmt << "|%8d|";
        for( auto it = loci.begin(); it !=loci.end(); ++it )
            fmt << "%" << longest[*it] << "s|";
        fmt << "%9.3f|%7.3f|";
        
        // Table header:
        cout << "|" << std::setw(8) << "Genotype"; 
        for( auto it = loci.begin(); it !=loci.end(); ++it )
            cout << "|" << std::setw(longest[*it]) << *it;
        cout << "|" << std::setw(9) << "init freq" << "|" << std::setw(7) << "fitness" << "|" << endl;
        
        cout << "|" << std::setw(8) << "--------"; 
        for( auto it = loci.begin(); it !=loci.end(); ++it )
        {
            if(longest[*it] > 0)
                cout << "|" << "---------";
            else
                cout << "|";
        }

        cout << "|" << std::setw(9) << "---------" << "|" << std::setw(7) << "-------" << "|" << endl;

        cout << std::fixed << setprecision(3);
        for( size_t i = 0; i < GT::genotypes.size(); ++i ){
            const Genotype& genotype = GT::genotypes[i];
            map<string,string> locus_allele;    // to find allele for each locus
            for( auto a = genotype.alleles.begin(); a != genotype.alleles.end(); ++a ){
                assert( *a < allele_codes.size() );
                locus_allele[allele_codes[*a].first] = allele_codes[*a].second;
            }

            cout << "|" << std::setw(8) << (i*100000);
            for( auto it = loci.begin(); it !=loci.end(); ++it ){
                auto la = locus_allele.find( *it );
                assert( la != locus_allele.end() );
                cout << "|" << std::setw(longest[*it]) << la->second << "|";

            }
            cout << std::setw(9) << genotype.init_freq << "|" << std::setw(7) << genotype.fitness << "|" << endl;
        }
    }
}

// ———  diagnostics  ———

Diagnostic::Diagnostic( const Parameters& parameters, const scnXml::Diagnostic& elt ){
    dens_lim = numeric_limits<double>::quiet_NaN();
    specificity = numeric_limits<double>::quiet_NaN();
    
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        return; // both parameters NaN (set above)
    }
    
    if( elt.getDeterministic().present() ){
        dens_lim = elt.getDeterministic().get().getMinDensity();
    }else if( elt.getStochastic().present() ){
        dens_lim = elt.getStochastic().get().getDens_50();
        if( dens_lim == 0.0 ){
            // The equation used for stochastic diagnostics breaks down when
            // dens=dens_lim=0 and for other cases the deterministic model is
            // the same when dens_lim=0.
        }else{
            specificity = elt.getStochastic().get().getSpecificity();
            if( specificity < 0.0 || specificity > 1.0 ){
                throw util::xml_scenario_error(
                    string("diagnostics/diagnostic(").append(elt.getName())
                    .append("): specificity must be in range [0,1]") );
            }
        }
    }else{
        // This should be impossible since according to schema one of these
        // elements must be present.
        throw SWITCH_DEFAULT_EXCEPTION;
    }
    if( dens_lim < 0.0 ){
        throw util::xml_scenario_error(
            string("diagnostics/diagnostic(").append(elt.getName())
            .append("): must have density ≥ 0") );
    }
    
    // We use a bias factor to adjust the "units" used to specify the density
    // of this diagnostic, since estimates from Garki and the standard
    // non-Garki sources are not equivalent to those from the Malariatherapy
    // data (which is used internally).
    if( !elt.getUnits().present() ){
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            // User must be explicit in this case, because presumably the Garki
            // bias is to be used for some diagnostics but likely not all
            // (e.g. neonatal mortality).
            throw util::xml_scenario_error( "diagnostics/diagnostic(*)/units: must specify this attribute when GARKI_DENSITY_BIAS is set" );
        }
        // otherwise we assume "Other"
        dens_lim *= parameters[Parameter::DENSITY_BIAS_NON_GARKI];
    }else if( elt.getUnits().get() == "Other" ){
        dens_lim *= parameters[Parameter::DENSITY_BIAS_NON_GARKI];
    }else if( elt.getUnits().get() == "Garki" ){
        dens_lim *= parameters[Parameter::DENSITY_BIAS_GARKI];
    }else{
        assert( elt.getUnits().get() == "Malariatherapy" );
        // in this case we don't need to use a bias factor
    }
    
    uses_hrp2 = elt.getMechanism() == "HRP2";
}

void diagnostics::init( const Parameters& parameters, const scnXml::Scenario& scenario ){
    diagnostic_set.clear(); // compatibility with unit tests
    
    if(scenario.getDiagnostics().present()){
        for( const scnXml::Diagnostic& diagnostic : scenario.getDiagnostics().get().getDiagnostic() ){
            string name = diagnostic.getName();     // conversion fails without this extra line
            bool inserted = diagnostic_set.insert( make_pair(std::move(name), Diagnostic(parameters, diagnostic)) ).second;
            if( !inserted ){
                throw util::xml_scenario_error( string("diagnostic with this name already set: ").append(diagnostic.getName()) );
            }
        }
    }

    const scnXml::Surveys& surveys = scenario.getMonitoring().getSurveys();
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        // So far the implemented Vivax code does not produce parasite
        // densities, thus this diagnostic model cannot be used.
        diagnostics::monitoring_diagnostic = &diagnostics::make_deterministic( numeric_limits<double>::quiet_NaN() );
    }else if( surveys.getDetectionLimit().present() ){
        if( surveys.getDiagnostic().present() ){
            throw util::xml_scenario_error( "monitoring/surveys: do not "
                "specify both detectionLimit and diagnostic" );
        }
        if( util::CommandLine::option( util::CommandLine::DEPRECATION_WARNINGS ) ){
            std::cerr << "Deprecation warning: monitoring/surveys: "
                "specification of \"diagnostic\" is suggested over \"detectionLimit\"" << std::endl;
        }

        // This controls whether the detection limit is specified relative to
        // the Garki or other methods.
        double densitybias = numeric_limits<double>::quiet_NaN();
        if (util::ModelOptions::option (util::GARKI_DENSITY_BIAS)) {
            densitybias = parameters[Parameter::DENSITY_BIAS_GARKI];
        } else {
            if( scenario.getAnalysisNo().present() ){
                int analysisNo = scenario.getAnalysisNo().get();
                if ((analysisNo >= 22) && (analysisNo <= 30)) {
                    cerr << "Warning: these analysis numbers used to mean "
                        "use Garki density bias. If you do want to use this, "
                        "specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
                }
            }
            densitybias = parameters[Parameter::DENSITY_BIAS_NON_GARKI];
        }
        double detectionLimit = surveys.getDetectionLimit().get() * densitybias;
        diagnostics::monitoring_diagnostic = &diagnostics::make_deterministic( detectionLimit );
    }else{
        if( !surveys.getDiagnostic().present() ){
            throw util::xml_scenario_error( "monitoring/surveys: require "
                "either detectionLimit or diagnostic" );
        }
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            throw util::xml_scenario_error( "Use of GARKI_DENSITY_BIAS is not "
                "appropriate when monitoring/surveys/diagnostic is used." );
        }
        diagnostics::monitoring_diagnostic = &diagnostics::get( surveys.getDiagnostic().get() );
    }
    
    if( diagnostics::monitoring_diagnostic->uses_hrp2 ){
        throw util::xml_scenario_error( "the diagnostic used for "
            "monitoring may not use HRP2 as its mechanism" );
    }
}
}
}
