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

#include "Clock.h"
#include "Context.h"
#include "Parameters.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "Host/WithinHost/Pathogenesis.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"
#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/ModelNameProvider.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "xml/DocumentLoader.h"
#include "schema/scenario.h"

#include <cassert>
#include <algorithm>
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
extern bool opt_vivax_simple, opt_dummy_whm, opt_empirical_whm, opt_molineaux_whm, opt_penny_whm, opt_common_whm;

// State defined in WHFalciparum.cpp (statics made external for this layer):
extern double sigma_i, immPenalty_22, asexImmRemain, immEffectorRemain, invCumulativeYstar, invCumulativeHstar, alpha_m,
    decayM;

// State defined in CommonWithinHost.cpp (external linkage):
extern double hetMassMultStdDev, minHetMassMult;
extern util::AgeGroupInterpolator massByAge;
extern bool reportInfectedOrPatentInfected;

// State defined in Diagnostic.cpp (external linkage):
extern map<string, Diagnostic> diagnostic_set;

void WHInterface::init(Context& ctx, const OM::Parameters& parameters, const scnXml::Scenario& scenario) {
    reportInfectionsByGenotype = mon::isUsed(ctx, mon::measure("nInfectByGenotype")) ||
                                 mon::isUsed(ctx, mon::measure("nPatentByGenotype")) ||
                                 mon::isUsed(ctx, mon::measure("logDensByGenotype"));

    if (util::ModelOptions::option(util::VIVAX_SIMPLE_MODEL)) {
        opt_vivax_simple = true;
        WHVivax::init(ctx, parameters, scenario.getModel());
    } else {
        WHFalciparum::init(ctx, parameters, scenario.getModel());

        if (util::ModelOptions::option(util::DUMMY_WITHIN_HOST_MODEL)) {
            opt_dummy_whm = true;
            DummyInfection::init(ctx);
        } else if (util::ModelOptions::option(util::EMPIRICAL_WITHIN_HOST_MODEL)) {
            opt_empirical_whm = true;
            EmpiricalInfection::init(ctx); // 1-day time step check
        } else if (util::ModelOptions::option(util::MOLINEAUX_WITHIN_HOST_MODEL)) {
            opt_molineaux_whm = true;
            MolineauxInfection::init(ctx, parameters);
        } else if (util::ModelOptions::option(util::PENNY_WITHIN_HOST_MODEL)) {
            opt_penny_whm = true;
            PennyInfection::init(ctx);
        } else {
            DescriptiveInfection::init(ctx, parameters); // 5-day time step check
        }
        opt_common_whm = opt_dummy_whm || opt_empirical_whm || opt_molineaux_whm || opt_penny_whm;

        if (opt_common_whm) {
            WHCommon::init(ctx, scenario);
        } else {
            WHDescriptive::init(ctx);
        }
    }
}

void WHFalciparum::init(Context& ctx, const OM::Parameters& parameters, const scnXml::Model& model) {
    sigma_i = sqrt(parameters[Parameter::SIGMA_I_SQ]);
    immPenalty_22 = 1 - exp(parameters[Parameter::IMMUNITY_PENALTY]);
    immEffectorRemain = exp(-parameters[Parameter::IMMUNE_EFFECTOR_DECAY]);
    asexImmRemain = exp(-parameters[Parameter::ASEXUAL_IMMUNITY_DECAY]);

    // calculate inverses here, so we can use multiplication later (faster):
    invCumulativeYstar = 1.0 / parameters[Parameter::CUMULATIVE_Y_STAR];
    invCumulativeHstar = 1.0 / parameters[Parameter::CUMULATIVE_H_STAR];
    alpha_m = 1.0 - exp(-parameters[Parameter::NEG_LOG_ONE_MINUS_ALPHA_M]);
    decayM = parameters[Parameter::DECAY_M];

    y_lag_len = ctx.clock.inSteps(sim::fromDays(20)) + 1;

    //NOTE: should also call cleanup() on the Model, but it only frees memory which the OS does anyway
    Pathogenesis::Model::init(ctx, parameters, model.getClinical(), false);

    // Default value to use if no <parameters> element is present in XML.
    SimTime latentP = read_latentP(ctx.clock, model);

    Infection::init(latentP);
}

void WHCommon::init(Context& ctx, const scnXml::Scenario& scenario) {
    const scnXml::Human human = scenario.getModel().getHuman();
    if (!human.getWeight().present()) {
        // Technically this is needed by the PK/PD and Molineaux models
        throw util::xml_error("model->human->weight element required by certain models");
    }
    massByAge.set(ctx.clock, human.getWeight().get(), "weight");
    hetMassMultStdDev = human.getWeight().get().getMultStdDev();
    // hetWeightMult must be large enough that birth weight is at least 0.5 kg:
    minHetMassMult = 0.5 / massByAge.eval(0.0);

    reportInfectedOrPatentInfected =
        mon::isUsed(ctx, mon::measure("totalInfs")) || mon::isUsed(ctx, mon::measure("totalPatentInf"));

    PkPd::LSTMModel::init(ctx, scenario);
}

// ———  parasite genetics  ———

namespace GT /*for genotype impl details*/ {
// Represent a set of loci: all possible combinations of alleles
// This is actually just machinery to calculate the list of all genotypes
struct LocusSet {
    LocusSet(Genotypes& state, const scnXml::ParasiteLocus& elt_l) : st(state) {
        alleles.reserve(elt_l.getAllele().size());
        double cum_p = 0.0;
        for (const scnXml::ParasiteAllele& elt_a : elt_l.getAllele()) {
            uint32_t alleleCode = (st.nextAlleleCode++);
            st.alleleCodes[elt_l.getName()][elt_a.getName()] = alleleCode;
            cum_p += elt_a.getInitialFrequency();
            alleles.push_back(Genotypes::Genotype(alleleCode, elt_a.getInitialFrequency(), elt_a.getFitness(),
                                                  elt_a.getHrp2_deletion()));
        }
        if (cum_p < 0.999 || cum_p > 1.001)
            throw util::xml_error("expected sum of initial probabilities of alleles to be 1, but for the " +
                                           to_string(alleles.size()) + " alleles under locus " +
                                           string(elt_l.getName()) + " this is " + to_string(cum_p));

        alleles[0].init_freq += 1.0 - cum_p; // account for any small errors by adjusting the first frequency
    }
    void include(LocusSet& that) {
        vector<Genotypes::Genotype> newAlleles;
        newAlleles.reserve(alleles.size() * that.alleles.size());
        for (size_t i = 0; i < alleles.size(); i += 1) {
            for (size_t j = 0; j < that.alleles.size(); j += 1) {
                // note that we generate new codes each time; we just waste the old ones
                newAlleles.push_back(alleles[i].cross(that.alleles[j]));
            }
        }
        alleles.swap(newAlleles);
    }

    vector<Genotypes::Genotype> alleles;
    Genotypes& st;
};

}

void Genotypes::init(Genotypes& state, const scnXml::Scenario& scenario) {
    if (scenario.getParasiteGenetics().present()) {
        const scnXml::ParasiteGenetics& genetics = scenario.getParasiteGenetics().get();

        state.current_mode = Genotypes::SAMPLE_INITIAL; // turn on sampling
        if (genetics.getSamplingMode() == "initial") {
            state.interv_mode = Genotypes::SAMPLE_INITIAL;
        } else if (genetics.getSamplingMode() == "tracking") {
            state.interv_mode = Genotypes::SAMPLE_TRACKING;
            if (!scenario.getEntomology().getVector().present()) {
                throw util::xml_error("incompatibility; either use "
                                               "entomology/vector (not nonVector) or set "
                                               "parasiteGenetics/samplingMode to \"initial\" (not \"tracking\")");
            }
            if (scenario.getEntomology().getMode() != "dynamic") {
                throw util::xml_error("incompatibility; either set "
                                               "entomology/mode to \"dynamic\" (not \"forced\") or set "
                                               "parasiteGenetics/samplingMode to \"initial\" (not \"tracking\")");
            }
        } else {
            throw util::xml_error("parasiteGenetics/samplingMode: expected \"initial\" or \"tracking\"");
        }

        // Build the list of all allele combinations by iterating over loci (plural of locus):
        GT::LocusSet loci(state, genetics.getLocus()[0]);
        for (size_t i = 1; i < genetics.getLocus().size(); i += 1) {
            GT::LocusSet newLocus(state, genetics.getLocus()[i]);
            loci.include(newLocus);
        }
        state.genotypes.swap(loci.alleles);
        state.N_genotypes = state.genotypes.size();

        double cum_p = 0.0;
        for (size_t i = 0; i < state.genotypes.size(); ++i) {
            cum_p += state.genotypes[i].init_freq;
            uint32_t genotype_id = static_cast<uint32_t>(i);
            assert(genotype_id == i);
            state.cum_initial_freqs.insert(make_pair(cum_p, genotype_id));
        }

        // Test cum_p is approx. 1.0 in case the input tree is wrong.
        if (cum_p < 0.999 || cum_p > 1.001)
            throw util::xml_error(
                "decision tree (random node): expected probability sum to be 1.0 but found " + to_string(cum_p));

        // last cum_p might be slightless less than 1 due to arithmetic errors; add a failsafe:
        state.cum_initial_freqs[1.0] = state.genotypes.size() - 1;
    } else {
        state.initSingle();
        state.cum_initial_freqs[1.0] = state.genotypes.size() - 1;
    }

    if (util::CommandLine::option(util::CommandLine::PRINT_GENOTYPES)) {
        // reorganise state.alleleCodes so that we can look up codes, not names
        vector<pair<string, string>> allele_codes(state.cum_initial_freqs.size());
        for (auto i = state.alleleCodes.begin(), iend = state.alleleCodes.end(); i != iend; ++i) {
            const string locus = i->first;
            for (auto j = i->second.begin(), jend = i->second.end(); j != jend; ++j) {
                uint32_t code = j->second;
                const string allele = j->first;
                assert(code < allele_codes.size());
                allele_codes[code] = make_pair(locus, allele);
            }
        }

        // determine our columns
        map<string, uint32_t> longest; // longest name in column; key is locus
        for (auto i = state.genotypes.begin(), iend = state.genotypes.end(); i != iend; ++i) {
            if (longest.size() == 0) {
                for (auto j = i->alleles.begin(); j != i->alleles.end(); ++j) {
                    const string& locus = allele_codes[*j].first;

                    longest[locus] = locus.length(); // locus name is included in column
                }
            } else
                assert(longest.size() == i->alleles.size());

            for (auto j = i->alleles.begin(); j != i->alleles.end(); ++j) {
                auto it = longest.find(allele_codes[*j].first);
                assert(it != longest.end());
                uint32_t len_allele = allele_codes[*j].second.length();
                if (len_allele > it->second)
                    it->second = len_allele;
            }
        }

        // find original loci order
        vector<string> loci;
        loci.reserve(longest.size());
        for (size_t i = 0; loci.size() < longest.size(); ++i) {
            assert(i < allele_codes.size());
            const string& locus = allele_codes[i].first;
            if (std::find(loci.begin(), loci.end(), locus) == loci.end()) {
                loci.push_back(locus);
            }
        }

        //         string fm = "|%1$-12d|%|14t|%2$-12d|";
        cout << endl;
        stringstream fmt;
        fmt << "|%8d|";
        for (auto it = loci.begin(); it != loci.end(); ++it)
            fmt << "%" << longest[*it] << "s|";
        fmt << "%9.3f|%7.3f|";

        // Table header:
        cout << "|" << std::setw(8) << "Genotype";
        for (auto it = loci.begin(); it != loci.end(); ++it)
            cout << "|" << std::setw(longest[*it]) << *it;
        cout << "|" << std::setw(9) << "init freq" << "|" << std::setw(7) << "fitness" << "|" << endl;

        cout << "|" << std::setw(8) << "--------";
        for (auto it = loci.begin(); it != loci.end(); ++it) {
            if (longest[*it] > 0)
                cout << "|" << "---------";
            else
                cout << "|";
        }

        cout << "|" << std::setw(9) << "---------" << "|" << std::setw(7) << "-------" << "|" << endl;

        cout << std::fixed << setprecision(3);
        for (size_t i = 0; i < state.genotypes.size(); ++i) {
            const Genotype& genotype = state.genotypes[i];
            map<string, string> locus_allele; // to find allele for each locus
            for (auto a = genotype.alleles.begin(); a != genotype.alleles.end(); ++a) {
                assert(*a < allele_codes.size());
                locus_allele[allele_codes[*a].first] = allele_codes[*a].second;
            }

            cout << "|" << std::setw(8) << (i * 100000);
            for (auto it = loci.begin(); it != loci.end(); ++it) {
                auto la = locus_allele.find(*it);
                assert(la != locus_allele.end());
                cout << "|" << std::setw(longest[*it]) << la->second << "|";
            }
            cout << std::setw(9) << genotype.init_freq << "|" << std::setw(7) << genotype.fitness << "|" << endl;
        }
    }
}

// ———  diagnostics  ———

Diagnostic::Diagnostic(const Parameters& parameters, const scnXml::Diagnostic& elt) {
    dens_lim = numeric_limits<double>::quiet_NaN();
    specificity = numeric_limits<double>::quiet_NaN();

    if (util::ModelOptions::option(util::VIVAX_SIMPLE_MODEL)) {
        return; // both parameters NaN (set above)
    }

    if (elt.getDeterministic().present()) {
        dens_lim = elt.getDeterministic().get().getMinDensity();
    } else if (elt.getStochastic().present()) {
        dens_lim = elt.getStochastic().get().getDens_50();
        if (dens_lim == 0.0) {
            // The equation used for stochastic diagnostics breaks down when
            // dens=dens_lim=0 and for other cases the deterministic model is
            // the same when dens_lim=0.
        } else {
            specificity = elt.getStochastic().get().getSpecificity();
            if (specificity < 0.0 || specificity > 1.0) {
                throw util::xml_error(string("diagnostics/diagnostic(")
                                                   .append(elt.getName())
                                                   .append("): specificity must be in range [0,1]"));
            }
        }
    } else {
        // This should be impossible since according to schema one of these
        // elements must be present.
        throw SWITCH_DEFAULT_EXCEPTION;
    }
    if (dens_lim < 0.0) {
        throw util::xml_error(
            string("diagnostics/diagnostic(").append(elt.getName()).append("): must have density ≥ 0"));
    }

    // We use a bias factor to adjust the "units" used to specify the density
    // of this diagnostic, since estimates from Garki and the standard
    // non-Garki sources are not equivalent to those from the Malariatherapy
    // data (which is used internally).
    if (!elt.getUnits().present()) {
        if (util::ModelOptions::option(util::GARKI_DENSITY_BIAS)) {
            // User must be explicit in this case, because presumably the Garki
            // bias is to be used for some diagnostics but likely not all
            // (e.g. neonatal mortality).
            throw util::xml_error(
                "diagnostics/diagnostic(*)/units: must specify this attribute when GARKI_DENSITY_BIAS is set");
        }
        // otherwise we assume "Other"
        dens_lim *= parameters[Parameter::DENSITY_BIAS_NON_GARKI];
    } else if (elt.getUnits().get() == "Other") {
        dens_lim *= parameters[Parameter::DENSITY_BIAS_NON_GARKI];
    } else if (elt.getUnits().get() == "Garki") {
        dens_lim *= parameters[Parameter::DENSITY_BIAS_GARKI];
    } else {
        assert(elt.getUnits().get() == "Malariatherapy");
        // in this case we don't need to use a bias factor
    }

    uses_hrp2 = elt.getMechanism() == "HRP2";
}

void diagnostics::init(Context& ctx, const Parameters& parameters, const scnXml::Scenario& scenario) {
    auto& diagnostic_set = ctx.diagnostic_set;
    diagnostic_set.clear(); // compatibility with unit tests

    if (scenario.getDiagnostics().present()) {
        for (const scnXml::Diagnostic& diagnostic : scenario.getDiagnostics().get().getDiagnostic()) {
            string name = diagnostic.getName(); // conversion fails without this extra line
            bool inserted =
                diagnostic_set.insert(make_pair(std::move(name), Diagnostic(parameters, diagnostic))).second;
            if (!inserted) {
                throw util::xml_error(
                    string("diagnostic with this name already set: ").append(diagnostic.getName()));
            }
        }
    }

    const scnXml::Surveys& surveys = scenario.getMonitoring().getSurveys();
    if (util::ModelOptions::option(util::VIVAX_SIMPLE_MODEL)) {
        // So far the implemented Vivax code does not produce parasite
        // densities, thus this diagnostic model cannot be used.
        ctx.monitoring_diagnostic = &make_deterministic(ctx, numeric_limits<double>::quiet_NaN());
    } else if (surveys.getDetectionLimit().present()) {
        if (surveys.getDiagnostic().present()) {
            throw util::xml_error("monitoring/surveys: do not "
                                           "specify both detectionLimit and diagnostic");
        }
        if (util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS)) {
            std::cerr << "Deprecation warning: monitoring/surveys: "
                         "specification of \"diagnostic\" is suggested over \"detectionLimit\""
                      << std::endl;
        }

        // This controls whether the detection limit is specified relative to
        // the Garki or other methods.
        double densitybias = numeric_limits<double>::quiet_NaN();
        if (util::ModelOptions::option(util::GARKI_DENSITY_BIAS)) {
            densitybias = parameters[Parameter::DENSITY_BIAS_GARKI];
        } else {
            if (scenario.getAnalysisNo().present()) {
                int analysisNo = scenario.getAnalysisNo().get();
                if ((analysisNo >= 22) && (analysisNo <= 30)) {
                    cerr << "Warning: these analysis numbers used to mean "
                            "use Garki density bias. If you do want to use this, "
                            "specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong."
                         << endl;
                }
            }
            densitybias = parameters[Parameter::DENSITY_BIAS_NON_GARKI];
        }
        double detectionLimit = surveys.getDetectionLimit().get() * densitybias;
        ctx.monitoring_diagnostic = &make_deterministic(ctx, detectionLimit);
    } else {
        if (!surveys.getDiagnostic().present()) {
            throw util::xml_error("monitoring/surveys: require "
                                           "either detectionLimit or diagnostic");
        }
        if (util::ModelOptions::option(util::GARKI_DENSITY_BIAS)) {
            throw util::xml_error("Use of GARKI_DENSITY_BIAS is not "
                                           "appropriate when monitoring/surveys/diagnostic is used.");
        }
        ctx.monitoring_diagnostic = &get(ctx, surveys.getDiagnostic().get());
    }

    if (ctx.monitoring_diagnostic->uses_hrp2) {
        throw util::xml_error("the diagnostic used for "
                                       "monitoring may not use HRP2 as its mechanism");
    }
}
}
}

namespace OM {

SimTime read_latentP(const SimClock& clock, const scnXml::Model& model) {
    SimTime latentP = UnitParse::readShortDuration(clock, "15d", UnitParse::DAYS);
    try {
        if (model.getParameters().present())
            latentP = UnitParse::readShortDuration(clock, model.getParameters().get().getLatentp(), UnitParse::NONE);
    } catch (const util::format_error& e) {
        throw util::xml_error(string("model/parameters/latentP: ").append(e.message()));
    }
    return latentP;
}

}

// ———  SimClock::init: scenario-level reading  ———

namespace OM {

void SimClock::init(Context& ctx, const scnXml::Scenario& scenario, util::ModelNameProvider mnp) {
    initInterval(scenario, mnp);
    steps_per_year = inSteps(sim::oneYear());
    years_per_step = 1.0 / steps_per_year;
    max_human_age = from_years_d(scenario.getDemography().getMaximumAgeYrs());

    start_date = sim::origin();
    auto mon = scenario.getMonitoring();
    if (mon.getStartDate().present()) {
        try {
            // on failure, this throws or returns sim::never()
            start_date = UnitParse::parseDate(*this, mon.getStartDate().get());
            if (start_date == sim::never()) {
                throw util::format_error("invalid format (expected YYYY-MM-DD)");
            }
        } catch (const util::format_error& e) {
            throw util::xml_error(string("monitoring/startDate: ").append(e.message()));
        }
    }

    interv = sim::never(); // large negative number
    end_date = mon::readSurveyDates(ctx, mon);
}

void SimClock::initInterval(const scnXml::Scenario& scenario, util::ModelNameProvider mnp) {
    // If scenario states to use a named model, set the relevant time step duration for
    // that named model here.  It may possibly get overwritten later in this method if
    // the scenario also states an explicit interval value - a behaviour that exists so
    // that users have the ability to override the time step duration if they need to.
    const util::ModelNames namedModelToUse = mnp.modelInUse;
    const bool useNamedModel = namedModelToUse != util::ModelNames::none;
    if (useNamedModel) {
        if (namedModelToUse == util::ModelNames::base) {
            // This is where we define the length of the default time step, in days,
            // for the base model.
            interval = 5;
        } else {
            throw util::xml_error("No pre-set interval value is available for the specified model name.");
        }
    }

    // Apply the interval duration explicitly written in the scenario, if one exists.
    const bool parameterBlockExistsInScenario = scenario.getModel().getParameters().present();
    if (parameterBlockExistsInScenario) {
        // Schema dictates that interval attribute is required if parameter element exists.
        // So no need to check for existsence of interval attribute in this case.
        interval = scenario.getModel().getParameters().get().getInterval();
    }
}

namespace util {

ModelNameProvider::ModelNameProvider(const scnXml::Model& model) {
    const bool useNamedModel = model.getModelName().present();
    if (useNamedModel) {
        const std::string nameFromXML = model.getModelName().get().getName();
        if (nameFromXML == "base") {
            modelInUse = ModelNames::base;
        } else {
            throw util::xml_error("Unrecognized model name: " + nameFromXML);
        }
    }
}

}

}

// ———  Parameters: model parameter block reading  ———

namespace OM {

Parameters::Parameters(const scnXml::Model::ParametersOptional& xmlParameterValues, util::ModelNameProvider mnp) {
    for (const std::pair<int, Parameter> pair : idCodeToNameMap) {
        const Parameter name = pair.second;
        nameToValueMap[name] = std::nullopt;
    }

    const util::ModelNames namedModelToUse = mnp.modelInUse;
    const bool useNamedModel = namedModelToUse != util::ModelNames::none;
    if (useNamedModel) {
        if (namedModelToUse == util::ModelNames::base) {
            initializeParamsBaseModel();
        } else {
            throw util::xml_error("No pre-set parameter values are available for the specified model name.");
        }
    }

    // Get parameters specified expicitly in input XML, if any.
    // In the event where a named model was used to set parameter values earlier, it is intended
    // that any parameter values written explicitly in XML will override any values set by the
    // named model.
    const bool useExplicitParamValues = xmlParameterValues.present();
    if (useExplicitParamValues) {
        initializeParamsFromXML(xmlParameterValues.get());
    }
}

void Parameters::initializeParamsFromXML(const scnXml::Parameters& parameters) {
    const scnXml::Parameters::ParameterSequence& paramSeq = parameters.getParameter();

    // It's okay if a user overwrites the value of a parameter that was set before this
    // method was called.  E.g. if they used a named model, and then chose to manually
    // override the value of a parameter set by it.  It's not okay if a user sets a value
    // for the same parameter twice.  That represents a mistake in the input XML which
    // we need to handle here.  This has the potential to save users time debugging inputs.
    std::set<int> paramIdsSetByUser;

    for (auto param : paramSeq) {
        const int paramId = param.getNumber();
        const double paramValueFromXML = param.getValue();

        // C++20 offers "contains" which would be more expressive than "count" here.
        const bool paramIdIsValid = idCodeToNameMap.count(paramId);
        if (!paramIdIsValid) {
            if (util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS)) {
                cerr << "Deprecation warning: <parameter> index" << to_string(paramId) << " is no longer used" << endl;
            }

            // If we throw here, old scenarios containing the deprecated parameter
            // will need more manual work to migrate to new OpenMalaria versions.
            // This isn't unacceptable.  But it's simpler here to just skip over the
            // deprecated parameter.  In particular, many tests specify values for
            // deprecated some parameters.
            continue;
        }
        const Parameter nameOfParamToSet = idCodeToNameMap.at(paramId);

        const bool paramValueAlreadySet = !paramIdsSetByUser.insert(paramId).second;
        if (paramValueAlreadySet) {
            throw util::xml_error("Parameter with index " + to_string(paramId) + " described twice in XML.");
        }

        nameToValueMap[nameOfParamToSet] = paramValueFromXML;
    }
}

}
