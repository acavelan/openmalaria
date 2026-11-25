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

// model/xml/ conversion layer: every scnXml -> plain-value function lives
// under model/xml/; only this directory includes generated schema headers
// (plan §5, §8 F7: XSD codegen + Xerces validation kept). Everything below was
// moved VERBATIM from the noted origin files; declarations stay in the normal
// model headers (scnXml types forward-declared there). Call sites unchanged.

#include "Global.h"
#include "core/RunContext.h"
#include "util/AgeGroupInterpolation.h"
#include "util/CommandLine.h"
#include "util/DecayFunction.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "sim.h"
#include "util/sampler.h"
#include <schema/healthSystem.h>
#include <schema/scenario.h>
#include <schema/util.h>

namespace OM { namespace util {

// ---- sampler readers (origins: util/sampler.{h,cpp}) ----

void NormalSampler::setParams(const scnXml::SampledValueN& elt){
    mu = elt.getMean();
    if( elt.getDistr() == "const" ){
        if( elt.getSD().present() && elt.getSD().get() != 0.0) {
            throw util::xml_scenario_error( "attribute SD must be zero or omitted when distr=\"const\" or is omitted" );
        }
        sigma = 0.0;
        return;
    }
    
    if( elt.getDistr() != "normal" ){
        throw util::xml_scenario_error( "expected distr to be one of \"const\", \"lognormal\" (note: not all distributions are supported here)" );
    }
    if( !elt.getSD().present() ){
        throw util::xml_scenario_error( "attribute \"SD\" required for sampled normal value when distr is not const" );
    }
    sigma = elt.getSD().get();
}

    /** Create a Sampler of the given kind from a SampledValueCV xml snippet. */
    inline unique_ptr<Sampler> createSampler(Sampler::Kind kind, double mean, const scnXml::SampledValueCV& elt)
    {
        const bool isLognormal = kind == Sampler::Kind::Lognormal;
        const char* expectedDistr = isLognormal ? "lognormal" : "gamma";

        const std::string& d = elt.getDistr();

        // Special handling for the "const" alias (lognormal with CV = 0)
        if( d == "const" )
        {
            if (!isLognormal) // only allowed for lognormal
                throw util::xml_scenario_error("\"distr=\"const\": expected \"lognormal\"");

            if( elt.getCV().present() && elt.getCV().get() != 0.0 )
                throw util::xml_scenario_error( "\"distr="+d+"\": attribute \"CV\" must be zero or omitted when distr=\"const\" or is omitted" );
            if( elt.getVariance().present() && elt.getVariance().get() != 0.0 )
                throw util::xml_scenario_error( "\"distr="+d+"\": attribute \"variance\" must be zero or omitted when distr=\"const\" or is omitted" );

            return Sampler::fromMeanCV(Sampler::Kind::Lognormal, mean, 0.0);
        }

        // Check that XML 'distr' matches the sampler kind
        if (d != expectedDistr)
            throw util::xml_scenario_error("\"distr=" + d + "\": expected \"" + expectedDistr + "\" or \"const\"");

        // Exactly one of CV / variance must be supplied
        const bool hasCV  = elt.getCV().present();
        const bool hasVar = elt.getVariance().present();

        if (hasCV == hasVar)   // both present or both absent
            throw util::xml_scenario_error("\"distr=" + d + "\": exactly one of attributes \"CV\" or \"variance\" must be specified");

        // The optional truncation parameter
        std::optional<double> truncate = elt.getTruncate().present()
            ? std::optional<double>{elt.getTruncate().get()}
            : std::nullopt;

        if(hasCV)
            return Sampler::fromMeanCV(kind, mean, elt.getCV().get(), truncate);
        else
            return Sampler::fromMeanVariance(kind, mean, elt.getVariance().get(), truncate);
    }

    /** Create a Sampler from a SampledValueCV xml snippet. */
    std::unique_ptr<util::Sampler> createSampler(double mean, const scnXml::SampledValueCV& elt)
    {
        const std::string& d = elt.getDistr();

        if (d == "const" || d == "lognormal")
            return createSampler(Sampler::Kind::Lognormal, mean, elt);
        else if (d == "gamma")
            return createSampler(Sampler::Kind::Gamma, mean, elt);
        else
            throw util::xml_scenario_error("\"distr=" + d + "\": expected \"const\", \"lognormal\" or \"gamma\"");
    }

    /** Create a lognormal Sampler from a SampledValueLN xml snippet. */
    std::unique_ptr<Sampler> createSampler(const scnXml::SampledValueLN& elt)
    {
        return createSampler(Sampler::Kind::Lognormal, elt.getMean(), static_cast<const scnXml::SampledValueCV&>(elt));
    }

// ---- BetaSampler / WeibullSampler XML readers (origin: util/sampler.h) ----
void BetaSampler::setParams( const scnXml::BetaMeanSample& elt ){
    setParamsMV( elt.getMean(), elt.getVariance() );
}

void WeibullSampler::setParams( const scnXml::WeibullSample& elt ){
    setScaleShape( elt.getScale(), elt.getShape() );
}

// ---- DecayFunction factory (origin: util/DecayFunction.cpp subclasses) ----
// Config construction only: Kind + parameters into the plain value type; the
// formulas live in util/DecayFunction.cpp. Error messages and their check
// order travel verbatim from the old subclass constructors.

inline double readLToDays( const SimClock& clock, const scnXml::DecayFunction& elt ){
    if( !elt.getL().present() ){
        throw xml_scenario_error( "decay function: attribute L required" );
    }
    return UnitParse::durationToDays(clock, elt.getL().get(), UnitParse::YEARS);
}

unique_ptr<DecayFunction> DecayFunction::makeObject(
    const SimClock& clock, const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    Kind kind;
    // parity: old factory dispatch order (error priority for unknown names)
    if( func == "constant" ) kind = Kind::Constant;
    else if( func == "step" ) kind = Kind::Step;
    else if( func == "linear" ) kind = Kind::Linear;
    else if( func == "exponential" ) kind = Kind::Exponential;
    else if( func == "weibull" ) kind = Kind::Weibull;
    else if( func == "hill" ) kind = Kind::Hill;
    else if( func == "smooth-compact" ) kind = Kind::SmoothCompact;
    else if( func == "plus" ) kind = Kind::Plus;
    else if( func == "minus" ) kind = Kind::Minus;
    else if( func == "divides" ) kind = Kind::Divides;
    else if( func == "multiplies" ) kind = Kind::Multiplies;
    else if( func == "emax" ) kind = Kind::Emax;
    else
        throw xml_scenario_error("decay function " + string(eltName) + " of type " + string(func) + " unrecognized");

    // base ctor first (het lognormal built for every kind; parity: old base
    // ctor ran before any subclass member init / body):
    auto o = make_unique<DecayFunction>(kind, elt.getIncreasing(),
                                        elt.getInitialEfficacy(), elt.getCV());
    o->simClock = &clock;
    switch(kind) {
    case Kind::Constant:
        break;
    case Kind::Step: case Kind::Linear:
        o->invL = 1.0 / readLToDays(clock, elt);
        break;
    case Kind::Exponential:
        o->invL = log(2.0) / readLToDays(clock, elt);
        break;
    case Kind::Weibull:
        o->k = elt.getK();
        o->invL = pow(log(2.0), 1.0 / elt.getK()) / readLToDays(clock, elt);
        break;
    case Kind::Hill: case Kind::SmoothCompact:
        o->invL = 1.0 / readLToDays(clock, elt);
        o->k = elt.getK();
        break;
    case Kind::Plus: case Kind::Minus: case Kind::Divides: case Kind::Multiplies: {
        const scnXml::DecayFunction::DecaySequence& decaySequence = elt.getDecay();
        if(decaySequence.size() != 2)
            throw xml_scenario_error("Operator decay function expects two decay functions, " + to_string(decaySequence.size()) +"  were given.");
        o->f1 = makeObject(clock, decaySequence[0], "Operator::f1");
        o->f2 = makeObject(clock, decaySequence[1], "Operator::f2");
        break;
    }
    case Kind::Emax: {
        // check order verbatim from the old EmaxFunction ctor (error priority):
        const scnXml::DecayFunction::DecaySequence& decaySequence = elt.getDecay();
        if(decaySequence.size() != 1)
            throw util::xml_scenario_error("The Emax function expects exactly one user function, but " + to_string(decaySequence.size()) +"  were given.");
        if(!elt.getEmax().present())
            throw util::xml_scenario_error("Emax function: the Emax parameter is not declared");
        if(!elt.getIC50().present())
            throw util::xml_scenario_error("Emax function: the IC50 parameter is not declared");
        if(!elt.getSlope().present())
            throw util::xml_scenario_error("Emax function: the Slope parameter is not declared");
        if(!elt.getInitialConcentration().present())
            throw util::xml_scenario_error("Emax function: the InitialConcentration parameter is not declared");
        o->f1 = makeObject(clock, decaySequence[0], "Emax::f");
        o->Emax = elt.getEmax().get();
        const double rawIC50 = elt.getIC50().get();
        o->slope = elt.getSlope().get();
        o->initialConcentration = elt.getInitialConcentration().get();
        if (o->Emax < 0.0 || o->Emax > 1.0) {
            throw util::xml_scenario_error(
                "Emax function: Emax must be in [0, 1], got " + std::to_string(o->Emax));
        }
        if (rawIC50 <= 0.0) {
            throw util::xml_scenario_error(
                "Emax function: IC50 must be > 0, got " + std::to_string(rawIC50));
        }
        if (o->slope <= 0.0) {
            throw util::xml_scenario_error(
                "Emax function: Slope must be > 0, got " + std::to_string(o->slope));
        }
        if (o->initialConcentration < 0.0) {
            throw util::xml_scenario_error(
                "Emax function: InitialConcentration must be >= 0, got " + std::to_string(o->initialConcentration));
        }
        // Store IC50^slope once since compute() evaluates this repeatedly.
        o->IC50 = std::pow(rawIC50, o->slope);
        break;
    }
    }
    return o;
}

// ---- AgeGroupInterpolator::set — config construction (origins: the old
// AgeGroupPiecewiseConstant / InterpolationPoints / PiecewiseLinear ctors in
// this file). Formulas + lookups live in util/AgeGroupInterpolation.{h,cpp}.

void AgeGroupInterpolator::set(const SimClock& clock,
    const scnXml::AgeGroupValues& ageGroups, const char* eltName
){
    reset();
    // Type mostly equivalent to a std::string:
    const scnXml::AgeGroupValues::InterpolationOptional& interp = ageGroups.getInterpolation();
    if( !interp.present() || interp.get() == "linear" ){
        // parity: AgeGroupInterpolationPoints + AgeGroupPiecewiseLinear ctors
        auto it = ageGroups.getGroup().begin();
        if( it == ageGroups.getGroup().end() )
            throw util::xml_scenario_error(string(eltName) + ": at least one age group required");
        // first lower-bound must be 0
        if( it->getLowerbound() != 0.0 )
            throw util::xml_scenario_error(string(eltName) + ": first lower-bound must be 0");

        // Our insert iterator (used for performance)
        auto pos = data.begin();
        double greatestLbound = 0.0;    // also last lbound: age groups must be ordered lowest to highest
        double lastValue = it->getValue();       // first value is repeated for const start
        for( ; it!=ageGroups.getGroup().end(); ++it ){
            double lbound = it->getLowerbound();
            if( lbound >= greatestLbound ){
                double insBound = 0.5*(greatestLbound + lbound);
                pos = data.insert( pos, make_pair( insBound, lastValue ) );
                greatestLbound = lbound;
                lastValue = it->getValue();     // always insert previous value
            } else
                throw util::xml_scenario_error(string(eltName) + ": lower bound " + to_string(lbound) + " less than previous " + to_string(greatestLbound));
        }
        // add a point in middle of last age group (taking upper bound as max-age-years:
        data[ 0.5*(greatestLbound + sim::inYears(clock.max_human_age)) ] = lastValue;
        // Add first point at zero
        data[ 0.0 ] = data.begin()->second;
        // add a point at infinity to catch remaining values
        data[ numeric_limits<double>::infinity() ] = data.rbegin()->second;

        kind = Kind::Linear;
        outputSamples( clock, eltName );
    }else if( interp.get() == "none" ){
        // parity: AgeGroupPiecewiseConstant ctor
        if( ageGroups.getGroup().size() == 0 ){
            throw util::xml_scenario_error(string(eltName) + ": at least one age group required" );
        }
        auto pos = data.begin();
        double greatestLbound = -1.0;
        for( const scnXml::Group& group : ageGroups.getGroup() ){
            double lbound = group.getLowerbound();
            if( lbound >= greatestLbound ){
                greatestLbound = lbound;
                pos = data.insert( pos, make_pair( lbound, group.getValue() ) );
            } else
                throw util::xml_scenario_error(string(eltName) + ": lower bound " + to_string(lbound) + " not greater than previous " + to_string(greatestLbound));
        }
        // first lower-bound must be 0
        if( data.begin()->first != 0.0 )
            throw util::xml_scenario_error(string(eltName) + ": first lower-bound must be 0");

        kind = Kind::PiecewiseConstant;
        outputSamples( clock, eltName );
    }else
        throw util::xml_scenario_error(string("age group interpolation ") + interp.get() + " not implemented" );
}

} }

// ———  ModelOptions: XML option reading  ———

namespace OM { namespace util {

    // Utility: converts option strings to codes and back
    class OptionCodeMap {
	// Lookup table to translate the strings used in the XML file to the internal enumerated values:
	map<string,OptionCodes> codeMap;
        set<string> ignoreOptions;
	
    public:
	OptionCodeMap () {
// 	    codeMap["PENALISATION_EPISODES"] = PENALISATION_EPISODES;
	    codeMap["NEGATIVE_BINOMIAL_MASS_ACTION"] = NEGATIVE_BINOMIAL_MASS_ACTION;
	    // codeMap["ATTENUATION_ASEXUAL_DENSITY"] = ATTENUATION_ASEXUAL_DENSITY;
	    codeMap["LOGNORMAL_MASS_ACTION"] = LOGNORMAL_MASS_ACTION;
	    codeMap["NO_PRE_ERYTHROCYTIC"] = NO_PRE_ERYTHROCYTIC;
	    codeMap["MAX_DENS_CORRECTION"] = MAX_DENS_CORRECTION;
	    codeMap["INNATE_MAX_DENS"] = INNATE_MAX_DENS;
	    // 	codeMap["MAX_DENS_RESET"] = MAX_DENS_RESET;
	    codeMap["DUMMY_WITHIN_HOST_MODEL"] = DUMMY_WITHIN_HOST_MODEL;
	    codeMap["PREDETERMINED_EPISODES"] = PREDETERMINED_EPISODES;
	    codeMap["NON_MALARIA_FEVERS"] = NON_MALARIA_FEVERS;
        ignoreOptions.insert("INCLUDES_PK_PD");
	    codeMap["CLINICAL_EVENT_SCHEDULER"] = CLINICAL_EVENT_SCHEDULER;
	    codeMap["MUELLER_PRESENTATION_MODEL"] = MUELLER_PRESENTATION_MODEL;
	    codeMap["TRANS_HET"] = TRANS_HET;
	    codeMap["COMORB_HET"] = COMORB_HET;
	    codeMap["TREAT_HET"] = TREAT_HET;
	    codeMap["COMORB_TRANS_HET"] = COMORB_TRANS_HET;
	    codeMap["TRANS_TREAT_HET"] = TRANS_TREAT_HET;
	    codeMap["COMORB_TREAT_HET"] = COMORB_TREAT_HET;
	    codeMap["TRIPLE_HET"] = TRIPLE_HET;
	    codeMap["EMPIRICAL_WITHIN_HOST_MODEL"] = EMPIRICAL_WITHIN_HOST_MODEL;
        codeMap["MOLINEAUX_WITHIN_HOST_MODEL"] = MOLINEAUX_WITHIN_HOST_MODEL;
        codeMap["PENNY_WITHIN_HOST_MODEL"] = PENNY_WITHIN_HOST_MODEL;
        codeMap["GARKI_DENSITY_BIAS"] = GARKI_DENSITY_BIAS;
//      codeMap["IPTI_SP_MODEL"] = IPTI_SP_MODEL;
//      codeMap["REPORT_ONLY_AT_RISK"] = REPORT_ONLY_AT_RISK; 
	    codeMap["MEAN_DURATION_GAMMA"] = MEAN_DURATION_GAMMA;
	    codeMap["FIRST_LOCAL_MAXIMUM_GAMMA"] = FIRST_LOCAL_MAXIMUM_GAMMA;
	    codeMap["IMMUNE_THRESHOLD_GAMMA"] = IMMUNE_THRESHOLD_GAMMA;
	    codeMap["UPDATE_DENSITY_GAMMA"] = UPDATE_DENSITY_GAMMA;
	    codeMap["PARASITE_REPLICATION_GAMMA"] = PARASITE_REPLICATION_GAMMA;
        codeMap["VECTOR_LIFE_CYCLE_MODEL"] = VECTOR_LIFE_CYCLE_MODEL;
        codeMap["VECTOR_SIMPLE_MPD_MODEL"] = VECTOR_SIMPLE_MPD_MODEL;
        codeMap["MOLINEAUX_PAIRWISE_SAMPLE"] = MOLINEAUX_PAIRWISE_SAMPLE;
        ignoreOptions.insert("PROPHYLACTIC_DRUG_ACTION_MODEL");
        codeMap["VIVAX_SIMPLE_MODEL"] = VIVAX_SIMPLE_MODEL;
        codeMap["INDIRECT_MORTALITY_FIX"] = INDIRECT_MORTALITY_FIX;
        codeMap["VACCINE_GENOTYPE"] = VACCINE_GENOTYPE;
        codeMap["CFR_PF_USE_HOSPITAL"] = CFR_PF_USE_HOSPITAL;
        codeMap["HEALTH_SYSTEM_MEMORY_FIX"] = HEALTH_SYSTEM_MEMORY_FIX;
        codeMap["USE_EXACT_NV0_SOLVER"] = USE_EXACT_NV0_SOLVER;
	}
	
	OptionCodes operator[] (const string s) {
	    auto codeIt = codeMap.find (s);
	    if( codeIt != codeMap.end() ) return codeIt->second;
            if( ignoreOptions.count(s) ){
                if( CommandLine::option(CommandLine::DEPRECATION_WARNINGS) ){
                    cerr << "Deprecation warning: model option " << s << " is no longer used" << endl;
                }
                return IGNORE;
            }
            
            ostringstream msg;
            if( s == "PENALISATION_EPISODES" || s == "ATTENUATION_ASEXUAL_DENSITY" ){
                msg << "Please use schema 31 or earlier to use option "
                    << s << "; it is not available in later versions.";
            }else if( s == "IPTI_SP_MODEL" ){
                msg << "The IPT model is no longer available. Use MDA instead.";
            }else if( s == "REPORT_ONLY_AT_RISK" ){
                msg << "Option " << s << " has been replaced by <SurveyOptions onlyNewEpisode=\"true\">.";
            }else{
                msg << "Unrecognised model option: " << s;
            }
            throw xml_scenario_error(msg.str());
	}
	// reverse-lookup in map; only used for error/debug printing so efficiency is unimportant
	// doesn't ensure code is unique in the map either
	string toString (const OptionCodes code) {
	    for(auto codeIt = codeMap.begin(); codeIt != codeMap.end(); ++codeIt) {
		if (codeIt->second == code)
		    return codeIt->first;
	    }
	    throw TRACED_EXCEPTION_DEFAULT ("toString called with unknown code");	// this is a code error
	}
    };

    void ModelOptions::init(const SimClock& clock, scnXml::Model::ModelOptionsOptional& optionsElt,
                            util::ModelNameProvider mnp)
    {
    OptionCodeMap codeMap;

	// Set options to legacy defaults, then override any given in the XML file:
	options = getLegacyDefaultModelOptions();

    const util::ModelNames namedModelToUse = mnp.GetModelName();
    const bool useNamedModel = namedModelToUse != util::ModelNames::none;
    if (useNamedModel)
    {
        if (namedModelToUse == util::ModelNames::base)
        {
            // This completely discards any legacy model options set above.
            options = getBaseModelOptions();
        }
        else
        {
            throw util::xml_scenario_error(
                "No pre-set model options are available for the specified model name.");
        }
    }

    // Apply any user-specified overrides.
    if (optionsElt.present())
    {
	    const scnXml::OptionSet::OptionSequence& optSeq = optionsElt.get().getOption();
	    for(auto it = optSeq.begin(); it != optSeq.end(); ++it) {
                OptionCodes opt = codeMap[it->getName()];
                if( opt != IGNORE ) options[opt] = it->getValue();
	    }
    }
	
	// Print non-default model options:
	if (CommandLine::option (CommandLine::PRINT_MODEL_OPTIONS)) {
        const std::bitset<NUM_OPTIONS> defaultOptSet = getLegacyDefaultModelOptions();
	    cout << "Non-default model options:";
	    for(int i = 0; i < NUM_OPTIONS; ++i) {
		if (options[i] != defaultOptSet[i])
		    cout << "\t" << codeMap.toString(OptionCodes(i)) << "=" << options[i];
	    }
	    cout << endl;
	}
	
	// Test for incompatible options
	
	// An incompatibility triangle, listing options incompatible with each option
	// Doesn't check required versions
	// Originally from "description of variables for interface" excel sheet.
	// Declarative rule table (same directed edges as the former .set() chains):
	static const struct { OptionCodes option; initializer_list<OptionCodes> with; }
	    INCOMPAT_RULES[] = {
	    { NEGATIVE_BINOMIAL_MASS_ACTION, { LOGNORMAL_MASS_ACTION, TRANS_HET, COMORB_TRANS_HET, TRANS_TREAT_HET, TRIPLE_HET } },
	    { LOGNORMAL_MASS_ACTION, { TRANS_HET, COMORB_TRANS_HET, TRANS_TREAT_HET, TRIPLE_HET } },
	    // Note: MAX_DENS_CORRECTION is irrelevant when using new
	    // within-host models, but we don't mark it incompatible so that we can
	    // leave MAX_DENS_CORRECTION on by default.
	    { DUMMY_WITHIN_HOST_MODEL, { PENNY_WITHIN_HOST_MODEL, EMPIRICAL_WITHIN_HOST_MODEL, MOLINEAUX_WITHIN_HOST_MODEL, MEAN_DURATION_GAMMA, FIRST_LOCAL_MAXIMUM_GAMMA, PARASITE_REPLICATION_GAMMA, IMMUNE_THRESHOLD_GAMMA, UPDATE_DENSITY_GAMMA } },
	    { EMPIRICAL_WITHIN_HOST_MODEL, { MOLINEAUX_WITHIN_HOST_MODEL, PENNY_WITHIN_HOST_MODEL, MEAN_DURATION_GAMMA, FIRST_LOCAL_MAXIMUM_GAMMA, PARASITE_REPLICATION_GAMMA, IMMUNE_THRESHOLD_GAMMA, UPDATE_DENSITY_GAMMA } },
	    { MOLINEAUX_WITHIN_HOST_MODEL, { PENNY_WITHIN_HOST_MODEL, IMMUNE_THRESHOLD_GAMMA, UPDATE_DENSITY_GAMMA } },
	    { PENNY_WITHIN_HOST_MODEL, { MEAN_DURATION_GAMMA, FIRST_LOCAL_MAXIMUM_GAMMA, PARASITE_REPLICATION_GAMMA } },
	    { NON_MALARIA_FEVERS, { MUELLER_PRESENTATION_MODEL } },
	    { TRANS_HET, { COMORB_TRANS_HET, TRANS_TREAT_HET, COMORB_TREAT_HET, TRIPLE_HET } },
	    { COMORB_HET, { COMORB_TRANS_HET, TRANS_TREAT_HET, COMORB_TREAT_HET, TRIPLE_HET } },
	    { TREAT_HET, { COMORB_TRANS_HET, TRANS_TREAT_HET, COMORB_TREAT_HET, TRIPLE_HET } },
	    { COMORB_TRANS_HET, { TRANS_TREAT_HET, COMORB_TREAT_HET, TRIPLE_HET } },
	    { TRANS_TREAT_HET, { COMORB_TREAT_HET, TRIPLE_HET } },
	    { COMORB_TREAT_HET, { TRIPLE_HET } },
	    { MOLINEAUX_PAIRWISE_SAMPLE, { FIRST_LOCAL_MAXIMUM_GAMMA, MEAN_DURATION_GAMMA } },
	    { VIVAX_SIMPLE_MODEL, { DUMMY_WITHIN_HOST_MODEL, EMPIRICAL_WITHIN_HOST_MODEL, MOLINEAUX_WITHIN_HOST_MODEL, PENNY_WITHIN_HOST_MODEL } },
	};
	bitset<NUM_OPTIONS> incompatibilities[NUM_OPTIONS];	// all default to 0
	for( const auto& rule : INCOMPAT_RULES ) for( OptionCodes w : rule.with ) incompatibilities[rule.option].set (w);
	
	for(size_t i = 0; i < NUM_OPTIONS; ++i) {
	    if (options [i] && (options & incompatibilities[i]).any()) {
		ostringstream msg;
		msg << "Incompatible model options: " << codeMap.toString(OptionCodes(i)) << "=" << options[i]
			<< " is incompatible with flags:";
		bitset<NUM_OPTIONS> incompat = (options & incompatibilities[i]);
		for(int j = 0; j < NUM_OPTIONS; ++j) {
		    if (incompat[j])
			msg << "\t" << codeMap.toString(OptionCodes(j)) << "=" << options[j];
		}
		throw xml_scenario_error (msg.str());
	    }
	}
	
        // Required options (above table can't check these):
        if (options[INNATE_MAX_DENS] && !options[MAX_DENS_CORRECTION])
            throw xml_scenario_error ("INNATE_MAX_DENS requires MAX_DENS_CORRECTION");
        if( !options[MOLINEAUX_WITHIN_HOST_MODEL] && (
            options[FIRST_LOCAL_MAXIMUM_GAMMA] ||
            options[MEAN_DURATION_GAMMA] ||
            options[PARASITE_REPLICATION_GAMMA] ) )
            throw xml_scenario_error( "Molineaux model option used without MOLINEAUX_WITHIN_HOST_MODEL option" );
        if( !options[PENNY_WITHIN_HOST_MODEL] && (
            options[IMMUNE_THRESHOLD_GAMMA] ||
            options[UPDATE_DENSITY_GAMMA] ) )
            throw xml_scenario_error( "Penny model option used without PENNY_WITHIN_HOST_MODEL option" );
        
        if( clock.one_TS() == sim::fromDays(5) ){
            // 5 day TS is okay; some tests specific to this TS:
            bitset<NUM_OPTIONS> require1DayTS;
            require1DayTS
                .set( CLINICAL_EVENT_SCHEDULER );
            
            for(size_t i = 0; i < NUM_OPTIONS; ++i) {
                if (options [i] && require1DayTS[i]) {
                    ostringstream msg;
                    msg << "Model option " << codeMap.toString(OptionCodes(i)) << " is only compatible with a 1-day time step.";
                    throw xml_scenario_error (msg.str());
                }
            }
        }else if( clock.one_TS() == sim::fromDays(1) ){
            // 1 day TS is also okay
        }else{
            ostringstream msg;
            msg << "Time step set to " << clock.one_TS() << " days but only 1 and 5 days are supported.";
            throw xml_scenario_error (msg.str());
        }
    }

} }
