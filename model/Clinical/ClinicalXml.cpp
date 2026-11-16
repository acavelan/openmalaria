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

#include "Clinical/ClinicalXml.h"
#include "Clinical/ClinicalParameters.h"
#include "ScenarioXml.h"
#include "Host/WithinHost/WHVivax.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"

namespace OM::xml {
using Clinical::ClinicalParameters;
using Clinical::CMDecisionTree;
using namespace util;
namespace {
WithinHost::Treatment::Stages stageFromString( const std::string& str ){
    if( str == "liver" ) return WithinHost::Treatment::LIVER;
    if( str == "blood" ) return WithinHost::Treatment::BLOOD;
    if( str == "both" ) return WithinHost::Treatment::BOTH;
    throw util::format_error( std::string("stage must be liver, blood or both, not ").append(str) );
}
}

WithinHost::Treatment readTreatment(const scnXml::TreatmentOption& elt,
        const interventions::Interventions& definitions) {
    WithinHost::Treatment result{sim::zero(), sim::zero(),
        interventions::TriggeredDeployments(elt, definitions)};
    for( auto it = elt.getClearInfections().begin(), end = elt.getClearInfections().end();
        it != end; ++it )
    {
        try{
            //NOTE: if changing XSD, this should not be called "timesteps" or have a default unit
            SimTime len = UnitParse::readShortDuration( it->getTimesteps(), UnitParse::STEPS );
            if( len < -sim::oneTS() || len == sim::zero() ){
                throw util::format_error( "timesteps must be ≥ 1 or have special value -1" );
            }
            auto stage = stageFromString( it->getStage() );
            if( stage & WithinHost::Treatment::LIVER ){
                if( result.timeLiver != sim::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of liver stage effect" );
                result.timeLiver = len;
            }
            if( stage & WithinHost::Treatment::BLOOD ){
                if( result.timeBlood != sim::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of blood stage effect" );
                result.timeBlood = len;
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string(".../clearInfections/timesteps: ").append(e.message()) );
        }
    }
    return result;
}

namespace {
enum CaseType { FirstLine, SecondLine };

double getHealthSystemACRByName (const scnXml::TreatmentDetails& td, const string& drug)
{
    double val;
    if (drug == "CQ")
        val = td.getCQ().present() ? td.getCQ().get().getValue() : 0;
    else if (drug == "SP")
        val = td.getSP().present() ? td.getSP().get().getValue() : 0;
    else if (drug == "AQ")
        val = td.getAQ().present() ? td.getAQ().get().getValue() : 0;
    else if (drug == "SPAQ")
        val = td.getSPAQ().present() ? td.getSPAQ().get().getValue() : 0;
    else if (drug == "ACT")
        val = td.getACT().present() ? td.getACT().get().getValue() : 0;
    else if (drug == "QN")
        val = td.getQN().present() ? td.getQN().get().getValue() : 0;
    else if (drug == "selfTreatment")
        val = td.getSelfTreatment().getValue();
    else {
        throw util::xml_scenario_error ("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
    }

    if( !(val >= 0.0 && val <= 1.0) ){
        throw util::xml_scenario_error (string("healthSystem initialACR/compliance/nonCompliersEffective: ").append(drug).append(" value must be in range [0,1]"));
    }
    return val;
}

// note: "CQOptional" is the same type as for other drugs
const scnXml::TreatmentActions::CQOptional& getHealthSystemTreatmentOptional(
    const scnXml::TreatmentActions& treatActs, const string& drug )
{
    if( drug == "CQ" ) return treatActs.getCQ();
    if( drug == "SP" ) return treatActs.getSP();
    if( drug == "AQ" ) return treatActs.getAQ();
    if( drug == "SPAQ" ) return treatActs.getSPAQ();
    if( drug == "ACT" ) return treatActs.getACT();
    if( drug == "QN" ) return treatActs.getQN();
    throw util::xml_scenario_error ("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
}
WithinHost::Treatment getHealthSystemTreatmentByName( const scnXml::TreatmentActions& treatActs, const string& drug, const interventions::Interventions& definitions )
{
    const scnXml::TreatmentActions::CQOptional& elt = getHealthSystemTreatmentOptional( treatActs, drug );
    if( !elt.present() )
        throw util::xml_scenario_error( string("healthSystem.treatmentActions: description required for ").append(drug) );
    return readTreatment( elt.get(), definitions );
}


// ———  static, init  ———

void setImmediateOutcomes (const scnXml::HSImmediateOutcomes& hsDescription, const interventions::Interventions& definitions, ClinicalParameters& config) {
    const string &firstLine = hsDescription.getDrugRegimen().getFirstLine(),
        &secondLine = hsDescription.getDrugRegimen().getSecondLine(),
        &inpatient = hsDescription.getDrugRegimen().getInpatient();

    // ———  calculate probability of getting treatment  ———

    double accessUCOfficial1 = hsDescription.getPSeekOfficialCareUncomplicated1().getValue();
    double accessUCOfficial2 = hsDescription.getPSeekOfficialCareUncomplicated2().getValue();
    // Note: this asymmetry is historical, and probably matters little:
    config.accessUCSelfTreat[FirstLine] = hsDescription.getPSelfTreatUncomplicated().getValue();
    config.accessUCSelfTreat[SecondLine] = 0.0;
    config.accessUCAny[FirstLine] = accessUCOfficial1 + config.accessUCSelfTreat[FirstLine];
    config.accessUCAny[SecondLine] = accessUCOfficial1 + config.accessUCSelfTreat[SecondLine];
    config.accessSevere = hsDescription.getPSeekOfficialCareSevere().getValue();

    if( accessUCOfficial1 < 0.0 || config.accessUCSelfTreat[FirstLine] < 0.0 || config.accessUCAny[FirstLine] > 1.0 ||
        accessUCOfficial2 < 0.0 || config.accessUCSelfTreat[SecondLine] < 0.0 || config.accessUCAny[SecondLine] > 1.0 ||
        config.accessSevere < 0.0 || config.accessSevere > 1.0 )
    {
        throw util::xml_scenario_error ("healthSystem: "
            "pSeekOfficialCareUncomplicated1 and pSelfTreatUncomplicated must be"
            " at least 0 and their sum must be at most 1, and"
            "pSeekOfficialCareUncomplicated2 and pSeekOfficialCareSevere must "
            "be in range [0,1]");
    }

    // ———  calculate probability of clearing parasites  ———

    const double complianceFirstLine = getHealthSystemACRByName (hsDescription.getCompliance(), firstLine);
    const double complianceSecondLine = getHealthSystemACRByName (hsDescription.getCompliance(), secondLine);

    const double cureRateFirstLine = getHealthSystemACRByName (hsDescription.getInitialACR(), firstLine);
    const double cureRateSecondLine = getHealthSystemACRByName (hsDescription.getInitialACR(), secondLine);

    const double nonCompliersEffectiveFirstLine = getHealthSystemACRByName (hsDescription.getNonCompliersEffective(), firstLine);
    const double nonCompliersEffectiveSecondLine = getHealthSystemACRByName (hsDescription.getNonCompliersEffective(), secondLine);

    const double complianceSelfTreatment = hsDescription.getCompliance().getSelfTreatment().getValue();
    const double cureRateSelfTreatment = hsDescription.getInitialACR().getSelfTreatment().getValue();
    if( !(
        complianceSelfTreatment >= 0.0 && complianceSelfTreatment <= 1.0
        && cureRateSelfTreatment >= 0.0 && cureRateSelfTreatment <= 1.0
    ) ){
        throw util::xml_scenario_error ("healthSystem initialACR/compliance/nonCompliersEffective: pSelfTreatment must be in range [0,1]");
    }


    config.cureRateUCOfficial[FirstLine] = complianceFirstLine * cureRateFirstLine
            + (1 - complianceFirstLine) * nonCompliersEffectiveFirstLine;
    config.cureRateUCOfficial[SecondLine] = complianceSecondLine * cureRateSecondLine
            + (1 - complianceSecondLine) * nonCompliersEffectiveSecondLine;
    config.cureRateUCSelfTreat[FirstLine] = complianceSelfTreatment * cureRateSelfTreatment
            + (1 - complianceSelfTreatment) * nonCompliersEffectiveFirstLine;
    config.cureRateUCSelfTreat[SecondLine] = 0.0;      // as with access, this is a historical asymmetry

    config.cureRateSevere = getHealthSystemACRByName (hsDescription.getInitialACR(), inpatient);

    config.treatmentUC[FirstLine] = getHealthSystemTreatmentByName(hsDescription.getTreatmentActions(), firstLine, definitions);

    if( secondLine == firstLine ) config.treatmentUC[SecondLine] = config.treatmentUC[FirstLine];
    else config.treatmentUC[SecondLine] = getHealthSystemTreatmentByName(hsDescription.getTreatmentActions(), secondLine, definitions);

    if( inpatient == firstLine ) config.treatmentSevere = config.treatmentUC[FirstLine];
    else if( inpatient == secondLine ) config.treatmentSevere = config.treatmentUC[SecondLine];
    else config.treatmentSevere = getHealthSystemTreatmentByName(hsDescription.getTreatmentActions(), inpatient, definitions);

    config.useDiagnosticUC = hsDescription.getUseDiagnosticUC();

    if( ModelOptions::option(util::VIVAX_SIMPLE_MODEL) ){
        WithinHost::WHVivax::setHSParameters(
            hsDescription.getLiverStageDrug().present() ?
            &hsDescription.getLiverStageDrug().get() : 0 );
    }else if( hsDescription.getLiverStageDrug().present() ){
        throw util::xml_scenario_error( "health-system's liverStageDrug element only supported by vivax" );
    }
}


void setDecisionTree5Day(const scnXml::HSDT5Day& hsDescription, const interventions::Interventions& definitions, ClinicalParameters& config){
    double accessUCOfficial1 = hsDescription.getPSeekOfficialCareUncomplicated1().getValue();
    double accessUCOfficial2 = hsDescription.getPSeekOfficialCareUncomplicated2().getValue();
    // Note: this asymmetry is historical, and probably matters little:
    config.accessUCSelfTreat[FirstLine] = hsDescription.getPSelfTreatUncomplicated().getValue();
    config.accessUCSelfTreat[SecondLine] = 0.0;
    config.accessUCAny[FirstLine] = accessUCOfficial1 + config.accessUCSelfTreat[FirstLine];
    config.accessUCAny[SecondLine] = accessUCOfficial1 + config.accessUCSelfTreat[SecondLine];
    config.accessSevere = hsDescription.getPSeekOfficialCareSevere().getValue();

    if( accessUCOfficial1 < 0.0 || config.accessUCSelfTreat[FirstLine] < 0.0 || config.accessUCAny[FirstLine] > 1.0 ||
        accessUCOfficial2 < 0.0 || config.accessUCSelfTreat[SecondLine] < 0.0 || config.accessUCAny[SecondLine] > 1.0 ||
        config.accessSevere < 0.0 || config.accessSevere > 1.0 )
    {
        throw util::xml_scenario_error ("healthSystem: "
            "pSeekOfficialCareUncomplicated1 and pSelfTreatUncomplicated must be"
            " at least 0 and their sum must be at most 1, and"
            "pSeekOfficialCareUncomplicated2 and pSeekOfficialCareSevere must "
            "be in range [0,1]");
    }

    config.treeUCOfficial = &CMDecisionTree::create( hsDescription.getTreeUCOfficial(), true, definitions, config);
    config.treeUCSelfTreat = &CMDecisionTree::create( hsDescription.getTreeUCSelfTreat(), true, definitions, config);

    config.cureRateSevere = hsDescription.getCureRateSevere().getValue();
    config.treatmentSevere = readTreatment( hsDescription.getTreatmentSevere(), definitions );

    if( ModelOptions::option(util::VIVAX_SIMPLE_MODEL) ){
        WithinHost::WHVivax::setHSParameters(
            hsDescription.getLiverStageDrug().present() ?
            &hsDescription.getLiverStageDrug().get() : 0 );
    }else if( hsDescription.getLiverStageDrug().present() ){
        throw util::xml_scenario_error( "health-system's liverStageDrug element only supported by vivax" );
    }
}


void initEventScheduler( const Parameters& parameters, const scnXml::Clinical& clinical, ClinicalParameters& config )
{
    if( sim::oneTS() != sim::oneDay() )
        throw util::xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day time step.");

    config.opt_non_malaria_fevers = util::ModelOptions::option( util::NON_MALARIA_FEVERS );

    config.alpha = exp( -parameters[Parameter::CFR_NEG_LOG_ALPHA] );
    if( !(0.0<=config.alpha && config.alpha<=1.0) ){
        throw util::xml_scenario_error(
            "Clinical outcomes: propDeathsFirstDay should be within range [0,1]");
    }
    config.neg_v = -parameters[Parameter::CFR_SCALE_FACTOR];

    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc2 =
            clinical.getNonMalariaFevers().get();
        if( !nmfDesc2.getPrNeedTreatmentMF().present() ||
            !nmfDesc2.getPrNeedTreatmentNMF().present() ){
            throw util::xml_scenario_error( "prNeedTreatmentMF and "
                "prNeedTreatmentNMF elements required in model->clinical->"
                "NonMalariaFevers" );
        }
        xml::setAgeGroups(config.NMF_need_antibiotic, nmfDesc2.getPrNeedTreatmentNMF().get(), "prNeedTreatmentNMF" );
        xml::setAgeGroups(config.MF_need_antibiotic, nmfDesc2.getPrNeedTreatmentMF().get(), "prNeedTreatmentMF" );
    }
}

void setEventScheduler(const scnXml::HSEventScheduler& esData, const interventions::Interventions& definitions, ClinicalParameters& config) {
    config.escm_uncomplicated = &Clinical::CMDecisionTree::create(esData.getUncomplicated(), true, definitions, config);
    config.escm_complicated = &Clinical::CMDecisionTree::create(esData.getComplicated(), false, definitions, config);
    const scnXml::ClinicalOutcomes& coData = esData.getClinicalOutcomes();

    config.maxUCSeekingMemory = sim::fromDays(coData.getMaxUCSeekingMemory());
    config.uncomplicatedCaseDuration = sim::fromDays(coData.getUncomplicatedCaseDuration());
    config.complicatedCaseDuration = sim::fromDays(coData.getComplicatedCaseDuration());
    config.extraDaysAtRisk = sim::fromDays(coData.getComplicatedRiskDuration()) - config.complicatedCaseDuration;
    if( config.uncomplicatedCaseDuration < sim::fromDays(1)
	|| config.complicatedCaseDuration < sim::fromDays(1)
	|| config.maxUCSeekingMemory < sim::zero()
	|| config.extraDaysAtRisk + config.complicatedCaseDuration < sim::fromDays(1) // at risk at least 1 day
	|| config.extraDaysAtRisk > sim::zero()        // at risk longer than case duration
    ){
	throw util::xml_scenario_error(
            "Clinical outcomes: constraints on case/risk/memory duration not met (see documentation)");
    }

    config.cumDailyPrImmUCTS.reserve( coData.getDailyPrImmUCTS().size() );
    double cumP = 0.0;
    for( auto it = coData.getDailyPrImmUCTS().begin(); it != coData.getDailyPrImmUCTS().end(); ++it ){
        cumP += *it;
        config.cumDailyPrImmUCTS.push_back( cumP );
    }
    if( cumP < 0.999 || cumP > 1.001 ){
        throw util::xml_scenario_error( "Event scheduler: dailyPrImmUCTS seq must add up to 1" );
    }
    config.cumDailyPrImmUCTS.back() = 1.0;

    config.caseFatalityRate.scale( config.alpha );

    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !esData.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of healthSystem->EventScheduler required");
        }
        const scnXml::HSESNMF& nmfDesc = esData.getNonMalariaFevers().get();

        double pT = nmfDesc.getPrTreatment();
        config.logOddsAbBase = log( pT / (1.0 - pT) );
        config.logOddsAbNegTest = log(nmfDesc.getEffectNegativeTest());
        config.logOddsAbPosTest = log(nmfDesc.getEffectPositiveTest());
        config.logOddsAbNeed = log(nmfDesc.getEffectNeed());
        config.logOddsAbInformal = log(nmfDesc.getEffectInformal());
        config.oneMinusEfficacyAb = 1.0 - nmfDesc.getTreatmentEfficacy();
        xml::setAgeGroups(config.severeNmfMortality, nmfDesc.getCFR(), "CFR" );
    }
}


}

void initClinical( const Parameters& parameters, const scnXml::Scenario& scenario, ClinicalParameters& config ) {
    const scnXml::Clinical& clinical = scenario.getModel().getClinical();
    try{
        config.indirectMortBugfix = util::ModelOptions::option (util::INDIRECT_MORTALITY_FIX);
        //NOTE: if changing XSD, this should not have a default unit:
        config.healthSystemMemory = UnitParse::readShortDuration(clinical.getHealthSystemMemory(), UnitParse::STEPS);
        config.oddsRatioThreshold = exp( parameters[Parameter::LOG_ODDS_RATIO_CF_COMMUNITY] );
    }catch( const util::format_error& e ){
        throw util::xml_scenario_error( string("model/clinical/healthSystemMemory: ").append(e.message()) );
    }

    if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER)){
        config.kind = ClinicalParameters::Kind::EventScheduler;
        initEventScheduler( parameters, clinical, config );
    }else{
        if( scenario.getHealthSystem().getImmediateOutcomes().present() ){
            config.kind = ClinicalParameters::Kind::ImmediateOutcomes;

            if( util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS) ){
                cerr << "Deprecation warning: healthSystem: use of ImmediateOutcomes can be replaced by the more flexible DecisionTree5Day (optional)" << endl;
            }
        }
        // else: decision tree 5 day

        config.cfr_pf_use_hospital = util::ModelOptions::option(util::CFR_PF_USE_HOSPITAL);
    }
}

void setHealthSystem( const scnXml::HealthSystem& healthSystem, const interventions::Interventions& definitions, ClinicalParameters& config){
    xml::setAgeGroups(config.caseFatalityRate, healthSystem.getCFR(), "CFR" );
    xml::setAgeGroups(config.pSequelaeInpatient, healthSystem.getPSequelaeInpatient(), "pSequelaeInpatient" );
    if( config.kind == ClinicalParameters::Kind::EventScheduler ){
        if( !healthSystem.getEventScheduler().present() ){
            throw util::xml_scenario_error ("Expected EventScheduler "
                "section in healthSystem data (initial or intervention)");
        }
        setEventScheduler(healthSystem.getEventScheduler().get(), definitions, config);
    }else if( config.kind == ClinicalParameters::Kind::ImmediateOutcomes ){
        if ( !healthSystem.getImmediateOutcomes().present() ){
            throw util::xml_scenario_error ("Expected ImmediateOutcomes "
                "section in healthSystem data (initial or intervention)");
        }
        setImmediateOutcomes(healthSystem.getImmediateOutcomes().get(), definitions, config);
    }else{
        if( !healthSystem.getDecisionTree5Day().present() ){
            throw util::xml_scenario_error ("Expected DecisionTree5Day "
                "section in healthSystem data (initial or intervention)");
        }
        setDecisionTree5Day(healthSystem.getDecisionTree5Day().get(), definitions, config);
    }
}

}
