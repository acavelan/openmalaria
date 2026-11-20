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

// model/xml/ conversion layer — pharmacology module: drug types/phenotypes,
// treatment schedules and dosage tables from scnXml. Moved verbatim from
// PkPd/LSTMDrugType.cpp, PkPd/LSTMTreatments.cpp and PkPd/LSTMModel.cpp
// (the *.load() members became file-local read_* functions; all their
// callers live in this file). Declarations for the member entry points stay
// in the PkPd headers (scnXml forward-declared).

#include "Clock.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "PkPd/LSTMDrugType.h"
#include "PkPd/LSTMDrug.h"
#include "PkPd/LSTMMedicate.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include <schema/pharmacology.h>
#include <schema/scenario.h>

#include <cmath>
#include <iomanip>

using namespace std;

namespace OM {
namespace PkPd {
using namespace OM::util;
using WithinHost::Genotypes;

LSTMDrugPD::LSTMDrugPD(const scnXml::Phenotype& phenotype) {
    n = phenotype.getSlope();
    IC50 = util::createSampler(phenotype.getIC50());
    V = phenotype.getMax_killing_rate();
}

void LSTMDrugType::init(RunContext& ctx, const scnXml::Pharmacology::DrugsType& drugData) {
    DrugRegistry& reg = ctx.drugs;
    vector<LSTMDrugType>& drugTypes = reg.drugTypes;
    map<string, size_t>& drugTypeNames = reg.drugTypeNames;
    drugTypes.reserve(drugData.getDrug().size());
    for (const scnXml::PKPDDrug& drug : drugData.getDrug()) {
        const string& abbrev = drug.getAbbrev();
        // Check drug doesn't already exist
        if (drugTypeNames.find(abbrev) != drugTypeNames.end())
            throw TRACED_EXCEPTION_DEFAULT(string("Drug added twice: ").append(abbrev));

        size_t i = drugTypes.size();
        drugTypes.emplace_back(ctx, i, drug);
        drugTypeNames[abbrev] = i;
    }
}

LSTMDrugType::LSTMDrugType(RunContext& ctx, size_t index, const scnXml::PKPDDrug& drugData)
    : index(index), metabolite(0), negligible_concentration(numeric_limits<double>::quiet_NaN()),
      neg_m_exp(numeric_limits<double>::quiet_NaN()), mwr(numeric_limits<double>::quiet_NaN()),
      ic50_log_corr(numeric_limits<double>::quiet_NaN()), ic50_corr_factor(numeric_limits<double>::quiet_NaN()) {
    // ———  PK parameters  ———
    const scnXml::PK& pk = drugData.getPK();
    negligible_concentration = pk.getNegligible_concentration();
    if (pk.getHalf_life().present()) {
        if (pk.getK().present() || pk.getM_exponent().present()) {
            throw util::xml_scenario_error(
                "PK data must specify one of half_life or (k, m_exponent); it specifies both");
        }
        elimination_rate = Sampler::fromMeanCV(Sampler::Kind::Lognormal, log(2.0) / pk.getHalf_life().get(), 0.0);
        neg_m_exp = 0.0; // no dependence on body mass
    } else {
        if (!(pk.getK().present() && pk.getM_exponent().present())) {
            throw util::xml_scenario_error("PK data must include either half_life or (k and m_exponent)");
        }
        elimination_rate = util::createSampler(pk.getK().get());
        neg_m_exp = -pk.getM_exponent().get();
    }
    vol_dist = util::createSampler(pk.getVol_dist());
    if (pk.getCompartment2().present()) {
        k12 = util::createSampler(pk.getCompartment2().get().getK12());
        k21 = util::createSampler(pk.getCompartment2().get().getK21());
        if (pk.getCompartment3().present()) {
            k13 = util::createSampler(pk.getCompartment3().get().getK13());
            k31 = util::createSampler(pk.getCompartment3().get().getK31());
        } else {
            // 2-compartment model: use 3-compartment code with these parameters set to zero
            // lognormal Sampler supports 0,0 as a special case
            k13 = Sampler::fromMeanCV(Sampler::Kind::Lognormal, 0.0, 0.0);
            k31 = Sampler::fromMeanCV(Sampler::Kind::Lognormal, 0.0, 0.0);
        }
    } else if (pk.getCompartment3().present()) {
        throw util::xml_scenario_error("PK model specifies parameters for "
                                       "compartment3 without compartment2");
    }
    if (pk.getConversion().present()) {
        if (k12) {
            throw util::xml_scenario_error("PK conversion model is incompatible with 2/3-compartment model");
        }
        const scnXml::Conversion& conv = pk.getConversion().get();
        try {
            metabolite = ctx.drugs.findDrug(conv.getMetabolite());
        } catch (util::xml_scenario_error& e) {
            throw util::xml_scenario_error(
                "PK: metabolite drug not found; metabolite must be defined *before* parent drug!");
        }
        conversion_rate = util::createSampler(conv.getRate());
        mwr = conv.getMolRatio();
        double correlation = conv.getIC50_log_correlation();
        if (correlation < 0.0 || correlation > 1.0) {
            throw util::xml_scenario_error("PK conversion model: IC50 correlation must be between 0 and 1");
        }
        ic50_log_corr = correlation;
        ic50_corr_factor = sqrt(1.0 - correlation * correlation);
    }
    if (pk.getK_a().present()) {
        if (!k12 && !conversion_rate) {
            throw util::xml_scenario_error("PK models only allow an "
                                           "absorption rate parameter (k_a) when compartment2 or "
                                           " conversion parameters are present");
        }
        absorption_rate = util::createSampler(pk.getK_a().get());
    } else {
        if (k12) {
            throw util::xml_scenario_error("PK models require an absorption "
                                           "rate parameter (k_a) when compartment2 is present");
        }
    }

    // ———  PD parameters  ———
    const scnXml::PD::PhenotypeSequence& pd = drugData.getPD().getPhenotype();
    assert(pd.size() > 0); // required by XSD

    set<string> loci_per_phenotype;
    string first_phenotype_name;
    size_t n_phenotypes = pd.size();
    PD.reserve(n_phenotypes);
    // per phenotype (first index), per locus-of-restriction (second list), a list of alleles
    vector<vector<vector<uint32_t>>> phenotype_restrictions;
    phenotype_restrictions.reserve(n_phenotypes);
    for (size_t i = 0; i < n_phenotypes; ++i) {
        const scnXml::Phenotype& phenotype = pd[i];
        // phenotype name for error messages: name attribute, else index
        auto pheno_name = [&phenotype, i]() {
            return phenotype.getName().present() ? string(phenotype.getName().get()) : to_string(i);
        };
        if (i == 0) {
            first_phenotype_name = pheno_name();
            for (const scnXml::Restriction& restriction : phenotype.getRestriction()) {
                loci_per_phenotype.insert(restriction.getOnLocus());
            }
        } else {
            set<string> expected_loci = loci_per_phenotype;
            for (const scnXml::Restriction& restriction : phenotype.getRestriction()) {
                size_t n = expected_loci.erase(restriction.getOnLocus());
                if (n == 0) {
                    throw util::xml_scenario_error("pharmacology/drugs/drug/pd/phenotype/restriction/onLocus:"
                                                   "locus " +
                                                   string(restriction.getOnLocus()) +
                                                   " included in restriction of phenotype " + string(pheno_name()) +
                                                   " but not " + string(first_phenotype_name));
                }
            }
            if (!expected_loci.empty()) {
                throw util::xml_scenario_error("pharmacology/drugs/drug/pd/phenotype/restriction/onLocus:"
                                               "locus " +
                                               string(*expected_loci.begin()) +
                                               " included in restriction of phenotype " + string(first_phenotype_name) +
                                               " but not " + string(pheno_name()));
            }
        }
        map<string, size_t> loci;
        vector<vector<uint32_t>> loc_alleles;
        for (const scnXml::Restriction& restriction : phenotype.getRestriction()) {
            uint32_t allele = ctx.genotypes.findAlleleCode(restriction.getOnLocus(), restriction.getToAllele());
            if (allele == numeric_limits<uint32_t>::max()) {
                throw util::xml_scenario_error("phenotype has "
                                               "restriction on locus " +
                                               string(restriction.getOnLocus()) + " allele " +
                                               string(restriction.getToAllele()) +
                                               " but this locus/allele "
                                               "has not been defined in parasiteGenetics section");
            }
            auto iter = loci.find(restriction.getOnLocus());
            if (iter == loci.end()) {
                loci[restriction.getOnLocus()] = loc_alleles.size();
                loc_alleles.push_back(vector<uint32_t>(1, allele));
            } else {
                loc_alleles[iter->second].push_back(allele);
            }
        }
        phenotype_restrictions.push_back(loc_alleles);
        PD.push_back(LSTMDrugPD(pd[i]));
    }

    if (loci_per_phenotype.size() == 0) {
        if (pd.size() > 1) {
            throw util::xml_scenario_error("pharmacology/drugs/drug/pd/phenotype:"
                                           " restrictions required when num. phenotypes > 1");
        } else {
            // All genotypes map to phenotype 0
            genotype_mapping.assign(ctx.genotypes.genotypes.size(), 0);
        }
    } else {
        const vector<Genotypes::Genotype>& genotypes = ctx.genotypes.genotypes;
        genotype_mapping.assign(genotypes.size(), 0);
        for (size_t j = 0; j < genotypes.size(); ++j) {
            uint32_t phenotype = numeric_limits<uint32_t>::max();
            for (size_t i = 0; i < n_phenotypes; ++i) {
                // genotype matches this phenotype when, for every locus among
                // the phenotype restriction rules, one allele matches a
                // genotype allele
                bool match = true;
                for (const vector<uint32_t>& loc_alleles : phenotype_restrictions[i]) {
                    bool match_for_locus = false;
                    for (uint32_t restrict_allele : loc_alleles) {
                        if (genotypes[j].alleles.count(restrict_allele) > 0) {
                            assert(!match_for_locus); // shouldn't be two matches
                            match_for_locus = true;
                        }
                    }
                    if (!match_for_locus) {
                        match = false;
                        break; // can skip rest of rules
                    }
                }
                if (match) {
                    if (phenotype == numeric_limits<uint32_t>::max()) {
                        phenotype = i;
                    } else {
                        // We could try to convert genotype's allele codes to
                        // locus/allele names, but requires several lookups
                        throw util::xml_scenario_error("phenotype restrictions "
                                                       "not restrictive enough: multiple phenotypes match "
                                                       "some genotype");
                    }
                }
            }
            if (phenotype == numeric_limits<uint32_t>::max()) {
                // We could try to convert genotype's allele codes to
                // locus/allele names, but requires several lookups
                throw util::xml_scenario_error("phenotype restrictions too"
                                               "restrictive: no phenotype matching some genotype");
            }
            genotype_mapping[j] = phenotype;
        }

        if (util::CommandLine::option(util::CommandLine::PRINT_GENOTYPES)) {
            // Debug output to see which genotypes correspond to phenotypes
            cout << endl;
            cout << "Phenotype mapping: " << endl << "----------" << endl;
            cout << "|" << setw(8) << "genotype" << "|" << setw(14) << "phenotype" << "|" << endl;
            cout << "|" << setw(8) << "--------" << "|" << setw(14) << "-------------" << "|" << endl;

            uint32_t genotype = 0;
            stringstream phenotypeName;
            for (auto phen = genotype_mapping.begin(); phen != genotype_mapping.end(); ++phen) {
                phenotypeName.str("");
                if (pd[*phen].getName().present()) {
                    phenotypeName << pd[*phen].getName().get();
                } else {
                    phenotypeName << "no: " << *phen;
                }
                cout << "|" << setw(8) << (genotype * 100000) << "|" << setw(14) << phenotypeName.str() << "|" << endl;
                genotype += 1;
            }
        }
    }
}

// ———  treatments (origins: LSTMTreatments.cpp; *.load members -> read_*)  ———

void MedicateData::load(DrugRegistry& drugs, const scnXml::PKPDMedication& med) {
    drug = drugs.findDrug(med.getDrug());
    qty = med.getMg();
    time = med.getHour() / 24.0; // convert from hours to days
}

void read_schedule(Schedule& out, DrugRegistry& drugs, const scnXml::PKPDSchedule::MedicateSequence& seq) {
    out.medications.resize(seq.size());
    size_t i = 0;
    for (const scnXml::PKPDMedication& med : seq) {
        out.medications[i].load(drugs, med);
        i += 1;
    }
}

void read_dosage_table(DosageTable& out, const xsd::cxx::tree::sequence<scnXml::PKPDDosageRange>& seq,
                       bool isBodyMass) {
    out.useMass = isBodyMass;
    out.multMassKg = false;
    double lastMult = 0.0, lastAge = numeric_limits<double>::quiet_NaN();
    for (const scnXml::PKPDDosageRange& age : seq) {
        if (lastAge != lastAge) {
            if (age.getLowerbound() != 0.0)
                throw util::xml_scenario_error("dosage table must have first lower bound equal 0");
        } else {
            if (age.getLowerbound() <= lastAge) {
                throw util::xml_scenario_error("dosage table must list age groups in increasing order");
            }
            out.table[age.getLowerbound()] = lastMult;
        }
        lastMult = age.getDose_mult();
        lastAge = age.getLowerbound();
    }
    out.table[numeric_limits<double>::infinity()] = lastMult;
}

void LSTMTreatments::init(RunContext& ctx, const scnXml::Treatments& data) {
    TreatmentRegistry& reg = ctx.treatments;
    vector<Schedule>& schedules = reg.schedules;
    map<string, size_t>& scheduleNames = reg.scheduleNames;
    vector<DosageTable>& dosages = reg.dosages;
    map<string, size_t>& dosagesNames = reg.dosagesNames;
    schedules.resize(data.getSchedule().size());
    size_t i = 0;
    for (const scnXml::PKPDSchedule& schedule : data.getSchedule()) {
        read_schedule(schedules[i], ctx.drugs, schedule.getMedicate());
        scheduleNames[schedule.getName()] = i;
        i += 1;
    }

    dosages.resize(data.getDosages().size());
    i = 0;
    for (const scnXml::PKPDDosages& elt : data.getDosages()) {
        if (elt.getAge().size())
            read_dosage_table(dosages[i], elt.getAge(), false);
        else if (elt.getBodymass().size())
            read_dosage_table(dosages[i], elt.getBodymass(), true);
        else {
            if (!elt.getMultiply().present()) {
                throw util::xml_scenario_error("dosages table does not "
                                               "include age, bodymass or multiply element");
            }
            if (elt.getMultiply().get().getBy() != "kg") {
                throw util::xml_scenario_error("dosages table \"multipy by\" "
                                               "option is not \"kg\"");
            }
            dosages[i].useMass = true;
            dosages[i].multMassKg = true;
        }
        dosagesNames[elt.getName()] = i;
        i += 1;
    }
}

// ———  LSTMModel::init (origin: LSTMModel.cpp)  ———

void LSTMModel::init(RunContext& ctx, const scnXml::Scenario& scenario) {
    if (scenario.getPharmacology().present()) {
        LSTMDrugType::init(ctx, scenario.getPharmacology().get().getDrugs());
        LSTMTreatments::init(ctx, scenario.getPharmacology().get().getTreatments());
    }
}

} // namespace PkPd
} // namespace OM
