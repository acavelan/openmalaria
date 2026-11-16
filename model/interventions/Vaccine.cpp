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

#include "ScenarioXml.h"
#include "interventions/Vaccine.h"
#include "interventions/HumanComponents.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "schema/interventions.h"
#include "util/StreamValidator.h"
#include "Host/WithinHost/Genotypes.h"

#include <limits>
#include <cmath>

namespace OM {
namespace interventions {
using namespace OM::util;

VaccineComponent::VaccineComponent(ComponentId component, const scnXml::VaccineDescription& vd,
        Vaccine::Types type, bool report) : VaccineComponent(component, type, report)
{
    decayFunc = xml::readDecay(vd.getDecay(), "decay");
    efficacyB = vd.getEfficacyB().getValue();
    if( type == Vaccine::BSV && ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) )
        throw util::unimplemented_exception( "blood stage vaccines (BSV) cannot be used with vivax model" );

    const bool opt_vaccine_genotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    const scnXml::VaccineDescription::InitialEfficacySequence ies = vd.getInitialEfficacy();
    const scnXml::VaccineDescription::PhenotypeSequence ps = vd.getPhenotype();

    if(opt_vaccine_genotype == false && ies.size() == 0)
        throw util::xml_scenario_error( "Vaccine: initialEfficacy is required but is missing" );

    if(opt_vaccine_genotype == true && ies.size() != 0)
        throw util::xml_scenario_error( "Vaccine: initialEfficacy outside of phenotype definition with VACCINE_GENOTYPE option enabled is not allowed" );

    if(opt_vaccine_genotype == true)
    {
        perGenotypeInitialMeanEfficacy.resize(WithinHost::Genotypes::N());

        for(size_t i = 0; i < ps.size(); ++i)
        {

            const scnXml::VaccinePhenotype::RestrictionSequence restrictions = ps[i].getRestriction();
            const scnXml::VaccinePhenotype::InitialEfficacySequence vpies = ps[i].getInitialEfficacy();

            for(size_t j = 0; j < restrictions.size(); ++j)
            {

                uint32_t allele_code = WithinHost::Genotypes::findAlleleCode(restrictions[j].getOnLocus(), restrictions[j].getToAllele());

                if( allele_code == numeric_limits<uint32_t>::max() ){
                    throw util::xml_scenario_error("phenotype has "
                    "restriction on locus " + string(restrictions[j].getOnLocus())
                    + " allele " + string(restrictions[j].getToAllele()) + " but this locus/allele "
                    "has not been defined in parasiteGenetics section");
                }

                if (allele_code >= WithinHost::Genotypes::N())
                    throw util::xml_scenario_error( string("Vaccine phenotype: Invalid Locus ") + restrictions[j].getOnLocus() + string(" and allele ") + restrictions[j].getToAllele());

                if (perGenotypeInitialMeanEfficacy[allele_code].size() != 0)
                    throw util::xml_scenario_error( string("Vaccine phenotype: there are multiple restrictions locus " + string(restrictions[j].getOnLocus() + string(" and allele ") + restrictions[j].getToAllele())));

                if(vpies.size() == 0)
                    throw util::xml_scenario_error( string("Vaccine phenotype: initialEfficacy is required but is missing for phenotype: " + string(restrictions[j].getOnLocus() + string(" ") + restrictions[j].getToAllele())));

                perGenotypeInitialMeanEfficacy[allele_code].resize(vpies.size());

                for(size_t k = 0; k < vpies.size(); ++k)
                    perGenotypeInitialMeanEfficacy[allele_code][k] = vpies[k].getValue();
            }
        }

        // Check that all genotypes have been adressed
        for(size_t i = 0; i < WithinHost::Genotypes::N(); ++i)
        {
            if(perGenotypeInitialMeanEfficacy[i].size() == 0)
                throw util::xml_scenario_error( "Vaccine phenotype: missing restriction for at least one locus and allele pair" );
        }
    }
    else
    {
        perGenotypeInitialMeanEfficacy.resize(1);
        perGenotypeInitialMeanEfficacy[0].resize(ies.size());
        for(size_t i = 0; i < ies.size(); ++i)
            perGenotypeInitialMeanEfficacy[0][i] = ies[i].getValue();
    }

}

void VaccineComponent::deploy(Host::Human& human, mon::Deploy::Method method, VaccineLimits vaccLimits) const
{
    bool administered = vaccinate(human.interventions.vaccines, human.rng, *this, vaccLimits,
        sim::nowOrTs1(), WithinHost::Genotypes::N());
    if( administered && reportDeployments ){
        mon::recordDeploy(mon::measure("nMassVaccinations"), mon::measure("nEPIVaccinations"), human, method);
    }
    assert(type == Vaccine::PEV || type == Vaccine::BSV || type == Vaccine::TBV);
    if (method == mon::Deploy::TREAT)
        mon::recordEvent(mon::measure("nTreatDeployments"), human);
}

Component::Type VaccineComponent::componentType() const
{
    if( type == Vaccine::PEV ) return Component::PEV;
    else if( type == Vaccine::BSV ) return Component::BSV;
    else if( type == Vaccine::TBV ) return Component::TBV;
    else throw SWITCH_DEFAULT_EXCEPTION;
}

void VaccineComponent::print_details(ostream& out) const
{
        out << id().id << "\t" << (type == Vaccine::PEV ? "PEV" : (type == Vaccine::BSV ? "BSV" : "TBV"));
}

double VaccineComponent::getInitialEfficacy (LocalRng& rng, size_t numPrevDoses, uint32_t genotype) const
{
    /* If initialMeanEfficacy.size or more doses have already been given, use
     * the last efficacy. */
    if (numPrevDoses >= perGenotypeInitialMeanEfficacy[genotype].size())
        numPrevDoses = perGenotypeInitialMeanEfficacy[genotype].size() - 1;
    double ime = perGenotypeInitialMeanEfficacy[genotype][numPrevDoses];
    //NOTE(validation): With extra valiadation in random, the first difference is noticed here:
    util::streamValidate(ime);
    util::streamValidate(efficacyB);
    if (ime == 0.0){
        return 0.0;
    } else if (ime < 1.0) {
        double result = rng.betaWithMean (ime, efficacyB);
        //NOTE(validation):: Without extra validation in random, the first difference is noticed here:
        util::streamValidate(result);
        //TODO(validation):: Why the difference? Bug in/limitation of StreamValidatior?
        // Extra memory allocation due to the extra logging causes some bad
        // memory usage to manifest differently? Am I forgetting to initialise
        // some variable?
        return result;
    } else {
        return 1.0;
    }
}

double vaccineFactor(const vector<VaccineEffect>& effects, Vaccine::Types type,
        SimTime effectTime, uint32_t genotype) {
    double factor = 1.0;
    for (const auto& effect : effects) {
        if (effect.parameters->type == type) {
            SimTime age = effectTime - effect.timeLastDeployment;
            double decayFactor = util::eval(effect.parameters->decayFunc, age, *effect.hetSample);
            factor *= 1.0 - effect.perGenotypeInitialEfficacy[genotype] * decayFactor;
        }
    }
    return factor;
}

bool vaccinate(vector<VaccineEffect>& effects, LocalRng& rng, const VaccineComponent& params,
        VaccineLimits limits, SimTime date, size_t numGenotypes) {
    VaccineEffect* effect = nullptr;
    for (auto& candidate : effects) {
        if (candidate.parameters == &params) {
            effect = &candidate;
            break;
        }
    }
    uint32_t previousDoses = effect ? effect->numDosesAdministered : 0;
    if (previousDoses < limits.minPrevDoses || previousDoses >= limits.maxCumDoses)
        return false;

    if (!effect) {
        VaccineEffect initial;
        initial.parameters = &params;
        initial.hetSample = params.decayFunc.heterogeneity.sample(rng);
        effects.push_back(std::move(initial));
        effect = &effects.back();
    }
    effect->perGenotypeInitialEfficacy.resize(numGenotypes);
    for (size_t g = 0; g < params.perGenotypeInitialMeanEfficacy.size(); ++g)
        effect->perGenotypeInitialEfficacy[g] = params.getInitialEfficacy(rng, previousDoses, g);
    util::streamValidate(effect->perGenotypeInitialEfficacy);
    effect->numDosesAdministered = previousDoses + 1;
    effect->timeLastDeployment = date;
    return true;
}

void checkpointVaccines(vector<VaccineEffect>& effects, ostream& stream) {
    effects.size() & stream;
    for (auto& effect : effects) {
        effect.parameters->id() & stream;
        effect.checkpointState(stream);
    }
}

void checkpointVaccines(vector<VaccineEffect>& effects, istream& stream,
        const vector<unique_ptr<HumanInterventionComponent>>& definitions) {
    size_t count;
    count & stream;
    util::checkpoint::validateListSize(count);
    effects.resize(count);
    for (auto& effect : effects) {
        ComponentId id(stream);
        if (id.id >= definitions.size() ||
            !(effect.parameters = dynamic_cast<const VaccineComponent*>(definitions[id.id].get())))
            throw util::checkpoint_error("invalid vaccine component in checkpoint");
        effect.checkpointState(stream);
    }
}

}
}
