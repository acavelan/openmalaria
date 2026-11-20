#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Pathogenesis.h"
#include "Host/WithinHost/Treatments.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/errors.h"

#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

bool reportInfectionsByGenotype = false;

bool opt_vivax_simple = false, opt_dummy_whm = false, opt_empirical_whm = false, opt_molineaux_whm = false,
     opt_penny_whm = false, opt_common_whm = false;

// -----  static functions  -----

TreatmentId WHInterface::addTreatment(Context& ctx, const scnXml::TreatmentOption& desc) {
    return Treatments::addTreatment(ctx, desc);
}

// Variant selection (ModelOptions flags set by WHInterface::init): one line
// per derived class; the Mock variant is created only by unittest/WHMock.cpp.
unique_ptr<WHInterface> WHInterface::createWithinHostModel(Context& ctx, LocalRng& rng, double comorbidityFactor) {
    if (opt_vivax_simple)
        return unique_ptr<WHInterface>(new WHVivax(ctx, rng, comorbidityFactor));
    if (opt_common_whm)
        return unique_ptr<WHInterface>(new WHCommon(ctx, rng, comorbidityFactor));
    return unique_ptr<WHInterface>(new WHDescriptive(ctx, rng, comorbidityFactor));
}

// -----  shared treatment (was treatment_falciparum; Vivax uses the same
// code — treatSimple dispatches to the variant)  -----

void WHInterface::treatment(Host::Human& human, TreatmentId treatId) {
    const Treatments& treat = Treatments::select(ctx, treatId);
    treatSimple(human, treat.timeLiver, treat.timeBlood);

    // triggered intervention deployments:
    treat.deploy(ctx, human, mon::Deploy::TREAT, interventions::VaccineLimits(/*default initialise: no limits*/));
}

}
}
