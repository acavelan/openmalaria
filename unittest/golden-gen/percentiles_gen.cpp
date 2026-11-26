// Compile ONLY against cfc3c8e0 (build-goldens.sh), never the rewritten model.
#include "Transmission/PerHost.h"
#include <cstdio>

int main() {
    scnXml::IntValue rest(3), eip(10);
    scnXml::DoubleValue laidEggs(0.95), seeking(0.15), survFeed(0.62);
    scnXml::SampledValueCV avail;
    avail.setDistr("gamma");
    avail.setCV(0.3);
    scnXml::BetaMeanSample biting(0.11, 0.0), restSite(0.88, 0.0), resting(0.89, 0.0);
    scnXml::DoubleValue oviposit(0.0), bloodIndex(0.9);
    scnXml::Mosq mosq(rest, eip, laidEggs, seeking, survFeed, avail, biting, restSite,
                     resting, oviposit, bloodIndex);
    using OM::Transmission::PerHostAnophParams;
    PerHostAnophParams::init(mosq);
    PerHostAnophParams::setEntoAvailabilityFactor(0, 1.0);
    PerHostAnophParams::calcAvailabilityPercentiles();
    for (int p = 0; p < 100; ++p)
        printf("%.17g\n", PerHostAnophParams::getEntoAvailabilityPercentile(p));
}
