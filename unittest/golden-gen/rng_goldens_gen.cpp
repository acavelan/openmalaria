/* Generator for unittest/goldens/rng_dist_v50.csv — run via build-goldens.sh,
 * do not build by hand (that script records the exact recipe).
 *
 * Emits one row per draw: dist, params, draw_index, value (%.17g → bit-exact
 * round-trip), and the 4 post-draw Xoshiro256P state words (hex).
 */
#include "RngGoldenAdapter.h"
#include <cstdio>
#include <gsl/gsl_version.h>

using namespace GoldenRng;

int main() {
    printf(
        "# rng_dist_v50.csv — T2 goldens for the OpenMalaria RNG primitive layer (Xoshiro256P + GSL-backed distributions).\n");
    printf(
        "# Seed: Xoshiro256P(%llu,%llu,%llu,%llu); plumbing/distributions mirror model/util/random.h (see golden-gen/RngGoldenAdapter.h, file:line cites there).\n",
        (unsigned long long)SEED[0], (unsigned long long)SEED[1], (unsigned long long)SEED[2],
        (unsigned long long)SEED[3]);
    printf(
        "# Each series restarts from the seed; value is the %%.17g draw; s0..s3 are the post-draw xoshiro state words (hex, via Xoshiro256P::binary_checkpoint).\n");
    printf(
        "# Compiler: g++ %s; build: see golden-gen/build-goldens.sh (-O2 -ffp-contract=off); GSL %s; x86-64 IEEE-754 binary64.\n",
        __VERSION__, GSL_VERSION);
    printf(
        "# gaussian: the model path is gsl_ran_gaussian = polar Box-Muller (random.h:155); gsl_ran_gaussian_ziggurat is on NO model code path, so no ziggurat tail seed applies.\n");
    printf(
        "# binomial/exponential: random.h exposes no wrapper for these; pinned via direct gsl_ran_binomial / gsl_ran_exponential on the same stream (gsl_ran_binomial picks Knuth vs BTPE from n,p internally).\n");
    printf(
        "# branch coverage: gamma a<1 small-shape vs a>=1 Marsaglia-Tsang; beta (0.5,0.7) vs (2.0,3.0) sampler branches; poisson mu=3.2 direct vs mu=25 rejection (branch at mu=10).\n");
    printf("# dist,params,draw_index,value,s0,s1,s2,s3\n");

    size_t nSeries;
    const DrawSpec* S = schedule(nSeries);
    for (size_t i = 0; i < nSeries; ++i) {
        const DrawSpec& s = S[i];
        Rng rng(SEED[0], SEED[1], SEED[2], SEED[3]);
        for (int d = 0; d < s.n; ++d) {
            double v = golden_draw(rng, s);
            uint64_t st[4];
            rng.state(st);
            printf("%s,%s,%d,%.17g,%016llx,%016llx,%016llx,%016llx\n", s.dist, s.params, d, v,
                   (unsigned long long)st[0], (unsigned long long)st[1], (unsigned long long)st[2],
                   (unsigned long long)st[3]);
        }
    }

    // Pinv spots: deterministic, RNG-independent; state words are the seed
    // state before the first spot and unchanged after each (Pinv consumes none).
    Rng rng(SEED[0], SEED[1], SEED[2], SEED[3]);
    size_t nSpots;
    const double* U = pinv_spots(nSpots);
    for (size_t i = 0; i < nSpots; ++i) {
        double v = gsl_cdf_ugaussian_Pinv(U[i]);
        uint64_t st[4];
        rng.state(st);
        printf("ugaussian_Pinv,u=%.17g,%d,%.17g,%016llx,%016llx,%016llx,%016llx\n", U[i], static_cast<int>(i), v,
               (unsigned long long)st[0], (unsigned long long)st[1], (unsigned long long)st[2],
               (unsigned long long)st[3]);
    }
    return 0;
}
