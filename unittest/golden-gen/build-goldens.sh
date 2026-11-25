#!/bin/sh
# Regenerate unittest/goldens/rng_dist_v50.csv.
#
# Exact compile command (the ONLY supported way to rebuild the generator):
#   g++ -O2 -std=c++20 -ffp-contract=off -I model -I unittest \
#       -I/usr/include/gsl unittest/golden-gen/rng_goldens_gen.cpp \
#       model/util/errors.cpp \
#       -o /tmp/opencode/rng_goldens_gen $(gsl-config --libs) -lbfd -ldl -lstdc++exp -lm
# The binary prints its compiler and GSL versions into the CSV header.
# Requires GSL headers (explicit -I for /usr/include/gsl: this box's
# gsl-config --cflags is wrong (-I/usr/include) and pkg-config --cflags gsl is
# empty; CMake's FindGSL cache value is GSL_INCLUDE_DIRS=/usr/include/gsl).
# errors.cpp is compiled in for the single
# symbol the mirrored random.h poisson wrapper (random.h:214) can throw.
set -e
cd "$(dirname "$0")/../.."
BIN=/tmp/opencode/rng_goldens_gen
g++ -O2 -std=c++20 -ffp-contract=off -I model -I unittest \
    -I/usr/include/gsl unittest/golden-gen/rng_goldens_gen.cpp \
    model/util/errors.cpp \
    -o "$BIN" $(gsl-config --libs) -lbfd -ldl -lstdc++exp -lm
mkdir -p unittest/goldens
"$BIN" > unittest/goldens/rng_dist_v50.csv
wc -l unittest/goldens/rng_dist_v50.csv
