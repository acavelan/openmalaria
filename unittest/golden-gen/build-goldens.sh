#!/usr/bin/env bash
# Opt-in fixture capture. Never run as part of the build/test gate.
# Usage: build-goldens.sh NEW_OUTPUT_DIRECTORY
# Requires the two pinned commits below in local git history.
set -euo pipefail
cd "$(dirname "$0")/../.."
out=$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve())' "${1:?provide a NEW output directory}")
mkdir "$out" # refuse to overwrite a previous capture
baseline=cfc3c8e0caee6fdcd6b488a5ab50fdbf2b7a5a90
generator=6e24e51eb576a12fd701ab54a1b6c5dd375ca081
mkdir "$out/source"
git archive "$baseline" | tar -x -C "$out/source"
mkdir -p "$out/source/golden-gen"
for file in RngGoldenAdapter.h rng_goldens_gen.cpp; do
    git show "$generator:unittest/golden-gen/$file" > "$out/source/golden-gen/$file"
done
cp unittest/golden-gen/percentiles_gen.cpp "$out/source/golden-gen/"
cat >> "$out/source/CMakeLists.txt" <<'CMAKE'

# Fixture drivers only; all linked model sources above are the pinned baseline.
foreach(driver rng_goldens percentiles)
    add_executable(${driver}_gen golden-gen/${driver}_gen.cpp)
    target_link_libraries(${driver}_gen model schema contrib ${GSL_LIBRARIES}
        ${XERCESC_LIBRARIES} ${Z_LIBRARIES} ${PTHREAD_LIBRARIES} ${OM_STD_LIBS})
endforeach()
CMAKE
cmake -S "$out/source" -B "$out/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS=-ffp-contract=off -DCMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE=OFF \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DOM_CXXTEST_ENABLE=OFF -DOM_BOXTEST_ENABLE=OFF
cmake --build "$out/build" --target rng_goldens_gen percentiles_gen -j "${JOBS:-4}"
"$out/build/rng_goldens_gen" > "$out/rng_dist_v50.csv"
"$out/build/percentiles_gen" > "$out/availability_percentiles.txt"
# Store only differing values; all platforms share the original schedule/states.
# Decimal strings are copied, not rounded or recomputed.
python3 - "$out" <<'PY'
import csv, pathlib, struct, sys
out = pathlib.Path(sys.argv[1])
def rows(path):
    return list(csv.reader(line for line in path.read_text().splitlines() if not line.startswith("#")))
original = rows(pathlib.Path("unittest/goldens/rng_dist_v50.csv"))
captured = rows(out / "rng_dist_v50.csv")
assert len(original) == len(captured) == 14012
with (out / "rng_values.csv").open("w", newline="") as f:
    writer = csv.writer(f, lineterminator="\n")
    for old, new in zip(original, captured):
        assert old[:3] + old[4:] == new[:3] + new[4:], "schedule/state drift"
        if struct.pack("d", float(old[3])) != struct.pack("d", float(new[3])):
            assert old[0] not in ("uniform_01", "poisson", "binomial"), "exact primitive drift"
            writer.writerow(new[:4])
values = (out / "availability_percentiles.txt").read_text().splitlines()
assert len(values) == 100
(out / "availability_percentiles.inc").write_text("".join(v + ",\n" for v in values))
PY
{
    echo "baseline=$baseline"
    echo "rng_generator=$generator"
    echo "flags=-O3 -DNDEBUG -ffp-contract=off; LTO/PGO disabled"
    echo "runner=${ImageOS:-local} ${ImageVersion:-}"
    uname -a
    "${CXX:-c++}" --version
    gsl-config --version
    if [[ $(uname -s) == Darwin ]]; then
        sw_vers
        otool -L "$out/build/rng_goldens_gen"
    else
        cat /etc/os-release
        ldd "$out/build/rng_goldens_gen"
        ldd "$out/build/rng_goldens_gen" | awk '/=> \// {print $3}' | xargs sha256sum
    fi
    sha256sum "$out"/{rng_dist_v50.csv,rng_values.csv,availability_percentiles.txt,availability_percentiles.inc} \
        "$out"/source/golden-gen/* "$out"/source/model/util/{random.h,xoshiro.h} \
        "$out"/source/contrib/chacha.h \
        "$(sed -n 's/^GSL_LIB_OPT:FILEPATH=//p' "$out/build/CMakeCache.txt")"
} > "$out/provenance.txt"
