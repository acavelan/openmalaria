#!/usr/bin/env python3
"""T1 bit-exact A/B parity runner (plan v3.3 §4).

Runs two openMalaria binaries (old reference vs new) on:
  tier a: every test/scenario*.xml, single continuous run
  tier b: util/example/example_scenario.xml with iseed patched to each seed
and compares output.bin bit-exactly (util/compareBin).

Hard-fails if the reference recipe (compiler/GSL/libc/kernel) drifted: run
with --recipe /path/recipe.json next to the reference binary.

Usage: parity.py --old OLDBIN --new NEWBIN [--recipe recipe.json]
                 [--tiers a,b] [--seeds 1-10] [-j N] [--work DIR]
Exit: 0 all identical, 1 any mismatch, 2 setup/recipe error.
"""
import argparse, glob, json, os, re, shutil, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path[0] = os.path.dirname(os.path.abspath(__file__))
import compareBin

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = os.path.join(REPO, "util", "example", "example_scenario.xml")

def schema_of(scenario):
    with open(scenario) as f:
        head = f.read(4096)
    m = re.search(r'schemaLocation="[^"]*\s+([^"\s]+)"|noNamespaceSchemaLocation="([^"]+)"', head)
    if not m:
        raise RuntimeError("no schema location in " + scenario)
    return m.group(1) or m.group(2)

def run_one(binary, scenario, seed, tag, work):
    d = tempfile.mkdtemp(prefix="%s-s%s-" % (tag, seed), dir=work)
    scen = scenario
    if seed is not None:  # patch iseed in a throwaway copy (binary-agnostic)
        src = open(scenario).read()
        new, n = re.subn(r'iseed="\d+"', 'iseed="%d"' % seed, src)
        if n == 0:
            raise RuntimeError("no iseed attribute in " + scenario)
        scen = os.path.join(d, os.path.basename(scenario))
        open(scen, "w").write(new)
    xsd = schema_of(scen)
    for cand in (os.path.join(REPO, "schema", xsd), os.path.join(REPO, "test", xsd),
                 os.path.join(REPO, "build", "schema", xsd)):
        if os.path.isfile(cand):
            os.symlink(cand, os.path.join(d, xsd))
            break
    out = os.path.join(d, "output.bin")
    r = subprocess.run([binary, "--resource-path", os.path.join(REPO, "test"),
                        "--scenario", os.path.abspath(scen), "--output-format", "bin",
                        "--output", "output.bin"],
                       cwd=d, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.isfile(out):
        return ("ERROR", "%s %s seed=%s rc=%d %s" % (os.path.basename(binary), tag, seed, r.returncode, r.stderr[:200]))
    return (out, d)

def check_recipe(recipe_path, old_binary):
    ref = json.load(open(recipe_path))
    def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout.strip()
    gsl = sh("ldd %s | awk '/libgsl.so/{print $3}'" % old_binary)
    cur = {
        "compiler": sh("gcc --version | head -1"),
        "gsl": sh("gsl-config --version"),
        "gsl_so_sha256": sh("sha256sum %s | cut -d' ' -f1" % gsl) if gsl else None,
        "machine": sh("uname -sr"),
    }
    bad = {k: (ref.get(k), v) for k, v in cur.items() if ref.get(k) not in (None, v)}
    if bad:
        print("RECIPE MISMATCH (bit-parity claims invalid):\n  " + json.dumps(bad, indent=2))
        sys.exit(2)
    print("recipe ok:", " ".join("%s=%s" % (k, v) for k, v in cur.items()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--recipe")
    ap.add_argument("--tiers", default="a,b")
    ap.add_argument("--seeds", default="1-10")
    ap.add_argument("-j", type=int, default=max(1, min(16, (os.cpu_count() or 4) - 6)))
    ap.add_argument("--work", default=None)
    a = ap.parse_args()

    if a.recipe:
        check_recipe(a.recipe, a.old)
    for b in (a.old, a.new):
        if not os.path.isfile(b):
            print("missing binary:", b)
            sys.exit(2)
    a.old, a.new = os.path.abspath(a.old), os.path.abspath(a.new)

    seeds = []
    for part in a.seeds.split(","):
        lo, _, hi = part.partition("-")
        seeds += list(range(int(lo), int(hi or lo) + 1))

    jobs = []  # (scenario, seed, label)
    if "a" in a.tiers:
        for p in sorted(glob.glob(os.path.join(REPO, "test", "scenario*.xml"))):
            jobs.append((p, None, os.path.basename(p)))
    if "b" in a.tiers:
        for s in seeds:
            jobs.append((EXAMPLE, s, "example#%d" % s))

    work = a.work or tempfile.mkdtemp(prefix="parity-")
    fails = []
    with ThreadPoolExecutor(max_workers=a.j) as ex:
        futs = {}
        for scen, seed, label in jobs:
            futs[ex.submit(run_one, a.old, scen, seed, "old", work)] = (scen, seed, label)
        done = {}
        for f, (scen, seed, label) in futs.items():
            res, info = f.result()
            if res == "ERROR":
                fails.append(info)
                continue
            done[label] = (info, res)
        futs2 = {ex.submit(run_one, a.new, scen, seed, "new", work): label
                 for scen, seed, label in jobs}
        for f, label in futs2.items():
            res, info = f.result()
            if res == "ERROR":
                fails.append(info)
                continue
            odir, obin = done[label]
            nbin = os.path.join(info, "output.bin")
            ok = compareBin.compare(obin, nbin)
            print(("PASS " if ok else "FAIL ") + label)
            if not ok:
                fails.append("%s: output.bin differs" % label)

    print("\n%d/%d identical; %d problems" % (len(jobs) - len(fails), len(jobs), len(fails)))
    for f in fails:
        print(" ", f)
    if not fails and a.work is None:
        shutil.rmtree(work, ignore_errors=True)
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
