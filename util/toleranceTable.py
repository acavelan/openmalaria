#!/usr/bin/env python3
"""T3 calibration: per-(survey,group,measure) CV table from an N-seed sweep of
the REFERENCE binary on the representative scenario (plan v3.3 §4).

The frozen table gives phase-8+ statistical gates their per-key tolerance:
    rel(key) = max(--floor, cv_mult * cv(key) * sqrt(2 / n_seeds))
i.e. the expected relative spread of the mean-of-N between two binaries that
differ only in sampling noise. Smaller floor (default 2%) never applies per
key -- blanket tolerances are what this table replaces.

Usage: toleranceTable.py BINARY [--scenario util/example/example_scenario.xml]
                          [--seeds 1-10] [-j N] [-o test/tolerance_example_p2000.csv.gz]
Output: CSV.gz (comments + header survey,group,measure,n,mean,sd,cv,rel).
"""
import argparse, json, math, os, statistics, sys

sys.path[0] = os.path.dirname(os.path.abspath(__file__))
import parity

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("--scenario", default=parity.EXAMPLE)
    ap.add_argument("--seeds", default="1-10")
    ap.add_argument("-j", type=int, default=10)
    ap.add_argument("-o", "--out", default=os.path.join(parity.REPO, "test", "tolerance_example_p2000.csv.gz"))
    a = ap.parse_args()
    seeds = []
    for part in a.seeds.split(","):
        lo, _, hi = part.partition("-")
        seeds += list(range(int(lo), int(hi or lo) + 1))

    import concurrent.futures as cf
    work = "/tmp/om-tol-%d" % os.getpid()
    os.makedirs(work, exist_ok=True)
    bins = {}
    with cf.ThreadPoolExecutor(max_workers=a.j) as ex:
        futs = {ex.submit(parity.run_one, os.path.abspath(a.binary), a.scenario, s, "s%d" % s, work): s for s in seeds}
        for f, s in futs.items():
            path, d = f.result()
            if path == "ERROR":
                print(d); sys.exit(2)
            bins[s] = path

    from compareBin import load
    series = {}  # key -> [values per seed]
    for s in seeds:
        vals = load(bins[s])
        for k, v in vals.items():
            series.setdefault(k, []).append(v)

    table, cvs = {}, []
    for k, v in sorted(series.items()):
        n = len(v)
        mean = sum(v) / n
        sd = statistics.stdev(v) if n > 1 else 0.0
        cv = sd / abs(mean) if mean != 0 else (float("inf") if sd != 0 else 0.0)
        rel = max(cv * math.sqrt(2.0 / n) * 4.0, 1e-12) if math.isfinite(cv) else float("inf")
        table["%d,%d,%d" % k] = {"n": n, "mean": mean, "sd": sd, "cv": cv, "rel": rel}
        if math.isfinite(cv):
            cvs.append(cv)
    import gzip
    rows = sorted(table.items(), key=lambda kv: tuple(map(int, kv[0].split(","))))
    out_path = a.out
    opener = gzip.open if out_path.endswith(".gz") else open
    with opener(out_path, "wt") as f:
        f.write("# T3 tolerance table (plan v3.3 §4). key=survey,group,measure.\n")
        f.write("# scenario=%s binary=%s seeds=%s\n" % (os.path.abspath(a.scenario), os.path.abspath(a.binary), a.seeds))
        f.write("# rule: rel(key)=max(cv*%.1f*sqrt(2/n),abs_floor); cv=inf means mean==0 with sd>0 (compare abs).\n" % 4.0)
        f.write("survey,group,measure,n,mean,sd,cv,rel\n")
        for k, d in rows:
            f.write("%s,%d,%.9g,%.9g,%s,%.9g\n" % (k, d["n"], d["mean"], d["sd"],
                    ("%.9g" % d["cv"]) if math.isfinite(d["cv"]) else "inf", d["rel"]))
    cvs.sort()
    sm = {"n_keys": len(table), "cv_median": cvs[len(cvs)//2] if cvs else None,
          "cv_p95": cvs[int(len(cvs)*0.95)] if cvs else None, "cv_max": cvs[-1] if cvs else None}
    print("wrote %s: %d keys, CV median %.4g p95 %.4g max %.4g" % (out_path, sm["n_keys"], sm["cv_median"], sm["cv_p95"], sm["cv_max"]))

if __name__ == "__main__":
    main()
