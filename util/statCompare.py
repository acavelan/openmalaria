#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file is part of OpenMalaria.
#
# OpenMalaria is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.

"""T3 statistical correctness gate: compare two sets of repeated scenario runs
(one file per seed) of the same scenario, e.g. old vs new binary.

Each directory holds the per-seed outputs (text or output.bin / output.bin.gz,
see compareBin). Values are aggregated as the mean over seeds per
(survey, group, measure) key and the two means are compared to rel/abs
tolerance (default 2% / 1e-6, reusing util/approxEqual.py).

Usage: statCompare.py DIR_A DIR_B [--rel 0.02] [--abs 1e-6]
Exit status: 0 within tolerance, 1 outside or key-set mismatch, 2 unreadable.
"""

import glob
import os
import sys

sys.path[0] = os.path.dirname(os.path.abspath(__file__))
from approxEqual import ApproxSame
from compareBin import load

def meanByFile(files):
    """files -> ({key: mean}, nFiles)"""
    sums = {}
    counts = {}
    n = 0
    for path in files:
        values = load(path)
        n += 1
        for key, value in values.items():
            sums[key] = sums.get(key, 0.0) + value
            counts[key] = counts.get(key, 0) + 1
    return {key: sums[key] / counts[key] for key in sums}, n

def compareDirs(dirA, dirB, relPrec, absPrec, report=sys.stderr):
    files = lambda d: sorted(glob.glob(os.path.join(d, "*.txt")) +
                             glob.glob(os.path.join(d, "*.bin")) +
                             glob.glob(os.path.join(d, "*.bin.gz")))
    fa, fb = files(dirA), files(dirB)
    if not fa or not fb:
        raise ValueError("no output files in %s or %s" % (dirA, dirB))
    meansA, na = meanByFile(fa)
    meansB, nb = meanByFile(fb)
    approx = ApproxSame(relPrec, absPrec)
    ok = True
    bad = 0
    for key in sorted(set(meansA) ^ set(meansB)):
        print("%s: key missing in %s" % (key, dirA if key in meansB else dirB), file=report)
        ok = False
    for key in sorted(set(meansA) & set(meansB)):
        if not approx(meansA[key], meansB[key]):
            if bad < 50: # don't flood; total count is what matters
                print("%s: %r vs %r" % (key, meansA[key], meansB[key]), file=report)
            bad += 1
            ok = False
    common = set(meansA) & set(meansB)
    print("%s vs %s: %d seeds vs %d, %d/%d common keys within rel=%g abs=%g"
          % (dirA, dirB, na, nb, len(common) - bad, len(common), relPrec, absPrec),
          file=report)
    if bad:
        print("%d keys outside tolerance" % bad, file=sys.stderr)
    return ok

if __name__ == "__main__":
    args = sys.argv[1:]
    def opt(name, default):
        if name in args:
            i = args.index(name)
            v = float(args[i + 1])
            del args[i:i + 2]
            return v
        return default
    rel, abs = opt("--rel", 2e-2), opt("--abs", 1e-6)
    if len(args) != 2:
        sys.exit(__doc__)
    try:
        sys.exit(0 if compareDirs(args[0], args[1], rel, abs) else 1)
    except (OSError, ValueError) as e:
        print(e, file=sys.stderr)
        sys.exit(2)
