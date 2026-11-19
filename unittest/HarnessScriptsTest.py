#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Self-checks for util/compareBin.py and util/statCompare.py (no fixtures)."""

import math
import os
import struct
import sys
import tempfile
import unittest

sys.path[0] = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "util")
import compareBin
import statCompare

def writeBin(path, rows):
    """rows: list of (survey, group, measure, value) all doubles."""
    n = len(rows)
    with open(path, "wb") as f:
        f.write(struct.pack("<8sIQQ", b"OMOUTBIN", 2, 0, n))
        for col in range(3):
            f.write(struct.pack("<%di" % n, *[r[col] for r in rows]))
        f.write(struct.pack("<%dd" % n, *[r[3] for r in rows]))

class HarnessTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.p = lambda name: os.path.join(self.dir.name, name)

    def test_bin_exact_roundtrip(self):
        rows = [(1, 0, 2, 3.5), (1, 0, 3, float("nan")), (2, 1, 3, 1e300)]
        a, b = self.p("a.bin"), self.p("b.bin")
        writeBin(a, rows)
        writeBin(b, rows)
        self.assertTrue(compareBin.compare(a, b))
        writeBin(b, rows[:2])
        self.assertFalse(compareBin.compare(a, b)) # missing row detected
        writeBin(b, [(1, 0, 2, 3.5000000001), (1, 0, 3, float("nan")), (2, 1, 3, 1e300)])
        self.assertFalse(compareBin.compare(a, b)) # 1e-10 drift fails exact check

    def test_text_vs_bin(self):
        rows = [(1, 0, 2, 3.5), (2, 1, 3, 7.0)]
        b = self.p("o.bin")
        writeBin(b, rows)
        t = self.p("o.txt")
        with open(t, "w") as f:
            for r in rows:
                f.write("%d\t%d\t%d\t%s\n" % r)
        self.assertEqual(compareBin.load(b), compareBin.load(t))

    def test_stat_means(self):
        # dir A: seeds 1,2,3; dir B: each value +0.1 -> ~5% drift on the mean
        keys = [(9, 0, 4)]
        for name, offset in (("a", 0.0), ("b", 0.1)):
            d = os.path.join(self.dir.name, name)
            os.mkdir(d)
            for seed in (1, 2, 3):
                writeBin(os.path.join(d, "s%d.bin" % seed),
                         [(k[0], k[1], k[2], seed + offset) for k in keys])
        ok = statCompare.compareDirs(os.path.join(self.dir.name, "a"),
                                     os.path.join(self.dir.name, "b"),
                                     2e-2, 1e-6, report=open(os.devnull, "w"))
        self.assertFalse(ok) # 5% mean drift fails at rel=2e-2
        ok = statCompare.compareDirs(os.path.join(self.dir.name, "a"),
                                     os.path.join(self.dir.name, "b"),
                                     1e-1, 1e-6, report=open(os.devnull, "w"))
        self.assertTrue(ok)

if __name__ == "__main__":
    unittest.main()
