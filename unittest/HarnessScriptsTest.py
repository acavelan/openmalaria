#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Self-checks for the output comparators and parity runner (no fixtures)."""

import contextlib
import gzip
import io
import math
import os
import struct
import sys
import tempfile
import unittest
from unittest import mock

sys.path[0] = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "util")
import compareBin
import statCompare
import parity

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
        self.addCleanup(self.dir.cleanup)
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

    def test_binary_bytes_not_just_numeric_values(self):
        a, b = self.p("a.bin"), self.p("b.bin")
        nan1, nan2 = [struct.unpack("<d", struct.pack("<Q", bits))[0]
                      for bits in (0x7ff8000000000001, 0x7ff8000000000002)]
        rows = [(1, 0, 2, 1.0), (2, 0, 2, 2.0)]
        for left, right in (([(1, 0, 2, 0.0)], [(1, 0, 2, -0.0)]),
                            ([(1, 0, 2, nan1)], [(1, 0, 2, nan2)]),
                            (rows, rows[::-1])):
            with self.subTest(left=left, right=right):
                writeBin(a, left)
                writeBin(b, right)
                self.assertFalse(compareBin.compare(a, b, report=io.StringIO()))

    def test_binary_gzip_compares_uncompressed_bytes(self):
        a, b = self.p("a.bin"), self.p("b.bin.gz")
        writeBin(a, [(1, 0, 2, -0.0)])
        with open(a, "rb") as source, gzip.open(b, "wb") as compressed:
            compressed.write(source.read())
        self.assertTrue(compareBin.compare(a, b))
        writeBin(a, [(1, 0, 2, 0.0)])
        self.assertFalse(compareBin.compare(a, b, report=io.StringIO()))

    def test_parity_checks_continuous_output(self):
        old, new = self.p("old"), self.p("new")
        for directory in (old, new):
            os.mkdir(directory)
            writeBin(os.path.join(directory, "output.bin"), [(1, 0, 2, 1.0)])
        def run_one(binary, scenario, seed, tag, work):
            directory = old if tag == "old" else new
            return os.path.join(directory, "output.bin"), directory
        for before, after, exit_code in ((None, None, 0), ("1\n2\n", "1\n2\n", 0),
                                         ("1\n2\n", "2\n1\n", 1), ("1\n", None, 1),
                                         (None, "1\n", 1)):
            with self.subTest(before=before, after=after):
                for directory, contents in ((old, before), (new, after)):
                    path = os.path.join(directory, "ctsout.txt")
                    if contents is None:
                        if os.path.exists(path):
                            os.remove(path)
                    else:
                        with open(path, "w") as stream:
                            stream.write(contents)
                args = ["parity.py", "--old", os.path.join(old, "output.bin"),
                        "--new", os.path.join(new, "output.bin"), "--tiers", "b", "--seeds", "1",
                        "--work", self.dir.name]
                with mock.patch.object(sys, "argv", args), mock.patch.object(parity, "run_one", run_one), \
                        contextlib.redirect_stdout(io.StringIO()) as output:
                    with self.assertRaises(SystemExit) as result:
                        parity.main()
                self.assertEqual(result.exception.code, exit_code)
                self.assertIn("FAIL " if exit_code else "PASS ", output.getvalue())
                if exit_code:
                    self.assertNotIn("PASS ", output.getvalue())

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
                                     2e-2, 1e-6, report=io.StringIO())
        self.assertFalse(ok) # 5% mean drift fails at rel=2e-2
        ok = statCompare.compareDirs(os.path.join(self.dir.name, "a"),
                                     os.path.join(self.dir.name, "b"),
                                     1e-1, 1e-6, report=io.StringIO())
        self.assertTrue(ok)

if __name__ == "__main__":
    unittest.main()
