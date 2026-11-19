#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file is part of OpenMalaria.
#
# OpenMalaria is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.

"""Exact comparison of two openMalaria output files.

Reads binary output (output.bin, format v2, optionally .gz) or text output
(survey<tab>group<tab>measure<tab>value rows), keyed on
(survey, group, measure). Binary values are compared bit-exactly (NaN==NaN).

Usage: compareBin.py FILE1 FILE2
Exit status: 0 identical, 1 differ, 2 unreadable/unparsable.
"""

import gzip
import struct
import sys

def readBin(path):
    """Parse output.bin v2 (see model/mon/Monitoring.cpp) -> {(survey, group, measure): value}."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as stream:
        magic, version, integerRows, doubleRows = struct.unpack("<8sIQQ", stream.read(28))
        if magic != b"OMOUTBIN" or version != 2:
            raise ValueError("unsupported binary output (need OMOUTBIN v2): " + path)
        rows = integerRows + doubleRows
        def readValues(code, count):
            size = struct.calcsize(code) * count
            data = stream.read(size)
            if len(data) != size:
                raise ValueError("truncated binary output: " + path)
            return struct.unpack("<%d%s" % (count, code), data)
        surveys = readValues("i", rows)
        groups = readValues("i", rows)
        measures = readValues("i", rows)
        values = readValues("i", integerRows) + readValues("d", doubleRows)
        if stream.read(1):
            raise ValueError("extra data in binary output: " + path)
    out = {}
    for key, value in zip(zip(surveys, groups, measures), values):
        if key in out:
            raise ValueError("duplicate (survey, group, measure) row in " + path)
        out[key] = value
    return out

def readText(path):
    """Parse text output rows -> {(survey, group, measure): float(value)}."""
    out = {}
    with open(path) as stream:
        for line in stream:
            survey, group, measure, value = line.split("\t")
            key = (int(survey), int(group), int(measure))
            if key in out:
                raise ValueError("duplicate (survey, group, measure) row in " + path)
            out[key] = float(value)
    return out

def load(path):
    if path.endswith(".bin") or path.endswith(".bin.gz"):
        return readBin(path)
    return readText(path)

def equal(a, b):
    return a == b or (a != a and b != b) # NaN equals itself

def compare(file1, file2, report=sys.stderr):
    """Return True iff every value is exactly equal."""
    d1, d2 = load(file1), load(file2)
    ok = True
    for key in sorted(set(d1) ^ set(d2)):
        print("%s: missing in %s" % (key, file1 if key in d2 else file2), file=report)
        ok = False
    for key in sorted(set(d1) & set(d2)):
        if not equal(d1[key], d2[key]):
            print("%s: %r != %r" % (key, d1[key], d2[key]), file=report)
            ok = False
    return ok

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    try:
        sys.exit(0 if compare(sys.argv[1], sys.argv[2]) else 1)
    except (OSError, ValueError) as e:
        print(e, file=sys.stderr)
        sys.exit(2)
