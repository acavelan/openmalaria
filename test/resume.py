#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Checkpoint resume-equivalence test (model/checkpoint/).
#
# For each scenario NAME (scenarioNAME.xml in this directory):
#   full:  run to completion
#   split: run with --save-to ck1 --save-at main (exits 0), then resume ck1
#          with --save-to ck2 --save-at <main start + 50 days> (exits 0),
#          then resume ck2 to completion.
# Both must produce byte-identical survey output (.bin) and ctsout files.
# Exit status: 0 all passed, 1 a failure (details printed).

import os
import struct
import subprocess
import sys
import tempfile

testSrcDir="@CMAKE_CURRENT_SOURCE_DIR@"
testBuildDir="@CMAKE_CURRENT_BINARY_DIR@"
exe=os.path.abspath(os.path.join(testBuildDir, "openMalaria"))
if not os.path.isfile(exe):
    exe=os.path.abspath(os.path.join(testBuildDir, "..", "openMalaria"))
schemaName=None

def run (args, cwd):
    with open (os.path.join (cwd, "err.log"), "ab") as err:
        ret=subprocess.call ([exe]+args, cwd=cwd,
            stdout=subprocess.DEVNULL, stderr=err)
    return ret

def checkpointClock (path):
    # header: u32 magic, u32 version, u8 phase, then t0, t1, interv, endTime
    # (all sim days)
    with open (path, "rb") as f:
        head=f.read (25)
    magic, version = struct.unpack ("<II", head[:8])
    assert magic==0x4B634D4F and version==3, "not a checkpoint"
    t0, t1, interv, end = struct.unpack ("<iiii", head[9:25])
    return t0, end

def check (name):
    scenario=os.path.join (testSrcDir, "scenario%s.xml" % name)
    assert os.path.isfile (scenario), scenario
    # The generated schema must sit next to the scenario copy's working dir.
    import re
    with open (scenario) as f:
        m=re.search (r'schemaLocation="[^"]* (\S+\.xsd)"', f.read())
    xsd=os.path.abspath (os.path.join (testBuildDir, "..", "schema", m.group(1)))
    assert os.path.isfile (xsd), xsd
    # scenario's own iseed (not overridden: some solvers are seed-sensitive)
    base=["--scenario", os.path.abspath (scenario), "--output-format", "bin"]
    src=["--resource-path", testSrcDir]

    simDir=tempfile.mkdtemp (prefix="resume-%s-" % name, dir=testBuildDir)
    os.link (xsd, os.path.join (simDir, m.group(1)))
    def om (extra):
        return run (base+src+extra, simDir)

    rc=om (["-o", "full.bin", "--ctsout", "full_cts.txt"])
    if rc: return "%s: full run failed (rc %d)" % (name, rc)

    rc=om (["-o", "out.bin", "--ctsout", "cts.txt",
            "--save-to", "ck1", "--save-at", "main"])
    if rc: return "%s: save@main failed (rc %d)" % (name, rc)
    day, end=checkpointClock (os.path.join (simDir, "ck1"))

    # resume to mid-main, checkpoint again, resume again (chain: the resumed
    # run must be able to save and be resumed itself)
    rc=om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck1",
            "--save-to", "ck2", "--save-at", str (day+50)])
    if rc: return "%s: resume+save@mid-main failed (rc %d)" % (name, rc)

    # late-main save day, chained after the mid-main checkpoint
    rc=om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck2",
            "--save-to", "ck4", "--save-at", str (end-25)])
    if rc: return "%s: resume+save@late-main failed (rc %d)" % (name, rc)
    if om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck4"]):
        return "%s: late resume failed" % name

    # mid-warmup save day, independent chain (outputs suffixed B)
    if om (["-o", "outB.bin", "--ctsout", "ctsB.txt",
            "--save-to", "ck3", "--save-at", str (day//2)]):
        return "%s: save@mid-warmup failed" % name
    if om (["-o", "outB.bin", "--ctsout", "ctsB.txt", "--resume", "ck3"]):
        return "%s: warmup resume failed" % name

    # truncated checkpoint must fail cleanly (nonzero exit, no crash)
    with open (os.path.join (simDir, "ck2"), "rb") as f:
        half=f.read (1000)
    with open (os.path.join (simDir, "ck2.trunc"), "wb") as f:
        f.write (half)
    if not om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck2.trunc"]):
        return "%s: truncated checkpoint accepted" % name

    # bumped format version must be rejected cleanly
    bad=open (os.path.join (simDir, "ck2"), "rb").read()
    open (os.path.join (simDir, "ck2.ver"), "wb").write (
        bad[:4] + struct.pack ("<I", 99) + bad[8:])
    if not om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck2.ver"]):
        return "%s: wrong format version accepted" % name

    rc=om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck2"])
    if rc: return "%s: final resume failed (rc %d)" % (name, rc)

    for a, b in (("full.bin", "out.bin"), ("full_cts.txt", "cts.txt"),
                 ("full.bin", "outB.bin"), ("full_cts.txt", "ctsB.txt")):
        fa=os.path.join (simDir, a)
        fb=os.path.join (simDir, b)
        if os.path.isfile (fa) != os.path.isfile (fb):
            return "%s: %s exists but %s doesn't" % (name, a, b)
        if not os.path.isfile (fa):
            continue    # scenario has no continuous output
        if open (fa, "rb").read() != open (fb, "rb").read():
            return "%s: %s != %s (split run diverged)" % (name, a, b)
    return None

def main():
    failed=False
    for name in sys.argv[1:]:
        err=check (name)
        if err:
            print ("FAIL %s" % err)
            failed=True
        else:
            print ("PASS %s" % name)
    return 1 if failed else 0

sys.exit (main())
