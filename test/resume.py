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

def checkpointDay (path):
    # header: u32 magic, u32 version, u8 phase, then t0 (sim days) at offset 9
    with open (path, "rb") as f:
        head=f.read (13)
    magic, version, phase = struct.unpack ("<IIB", head[:9])
    assert magic==0x4B634D4F and version==1, "not a checkpoint"
    return struct.unpack ("<i", head[9:13])[0]

def check (name):
    scenario=os.path.join (testSrcDir, "scenario%s.xml" % name)
    assert os.path.isfile (scenario), scenario
    # The generated schema must sit next to the scenario copy's working dir.
    import re
    with open (scenario) as f:
        m=re.search (r'schemaLocation="[^"]* (\S+\.xsd)"', f.read())
    xsd=os.path.abspath (os.path.join (testBuildDir, "..", "schema", m.group(1)))
    assert os.path.isfile (xsd), xsd
    base=["--scenario", os.path.abspath (scenario), "--seed", "112",
          "--output-format", "bin"]
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
    day=checkpointDay (os.path.join (simDir, "ck1"))

    # resume to mid-main, checkpoint again, resume again (chain: the resumed
    # run must be able to save and be resumed itself)
    rc=om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck1",
            "--save-to", "ck2", "--save-at", str (day+50)])
    if rc: return "%s: resume+save@mid-main failed (rc %d)" % (name, rc)

    # truncated checkpoint must fail cleanly (nonzero exit, no crash)
    with open (os.path.join (simDir, "ck2"), "rb") as f:
        half=f.read (1000)
    with open (os.path.join (simDir, "ck2.trunc"), "wb") as f:
        f.write (half)
    if not om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck2.trunc"]):
        return "%s: truncated checkpoint accepted" % name

    rc=om (["-o", "out.bin", "--ctsout", "cts.txt", "--resume", "ck2"])
    if rc: return "%s: final resume failed (rc %d)" % (name, rc)

    for a, b in (("full.bin", "out.bin"), ("full_cts.txt", "cts.txt")):
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
