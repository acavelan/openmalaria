#!/usr/bin/env python3
"""Capture or compare exact scenario results using a fixed executable/toolchain.

Example (run capture before rebuilding):
  python3 test/exact_regression.py --build build --output /tmp/om-before
  python3 test/exact_regression.py --build build --output /tmp/om-after \
      --baseline /tmp/om-before

Runs both uninterrupted and checkpoint-resumed simulations. Keeps logs, outputs,
and a manifest in a new directory; never updates the repository's fixtures.
Each execution mode is compared with its own baseline, since existing scenarios
may differ between uninterrupted and resumed execution at binary precision.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def digest(data):
    return hashlib.sha256(data).hexdigest()


def run_scenario(args, name):
    source = Path(__file__).resolve().parent
    scenario = Path(name).resolve() if Path(name).is_file() else source / f"scenario{name}.xml"
    label = scenario.stem if Path(name).is_file() else name
    location = ET.parse(scenario).getroot().attrib[
        "{http://www.w3.org/2001/XMLSchema-instance}schemaLocation"]
    schema_name = location.split()[-1]
    schema = source.parent / "schema" / schema_name
    if not schema.is_file():
        schema = args.build / "schema" / schema_name
    results = {}
    for mode in ("direct", "resumed"):
        directory = args.output / label / mode
        directory.mkdir(parents=True)
        (directory / schema_name).symlink_to(schema)
        command = [str(args.executable), "--resource-path", str(source),
                   "--scenario", str(scenario), "--output-format",
                   "bin" if args.format == "bin" else "txt"]
        if mode == "resumed":
            command.append("--checkpoint-stop")
        with (directory / "stdout.log").open("wb") as out, \
                (directory / "stderr.log").open("wb") as err:
            for _ in range(2 if mode == "resumed" else 1):
                subprocess.run(command, cwd=directory, stdout=out, stderr=err,
                               check=True, timeout=args.timeout)
        files = {}
        for stem in ("output.bin" if args.format == "bin" else "output.txt", "ctsout.txt"):
            path = directory / stem
            if path.is_file():
                data = path.read_bytes()
            elif path.with_suffix(path.suffix + ".gz").is_file():
                data = gzip.decompress(path.with_suffix(path.suffix + ".gz").read_bytes())
                path.write_bytes(data)
            elif stem.startswith("output"):
                raise RuntimeError(f"{name}/{mode}: missing {stem}")
            else:
                continue
            files[stem] = digest(data)
        results[mode] = files
        if mode == "resumed":
            checkpoint = (directory / "checkpoint").read_text().strip()
            payload = gzip.decompress((directory / f"checkpoint{checkpoint}.gz").read_bytes())
            results["checkpoint"] = digest(payload)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--format", choices=("text", "bin"), default="text")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("scenarios", nargs="*")
    args = parser.parse_args()
    args.build = args.build.resolve()
    args.executable = (args.executable or args.build / "openMalaria").resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    tests = json.loads(subprocess.check_output(
        ["ctest", "--test-dir", str(args.build), "--show-only=json-v1"]))
    names = args.scenarios or [t["name"] for t in tests["tests"]
                              if any(Path(arg).name == "run.py" for arg in t.get("command", []))]
    if not names:
        parser.error("no scenario tests found; enable OM_BOXTEST_ENABLE or name scenarios explicitly")
    baseline = None
    if args.baseline:
        baseline = json.loads((args.baseline / "manifest.json").read_text())
        if baseline["format"] != args.format:
            parser.error("baseline output format differs")
    manifest = {"executable_sha256": digest(args.executable.read_bytes()),
                "format": args.format, "results": {}, "errors": {}}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {pool.submit(run_scenario, args, name): name for name in names}
        for future in as_completed(pending):
            name = pending[future]
            try:
                result = future.result()
                manifest["results"][name] = result
                if baseline and baseline["results"].get(name) != result:
                    raise RuntimeError("outputs or checkpoint payload differ from baseline")
                print(f"PASS {name}", flush=True)
                if result["direct"] != result["resumed"]:
                    print(f"NOTE {name}: uninterrupted and resumed outputs differ", flush=True)
            except Exception as error:
                manifest["errors"][name] = str(error)
                print(f"FAIL {name}: {error}", flush=True)
            (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return bool(manifest["errors"])


if __name__ == "__main__":
    raise SystemExit(main())
