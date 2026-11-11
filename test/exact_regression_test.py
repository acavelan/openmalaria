"""Exercise the exact regression gate without running scientific simulations."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import exact_regression


class ExactRegressionTest(unittest.TestCase):
    def test_capture_accepts_equal_bytes_and_rejects_changed_outputs_or_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "openMalaria"
            executable.write_text("test executable")
            results = {"direct": {"output.txt": "original", "ctsout.txt": "continuous"},
                       "resumed": {"output.txt": "original", "ctsout.txt": "continuous"},
                       "checkpoint": "checkpoint"}
            tests = {"tests": [{"name": "unittest", "command": ["unittest"]},
                               {"name": "1", "command": ["python3", "/build/test/run.py", "1"]}]}

            def run(directory, baseline=None):
                argv = ["exact_regression.py", "--build", str(root),
                        "--output", str(root / directory)]
                if baseline:
                    argv += ["--baseline", str(root / baseline)]
                with patch.object(sys, "argv", argv), \
                        patch.object(exact_regression.subprocess, "check_output", return_value=json.dumps(tests)), \
                        patch.object(exact_regression, "run_scenario", return_value=results) as simulation, \
                        contextlib.redirect_stdout(io.StringIO()):
                    status = exact_regression.main()
                    self.assertEqual(simulation.call_count, 1)
                    self.assertEqual(simulation.call_args.args[1], "1")
                    return status

            self.assertFalse(run("baseline"))
            self.assertFalse(run("equal", "baseline"))
            results["direct"]["output.txt"] = "one changed byte"
            self.assertTrue(run("changed-output", "baseline"))
            results["direct"]["output.txt"] = "original"
            results["resumed"]["output.txt"] = "changed after resume"
            self.assertTrue(run("changed-resumed-output", "baseline"))
            results["resumed"]["output.txt"] = "original"
            results["checkpoint"] = "changed RNG state"
            self.assertTrue(run("changed-checkpoint", "baseline"))

    def test_missing_output_is_a_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = type("Args", (), dict(build=root, output=root, executable=root / "unused",
                                        format="text", timeout=1))()
            with patch.object(exact_regression.subprocess, "run"):
                with self.assertRaisesRegex(RuntimeError, "missing output.txt"):
                    exact_regression.run_scenario(args, "1")


if __name__ == "__main__":
    unittest.main()
