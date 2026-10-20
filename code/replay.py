#!/usr/bin/env python3
"""Run ordinary and optimized exact replays plus fail-closed controls."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
VERIFY = ROOT / "code" / "verify.py"
RESULTS = ROOT / "results"


def run(args: list[str], expected_success: bool) -> subprocess.CompletedProcess[str]:
    proc = subprocess.run(args, text=True, capture_output=True)
    if (proc.returncode == 0) != expected_success:
        raise RuntimeError(
            f"unexpected exit {proc.returncode} for {' '.join(args)}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    if not expected_success and "VerificationError" not in proc.stderr:
        raise RuntimeError(
            f"expected a VerificationError for {' '.join(args)}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc


def normalized(path: Path) -> dict:
    value = json.loads(path.read_text())
    value.pop("python_optimization_flag", None)
    return value


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="dys-replay-") as tmp:
        tmpdir = Path(tmp)
        normal = tmpdir / "normal.json"
        optimized = tmpdir / "optimized.json"
        run([sys.executable, "-B", str(VERIFY), "--output", str(normal)], True)
        run([sys.executable, "-B", "-O", str(VERIFY), "--output", str(optimized)], True)
        normal_report = json.loads(normal.read_text())
        optimized_report = json.loads(optimized.read_text())
        if normal_report.get("python_optimization_flag") != 0:
            raise RuntimeError("ordinary child did not report optimization flag 0")
        if optimized_report.get("python_optimization_flag") != 1:
            raise RuntimeError("optimized child did not report optimization flag 1")
        for control in ("trajectory", "constant"):
            run([sys.executable, "-B", str(VERIFY), "--negative-control", control], False)
            run([sys.executable, "-B", "-O", str(VERIFY), "--negative-control", control], False)
        if normalized(normal) != normalized(optimized):
            raise RuntimeError("ordinary and optimized semantic receipts differ")
        report = normal_report
        (RESULTS / "davis_yin_exact.json").write_text(json.dumps(report, indent=2) + "\n")

    receipt = {
        "ordinary_replay": "passed",
        "optimized_replay": "passed",
        "ordinary_negative_controls": ["trajectory_failed_as_required", "constant_failed_as_required"],
        "optimized_negative_controls": ["trajectory_failed_as_required", "constant_failed_as_required"],
        "semantic_receipts_equal_ignoring_optimization_flag": True,
        "scope": "exact arithmetic checks of formulas and explicit trajectories; analytic proof remains controlling",
    }
    (RESULTS / "replay_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
