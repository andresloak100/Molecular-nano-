"""Run E4 synthetic checks and preserve exact output plus inspected source hashes."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New receipt path; never overwrites")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output already exists; preserve previous review evidence")
    root = Path(__file__).resolve().parents[2]
    source_paths = [root / "nanodesign/electronic_state.py", root / "nanodesign/quantum.py"]
    reviewed_paths = sorted(Path(__file__).resolve().parent.glob("*.py"))

    def hashes(paths):
        return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}

    before = hashes(source_paths)
    command = [sys.executable, "-m", "pytest", "-q", "research/electronic-state-review"]
    result = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
    after = hashes(source_paths)
    receipt = {
        "schema_version": 1,
        "reviewed_utc": datetime.now(timezone.utc).isoformat(),
        "reviewer": "codex-712e / E4",
        "scope": "synthetic electronic-state math, capture conventions and evidence contracts",
        "command": command,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "source_sha256_before": before,
        "source_sha256_after": after,
        "source_unchanged": before == after,
        "test_sha256": hashes(reviewed_paths),
        "new_scf_calls": 0,
        "new_gradient_calls": 0,
        "new_ao_integral_calls": 0,
        "solver_objects": "synthetic fake mean fields; capture kernel is forbidden by fixture",
        "electronic_state_identity_verified": False,
        "ground_state_verified": False,
        "chemical_accuracy_validated": False,
        "limitation": "Snapshot arithmetic/context and synthetic producer behavior only; actual solver lifecycle and cross-geometry integral reconstruction are separately assigned reviews.",
    }
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    print(json.dumps({"receipt": str(args.output), "source_unchanged": before == after}))
    return result.returncode or (0 if before == after else 3)


if __name__ == "__main__":
    raise SystemExit(main())
