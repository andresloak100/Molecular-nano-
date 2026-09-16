"""Record the bounded mock review with source hashes; starts no quantum jobs."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
REVIEW = Path(__file__).resolve().parent
SOURCES = [ROOT / "nanodesign/state_scan.py", REVIEW / "mock_backend.py",
           REVIEW / "test_execution_contract.py", Path(__file__).resolve()]


def hashes():
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in SOURCES}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New receipt path; never overwrites.")
    args = parser.parse_args()
    # Reserve the new receipt before starting, rather than running and then
    # discovering that its destination would overwrite an earlier review.
    with args.output.open("x") as destination:
        before = hashes()
        started = time.monotonic()
        command = [sys.executable, "-m", "pytest", str(REVIEW / "test_execution_contract.py"),
                   "-q", "--tb=short", "--disable-warnings"]
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=60)
        after = hashes()
        report = {
            "schema_version": 1, "checked_utc": datetime.now(timezone.utc).isoformat(),
            "scope": "S3 independent mock software-contract checks; no chemical validation",
            "quantum_jobs_started": 0, "command": command, "exit_code": result.returncode,
            "elapsed_seconds": time.monotonic() - started,
            "stdout": result.stdout, "stderr": result.stderr,
            "source_sha256_before": before, "source_sha256_after": after,
            "sources_changed_during_check": before != after,
        }
        json.dump(report, destination, indent=2, allow_nan=False)
        destination.write("\n")
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    return result.returncode if before == after else 2


if __name__ == "__main__":
    raise SystemExit(main())
