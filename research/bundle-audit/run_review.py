"""Run the bounded B2 review and optionally preserve one new JSON receipt.

All checked bytes are synthetic, created in temporary directories. This never
changes source evidence or imports an electronic solver.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from nanodesign.bundle import export_bundle, verify_bundle
from fixtures import file_bytes, fingerprints, make_source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional new receipt file; existing files are never overwritten")
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error("Receipt already exists; use a new path to preserve earlier evidence")
    reviewed = [ROOT / "nanodesign/bundle.py", *sorted(LANE.glob("test_*.py")), LANE / "fixtures.py"]
    hashes = lambda: {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in reviewed}
    before = hashes()
    command = [sys.executable, "-m", "unittest", "discover", "-s", "research/bundle-audit", "-p", "test_*.py", "-v"]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=30)
    with tempfile.TemporaryDirectory(prefix="nanodesign-b2-") as temp:
        sample = Path(temp) / "synthetic-source"
        make_source(sample)
        source_before = fingerprints(sample)
        dest = Path(temp) / "portable"
        manifest = export_bundle(sample, dest)
        payload_equal = file_bytes(sample) == file_bytes(dest / "payload")
        verification = verify_bundle(dest)
        manifest_digest = hashlib.sha256((dest / "manifest.json").read_bytes()).hexdigest()
        source_after = fingerprints(sample)
    after = hashes()
    smoke_ok = payload_equal and verification["integrity_verified"] and source_before == source_after
    ok = result.returncode == 0 and before == after and smoke_ok
    report = {
        "review": "B2", "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": command, "exit_code": result.returncode,
        "stdout": result.stdout, "stderr": result.stderr,
        "source_sha256_before": before, "source_sha256_after": after,
        "review_source_unchanged": before == after,
        "synthetic_roundtrip": {
            "source": "nine literal B2 fixture files; no scientific archive",
            "file_count": manifest["file_count"], "total_bytes": manifest["total_bytes"],
            "manifest_sha256": manifest_digest,
            "source_before": source_before, "source_after": source_after,
            "exact_payload_bytes_preserved": payload_equal,
            "verification": verification,
            "export_location": "temporary directory, removed after checking",
        },
        "passed": ok, "solver_jobs_launched": 0,
        "scientific_validation": "not_assessed",
    }
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(encoded)
    print(result.stderr, end="")
    print(f"Synthetic roundtrip: {manifest['file_count']} files / {manifest['total_bytes']} bytes; passed={smoke_ok}")
    print(f"Reviewed source unchanged: {before == after}; combined passed={ok}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
