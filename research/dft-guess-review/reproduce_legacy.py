"""Reproduce default-field compatibility risk without changing production files.

Run from repository root using .venv/bin/python. Before the feature lands,
--simulate-added-default exercises the future settings serialization through a
temporary process-local patch. After it lands omit that option.
"""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import nanodesign.campaign as campaign
from nanodesign.design import sha256
from compatibility_fixtures import fixture, file_hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simulate-added-default", action="store_true")
    args = parser.parse_args()
    report = {
        "review": "DFT initial-guess legacy campaign compatibility",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_evidence": False, "quantum_evaluations": 0,
        "simulated_default_extension": args.simulate_added_default,
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": {name: sha256(ROOT / name) for name in
                          ("nanodesign/quantum.py", "nanodesign/campaign.py", "nanodesign/design.py")},
        "cases": [],
    }
    def expanded_asdict(value):
        result = asdict(value)
        result.setdefault("scf_initial_guess", "minao")
        return result
    for saved_status in (None, "completed", "failed", "not_converged"):
        with tempfile.TemporaryDirectory(prefix="dft-guess-review-") as temporary:
            root = fixture(temporary, saved_status=saved_status)
            before = file_hashes(root)
            case = {"saved_status": saved_status or "pending"}
            try:
                if args.simulate_added_default:
                    with patch.object(campaign, "asdict", expanded_asdict):
                        result = campaign.campaign_report(root)
                else:
                    result = campaign.campaign_report(root)
                case.update(readable=True, counts=result["counts"])
            except Exception as error:
                case.update(readable=False, error_type=type(error).__name__, error=str(error))
            case["all_fixture_files_unchanged"] = before == file_hashes(root)
            report["cases"].append(case)
    report["all_legacy_cases_readable"] = all(case["readable"] for case in report["cases"])
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
