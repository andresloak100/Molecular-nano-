"""Exercise documented user entrypoints without launching quantum chemistry.

Uses an already installed environment, isolated temporary working directories,
and a strict command allowlist. It does not test a fresh dependency installation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
COMMANDS = ("candidate", "check", "calculate", "audit", "benchmark",
            "compare-methods", "campaign-create", "campaign-run",
            "campaign-report", "characterize")
READ_OR_CREATE = {"candidate", "check", "audit", "campaign-create", "campaign-report"}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="New evidence directory")
    parser.add_argument("--executable", type=Path, default=Path(sys.executable).parent / "nanodesign")
    args = parser.parse_args()
    executable = args.executable.absolute()
    if not executable.is_file():
        parser.error(f"Installed entrypoint not found: {executable}")
    args.out.mkdir(parents=True, exist_ok=False)
    source_paths = [ROOT / name for name in (
        "README.md", "pyproject.toml", "docs/CAMPAIGNS.md", "nanodesign/cli.py",
        "nanodesign/candidates.py", "nanodesign/campaign.py", "nanodesign/quantum.py",
        "nanodesign/workflow.py", "examples/pose-campaign/README.md",
    )]
    before = {str(p.relative_to(ROOT)): digest(p) for p in source_paths}
    report = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "executable": str(executable),
        "source_sha256_before": before,
        "scope": "Installed environment; CLI help/create/check/audit/report only. No fresh install, solver, path, numerical-accuracy or UI verification.",
        "packages": {name: importlib.metadata.version(name) for name in (
            "molecular-nanodesign", "numpy", "scipy", "ase", "pyscf", "dftd3")},
        "cases": [], "artifact_checks": [],
    }
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")

    def run(name, argv, cwd, *, expected=0, module=False):
        if argv and "--help" not in argv and argv[0] not in READ_OR_CREATE:
            raise ValueError(f"Non-help command not permitted: {argv}")
        command = ([sys.executable, "-m", "nanodesign"] if module else [str(executable)]) + list(map(str, argv))
        case = {"name": name, "argv": command, "cwd": str(cwd), "expected_exit": expected}
        try:
            result = subprocess.run(command, cwd=cwd, env=env, text=True, capture_output=True, timeout=30)
            case.update(exit_code=result.returncode, stdout=result.stdout, stderr=result.stderr,
                        passed=result.returncode == expected)
        except subprocess.TimeoutExpired as exc:
            case.update(exit_code=None, stdout=str(exc.stdout or ""), stderr=str(exc.stderr or ""), passed=False)
        report["cases"].append(case)
        return case

    def check(name, passed, observed):
        report["artifact_checks"].append({"name": name, "passed": bool(passed), "observed": observed})

    def parsed(case):
        try:
            return json.loads(case["stdout"])
        except (KeyError, ValueError):
            return {}

    with tempfile.TemporaryDirectory(prefix="nanodesign-r1-") as directory:
        work = Path(directory)
        run("top_level_help", ["--help"], work)
        run("module_help_outside_checkout", ["--help"], work, module=True)
        for command in COMMANDS:
            run(f"help_{command}", [command, "--help"], work)
        run("missing_command_usage", [], work, expected=2)
        run("candidate_readme", ["candidate", "--out", "designs/h-abstraction", "--separation", "3.6", "--offset", "0.0"], work)
        design_dir = work / "designs/h-abstraction"
        candidate_files = sorted(p.name for p in design_dir.iterdir()) if design_dir.exists() else []
        check("candidate_files", candidate_files == ["design.json", "final.xyz", "initial.xyz"], candidate_files)
        candidate = parsed(run("check_created_design", ["check", "designs/h-abstraction/design.json"], work))
        check("candidate_check_summary", all(candidate.get(k) == v for k, v in {
            "input_valid": True, "atoms": 53, "formula": "C22H31", "electrons": 163,
            "spin_2S": 1, "fixed_atoms": 6, "design_validated": False}.items()), candidate)
        run("check_supplied_example", ["check", ROOT / "examples/h-abstraction/design.json"], work)
        run("existing_candidate_refused", ["candidate", "--out", "designs/h-abstraction"], work, expected=2)
        run("missing_design_clear_error", ["check", "missing/design.json"], work, expected=2)
        run("missing_audit_clear_error", ["audit", "missing/result.json"], work, expected=2)
        audit = parsed(run("audit_archived_singlepoint", ["audit", ROOT / "data/validation/h-abstraction-df-initial/result.json"], work))
        check("audit_does_not_validate_design", audit.get("design_validated") is False and bool(audit.get("missing_evidence")), audit)
        create = ["campaign-create", "--out", "campaigns/pose-study", "--separations", "3.4", "3.6", "3.8", "--offsets", "-0.2", "0.0", "0.2", "--stage", "singlepoint"]
        run("campaign_create_readme_grid", create, work)
        generated = parsed(run("campaign_report_created_grid", ["campaign-report", "campaigns/pose-study"], work))
        check("generated_grid_pending", generated.get("counts", {}).get("pending") == 9, generated.get("counts"))
        campaign_dir = work / "campaigns/pose-study"
        required = ["plan.json", "campaign.json", "designs/pose-0001/design.json", "designs/pose-0001/source-design.json", "designs/pose-0009/design.json"]
        check("documented_campaign_artifacts", all((campaign_dir / p).is_file() for p in required), {p: (campaign_dir / p).is_file() for p in required})
        run("existing_campaign_refused", create, work, expected=2)
        run("campaign_create_existing_design", ["campaign-create", "designs/h-abstraction/design.json", "--out", "campaigns/from-design"], work)
        existing = parsed(run("campaign_report_existing_design", ["campaign-report", "campaigns/from-design"], work))
        check("existing_design_pending", existing.get("counts", {}).get("pending") == 1, existing.get("counts"))
        copied = work / "campaigns/copied-example"
        shutil.copytree(ROOT / "examples/pose-campaign", copied)
        copied_report = parsed(run("campaign_report_copied_example", ["campaign-report", copied], work))
        check("copied_example_pending", copied_report.get("counts", {}).get("pending") == 9, copied_report.get("counts"))
        created_quantum_outputs = sorted(str(p.relative_to(work)) for pattern in ("electronic.jsonl", "result.json") for p in work.rglob(pattern))
        check("no_quantum_output_created", not created_quantum_outputs, created_quantum_outputs)

    report["source_sha256_after"] = {str(p.relative_to(ROOT)): digest(p) for p in source_paths}
    report["source_changed_during_audit"] = before != report["source_sha256_after"]
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    all_checks = report["cases"] + report["artifact_checks"]
    report["summary"] = {"passed": sum(c["passed"] for c in all_checks), "failed": sum(not c["passed"] for c in all_checks), "command_count": len(report["cases"]), "artifact_check_count": len(report["artifact_checks"])}
    (args.out / "results.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str((args.out / "results.json").absolute()), **report["summary"], "source_changed_during_audit": report["source_changed_during_audit"]}, indent=2))
    return 1 if report["summary"]["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
