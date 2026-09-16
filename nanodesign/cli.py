"""Local, reproducible research workflow; computation requires explicit stages."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import sys

from .candidates import create_design
from .design import load_design
from .quantum import QuantumSettings
from .workflow import run, audit_result, run_characterization


def _snapshot_file(path):
    """Bound and parse the same regular-file bytes that the report identifies."""
    from .electronic_state import DEFAULT_MAX_BYTES, read_snapshot_bytes
    path = Path(path).resolve(strict=True)
    descriptor = os.open(path, os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW)
    with os.fdopen(descriptor, "rb") as stream:
        if not stat.S_ISREG(os.fstat(stream.fileno()).st_mode):
            raise ValueError("Electronic snapshots must be regular files.")
        raw = stream.read(DEFAULT_MAX_BYTES + 1)
    snapshot = read_snapshot_bytes(raw)
    return snapshot, {"path": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
                      "size_bytes": len(raw)}


def main(argv=None):
    parser = argparse.ArgumentParser(description="Quantum-chemistry research for positional diamondoid assembly.")
    sub = parser.add_subparsers(dest="command", required=True)
    init = sub.add_parser("candidate", help="Generate an explicit diamondoid H-abstraction structure.")
    init.add_argument("--out", required=True)
    init.add_argument("--separation", type=float, default=3.6, help="Target carbon to tip-apex axial distance in Å")
    init.add_argument("--offset", type=float, default=0.0, help="Lateral tool offset in Å")
    check = sub.add_parser("check", help="Check geometry, charge/spin and mechanical boundary input.")
    check.add_argument("design")
    calc = sub.add_parser("calculate", help="Run actual DFT; outputs never overwrite a prior run.")
    calc.add_argument("design")
    calc.add_argument("--out", required=True)
    calc.add_argument("--stage", choices=["singlepoint", "relax", "path"], default="singlepoint")
    calc.add_argument("--state", choices=["initial", "final"], default="initial")
    calc.add_argument("--fmax", type=float, default=0.03, help="Convergence force in eV/Å")
    calc.add_argument("--steps", type=int, default=200)
    calc.add_argument("--images", type=int, default=7, help="Total NEB images including endpoints")
    audit = sub.add_parser("audit", help="Show numerical status and missing physical validation.")
    audit.add_argument("result")
    benchmark = sub.add_parser("benchmark", help="Compare energies at published reference geometries; this is not a calibrated tool prediction.")
    benchmark.add_argument("--out", required=True)
    benchmark.add_argument("--settings", help="JSON quantum settings, or a design JSON containing quantum settings")
    benchmark.add_argument("--reference-dir", help="Directory with the attributed published geometry package")
    paired = sub.add_parser("compare-methods", help="Compare DFT and CCSD(T) at identical small reference geometries and basis.")
    paired.add_argument("--out", required=True)
    paired.add_argument("--reference-dir")
    paired.add_argument("--basis", default="cc-pvdz", help="Currently cc-pVDZ only; coupled-cluster calculations are expensive")
    paired.add_argument("--cc-initial-guess", choices=["minao", "atom", "1e", "huckel"], default="minao",
                        help="Explicit HF starting guess; different converged solutions can give different reference energies")
    campaign_create = sub.add_parser("campaign-create", help="Snapshot a finite pose grid or existing designs; no quantum jobs are launched.")
    campaign_create.add_argument("designs", nargs="*", help="Existing design JSON files instead of a generated grid")
    campaign_create.add_argument("--out", required=True)
    campaign_create.add_argument("--separations", type=float, nargs="+")
    campaign_create.add_argument("--offsets", type=float, nargs="+", default=[0.0])
    campaign_create.add_argument("--settings", help="Quantum settings JSON for generated poses")
    campaign_create.add_argument("--stage", choices=["singlepoint", "relax", "path"], default="singlepoint")
    campaign_create.add_argument("--state", choices=["initial", "final"], default="initial")
    campaign_create.add_argument("--fmax", type=float, default=0.03)
    campaign_create.add_argument("--steps", type=int, default=200)
    campaign_create.add_argument("--images", type=int, default=7)
    campaign_run = sub.add_parser("campaign-run", help="Run a bounded number of jobs, preserving every attempt.")
    campaign_run.add_argument("directory")
    campaign_run.add_argument("--max-jobs", type=int, default=1)
    campaign_run.add_argument("--retry-incomplete", action="store_true", help="Explicitly restart failed/interrupted attempts into new output directories")
    campaign_report = sub.add_parser("campaign-report", help="Read campaign evidence without choosing an unvalidated winner.")
    campaign_report.add_argument("directory")
    scan_create = sub.add_parser("state-scan-create", help="Freeze a structure and explicit settings for a DFT starting-guess survey; no calculation.")
    scan_create.add_argument("structure")
    scan_create.add_argument("--out", required=True)
    scan_create.add_argument("--settings", required=True, help="Quantum settings JSON or design JSON; charge and spin come from these settings")
    scan_create.add_argument("--image", type=int, default=-1)
    scan_create.add_argument("--capture-electronic-state", action="store_true",
                             help="Save bounded converged-SCF orbital snapshots for each attempt; may require substantial disk space")
    scan_create.add_argument("--guesses", nargs="+", choices=["minao", "atom", "1e", "huckel"],
                             default=["minao", "atom", "1e", "huckel"])
    scan_run = sub.add_parser("state-scan-run", help="Evaluate a bounded number of pending DFT guesses, preserving every attempt.")
    scan_run.add_argument("directory")
    scan_run.add_argument("--max-jobs", type=int, default=1)
    scan_report = sub.add_parser("state-scan-report", help="Inspect saved guess dependence without selecting an electronic ground state.")
    scan_report.add_argument("directory")
    state_compare = sub.add_parser("state-compare", help="Compare saved occupied-orbital evidence at identical geometry; no state certification or calculations.")
    state_compare.add_argument("left")
    state_compare.add_argument("right")
    bundle_create = sub.add_parser("bundle-create", help="Copy evidence bytes to a new portable bundle; does not assess scientific validity.")
    bundle_create.add_argument("source")
    bundle_create.add_argument("--out", required=True)
    bundle_verify = sub.add_parser("bundle-verify", help="Check bundle contents and checksums without modifying evidence.")
    bundle_verify.add_argument("directory")
    for bundle_parser in (bundle_create, bundle_verify):
        for name in ("max_files", "max_file_bytes", "max_total_bytes", "max_entries", "max_depth"):
            bundle_parser.add_argument("--" + name.replace("_", "-"), type=int,
                                       help="Positive integer bound; uses the documented bundle default when omitted")
    characterize = sub.add_parser("characterize", help="Compute free-coordinate vibrational modes at a proposed stationary structure.")
    characterize.add_argument("design")
    characterize.add_argument("--structure", required=True, help="Optimized structure or trajectory; coordinates in Å")
    characterize.add_argument("--image", type=int, default=-1, help="Trajectory frame, default final frame")
    characterize.add_argument("--out", required=True)
    characterize.add_argument("--step", type=float, default=0.005, help="Central-difference displacement in Å")
    characterize.add_argument("--fmax", type=float, default=0.03)
    characterize.add_argument("--frequency-tolerance", type=float, default=20.0, help="Unresolved-mode threshold in cm^-1")
    characterize.add_argument("--max-free-coordinates", type=int, default=120, help="Cost guard: each free coordinate requires two quantum force evaluations")
    args = parser.parse_args(argv)
    try:
        if args.command == "candidate":
            print(create_design(args.out, args.separation, args.offset))
        elif args.command == "check":
            data, initial, _, settings, hashes = load_design(args.design)
            electrons = int(initial.numbers.sum()) - settings.charge
            if electrons < 1 or settings.spin > electrons or (electrons - settings.spin) % 2:
                raise ValueError("Charge and spin are inconsistent with the electron count.")
            print(json.dumps({"input_valid": True, "atoms": len(initial), "formula": initial.get_chemical_formula(),
                              "electrons": electrons, "spin_2S": settings.spin,
                              "fixed_atoms": len(data.get("fixed_indices", [])), "input_hashes": hashes,
                              "design_validated": False}, indent=2))
        elif args.command == "calculate":
            result = run(args.design, args.out, args.stage, args.state, args.fmax, args.steps, args.images)
            print(json.dumps({"status": result["status"], "hydrogen_basin_preserved": result.get("structure", {}).get("endpoint_identity_ok"), "topology_screen": result.get("structure", {}).get("topology_screen"), "result": str(Path(args.out).resolve() / "result.json"), "validation": result["validation"]}, indent=2))
            if result["status"] != "completed":
                return 2
        elif args.command == "audit":
            result = json.loads(Path(args.result).read_text())
            print(json.dumps(audit_result(result), indent=2))
        elif args.command == "benchmark":
            from .benchmark import run_benchmark
            settings_data = json.loads(Path(args.settings).read_text()) if args.settings else {}
            settings = QuantumSettings(**settings_data.get("quantum", settings_data))
            result = run_benchmark(settings, args.reference_dir, args.out)
            print(json.dumps({"result": str(Path(args.out).resolve() / "benchmark.json"),
                              "computed": result.get("computed"), "comparison": result.get("comparison"),
                              "method_validated": False}, indent=2))
        elif args.command == "characterize":
            result = run_characterization(args.design, args.structure, args.out, image=args.image,
                step=args.step, fmax=args.fmax, frequency_tolerance=args.frequency_tolerance,
                max_free_coordinates=args.max_free_coordinates)
            print(json.dumps({"result": str(Path(args.out).resolve() / "result.json"),
                              "status": result["status"], "classification": result["stationary"]["classification"],
                              "negative_mode_count": result["stationary"]["negative_mode_count"],
                              "unresolved_near_zero_mode_count": result["stationary"]["unresolved_near_zero_mode_count"],
                              "frequencies_cm1": result["stationary"]["frequencies_cm1"],
                              "design_validated": False}, indent=2))
        elif args.command == "compare-methods":
            from .method_comparison import run_method_comparison
            result = run_method_comparison(args.out, args.reference_dir, args.basis,
                                           cc_initial_guess=args.cc_initial_guess)
            print(json.dumps({"result": str(Path(args.out).resolve() / "method_comparison.json"),
                              "status": result["status"], "computed": result.get("computed"),
                              "method_validated": False}, indent=2))
        elif args.command == "campaign-create":
            from .campaign import create_campaign, create_pose_campaign
            controls = {name: getattr(args, name) for name in ("stage", "state", "fmax", "steps", "images")}
            if args.designs:
                if args.separations is not None or args.settings or args.offsets != [0.0]:
                    raise ValueError("Existing designs cannot be combined with generated-grid options.")
                result = create_campaign(args.designs, args.out, **controls)
            else:
                if args.separations is None:
                    raise ValueError("Supply existing designs or explicit --separations for a generated grid.")
                settings_data = json.loads(Path(args.settings).read_text()) if args.settings else {}
                settings = QuantumSettings(**settings_data.get("quantum", settings_data))
                result = create_pose_campaign(args.out, args.separations, args.offsets, settings=settings, **controls)
            print(json.dumps(result, indent=2))
        elif args.command == "campaign-run":
            from .campaign import run_campaign
            result = run_campaign(args.directory, max_jobs=args.max_jobs, retry_incomplete=args.retry_incomplete)
            print(json.dumps(result, indent=2))
            if any(result["counts"].get(status, 0) for status in ("failed", "not_converged", "abandoned", "interrupted")):
                return 2
        elif args.command == "campaign-report":
            from .campaign import campaign_report
            print(json.dumps(campaign_report(args.directory), indent=2))
        elif args.command == "state-scan-create":
            from .state_scan import create_state_scan, state_scan_report
            settings_data = json.loads(Path(args.settings).read_text())
            if not isinstance(settings_data, dict):
                raise ValueError("Settings must be a JSON object.")
            settings_data = settings_data.get("quantum", settings_data)
            if not isinstance(settings_data, dict) or not {"charge", "spin"}.issubset(settings_data):
                raise ValueError("State-scan settings must explicitly declare charge and spin.")
            settings = QuantumSettings(**settings_data)
            create_state_scan(args.structure, args.out, settings, image=args.image, guesses=args.guesses,
                              capture_electronic_state=args.capture_electronic_state)
            print(json.dumps(state_scan_report(args.out), indent=2))
        elif args.command == "state-scan-run":
            from .state_scan import run_state_scan
            result = run_state_scan(args.directory, max_jobs=args.max_jobs)
            print(json.dumps(result, indent=2))
            if any(result.get(name) for name in ("failed_guesses", "interrupted_guesses", "abandoned_guesses")):
                return 2
        elif args.command == "state-scan-report":
            from .state_scan import state_scan_report
            print(json.dumps(state_scan_report(args.directory), indent=2))
        elif args.command == "state-compare":
            from .electronic_state import compare_snapshots
            left, left_source = _snapshot_file(args.left)
            right, right_source = _snapshot_file(args.right)
            result = compare_snapshots(left, right)
            result["source_files"] = {"left": left_source, "right": right_source}
            print(json.dumps(result, indent=2, allow_nan=False))
            if result["comparison_status"] != "compared":
                return 2
        elif args.command in ("bundle-create", "bundle-verify"):
            from .bundle import BundleLimits, export_bundle, verify_bundle
            destination = args.out if args.command == "bundle-create" else args.directory
            try:
                limits = BundleLimits(**{name: getattr(args, name) for name in
                    ("max_files", "max_file_bytes", "max_total_bytes", "max_entries", "max_depth")
                    if getattr(args, name) is not None})
                if args.command == "bundle-create":
                    export_bundle(args.source, destination, limits=limits)
                result = verify_bundle(destination, limits=limits)
            except (ValueError, OSError, RuntimeError) as exc:
                result = {"integrity_verified": False, "scientific_validation": "not_assessed",
                          "reference_closure_verified": False, "file_count": 0, "total_bytes": 0,
                          "issues": [{"code": "command_failed", "path": str(destination), "message": str(exc)}]}
            print(json.dumps({"bundle": str(Path(destination).resolve()), "manifest_sha256": None, **result}, indent=2))
            if not result["integrity_verified"]:
                return 2
    except KeyboardInterrupt:
        print("Command interrupted; inspect saved evidence and incomplete output before continuing.", file=sys.stderr)
        return 130
    except (ValueError, OSError, RuntimeError, KeyError, TypeError) as exc:
        print(f"Calculation stopped: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
