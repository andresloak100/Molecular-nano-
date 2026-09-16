"""Local, reproducible research workflow; computation requires explicit stages."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .candidates import create_design
from .design import load_design
from .quantum import QuantumSettings
from .workflow import run, audit_result, run_characterization


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
    except KeyboardInterrupt:
        print("Calculation interrupted; inspect its result.json and saved trajectories.", file=sys.stderr)
        return 130
    except (ValueError, OSError, RuntimeError, KeyError, TypeError) as exc:
        print(f"Calculation stopped: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
