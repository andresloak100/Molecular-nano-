"""Local, reproducible research workflow; computation requires explicit stages."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

from ase.io import write

from .candidates import make_h_abstraction
from .design import load_design
from .quantum import QuantumSettings
from .workflow import json_write, run, audit_result, run_characterization


def create_design(output, separation, offset):
    initial, final, metadata = make_h_abstraction(separation, offset)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    write(out / "initial.xyz", initial)
    write(out / "final.xyz", final)
    design = {
        "schema_version": 1,
        "length_unit": "angstrom",
        "name": "Adamantane-supported ethynyl H-abstraction candidate",
        "scope": "Finite diamondoid cluster in vacuum with fixed distal carbon anchors. Unrelaxed candidate, not a diamond surface or validated assembly tool.",
        "initial": "initial.xyz", "final": "final.xyz",
        "fixed_indices": metadata["fixed_indices"],
        "hydrogen_transfer": {"donor": metadata["target_carbon"], "hydrogen": metadata["transferred_hydrogen"], "acceptor": metadata["tip_apex"]},
        "quantum": asdict(QuantumSettings()),
        "metadata": metadata,
    }
    json_write(out / "design.json", design)
    return out / "design.json"


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
    except KeyboardInterrupt:
        print("Calculation interrupted; inspect its result.json and saved trajectories.", file=sys.stderr)
        return 130
    except (ValueError, OSError, RuntimeError, KeyError, TypeError) as exc:
        print(f"Calculation stopped: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
