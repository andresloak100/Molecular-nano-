"""Per-evaluation cost of the candidate and its reduced variants, in CPU time.

The archived 53-atom timings are wall clock taken under contention nobody
recorded. They are upper bounds of unknown tightness. This script measures the
same quantity in *process CPU time*, which for a single-threaded PySCF build is
insensitive to how many other jobs share the host, so the numbers remain
comparable to each other and to a future uncontended run.

One energy-and-gradient evaluation per system, which is exactly the unit a NEB
optimizer step consumes. Nothing is relaxed. The point is cost, not chemistry:
a cheaper model is not a better model, and this file makes no claim about
whether any of these systems describes the reaction correctly.

Systems, all at the same pose, same method, same anchors where applicable:

    hydrogen handle    HCC* + adamantane        29 atoms
    methyl handle      CH3-CC* + adamantane     32 atoms
    adamantyl handle   the 53-atom candidate    53 atoms

Bounded: one system at a time, results written after each, and the 53-atom
system only when explicitly requested since it is the expensive one.

Run from the repository root:

    python research/candidate-feasibility/gradient_cost.py --handles hydrogen methyl
    python research/candidate-feasibility/gradient_cost.py --handles adamantyl
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

from nanodesign.quantum import PySCFCalculator, QuantumSettings  # noqa: E402
from reduced_models import HANDLES, reduced_candidate  # noqa: E402
from timing import load_snapshot, timed  # noqa: E402

EVIDENCE = LANE / "evidence"


def measure(handle: str, settings: QuantumSettings, event_log: Path) -> dict:
    atoms, info = reduced_candidate(handle)
    print(f"{handle}: {atoms.get_chemical_formula()}, {len(atoms)} atoms")
    record: dict = {
        "handle": handle,
        "formula": atoms.get_chemical_formula(),
        "n_atoms": len(atoms),
        "fixed_indices": info["fixed_indices"],
        "anchor_scheme": info["anchor_scheme"],
        "settings": settings.to_dict(),
        "geometry_status": info["geometry_status"],
    }
    atoms.calc = PySCFCalculator(settings, event_log=event_log)
    with timed(record, "timing"):
        energy = float(atoms.get_potential_energy())
        forces = atoms.get_forces()
    diagnostics = atoms.calc.diagnostics
    record.update(
        energy_ev=energy,
        max_free_force_ev_per_angstrom=float(
            np.linalg.norm(atoms.get_forces(apply_constraint=True), axis=1).max()
        ),
        max_raw_force_ev_per_angstrom=float(np.linalg.norm(forces, axis=1).max()),
        basis_functions=diagnostics["basis_functions"],
        scf_cycles=diagnostics["scf_cycles"],
        s2=diagnostics["s2"],
        effective_threads=diagnostics.get("effective_pyscf_threads"),
        threads_honored=diagnostics.get("threads_honored"),
        scf_converged=diagnostics["scf_converged"],
    )
    timing = record["timing"]
    record["cpu_seconds_per_energy_and_gradient"] = timing["cpu_seconds"]
    print(f"    {record['basis_functions']} basis functions, {record['scf_cycles']} SCF cycles")
    print(f"    cpu {timing['cpu_seconds']:.1f} s   wall {timing['wall_seconds']:.1f} s"
          f"   contention x{timing['contention_factor_wall_over_cpu']:.1f}"
          f"   load {timing['load_at_start']['load_average_1min']:.0f}")
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handles", nargs="+", default=["hydrogen", "methyl"], choices=list(HANDLES))
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--xc", default="pbe0")
    parser.add_argument("--no-density-fit", action="store_true")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    EVIDENCE.mkdir(parents=True, exist_ok=True)
    out = EVIDENCE / (args.out or f"gradient-cost-{'-'.join(args.handles)}.json")
    if out.exists():
        raise SystemExit(f"{out} exists; evidence is never overwritten. Choose --out.")

    settings = QuantumSettings(
        xc=args.xc, basis=args.basis, dispersion="d3bj",
        density_fit=not args.no_density_fit, threads=1,
    )
    report: dict = {
        "schema_version": 1,
        "question": "What does one energy+gradient evaluation cost for the candidate and its reduced variants?",
        "unit": "one energy-and-gradient evaluation, the unit a NEB optimizer step consumes",
        "method": f"{args.xc.upper()}-D3(BJ)/{args.basis}"
                  + (", density fitting" if not args.no_density_fit else ", no density fitting"),
        "timing_basis": (
            "Process CPU time is the transferable cost on this single-threaded "
            "build. Wall clock on this host today is inflated by contention "
            "from other lanes and is recorded alongside, not instead."
        ),
        "measures": "cost only; says nothing about whether any model describes the reaction correctly",
        "host_at_start": load_snapshot(),
        "systems": {},
        "status": "running",
    }
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    print(f"method: {report['method']}   host load {report['host_at_start']['load_average_1min']:.1f} "
          f"on {report['host_at_start']['logical_cores']} cores")
    for handle in args.handles:
        report["systems"][handle] = measure(handle, settings, EVIDENCE / f"gradient-cost-{handle}.jsonl")
        out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    if len(report["systems"]) > 1:
        costs = {h: s["cpu_seconds_per_energy_and_gradient"] for h, s in report["systems"].items()}
        sizes = {h: s["basis_functions"] for h, s in report["systems"].items()}
        report["scaling"] = {
            "cpu_seconds_per_evaluation": costs,
            "basis_functions": sizes,
            "note": (
                "Two or three points cannot establish a scaling exponent. These "
                "are the measured costs of these specific systems, not a law to "
                "extrapolate with."
            ),
        }
    report["status"] = "completed"
    report["host_at_end"] = load_snapshot()
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
