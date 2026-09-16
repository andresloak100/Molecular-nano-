"""What the cheap speed knobs actually buy, and what they cost.

A2 item 4. Three knobs are available to make the path affordable: density
fitting, a smaller basis, and a looser force tolerance. Each is usually
defended by assertion. This measures them on a system small enough to run both
ways, so the trade is a number instead of an opinion.

The test system is the hydrogen-handle reduced candidate (29 atoms): the real
substrate and the real tip at the real pose, with the tool handle replaced by a
hydrogen. It is the cheapest system that still contains the bond being broken.

Knob 1, density fitting. The archived 53-atom comparison put direct at 2044.8 s
against 701.5 s density-fitted, but those two ran concurrently, so the ratio
mixes the approximation with the scheduler. Here both run in the same process,
one after the other, timed in CPU seconds, which is a controlled comparison.

Knob 2, basis. def2-SVP against def2-TZVP: cost ratio, and how far the SVP
forces and relative energy sit from the larger basis at the same geometry.

Knob 3, force tolerance. Relaxation step count at fmax 0.03, 0.05 and 0.10
eV/A, from the same starting structure. Steps are what the NEB cost model
multiplies by, so this is the knob the projection is most sensitive to. The
looser tolerances are run first and each tighter run continues from the
previous structure, which is how a person would actually do it; the step counts
are therefore cumulative-to-reach, and reported that way.

None of this makes any statement about whether the reaction is real.

Run from the repository root:

    python research/candidate-feasibility/knobs.py --knob density-fitting
    python research/candidate-feasibility/knobs.py --knob basis
    python research/candidate-feasibility/knobs.py --knob tolerance
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
from ase.optimize import FIRE

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

from nanodesign.quantum import PySCFCalculator, QuantumSettings  # noqa: E402
from reduced_models import reduced_candidate  # noqa: E402
from timing import load_snapshot, timed  # noqa: E402

EVIDENCE = LANE / "evidence"
EV_KCAL = 23.060548
TEST_HANDLE = "hydrogen"


def _evaluate(settings: QuantumSettings, label: str) -> dict:
    """One energy+gradient on the test system with the given settings."""
    atoms, _ = reduced_candidate(TEST_HANDLE)
    atoms.calc = PySCFCalculator(settings, event_log=EVIDENCE / f"knobs-{label}.jsonl")
    record: dict = {"label": label, "settings": settings.to_dict(), "n_atoms": len(atoms)}
    with timed(record, "timing"):
        energy = float(atoms.get_potential_energy())
        forces = atoms.get_forces()
    diagnostics = atoms.calc.diagnostics
    record.update(
        energy_ev=energy,
        forces_ev_per_angstrom=forces.tolist(),
        max_force_ev_per_angstrom=float(np.linalg.norm(forces, axis=1).max()),
        basis_functions=diagnostics["basis_functions"],
        scf_cycles=diagnostics["scf_cycles"],
        s2=diagnostics["s2"],
        cpu_seconds=record["timing"]["cpu_seconds"],
    )
    print(f"  {label:26s} {record['basis_functions']:4d} bf  cpu {record['cpu_seconds']:8.1f} s"
          f"  E = {energy:.6f} eV")
    return record


def _force_deviation(a: dict, b: dict) -> dict:
    fa = np.array(a["forces_ev_per_angstrom"])
    fb = np.array(b["forces_ev_per_angstrom"])
    difference = fb - fa
    return {
        "max_component_ev_per_angstrom": float(np.abs(difference).max()),
        "rms_component_ev_per_angstrom": float(np.sqrt((difference ** 2).mean())),
        "max_vector_norm_ev_per_angstrom": float(np.linalg.norm(difference, axis=1).max()),
        "optimization_tolerance_ev_per_angstrom": 0.03,
    }


def knob_density_fitting(base: QuantumSettings) -> dict:
    print("density fitting, both arms in one process (controlled):")
    direct = _evaluate(replace(base, density_fit=False), "direct")
    fitted = _evaluate(replace(base, density_fit=True), "density-fitted")
    deviation = _force_deviation(direct, fitted)
    speedup = direct["cpu_seconds"] / fitted["cpu_seconds"]
    return {
        "knob": "density_fitting",
        "controlled": True,
        "direct": direct,
        "density_fitted": fitted,
        "cpu_speedup_direct_over_fitted": speedup,
        "force_deviation": deviation,
        "energy_difference_ev": fitted["energy_ev"] - direct["energy_ev"],
        "energy_difference_kcal_per_mol": (fitted["energy_ev"] - direct["energy_ev"]) * EV_KCAL,
        "verdict_basis": (
            "Buys a CPU-time factor of "
            f"{speedup:.2f} at a maximum force-component deviation of "
            f"{deviation['max_component_ev_per_angstrom']:.2e} eV/A against a 0.03 "
            "tolerance. Measured at one geometry on a 29-atom model; the "
            "absolute energy shift does not cancel between different systems "
            "and is only safe in energy differences at matched settings."
        ),
    }


def knob_basis(base: QuantumSettings) -> dict:
    print("basis, same geometry:")
    small = _evaluate(replace(base, basis="def2-svp", density_fit=True), "def2-svp")
    large = _evaluate(replace(base, basis="def2-tzvp", density_fit=True), "def2-tzvp")
    deviation = _force_deviation(small, large)
    return {
        "knob": "basis",
        "def2_svp": small,
        "def2_tzvp": large,
        "cpu_cost_ratio_tzvp_over_svp": large["cpu_seconds"] / small["cpu_seconds"],
        "force_deviation": deviation,
        "verdict_basis": (
            "Absolute energies at different bases are not comparable; only the "
            "forces and the cost ratio are read here. A force deviation of "
            f"{deviation['max_component_ev_per_angstrom']:.3f} eV/A against a 0.03 "
            "tolerance indicates whether def2-SVP geometries can be trusted "
            "as starting points for def2-TZVP energies."
        ),
    }


def knob_tolerance(base: QuantumSettings, tolerances=(0.10, 0.05, 0.03)) -> dict:
    """Relaxation step counts to reach successively tighter force tolerances."""
    print("force tolerance, cumulative relaxation from one starting structure:")
    settings = replace(base, density_fit=True)
    atoms, info = reduced_candidate(TEST_HANDLE)
    atoms.calc = PySCFCalculator(settings, event_log=EVIDENCE / "knobs-tolerance.jsonl")
    stages = []
    cumulative = 0
    for tolerance in tolerances:
        record: dict = {"fmax_ev_per_angstrom": tolerance}
        optimizer = FIRE(
            atoms,
            trajectory=str(EVIDENCE / f"knobs-tolerance-{tolerance:.2f}.traj"),
            logfile=str(EVIDENCE / f"knobs-tolerance-{tolerance:.2f}.log"),
            maxstep=0.08,
        )
        with timed(record, "timing"):
            converged = bool(optimizer.run(fmax=tolerance, steps=200))
        steps = int(optimizer.get_number_of_steps())
        cumulative += steps
        record.update(
            converged=converged,
            steps_this_stage=steps,
            cumulative_steps_from_start=cumulative,
            energy_ev=float(atoms.get_potential_energy()),
            max_free_force_ev_per_angstrom=float(
                np.linalg.norm(atoms.get_forces(apply_constraint=True), axis=1).max()
            ),
            cpu_seconds=record["timing"]["cpu_seconds"],
        )
        print(f"  fmax {tolerance:.2f}: {steps:3d} steps this stage, {cumulative:3d} cumulative,"
              f" converged={converged}, cpu {record['cpu_seconds']:.1f} s")
        stages.append(record)

    return {
        "knob": "force_tolerance",
        "system": {"handle": TEST_HANDLE, "n_atoms": len(atoms), "anchor_scheme": info["anchor_scheme"]},
        "stages": stages,
        "reading": (
            "Steps are cumulative from the same starting structure, each stage "
            "continuing from the previous one. A single endpoint relaxation on "
            "this reduced model is what these counts describe; the 53-atom "
            "candidate has more degrees of freedom and will generally need more."
        ),
        "why_it_matters": (
            "The NEB projection multiplies per-evaluation cost by step counts, "
            "so the tolerance is the knob the projection is most sensitive to. "
            "A looser tolerance leaves residual forces on the structure, which "
            "for a saddle search means a worse starting band, not just a less "
            "converged endpoint."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--knob", required=True, choices=["density-fitting", "basis", "tolerance"])
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    EVIDENCE.mkdir(parents=True, exist_ok=True)
    out = EVIDENCE / (args.out or f"knob-{args.knob}.json")
    if out.exists():
        raise SystemExit(f"{out} exists; evidence is never overwritten. Choose --out.")

    base = QuantumSettings(xc="pbe0", basis="def2-svp", dispersion="d3bj", threads=1)
    report: dict = {
        "schema_version": 1,
        "question": "What does each speed knob buy, and what does it cost?",
        "test_system": f"{TEST_HANDLE}-handle reduced candidate, 29 atoms, real substrate and tip at the real pose",
        "timing_basis": "process CPU time; this host is heavily contended and wall clock is not comparable",
        "host_at_start": load_snapshot(),
        "status": "running",
    }
    print(f"host load {report['host_at_start']['load_average_1min']:.1f} on "
          f"{report['host_at_start']['logical_cores']} cores\n")

    if args.knob == "density-fitting":
        report["result"] = knob_density_fitting(base)
    elif args.knob == "basis":
        report["result"] = knob_basis(base)
    else:
        report["result"] = knob_tolerance(base)

    report["status"] = "completed"
    report["host_at_end"] = load_snapshot()
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
