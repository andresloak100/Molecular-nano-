"""What shrinking the tool handle costs, measured rather than asserted.

A2 item 3 asks for the smallest faithful substitute for the 53-atom candidate.
The obvious reduction is the handle: the second adamantane cage is 27 of the
53 atoms and most of the cost. Whether that reduction is faithful is an
empirical question with at least two separable parts.

**Measured here: the electronic substituent effect.** How much the handle
changes the apex radical's appetite for a hydrogen, as the tip H-affinity

    A(R) = E(R-CC-H) - E(R-CC*) - E(H)

for R = H, methyl, adamantyl. Only differences across R matter, and E(H)
cancels in those differences; it is computed anyway so the absolute numbers
are interpretable. A small spread across R is evidence that the handle's
electronic role is minor. It is not evidence that the reduction is safe.

**Not measured here: mechanical and steric fidelity.** The handle also sets
the mount's stiffness, its mass, its sterics against the substrate, and where
the anchors sit. A single fixed atom on a methyl handle is a different
boundary condition from three fixed carbons in a cage, not a scaled-down one.
Nothing in this file addresses that, and the report says so.

Every open-shell species gets a four-guess SCF scan before its production
calculation, because this project has already published a wrong number from
the default guess landing on a higher stable solution. The scan is done at
SCF level only, then the production calculator is run with the explicitly
selected guess so the archived record carries full provenance.

Bounded by design: one process, one species at a time, results written after
each species so an interrupted run still leaves usable evidence.

Run from the repository root:

    python research/candidate-feasibility/handle_fidelity.py --handles hydrogen methyl
    python research/candidate-feasibility/handle_fidelity.py --handles adamantyl
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

from nanodesign.quantum import PySCFCalculator, QuantumSettings  # noqa: E402
from reduced_models import HANDLES, tool_fragment  # noqa: E402
from timing import load_snapshot, timed  # noqa: E402

EVIDENCE = LANE / "evidence"
GUESSES = ("minao", "atom", "huckel", "1e")
HARTREE_KCAL = 627.50947
EV_KCAL = 23.060548


def _pyscf_molecule(atoms, settings: QuantumSettings):
    from pyscf import gto

    molecule = gto.Mole()
    molecule.atom = [(symbol, tuple(position)) for symbol, position in zip(atoms.get_chemical_symbols(), atoms.positions)]
    molecule.basis = settings.basis
    molecule.charge = settings.charge
    molecule.spin = settings.spin
    molecule.unit = "Angstrom"
    molecule.verbose = 0
    molecule.max_memory = settings.memory_mb
    molecule.build()
    return molecule


def guess_scan(atoms, settings: QuantumSettings) -> dict:
    """SCF-only scan over starting guesses. Returns every attempt, converged or not.

    Agreement among the guesses that converged says nothing about the ones
    that did not, so the missing guesses are reported explicitly.
    """
    from pyscf import dft

    molecule = _pyscf_molecule(atoms, settings)
    attempts = []
    for guess in GUESSES:
        record: dict = {"initial_guess": guess}
        try:
            method = dft.UKS(molecule) if settings.spin else dft.RKS(molecule)
            method.xc = settings.xc
            method.grids.level = settings.grid_level
            method.conv_tol = settings.conv_tol
            method.max_cycle = settings.max_cycle
            if settings.density_fit:
                method = method.density_fit()
            method.init_guess = guess
            with timed(record, "timing"):
                energy = float(method.kernel())
            record.update(
                converged=bool(method.converged),
                scf_energy_hartree=energy,
                s2=float(method.spin_square()[0]) if settings.spin else 0.0,
            )
        except Exception as exc:  # noqa: BLE001 - a failed guess is a result
            record.update(converged=False, error=f"{type(exc).__name__}: {exc}")
        attempts.append(record)
        print(f"      guess {guess:7s} "
              + (f"E = {record['scf_energy_hartree']:.8f} Ha  S^2 = {record.get('s2', 0.0):.4f}"
                 if record.get("converged") else f"not converged ({record.get('error', 'no convergence')[:60]})"))

    converged = [a for a in attempts if a.get("converged")]
    if not converged:
        raise RuntimeError("No SCF starting guess converged for this species.")
    best = min(converged, key=lambda a: a["scf_energy_hartree"])
    spread = max(a["scf_energy_hartree"] for a in converged) - best["scf_energy_hartree"]
    return {
        "attempts": attempts,
        "guesses_requested": list(GUESSES),
        "guesses_converged": [a["initial_guess"] for a in converged],
        "guesses_failed": [a["initial_guess"] for a in attempts if not a.get("converged")],
        "scan_complete": len(converged) == len(GUESSES),
        "selected_initial_guess": best["initial_guess"],
        "selected_scf_energy_hartree": best["scf_energy_hartree"],
        "solution_spread_kcal_per_mol": spread * HARTREE_KCAL,
        "interpretation": (
            "Lowest converged SCF solution among the guesses that converged. "
            "This does not prove it is the ground state, nor that no other "
            "self-consistent solution exists. Agreement applies only to the "
            "converged subset."
        ),
    }


def species_energy(atoms, settings: QuantumSettings, label: str) -> dict:
    """Guess scan, then one production energy+gradient with the selected guess."""
    print(f"  {label}: {atoms.get_chemical_formula()}, {len(atoms)} atoms, spin={settings.spin}")
    record: dict = {
        "label": label,
        "formula": atoms.get_chemical_formula(),
        "n_atoms": len(atoms),
        "charge": settings.charge,
        "spin": settings.spin,
        "geometry_status": atoms.info.get("geometry_status", "unrelaxed_starting_guess"),
    }
    if len(atoms) > 1:
        record["guess_scan"] = guess_scan(atoms, settings)
        chosen = record["guess_scan"]["selected_initial_guess"]
    else:
        # A single atom has one solution; a scan would be noise.
        chosen = settings.scf_initial_guess
        record["guess_scan"] = {"performed": False, "reason": "single atom"}

    production = replace(settings, scf_initial_guess=chosen)
    atoms.calc = PySCFCalculator(production, event_log=EVIDENCE / f"{label}.jsonl")
    with timed(record, "production_timing"):
        energy_ev = float(atoms.get_potential_energy())
        forces = atoms.get_forces()
    record.update(
        energy_ev=energy_ev,
        energy_hartree=energy_ev / 27.211386024367243,
        max_force_ev_per_angstrom=float(np.linalg.norm(forces, axis=1).max()),
        production_initial_guess=chosen,
        quantum_diagnostics=atoms.calc.diagnostics,
        settings=production.to_dict(),
    )
    timing = record["production_timing"]
    print(f"      E = {energy_ev:.6f} eV   cpu {timing['cpu_seconds']:.1f} s"
          f"   wall {timing['wall_seconds']:.1f} s   contention x{timing['contention_factor_wall_over_cpu']:.1f}")
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handles", nargs="+", default=["hydrogen", "methyl"], choices=list(HANDLES))
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--xc", default="pbe0")
    parser.add_argument("--no-density-fit", action="store_true")
    parser.add_argument("--out", default=None, help="Evidence file name (default handle-fidelity-<handles>.json)")
    args = parser.parse_args()

    EVIDENCE.mkdir(parents=True, exist_ok=True)
    name = args.out or f"handle-fidelity-{'-'.join(args.handles)}.json"
    out = EVIDENCE / name

    base = QuantumSettings(
        xc=args.xc, basis=args.basis, dispersion="d3bj",
        density_fit=not args.no_density_fit, threads=1,
    )

    report: dict = {
        "schema_version": 1,
        "question": "Does replacing the adamantane tool handle change the tip's electronic appetite for a hydrogen?",
        "quantity": "tip H-affinity A(R) = E(R-CC-H) - E(R-CC*) - E(H); only differences across R are interpreted",
        "method": f"{args.xc.upper()}-D3(BJ)/{args.basis}"
                  + (", density fitting" if not args.no_density_fit else ", no density fitting"),
        "geometry_status": "rigid, taken from the 53-atom candidate pose; no relaxation",
        "measures": "electronic substituent effect only",
        "does_not_measure": (
            "mount stiffness, mass, sterics against the substrate, or anchor placement; "
            "a reduced handle is a different mechanical boundary condition, not a scaled one"
        ),
        "host_at_start": load_snapshot(),
        "species": {},
        "handles": {},
        "status": "running",
    }
    if out.exists():
        raise SystemExit(f"{out} exists; evidence is never overwritten. Choose --out.")
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    def save() -> None:
        out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    from ase import Atoms

    print(f"method: {report['method']}   host load {report['host_at_start']['load_average_1min']:.1f} "
          f"on {report['host_at_start']['logical_cores']} cores")

    print("\nhydrogen atom (reference; cancels in differences across handles)")
    hydrogen = Atoms("H", positions=[[0.0, 0.0, 0.0]])
    hydrogen.set_initial_magnetic_moments([1.0])
    report["species"]["H_atom"] = species_energy(hydrogen, replace(base, spin=1), "H_atom")
    save()

    for handle in args.handles:
        print(f"\nhandle: {handle}")
        radical, _ = tool_fragment(handle, hydrogenated=False)
        hydrogenated, _ = tool_fragment(handle, hydrogenated=True)
        radical_record = species_energy(radical, replace(base, spin=1), f"{handle}_radical")
        report["species"][f"{handle}_radical"] = radical_record
        save()
        hydrogenated_record = species_energy(hydrogenated, replace(base, spin=0), f"{handle}_hydrogenated")
        report["species"][f"{handle}_hydrogenated"] = hydrogenated_record
        save()

        affinity_ev = (
            hydrogenated_record["energy_ev"]
            - radical_record["energy_ev"]
            - report["species"]["H_atom"]["energy_ev"]
        )
        report["handles"][handle] = {
            "h_affinity_ev": affinity_ev,
            "h_affinity_kcal_per_mol": affinity_ev * EV_KCAL,
            "radical_atoms": radical_record["n_atoms"],
            "hydrogenated_atoms": hydrogenated_record["n_atoms"],
            "cpu_seconds_total": sum(
                record[key]["cpu_seconds"]
                for record in (radical_record, hydrogenated_record)
                for key in ("production_timing",)
            ),
        }
        print(f"    A({handle}) = {report['handles'][handle]['h_affinity_kcal_per_mol']:.2f} kcal/mol")
        save()

    if len(report["handles"]) > 1:
        values = {h: v["h_affinity_kcal_per_mol"] for h, v in report["handles"].items()}
        reference = values.get("adamantyl")
        report["comparison"] = {
            "h_affinity_kcal_per_mol": values,
            "spread_kcal_per_mol": max(values.values()) - min(values.values()),
            "shift_versus_adamantyl_kcal_per_mol": (
                {h: v - reference for h, v in values.items()} if reference is not None else None
            ),
            "interpretation": (
                "A spread small against the ~7.4 kcal/mol tertiary-versus-primary "
                "abstraction difference measured in this repository is evidence that "
                "the handle's electronic role is minor at fixed geometry. It is not "
                "evidence that the reduced model is mechanically or sterically faithful, "
                "and these are rigid-geometry energies, not barriers."
            ),
        }
    report["status"] = "completed"
    report["host_at_end"] = load_snapshot()
    save()
    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
