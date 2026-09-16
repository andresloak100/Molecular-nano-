"""Verify that the open-shell DFT energies used here are initial-guess independent.

The shared coordination notes record a trap that already produced one wrong
published result in this project: open-shell Hartree-Fock in this system has
several converged, *stable* SCF solutions, and PySCF's default ``minao`` guess
landed 8.7 kcal/mol above the solution reached from ``atom`` or ``1e`` for the
ethynyl radical, with ``stability()`` calling both stable.

DFT was reported guess-independent at the geometries tested previously.  This
script does not assume that carries over to the new radical geometries generated
in this lane; it re-tests each open-shell species at its own relaxed geometry.

Only the Kohn-Sham part is scanned.  The D3(BJ) dispersion correction is a
function of the nuclear geometry alone, so it is identical across guesses and
cannot affect a guess comparison; it is therefore omitted here and the omission
does not weaken the test.

A converged SCF is not a correct one, and agreement across four guesses is
evidence against a low-lying alternative solution rather than proof that none
exists.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from ase.io import read
from scipy.constants import Avogadro, calorie, electron_volt

from species import OPEN_SHELL, build_species

HARTREE_TO_KCAL_PER_MOL = (
    electron_volt * Avogadro / (1000.0 * calorie) * 27.211386245988
)
GUESSES = ("minao", "atom", "huckel", "1e")


def scan_species(name: str, geometry_path: Path, spin: int, *, xc: str, basis: str,
                 grid_level: int, conv_tol: float, stability: bool) -> dict:
    from pyscf import dft, gto, lib

    atoms = read(geometry_path)
    molecule = gto.M(
        atom=list(zip(atoms.get_chemical_symbols(), atoms.positions.tolist())),
        unit="Angstrom", basis=basis, charge=0, spin=spin,
        symmetry=False, verbose=0, max_memory=3000,
    )
    solutions = []
    for guess in GUESSES:
        started = time.monotonic()
        record: dict = {"initial_guess": guess}
        try:
            mean_field = dft.UKS(molecule)
            mean_field.xc = xc
            mean_field.disp = False
            mean_field.conv_tol = conv_tol
            mean_field.max_cycle = 200
            mean_field.grids.level = grid_level
            mean_field.init_guess = guess
            energy = float(mean_field.kernel())
            s2, multiplicity = mean_field.spin_square()
            record.update(
                converged=bool(mean_field.converged),
                energy_hartree=energy,
                s2=float(s2),
                effective_multiplicity=float(multiplicity),
            )
            if stability and mean_field.converged:
                internal, _external = mean_field.stability()
                # Stability analysis returned both solutions as stable in the
                # documented failure, so this is recorded and not relied upon.
                record["stability_internal_rotation_changed"] = bool(
                    not lib.finger(internal) == lib.finger(mean_field.mo_coeff)
                )
        except Exception as error:
            record.update(converged=False, error_type=type(error).__name__, error=str(error))
        record["wall_seconds"] = time.monotonic() - started
        solutions.append(record)
        print(
            f"  {name:<24} {guess:<7} "
            + (
                f"E = {record['energy_hartree']:.9f} Ha  S2 = {record['s2']:.4f}"
                if record.get("converged")
                else f"FAILED: {record.get('error_type')}"
            ),
            flush=True,
        )

    converged = [item for item in solutions if item.get("converged")]
    result = {
        "species": name,
        "geometry": str(geometry_path),
        "spin_n_alpha_minus_n_beta": spin,
        "method": f"{xc.upper()}/{basis}, grid level {grid_level}, no dispersion",
        "basis_functions": int(molecule.nao_nr()),
        "guesses_tested": list(GUESSES),
        "solutions": solutions,
        "converged_count": len(converged),
    }
    if converged:
        energies = [item["energy_hartree"] for item in converged]
        lowest = min(converged, key=lambda item: item["energy_hartree"])
        spread = (max(energies) - min(energies)) * HARTREE_TO_KCAL_PER_MOL
        result.update(
            lowest_energy_hartree=lowest["energy_hartree"],
            lowest_energy_guess=lowest["initial_guess"],
            spread_kcal_per_mol=spread,
            guess_independent=bool(spread < 1e-3),
            s2_values={item["initial_guess"]: item["s2"] for item in converged},
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometries", required=True,
                        help="Stage 1/2 evidence directory holding <species>/relaxed.xyz")
    parser.add_argument("--output", required=True)
    parser.add_argument("--xc", default="pbe0")
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--grid-level", type=int, default=3)
    parser.add_argument("--conv-tol", type=float, default=1e-9)
    parser.add_argument("--stability", action="store_true",
                        help="Also run PySCF stability analysis; recorded, not relied upon")
    arguments = parser.parse_args()

    geometries = Path(arguments.geometries).resolve()
    output = Path(arguments.output).resolve()
    if output.name == "runs":
        raise SystemExit("Refusing to write to a directory named 'runs'; .gitignore drops it")
    output.mkdir(parents=True, exist_ok=False)

    species = build_species()
    results = []
    for name in OPEN_SHELL:
        relaxed = geometries / name / "relaxed.xyz"
        if not relaxed.exists():
            # A single atom has no relaxed geometry written; use its input.
            relaxed = geometries / name / "input.xyz"
        if not relaxed.exists():
            print(f"  {name}: no geometry at {relaxed}; skipped", flush=True)
            results.append({"species": name, "skipped": "geometry not available"})
            continue
        results.append(
            scan_species(
                name, relaxed, species[name].spin,
                xc=arguments.xc, basis=arguments.basis,
                grid_level=arguments.grid_level, conv_tol=arguments.conv_tol,
                stability=arguments.stability,
            )
        )

    scanned = [item for item in results if "spread_kcal_per_mol" in item]
    summary = {
        "purpose": "Guess-independence test for the open-shell DFT energies used in stage 1 and 2",
        "guesses": list(GUESSES),
        "species": results,
        "max_spread_kcal_per_mol": max(
            (item["spread_kcal_per_mol"] for item in scanned), default=None
        ),
        "all_guess_independent": bool(scanned) and all(
            item.get("guess_independent") for item in scanned
        ),
        "limitations": [
            "Agreement across four initial guesses is evidence against a low-lying alternative SCF solution, not proof that none exists.",
            "PySCF stability analysis called both solutions stable in the documented Hartree-Fock failure, so it is recorded here and not treated as decisive.",
            "This tests the self-consistent-field solution only. It says nothing about functional accuracy, basis adequacy, or multireference character.",
        ],
    }
    (output / "guess_scan.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(
        "max spread across guesses:",
        f"{summary['max_spread_kcal_per_mol']:.2e} kcal/mol"
        if summary["max_spread_kcal_per_mol"] is not None else "n/a",
        "| all guess-independent:", summary["all_guess_independent"],
        flush=True,
    )
    print("wrote", output / "guess_scan.json", flush=True)


if __name__ == "__main__":
    main()
