"""Was the headline 53-atom calculation on the right electronic state?

The archived 53-atom records under ``data/validation/`` are the repository's
central evidence. They ran on PySCF's default ``minao`` starting guess, and
both record ``scf_initial_guess`` and ``initial_guess_scan_performed`` as null:
nobody scanned them.

That matters because the default guess has now been caught converging to a
higher, wrong solution for four species in this project, including the
candidate's own tool tip. The adamantyl-supported ethynyl radical, which is
26 of the 53 atoms, comes out 10.42 kcal/mol too high from ``minao``.

The reassuring argument, and why it is not enough. Both archived runs report
S^2 = 0.78462, near the correct tip solution's 0.7863 rather than the trap's
0.7521 - and since the added adamantane is closed-shell and contributes nothing
to S^2, that comparison is closer to like-for-like than a generic cross-system
one would be. So the records are probably fine. But this is using S^2 to select
a solution, and the single clearest lesson of this project is that S^2 selects
wrongly: in all four cases measured, the WRONG higher solution carried the
cleaner S^2.

So: probably fine on a physically motivated argument, unverified. This settles
it.

A NULL RESULT IS THE INTERESTING ONE. If ``minao`` agrees with the other three
guesses here, that is not a boring confirmation - it would mean a closed-shell
spectator changes which solution the default guess finds, since ``minao`` fails
for the bare 26-atom tip and would have succeeded for the 53-atom system
containing it. None of the four cases so far could show that. The alternative,
that ``minao`` sits ~10 kcal/mol high here too, would mean the archived
energies are wrong and the S^2 argument misled everyone.

Either way this reports ``minao``'s energy explicitly, agreement or not.

SCOPE. SCF only, no gradient: the question is which solution, not what force.
Writes to this lane's evidence directory and touches nothing under
``data/validation/``. Reads the archived geometry rather than rebuilding it, so
the comparison is against the same coordinates the archive used.

Run from the repository root:

    python research/candidate-feasibility/scan_candidate_scf.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

from timing import load_snapshot, timed  # noqa: E402

EVIDENCE = LANE / "evidence"
ARCHIVE = ROOT / "data" / "validation" / "h-abstraction-df-initial"
GUESSES = ("minao", "atom", "huckel", "1e")
HARTREE_KCAL = 627.50947


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="candidate-scf-scan.json")
    args = parser.parse_args()

    EVIDENCE.mkdir(parents=True, exist_ok=True)
    out = EVIDENCE / args.out
    if out.exists():
        raise SystemExit(f"{out} exists; evidence is never overwritten. Choose --out.")

    archived = json.loads((ARCHIVE / "result.json").read_text())
    diagnostics = archived["structure"]["quantum_diagnostics"]
    settings = diagnostics["settings"]
    atoms = read(ARCHIVE / "input-initial.extxyz")

    report: dict = {
        "schema_version": 1,
        "question": "Did the archived 53-atom calculation converge to the lowest SCF solution?",
        "geometry_source": str((ARCHIVE / "input-initial.extxyz").relative_to(ROOT)),
        "archived_record": str((ARCHIVE / "result.json").relative_to(ROOT)),
        "archived_energy_hartree": diagnostics["total_energy_hartree"],
        "archived_dft_energy_hartree": diagnostics["dft_energy_hartree"],
        "archived_s2": diagnostics["s2"],
        "archived_scf_initial_guess": settings.get("scf_initial_guess"),
        "archived_scan_performed": diagnostics.get("initial_guess_scan_performed"),
        "settings_used": settings,
        "scope": "SCF only, no gradient; dispersion excluded so energies compare to dft_energy_hartree",
        "touches_archive": False,
        "host_at_start": load_snapshot(),
        "attempts": [],
        "status": "running",
    }
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    from pyscf import dft, gto

    molecule = gto.Mole()
    molecule.atom = [
        (symbol, tuple(position))
        for symbol, position in zip(atoms.get_chemical_symbols(), atoms.positions)
    ]
    molecule.basis = settings["basis"]
    molecule.charge = settings["charge"]
    molecule.spin = settings["spin"]
    molecule.unit = "Angstrom"
    molecule.verbose = 0
    molecule.max_memory = settings["memory_mb"]
    molecule.build()

    print(f"53-atom candidate, {len(atoms)} atoms, {settings['xc']}/{settings['basis']}, "
          f"spin={settings['spin']}")
    print(f"archived DFT energy {diagnostics['dft_energy_hartree']:.8f} Ha, "
          f"S^2 {diagnostics['s2']:.5f}, guess recorded: {settings.get('scf_initial_guess')}\n")

    for guess in GUESSES:
        record: dict = {"initial_guess": guess}
        try:
            method = dft.UKS(molecule)
            method.xc = settings["xc"]
            method.grids.level = settings["grid_level"]
            method.conv_tol = settings["conv_tol"]
            method.max_cycle = settings["max_cycle"]
            if settings.get("density_fit"):
                method = method.density_fit()
            method.init_guess = guess
            with timed(record, "timing"):
                energy = float(method.kernel())
            record.update(
                converged=bool(method.converged),
                scf_energy_hartree=energy,
                s2=float(method.spin_square()[0]),
                difference_from_archive_kcal_per_mol=(
                    (energy - diagnostics["dft_energy_hartree"]) * HARTREE_KCAL
                ),
            )
            print(f"  {guess:7s} E = {energy:.8f} Ha  S^2 = {record['s2']:.5f}  "
                  f"vs archive {record['difference_from_archive_kcal_per_mol']:+.4f} kcal/mol  "
                  f"({record['timing']['user_seconds']:.0f} s user)")
        except Exception as exc:  # noqa: BLE001 - a failed guess is a result
            record.update(converged=False, error=f"{type(exc).__name__}: {exc}")
            print(f"  {guess:7s} FAILED: {record['error'][:70]}")
        report["attempts"].append(record)
        out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    converged = [a for a in report["attempts"] if a.get("converged")]
    if converged:
        best = min(converged, key=lambda a: a["scf_energy_hartree"])
        spread = max(a["scf_energy_hartree"] for a in converged) - best["scf_energy_hartree"]
        minao = next((a for a in converged if a["initial_guess"] == "minao"), None)
        report["summary"] = {
            "guesses_converged": [a["initial_guess"] for a in converged],
            "guesses_failed": [
                a["initial_guess"] for a in report["attempts"] if not a.get("converged")
            ],
            "scan_complete": len(converged) == len(GUESSES),
            "lowest_guess": best["initial_guess"],
            "lowest_energy_hartree": best["scf_energy_hartree"],
            "solution_spread_kcal_per_mol": spread * HARTREE_KCAL,
            "minao_energy_hartree": minao["scf_energy_hartree"] if minao else None,
            "minao_above_lowest_kcal_per_mol": (
                (minao["scf_energy_hartree"] - best["scf_energy_hartree"]) * HARTREE_KCAL
                if minao else None
            ),
            "archive_matches_lowest": (
                abs(diagnostics["dft_energy_hartree"] - best["scf_energy_hartree"]) * HARTREE_KCAL < 0.01
            ),
            "interpretation": (
                "Lowest converged solution among the guesses that converged. Does "
                "not prove it is the ground state or that no other solution exists."
            ),
        }
        summary = report["summary"]
        print(f"\nlowest: {summary['lowest_guess']}, spread {summary['solution_spread_kcal_per_mol']:.4f} kcal/mol")
        print(f"minao above lowest: {summary['minao_above_lowest_kcal_per_mol']:.4f} kcal/mol"
              if summary["minao_above_lowest_kcal_per_mol"] is not None else "minao did not converge")
        print(f"ARCHIVE ON LOWEST SOLUTION: {summary['archive_matches_lowest']}")

    report["status"] = "completed"
    report["host_at_end"] = load_snapshot()
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
