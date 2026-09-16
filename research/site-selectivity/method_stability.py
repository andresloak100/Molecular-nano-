"""Is the site difference method-stable, or is it an artifact of one functional?

A peer lane measured, at this lane's own level of theory, that the PBE0-D3(BJ)/
def2-SVP error against CCSD(T)/cc-pVDZ is *not constant across site type*:

    methane   primary  C-H   DFT -26.80   CCSD(T) -24.79   error -2.01
    isobutane tertiary C-H   DFT -38.38   CCSD(T) -32.23   error -6.15

a 4.1 kcal/mol swing, growing with substitution.  The usual defence of a DFT
energy *difference* is that method error cancels between similar species, and
that measurement says the cancellation is only partial when the two sites differ
in substitution -- which is exactly the bridgehead-versus-methylene case.  So a
site preference of a couple of kcal/mol could be substantially method error.

This script cannot fix that; only a high-level reference can, and CCSD(T) on a
215-basis-function C10H15 radical is not affordable here.  What it can do cheaply
is establish whether the difference is *method-stable*.  A difference that moves
by kcal/mol between functionals or basis sets must not be quoted as a value at
all, whatever its central number.

It exploits the cancellation identity: the site difference is
E(adamantyl_methylene) - E(adamantyl_bridgehead) and nothing else, so each
additional level of theory costs exactly **two** single-point SCF calculations,
with no relaxation and no other species.  Geometries stay fixed at the
PBE0-D3/def2-SVP relaxed structures, because the effect being probed is an
electronic-method effect rather than a geometry effect; that also means these are
*not* independent adiabatic differences, and they are labelled accordingly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from ase.io import read as ase_read
from scipy.constants import Avogadro, calorie, electron_volt

from nanodesign.quantum import PySCFCalculator, QuantumSettings

EV_TO_KCAL_PER_MOL = electron_volt * Avogadro / (1000.0 * calorie)

# (label, xc, basis, dispersion).  M06-2X has no D3(BJ) damping entry, so it
# needs d3zero or no dispersion; that is a documented finding in this repository.
DEFAULT_LEVELS = (
    ("b3lyp-d3bj-svp", "b3lyp", "def2-svp", "d3bj"),
    ("pbe0-d3bj-tzvp", "pbe0", "def2-tzvp", "d3bj"),
)

RADICALS = ("adamantyl_bridgehead", "adamantyl_methylene")


def single_point(geometry: Path, spin: int, *, xc: str, basis: str, dispersion: str | None,
                 density_fit: bool, log: Path) -> dict:
    atoms = ase_read(geometry)
    settings = QuantumSettings(
        charge=0, spin=spin, xc=xc, basis=basis, dispersion=dispersion,
        grid_level=3, conv_tol=1e-9, max_cycle=200, threads=1,
        memory_mb=3000, density_fit=density_fit,
    )
    calculator = PySCFCalculator(settings, event_log=log)
    atoms.calc = calculator
    started = time.monotonic()
    record: dict = {
        "geometry": str(geometry),
        "settings": settings.to_dict(),
        "geometry_source": "relaxed at PBE0-D3(BJ)/def2-SVP; not re-relaxed at this level",
        "scf_initial_guess": "pyscf default (minao); not configurable through PySCFCalculator",
    }
    try:
        record["energy_ev"] = float(atoms.get_potential_energy())
        record["completed"] = True
    except Exception as error:
        record.update(completed=False, error_type=type(error).__name__, error=str(error))
    record["wall_seconds"] = time.monotonic() - started
    diagnostics = dict(calculator.diagnostics)
    record["diagnostics"] = {
        key: diagnostics.get(key)
        for key in ("scf_converged", "basis_functions", "reference", "s2", "expected_s2",
                    "spin_contamination", "effective_pyscf_threads", "threads_honored",
                    "elapsed_seconds")
    }
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometries", required=True,
                        help="Stage 1/2 evidence directory holding <radical>/relaxed.xyz")
    parser.add_argument("--output", required=True)
    parser.add_argument("--density-fit", action="store_true")
    parser.add_argument("--levels", nargs="*", default=None,
                        help="Subset of level labels to run; default runs all")
    arguments = parser.parse_args()

    geometries = Path(arguments.geometries).resolve()
    output = Path(arguments.output).resolve()
    if output.name == "runs":
        raise SystemExit("Refusing to write to a directory named 'runs'")
    output.mkdir(parents=True, exist_ok=False)

    levels = [
        level for level in DEFAULT_LEVELS
        if arguments.levels is None or level[0] in arguments.levels
    ]
    if not levels:
        raise SystemExit(f"No levels selected; known: {[item[0] for item in DEFAULT_LEVELS]}")

    results = []
    for label, xc, basis, dispersion in levels:
        energies: dict[str, dict] = {}
        for radical in RADICALS:
            relaxed = geometries / radical / "relaxed.xyz"
            if not relaxed.exists():
                print(f"[skip] {label}/{radical}: no relaxed geometry at {relaxed}", flush=True)
                energies[radical] = {"completed": False, "error": "missing relaxed geometry"}
                continue
            print(f"[run ] {label} / {radical}", flush=True)
            record = single_point(
                relaxed, spin=1, xc=xc, basis=basis, dispersion=dispersion,
                density_fit=bool(arguments.density_fit),
                log=output / f"{label}-{radical}.jsonl",
            )
            energies[radical] = record
            if record["completed"]:
                print(
                    f"[done] {label} / {radical}: {record['energy_ev']:.6f} eV, "
                    f"nbf {record['diagnostics'].get('basis_functions')}, "
                    f"S2 {record['diagnostics'].get('s2')}, "
                    f"{record['wall_seconds']:.0f} s", flush=True
                )
            else:
                print(f"[FAIL] {label} / {radical}: {record.get('error')}", flush=True)

        entry = {
            "level": label,
            "xc": xc, "basis": basis, "dispersion": dispersion,
            "species": energies,
        }
        if all(energies[radical].get("completed") for radical in RADICALS):
            entry["site_difference_kcal_per_mol"] = (
                energies["adamantyl_methylene"]["energy_ev"]
                - energies["adamantyl_bridgehead"]["energy_ev"]
            ) * EV_TO_KCAL_PER_MOL
            entry["sign_convention"] = (
                "E(methylene) - E(bridgehead); negative means the methylene site is the "
                "thermodynamically easier abstraction"
            )
        results.append(entry)

    differences = [
        entry["site_difference_kcal_per_mol"] for entry in results
        if "site_difference_kcal_per_mol" in entry
    ]
    summary = {
        "purpose": "Method stability of the bridgehead-versus-methylene site difference",
        "identity_used": (
            "The site difference equals E(adamantyl_methylene) - E(adamantyl_bridgehead); "
            "every other term in both thermochemical cycles cancels exactly, so each level "
            "of theory costs two single-point calculations."
        ),
        "levels": results,
        "spread_kcal_per_mol": (max(differences) - min(differences)) if len(differences) > 1 else None,
        "peer_measured_site_dependent_dft_error": {
            "methane_primary": {"dft": -26.80, "ccsd_t": -24.79, "error": -2.01},
            "isobutane_tertiary": {"dft": -38.38, "ccsd_t": -32.23, "error": -6.15},
            "swing_kcal_per_mol": 4.14,
            "note": (
                "Measured by a peer lane at PBE0-D3(BJ)/def2-SVP against CCSD(T)/cc-pVDZ. "
                "The comparison conflates method and basis and is two points on one "
                "functional; the trend with substitution is the trustworthy part. It "
                "concerns primary-versus-tertiary, while this lane compares "
                "secondary-versus-tertiary, a smaller substitution gap."
            ),
        },
        "limitations": [
            "Method stability is not accuracy. Agreement across functionals would show the difference is robust to this choice, not that it is correct.",
            "Geometries are fixed at the PBE0-D3(BJ)/def2-SVP minima, so these are not independent adiabatic differences at each level.",
            "No high-level reference is computed. CCSD(T)/cc-pVDZ on a 215-basis-function C10H15 radical is not affordable on this host.",
            "Single-point energies inherit the open-shell initial-guess hazard; run the guess scan alongside.",
        ],
    }
    (output / "method_stability.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    for entry in results:
        if "site_difference_kcal_per_mol" in entry:
            print(f"{entry['level']:<20} {entry['site_difference_kcal_per_mol']:+8.3f} kcal/mol", flush=True)
    if summary["spread_kcal_per_mol"] is not None:
        print(f"spread across levels: {summary['spread_kcal_per_mol']:.3f} kcal/mol", flush=True)
    print("wrote", output / "method_stability.json", flush=True)


if __name__ == "__main__":
    main()
