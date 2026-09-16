"""Extend the supporting-information reproduction to the isobutane reaction.

The methane check (``reproduce.py``) traced the published energy inconsistency
to a single erroneous CH4 entry. This script applies the same test to the
second reaction in the source, ``C2H + iso-C4H10 -> C2H2 + t-C4H9``.

That reaction matters more than methane for this repository. Its abstraction
site is a *tertiary* C-H, which is the closest small-molecule proxy available
for the adamantane bridgehead C-H that the 53-atom candidate targets. A
calibration that only covers methane's primary C-H does not cover the bond the
tool is actually designed to break.

Species are computed cheapest first and results are written after each one, so
an interrupted run still leaves usable evidence. The transition structure is
by far the most expensive: 139 basis functions against the methane transition
structure's 67, and CCSD(T) triples scale steeply, so it may not complete in a
short session. That is recorded rather than worked around.

The SCF initial-guess scan is carried over from ``reproduce.py`` and is
load-bearing, not defensive padding; see that file for the ethynyl radical
case where the default guess silently gives a wrong answer.

Run from the repository root:

    python data/validation/si-energy-reproduction/reproduce_isobutane.py
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

from ase.io import read

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pyscf import cc, gto, scf  # noqa: E402

HARTREE_PER_KCAL_PER_MOL = 1.0 / 627.50947
REFERENCE = ROOT / "data" / "reference"
INITIAL_GUESSES = ("minao", "atom", "huckel", "1e")
AGREEMENT_HARTREE = 1e-7

# name, file, charge, spin (N_alpha - N_beta), SI label. Cheapest first.
SPECIES = [
    ("ethynyl_radical", "ethynyl_radical.xyz", 0, 1, "CCH"),
    ("acetylene", "acetylene.xyz", 0, 0, "HCCH"),
    ("tert_butyl_radical", "tert_butyl_radical.xyz", 0, 1, "T-butyl"),
    ("isobutane", "isobutane.xyz", 0, 0, "Isobutane"),
    ("isobutane_ethynyl_ts", "isobutane_ethynyl_ts.xyz", 0, 1, "C(CH3)3--H-CCH"),
]


def lowest_scf(filename: str, charge: int, spin: int, basis: str = "cc-pvdz") -> dict:
    """Return the lowest converged HF solution found over the guess scan."""
    atoms = read(REFERENCE / filename)
    attempts, best = [], None
    for guess in INITIAL_GUESSES:
        molecule = gto.M(
            atom=list(zip(atoms.get_chemical_symbols(), atoms.positions.tolist())),
            unit="Angstrom", basis=basis, charge=charge, spin=spin,
            symmetry=False, verbose=0, max_memory=8000,
        )
        mean_field = scf.UHF(molecule) if spin else scf.RHF(molecule)
        mean_field.conv_tol = 1e-10
        mean_field.max_cycle = 300
        mean_field.init_guess = guess
        energy = float(mean_field.kernel())
        if not mean_field.converged:
            attempts.append({"initial_guess": guess, "converged": False})
            continue
        attempts.append({"initial_guess": guess, "converged": True,
                         "hartree_fock_energy_hartree": energy,
                         "hartree_fock_s2": float(mean_field.spin_square()[0])})
        if best is None or energy < best["hartree_fock_energy_hartree"] - 1e-10:
            best = {"initial_guess": guess, "hartree_fock_energy_hartree": energy,
                    "hartree_fock_s2": float(mean_field.spin_square()[0]),
                    "mean_field": mean_field, "basis_functions": int(molecule.nao_nr())}
    if best is None:
        raise RuntimeError(f"No converged SCF solution for {filename}")
    converged = [a["hartree_fock_energy_hartree"] for a in attempts if a.get("converged")]
    best["scf_attempts"] = attempts
    best["scf_solution_spread_kcal_per_mol"] = (max(converged) - min(converged)) / HARTREE_PER_KCAL_PER_MOL
    return best


def main() -> int:
    provenance = json.loads((REFERENCE / "provenance.json").read_text())
    published = {e["source_label"]: e.get("source_energy_hartree") for e in provenance["structures"]}
    output = Path(__file__).resolve().parent / "si-energy-reproduction-isobutane.json"

    records: list[dict] = []
    result = {
        "schema_version": 1,
        "reaction": "C2H + iso-C4H10 -> C2H2 + t-C4H9",
        "why_this_reaction": (
            "The abstraction site is a tertiary C-H, the closest available small-molecule proxy "
            "for the adamantane bridgehead C-H targeted by the 53-atom candidate. The methane "
            "calibration covers a primary C-H and does not cover the bond the tool must break."
        ),
        "level_of_theory": "all-electron UHF/UCCSD(T) and RHF/RCCSD(T), cc-pVDZ, no frozen orbitals",
        "agreement_threshold_hartree": AGREEMENT_HARTREE,
        "initial_guesses_scanned": list(INITIAL_GUESSES),
        "status": "running",
        "species": records,
    }
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

    for name, filename, charge, spin, label in SPECIES:
        started = time.monotonic()
        solution = lowest_scf(filename, charge, spin)
        mean_field = solution["mean_field"]
        solver = cc.UCCSD(mean_field, frozen=0) if spin else cc.RCCSD(mean_field, frozen=0)
        solver.conv_tol = 1e-9
        solver.max_cycle = 300
        solver.kernel()
        if not solver.converged:
            raise RuntimeError(f"CCSD did not converge for {name}")
        computed = float(solver.e_tot + solver.ccsd_t())
        difference = computed - published[label]
        records.append({
            "species": name,
            "source_label": label,
            "variant": "UCCSD(T)" if spin else "RCCSD(T)",
            "basis_functions": solution["basis_functions"],
            "published_energy_hartree": published[label],
            "recomputed_energy_hartree": computed,
            "difference_kcal_per_mol": difference / HARTREE_PER_KCAL_PER_MOL,
            "reproduces_published_value": abs(difference) < AGREEMENT_HARTREE,
            "selected_initial_guess": solution["initial_guess"],
            "hartree_fock_s2": solution["hartree_fock_s2"],
            "scf_solution_spread_kcal_per_mol": solution["scf_solution_spread_kcal_per_mol"],
            "scf_attempts": solution["scf_attempts"],
            "elapsed_seconds": time.monotonic() - started,
        })
        output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        print(f"{name:<22} {solution['basis_functions']:>4} bf  published {published[label]:>17.9f}  "
              f"recomputed {computed:>17.9f}  {difference / HARTREE_PER_KCAL_PER_MOL:>+8.2f} kcal/mol  "
              f"(guess {solution['selected_initial_guess'] if 'selected_initial_guess' in solution else solution['initial_guess']}, "
              f"spread {solution['scf_solution_spread_kcal_per_mol']:.2f}, "
              f"{time.monotonic() - started:.0f}s)", flush=True)

    ours = {r["species"]: r["recomputed_energy_hartree"] for r in records}
    theirs = {r["species"]: r["published_energy_hartree"] for r in records}

    def barrier(source):
        return (source["isobutane_ethynyl_ts"] - source["isobutane"]
                - source["ethynyl_radical"]) / HARTREE_PER_KCAL_PER_MOL

    def reaction(source):
        return (source["acetylene"] + source["tert_butyl_radical"] - source["isobutane"]
                - source["ethynyl_radical"]) / HARTREE_PER_KCAL_PER_MOL

    not_reproduced = [r["species"] for r in records if not r["reproduces_published_value"]]
    result.update(
        status="completed",
        reproduced=[r["species"] for r in records if r["reproduces_published_value"]],
        not_reproduced=not_reproduced,
        barrier_kcal_per_mol={
            "from_published_absolute_energies": barrier(theirs),
            "from_recomputed_absolute_energies": barrier(ours),
        },
        reaction_energy_kcal_per_mol={
            "from_published_absolute_energies": reaction(theirs),
            "from_recomputed_absolute_energies": reaction(ours),
        },
        interpretation=(
            "A tertiary C-H barrier below the primary C-H barrier is the expected ordering and is "
            "not by itself a validation. The supplied transition structure is a nominal geometry "
            "from the source, not a verified first-order saddle, and no DFT saddle has been located "
            "for either reaction."
        ),
        method_validated=False,
    )
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

    print()
    print(f"isobutane barrier, SI absolutes   : {barrier(theirs):+.3f} kcal/mol")
    print(f"isobutane barrier, recomputed     : {barrier(ours):+.3f} kcal/mol")
    print(f"isobutane reaction, recomputed    : {reaction(ours):+.3f} kcal/mol")
    print(f"entries not reproduced            : {not_reproduced or 'none'}")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
