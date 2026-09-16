"""Measure the DFT method error on a SECONDARY-versus-TERTIARY site difference.

A1's adamantane site preference is a difference of two DFT reaction energies,
bridgehead (tertiary) against methylene (secondary). The usual defence of such a
difference is that method error cancels between similar species. My measured
errors say the cancellation is only partial when substitution differs:
PBE0-D3/def2-SVP is 2.01 kcal/mol too exothermic for a primary C-H abstraction
and 6.15 for a tertiary one, a 4.1 kcal/mol swing.

But that swing is primary-versus-tertiary, and A1's case is
secondary-versus-tertiary, a smaller substitution gap. Extrapolating is guesswork
in exactly the direction that decides whether A1's number is chemistry or method
error. A1 proposed measuring it directly, which is both cheaper and sharper than
the CCSD(T)-on-adamantyl calculation I had suggested and which A1 correctly showed
is unaffordable here (215 basis functions, roughly 22x the cost of a transition
structure already flagged as marginal, and past this machine's memory).

The acyclic analogue of A1's exact comparison:

    D(propane secondary) - D(isobutane tertiary)
      = E(isopropyl) - E(propane) - E(tert-butyl) + E(isobutane)

The hydrogen atom cancels, so only four species are needed and no abstracting
species enters at all. Evaluating that quantity at PBE0-D3/def2-SVP and at
CCSD(T)/cc-pVDZ, **at identical geometries**, makes the residual a direct
measurement of the DFT error on the quantity A1's result depends on.

Interpretation agreed with A1 in advance, so the reading is not chosen after
seeing the number: a residual near zero means A1's adamantane difference stands
as chemistry; a residual of 1-2 kcal/mol means it is not separable from method
error and the honest output is a bound and a direction, not a value.

Geometry provenance is mixed and that is deliberate. Isobutane and tert-butyl
use the published UCCSD(T)/cc-pVDZ structures already verified in this lane;
propane and isopropyl have no published counterpart and are relaxed here at
PBE0-D3/def2-SVP. This does not affect the measurement, because each species is
evaluated by *both* methods at one geometry, so the residual remains a pure
electronic-method comparison. It does mean the absolute bond dissociation
energies are not all of equal geometric quality, and they are not the point.

Run from the repository root:

    python data/validation/si-energy-reproduction/secondary_tertiary_calibration.py
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from ase.build import molecule
from ase.io import read, write
from ase.optimize import BFGS

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from nanodesign.quantum import PySCFCalculator, QuantumSettings  # noqa: E402
from pyscf import cc, gto, scf  # noqa: E402

KCAL_PER_EV = 1.0 / 0.0433641153087705
HARTREE_PER_KCAL = 1.0 / 627.50947
REFERENCE = ROOT / "data" / "reference"
OUT = Path(__file__).resolve().parent
GUESSES = ("minao", "atom", "huckel", "1e")


def build_propane_and_isopropyl():
    """Propane from the G2 set; isopropyl by removing one secondary hydrogen."""
    propane = molecule("C3H8")
    symbols = propane.get_chemical_symbols()
    distance = propane.get_all_distances()
    carbons = [i for i, s in enumerate(symbols) if s == "C"]
    central = [c for c in carbons if sum(1 for o in carbons if o != c and distance[c, o] < 1.8) == 2]
    if len(central) != 1:
        raise RuntimeError("Could not identify a unique secondary carbon in propane")
    secondary_h = [i for i, s in enumerate(symbols) if s == "H" and distance[central[0], i] < 1.3]
    if len(secondary_h) != 2:
        raise RuntimeError(f"Expected two secondary hydrogens, found {len(secondary_h)}")
    isopropyl = propane.copy()
    del isopropyl[secondary_h[0]]
    return propane, isopropyl, central[0]


def relax(atoms, spin, label, fmax=0.02, steps=200):
    settings = QuantumSettings(charge=0, spin=spin, xc="pbe0", basis="def2-svp",
                               dispersion="d3bj", grid_level=4, threads=1,
                               memory_mb=4000, density_fit=True)
    atoms = atoms.copy()
    atoms.calc = PySCFCalculator(settings)
    started, load0 = time.monotonic(), os.getloadavg()
    optimiser = BFGS(atoms, logfile=str(OUT / f"relax-{label}.log"))
    converged = bool(optimiser.run(fmax=fmax, steps=steps))
    return atoms, {
        "label": label, "geometry_converged": converged,
        "final_fmax_ev_per_angstrom": float(np.linalg.norm(atoms.get_forces(), axis=1).max()),
        "relax_seconds": time.monotonic() - started,
        "loadavg_at_start": list(load0),
        "relax_method": "pbe0-d3bj/def2-svp with density fitting",
    }


def dft_energy_ev(atoms, spin):
    settings = QuantumSettings(charge=0, spin=spin, xc="pbe0", basis="def2-svp",
                               dispersion="d3bj", grid_level=4, threads=1, memory_mb=4000)
    atoms = atoms.copy()
    atoms.calc = PySCFCalculator(settings)
    return float(atoms.get_potential_energy()), atoms.calc.diagnostics


def coupled_cluster_hartree(atoms, spin):
    """CCSD(T)/cc-pVDZ at the lowest SCF solution found over a guess scan."""
    best, attempts = None, []
    for guess in GUESSES:
        mol = gto.M(atom=list(zip(atoms.get_chemical_symbols(), atoms.positions.tolist())),
                    unit="Angstrom", basis="cc-pvdz", charge=0, spin=spin,
                    symmetry=False, verbose=0, max_memory=8000)
        mean_field = scf.UHF(mol) if spin else scf.RHF(mol)
        mean_field.conv_tol, mean_field.max_cycle, mean_field.init_guess = 1e-10, 300, guess
        energy = float(mean_field.kernel())
        if not mean_field.converged:
            attempts.append({"guess": guess, "converged": False})
            continue
        attempts.append({"guess": guess, "converged": True, "hf_energy_hartree": energy,
                         "s2": float(mean_field.spin_square()[0])})
        if best is None or energy < best["hf"] - 1e-10:
            best = {"hf": energy, "guess": guess, "mf": mean_field,
                    "s2": float(mean_field.spin_square()[0]), "nao": int(mol.nao_nr())}
        if spin == 0:
            break  # closed shell: one solution, no scan needed
    solver = cc.UCCSD(best["mf"], frozen=0) if spin else cc.RCCSD(best["mf"], frozen=0)
    solver.conv_tol, solver.max_cycle = 1e-9, 300
    solver.kernel()
    if not solver.converged:
        raise RuntimeError("CCSD did not converge")
    converged = [a["hf_energy_hartree"] for a in attempts if a.get("converged")]
    minao = [a for a in attempts if a["guess"] == "minao" and a.get("converged")]
    return float(solver.e_tot + solver.ccsd_t()), {
        "selected_guess": best["guess"], "hf_s2": best["s2"], "basis_functions": best["nao"],
        "scf_attempts": attempts,
        "minao_above_lowest_kcal_per_mol": ((minao[0]["hf_energy_hartree"] - best["hf"]) / HARTREE_PER_KCAL) if minao else None,
    }


def main() -> int:
    propane, isopropyl, _ = build_propane_and_isopropyl()
    species = [
        ("propane", propane, 0, "relax", "G2 set, relaxed here"),
        ("isopropyl", isopropyl, 1, "relax", "propane minus one secondary H, relaxed here"),
        ("isobutane", read(REFERENCE / "isobutane.xyz"), 0, "published", "published UCCSD(T)/cc-pVDZ"),
        ("tert_butyl", read(REFERENCE / "tert_butyl_radical.xyz"), 1, "published", "published UCCSD(T)/cc-pVDZ"),
    ]

    result = {
        "schema_version": 1,
        "question": "How large is the PBE0-D3/def2-SVP error on a SECONDARY-vs-TERTIARY C-H site difference?",
        "why": ("A1's adamantane site preference is a difference of two DFT reaction energies, "
                "bridgehead against methylene. Measured DFT error swings 4.1 kcal/mol between "
                "primary and tertiary sites, so cancellation across differing substitution is "
                "only partial. This measures the residual on A1's actual comparison."),
        "quantity": "D(propane secondary) - D(isobutane tertiary) = E(isopropyl) - E(propane) - E(tert_butyl) + E(isobutane)",
        "hydrogen_cancels": True,
        "agreed_interpretation": (
            "Residual near zero: A1's adamantane difference stands as chemistry. Residual 1-2 "
            "kcal/mol: not separable from method error, report a bound and a direction, not a value. "
            "Agreed with A1 before the number was known."),
        "geometry_note": (
            "Mixed provenance by design. Each species is evaluated by BOTH methods at ONE geometry, "
            "so the residual is a pure electronic-method comparison. Absolute BDEs are not the point "
            "and are not all of equal geometric quality."),
        "proposed_by": "A1 (andresarriaga-f2); executed in the forensics lane at its request",
        "species": [], "status": "running",
    }
    path = OUT / "secondary-tertiary-calibration.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

    geometries, records = {}, {}
    for name, atoms, spin, source, provenance in species:
        entry = {"species": name, "spin_2s": spin, "geometry_source": source,
                 "geometry_provenance": provenance}
        if source == "relax":
            atoms, info = relax(atoms, spin, name)
            entry.update(info)
            if not info["geometry_converged"]:
                raise RuntimeError(f"{name} did not relax")
        geometries[name] = atoms
        write(OUT / f"geom-{name}.xyz", atoms)

        started = time.monotonic()
        ev, diagnostics = dft_energy_ev(atoms, spin)
        entry.update(dft_energy_ev=ev, dft_seconds=time.monotonic() - started,
                     dft_s2=diagnostics.get("s2"), dft_basis_functions=diagnostics.get("basis_functions"))

        started = time.monotonic()
        hartree, cc_info = coupled_cluster_hartree(atoms, spin)
        entry.update(ccsd_t_energy_hartree=hartree, ccsd_t_seconds=time.monotonic() - started, **{
            "ccsd_t_" + k: v for k, v in cc_info.items()})
        records[name] = entry
        result["species"].append(entry)
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        print(f"{name:<12} DFT {ev:14.6f} eV   CCSD(T) {hartree:16.9f} Ha   "
              f"(guess {cc_info['selected_guess']}, {cc_info['basis_functions']} bf)", flush=True)

    dft = {k: records[k]["dft_energy_ev"] for k in records}
    ccc = {k: records[k]["ccsd_t_energy_hartree"] for k in records}
    dft_diff = ((dft["isopropyl"] - dft["propane"]) - (dft["tert_butyl"] - dft["isobutane"])) * KCAL_PER_EV
    cc_diff = ((ccc["isopropyl"] - ccc["propane"]) - (ccc["tert_butyl"] - ccc["isobutane"])) / HARTREE_PER_KCAL
    residual = dft_diff - cc_diff

    result.update(
        status="completed",
        secondary_minus_tertiary_kcal_per_mol={"pbe0_d3bj_def2svp": dft_diff, "ccsd_t_ccpvdz": cc_diff},
        dft_method_error_on_site_difference_kcal_per_mol=residual,
        comparison_to_primary_vs_tertiary_swing_kcal_per_mol=4.14,
        verdict=("A1 difference not separable from method error" if abs(residual) >= 1.0
                 else "DFT site difference survives this check"),
        caveats=[
            "One functional and one basis pair; def2-SVP against cc-pVDZ, so a small residual basis "
            "mismatch is folded in along with the method difference.",
            "Acyclic analogue. Adamantane's bridgehead radical is held pyramidal by the cage and does "
            "not behave like an unconstrained tertiary radical, so this bounds rather than reproduces "
            "A1's case.",
            "A small residual shows error cancellation for this pair; it does not validate PBE0-D3 for "
            "barriers, for the candidate, or for anything else.",
        ],
    )
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print()
    print(f"  secondary - tertiary, PBE0-D3/def2-SVP : {dft_diff:+.3f} kcal/mol")
    print(f"  secondary - tertiary, CCSD(T)/cc-pVDZ  : {cc_diff:+.3f} kcal/mol")
    print(f"  DFT METHOD ERROR ON THE DIFFERENCE     : {residual:+.3f} kcal/mol")
    print(f"  (primary-vs-tertiary swing for scale   : 4.14 kcal/mol)")
    print(f"  verdict: {result['verdict']}")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
