"""DFT barrier survey for the isobutane (tertiary C-H) abstraction reaction.

The existing functional survey covers only the methane reaction, whose site is a
primary C-H. This repository's candidate targets an adamantane bridgehead, which
is tertiary, so the question that matters is whether DFT's failure on the
methane barrier also appears at a tertiary site, or whether it is specific to
methane.

Evaluated at the published UCCSD(T)/cc-pVDZ geometries, exactly as the methane
survey was, so the two are directly comparable. This is a fixed-geometry
comparison, not a DFT barrier: the supplied structure is a coupled-cluster
stationary point and is not stationary on any DFT surface. Under PBE0 the
methane structure carries a 2.38 eV/Angstrom residual force, about eighty times
this repository's convergence threshold. Treat every number here as a single
point on a hillside.

Runs one species at a time in a single process. Records `os.getloadavg()` with
every timing, because this host has been carrying 200+ load and a wall-clock
figure taken under that describes the scheduler rather than the calculation.

Run from the repository root:

    python data/validation/si-energy-reproduction/dft_survey_isobutane.py
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import time

from ase.io import read

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from nanodesign.quantum import PySCFCalculator, QuantumSettings  # noqa: E402

KCAL_PER_EV = 1.0 / 0.0433641153087705
REFERENCE = ROOT / "data" / "reference"

# name, file, charge, spin (N_alpha - N_beta)
SPECIES = [
    ("ethynyl_radical", "ethynyl_radical.xyz", 0, 1),
    ("acetylene", "acetylene.xyz", 0, 0),
    ("isobutane", "isobutane.xyz", 0, 0),
    ("tert_butyl_radical", "tert_butyl_radical.xyz", 0, 1),
    ("isobutane_ethynyl_ts", "isobutane_ethynyl_ts.xyz", 0, 1),
]

METHODS = [
    ("pbe0", "def2-svp", "d3bj"),
    ("pbe0", "def2-tzvp", "d3bj"),
    ("b3lyp", "def2-tzvp", "d3bj"),
    ("m06-2x", "def2-tzvp", None),
]

# All-electron CCSD(T)/cc-pVDZ at these same geometries, lowest SCF solution,
# from reproduce.py and reproduce_isobutane.py. The isobutane barrier is filled
# in once its transition structure completes; the reaction energies are final.
COUPLED_CLUSTER = {
    "methane_barrier_kcal_per_mol": 2.3986,
    "methane_reaction_kcal_per_mol": -24.788,
    "isobutane_reaction_kcal_per_mol": -32.232,
    "isobutane_barrier_kcal_per_mol": None,
    "note": (
        "Isobutane barrier pending its 139-basis-function transition structure. "
        "The published SI absolute energies imply -0.627 kcal/mol for it, and all four "
        "non-transition isobutane species reproduce those energies to under 1e-7 Hartree."
    ),
}


def energy_ev(name: str, filename: str, charge: int, spin: int, xc: str, basis: str, dispersion):
    atoms = read(REFERENCE / filename)
    settings = QuantumSettings(charge=charge, spin=spin, xc=xc, basis=basis,
                               dispersion=dispersion, grid_level=4, threads=1, memory_mb=6000)
    atoms.calc = PySCFCalculator(settings)
    started, load_before = time.monotonic(), os.getloadavg()
    value = float(atoms.get_potential_energy())
    diagnostics = atoms.calc.diagnostics
    return value, {
        "species": name,
        "energy_ev": value,
        "basis_functions": diagnostics.get("basis_functions"),
        "s2": diagnostics.get("s2"),
        "expected_s2": diagnostics.get("expected_s2"),
        "elapsed_seconds": time.monotonic() - started,
        "effective_threads": diagnostics.get("effective_pyscf_threads"),
        "threads_honored": diagnostics.get("threads_honored"),
        "loadavg_before": list(load_before),
        "loadavg_after": list(os.getloadavg()),
    }


def main() -> int:
    output = Path(__file__).resolve().parent / "dft-survey-isobutane.json"
    result = {
        "schema_version": 1,
        "reaction": "C2H + iso-C4H10 -> C2H2 + t-C4H9",
        "site": "tertiary C-H, the closest small-molecule proxy for the adamantane bridgehead",
        "geometry_source": "published UCCSD(T)/cc-pVDZ structures, unmodified",
        "is_a_dft_barrier": False,
        "why_not": (
            "The supplied structure is a coupled-cluster stationary point and is not stationary "
            "on any DFT surface; under PBE0 the analogous methane structure shows a 2.38 eV/A "
            "residual force. These are single points at a fixed geometry, not barriers."
        ),
        "coupled_cluster_reference": COUPLED_CLUSTER,
        "timing_caveat": (
            "Host carried 200+ load average during these runs; every elapsed time describes a "
            "job receiving a fraction of one core, not the cost of the calculation."
        ),
        "methods": [],
        "status": "running",
    }
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

    for xc, basis, dispersion in METHODS:
        label = f"{xc}/{basis}" + (f"+{dispersion}" if dispersion else "")
        try:
            energies, records = {}, []
            for name, filename, charge, spin in SPECIES:
                value, record = energy_ev(name, filename, charge, spin, xc, basis, dispersion)
                energies[name] = value
                records.append(record)
            reactants = energies["isobutane"] + energies["ethynyl_radical"]
            barrier = (energies["isobutane_ethynyl_ts"] - reactants) * KCAL_PER_EV
            reaction = (energies["acetylene"] + energies["tert_butyl_radical"] - reactants) * KCAL_PER_EV
            entry = {
                "method": label, "xc": xc, "basis": basis, "dispersion": dispersion,
                "reference_geometry_energy_kcal_per_mol": barrier,
                "reaction_energy_kcal_per_mol": reaction,
                "species": records, "status": "completed",
            }
            print(f"{label:<24} fixed-geometry {barrier:+8.2f}   reaction {reaction:+8.2f} kcal/mol", flush=True)
        except Exception as error:
            entry = {"method": label, "status": "failed", "error": str(error)}
            print(f"{label:<24} FAILED: {str(error)[:90]}", flush=True)
        result["methods"].append(entry)
        output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

    result["status"] = "completed"
    completed = [m for m in result["methods"] if m["status"] == "completed"]
    if completed:
        result["summary"] = {
            "all_below_reactants": all(m["reference_geometry_energy_kcal_per_mol"] < 0 for m in completed),
            "range_kcal_per_mol": [
                min(m["reference_geometry_energy_kcal_per_mol"] for m in completed),
                max(m["reference_geometry_energy_kcal_per_mol"] for m in completed),
            ],
            "interpretation": (
                "If these sit below the reactants as they do for methane, the behaviour is not "
                "specific to methane's primary C-H and carries to the tertiary site the candidate "
                "targets. That still reflects evaluating DFT at a coupled-cluster geometry rather "
                "than a demonstrated functional failure; locating DFT saddles is what would settle it."
            ),
        }
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"\nwrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
