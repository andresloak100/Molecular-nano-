"""Stages 1 and 2: relaxed site energetics for bridgehead vs methylene C-H.

Stage 1 is the adiabatic C-H dissociation energy at each site,

    D(site) = E(adamantyl radical) + E(H atom) - E(adamantane)

with every species independently relaxed at one consistent method, so each D is
an *adiabatic electronic* dissociation energy.  No zero-point, thermal,
entropic, tunnelling or solvent term is included, and none is estimated.

Stage 2 is the same comparison with the real abstracting species,

    dE(site) = E(acetylene) + E(adamantyl radical) - E(ethynyl) - E(adamantane)

Both cycles share the two adamantyl radical energies, so the *difference between
sites* is algebraically identical in the two cycles:

    D(methylene) - D(bridgehead) = dE(methylene) - dE(bridgehead)
                                 = E(adamantyl_methylene) - E(adamantyl_bridgehead)

The hydrogen atom, ethynyl and acetylene terms cancel exactly.  This script
computes both cycles independently and asserts that identity numerically, which
is the arithmetic self-check Codex asked for.  It also means stage 2 supplies no
new *thermodynamic* site information: any abstractor-specific site preference
must come from activation barriers or accessibility, not from reaction energies.

All quantities here are thermodynamic.  Nothing in this file is a barrier, a
rate, or a selectivity.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import time

import numpy as np
from ase.io import read as ase_read, write
from ase.optimize import BFGS
from scipy.constants import Avogadro, calorie, electron_volt

from nanodesign.quantum import PySCFCalculator, QuantumSettings

from species import OPEN_SHELL, build_species

# 1 eV expressed in kcal/mol, derived rather than hard-coded.
EV_TO_KCAL_PER_MOL = electron_volt * Avogadro / (1000.0 * calorie)


def _settings(spin: int, *, xc: str, basis: str, density_fit: bool) -> QuantumSettings:
    return QuantumSettings(
        charge=0,
        spin=spin,
        xc=xc,
        basis=basis,
        dispersion="d3bj",
        grid_level=3,
        conv_tol=1e-9,
        max_cycle=150,
        threads=1,
        memory_mb=3000,
        density_fit=density_fit,
    )


def relax_species(species, directory: Path, *, xc: str, basis: str, density_fit: bool,
                  fmax: float, max_steps: int) -> dict:
    """Relax one species and return its record, including any failure."""
    directory.mkdir(parents=True, exist_ok=True)
    atoms = species.atoms.copy()
    settings = _settings(species.spin, xc=xc, basis=basis, density_fit=density_fit)
    calculator = PySCFCalculator(settings, event_log=directory / "electronic.jsonl")
    atoms.calc = calculator

    record = {
        "name": species.name,
        "description": species.description,
        "site_type": species.site_type,
        "formula": species.formula,
        "atom_count": len(atoms),
        "spin_n_alpha_minus_n_beta": species.spin,
        "charge": species.charge,
        "settings": settings.to_dict(),
        "relaxation": {
            "optimizer": "ASE BFGS",
            "fmax_ev_per_angstrom": fmax,
            "max_steps": max_steps,
        },
        "geometry_treatment": "relaxed (adiabatic)",
        # Provenance the shared evidence audit asked for explicitly: the one
        # corrupted committed artifact in this repository exists because this
        # field was missing.  PySCFCalculator does not expose an initial-guess
        # setting, so every DFT energy here uses PySCF's default, and the
        # guess-independence test lives in scf_guess_scan.py.
        "scf_initial_guess": "pyscf default (minao); not configurable through PySCFCalculator",
        "scf_initial_guess_verified_by": "scf_guess_scan.py over minao/atom/huckel/1e",
        "completed": False,
    }
    write(directory / "input.xyz", atoms)
    started = time.monotonic()
    try:
        if len(atoms) == 1:
            # A single atom has no geometry to relax; evaluate it directly so
            # the record still carries a real force residual of exactly zero.
            energy = float(atoms.get_potential_energy())
            steps = 0
            converged = True
        else:
            optimizer = BFGS(
                atoms,
                trajectory=str(directory / "relax.traj"),
                logfile=str(directory / "relax.log"),
            )
            converged = bool(optimizer.run(fmax=fmax, steps=max_steps))
            steps = int(optimizer.get_number_of_steps())
            energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=float)
        force_max = float(np.linalg.norm(forces, axis=1).max())
        record.update(
            completed=True,
            relaxation_converged=converged,
            optimizer_steps=steps,
            energy_ev=energy,
            energy_hartree=calculator.diagnostics.get("total_energy_hartree"),
            max_force_ev_per_angstrom=force_max,
            force_residual_within_fmax=bool(force_max <= fmax),
        )
        write(directory / "relaxed.xyz", atoms)
    except Exception as error:  # Preserve the failure as a result.
        record.update(
            completed=False,
            error_type=type(error).__name__,
            error=str(error),
        )
    finally:
        diagnostics = dict(calculator.diagnostics)
        record["wall_seconds"] = time.monotonic() - started
        record["diagnostics"] = {
            key: diagnostics.get(key)
            for key in (
                "scf_converged", "gradient_completed", "basis_functions", "electron_count",
                "alpha_electrons", "beta_electrons", "reference", "s2", "expected_s2",
                "spin_contamination", "effective_multiplicity", "requested_pyscf_threads",
                "effective_pyscf_threads", "threads_honored", "elapsed_seconds",
                "dispersion_energy_hartree", "versions",
            )
        }
        (directory / "result.json").write_text(
            json.dumps(record, indent=2, allow_nan=False) + "\n"
        )
    return record


def _kcal(value_ev: float) -> float:
    return value_ev * EV_TO_KCAL_PER_MOL


def analyse(records: dict[str, dict]) -> dict:
    """Both thermochemical cycles, plus the exact cancellation identity."""
    missing = [name for name, record in records.items() if not record.get("completed")]
    if missing:
        return {
            "status": "incomplete",
            "missing_or_failed": missing,
            "note": "No site comparison is reported while any species is missing or failed.",
        }
    energy = {name: record["energy_ev"] for name, record in records.items()}

    dissociation = {}
    abstraction = {}
    for site, radical in (
        ("bridgehead_tertiary", "adamantyl_bridgehead"),
        ("methylene_secondary", "adamantyl_methylene"),
    ):
        dissociation[site] = _kcal(
            energy[radical] + energy["hydrogen_atom"] - energy["adamantane"]
        )
        abstraction[site] = _kcal(
            energy["acetylene"] + energy[radical] - energy["ethynyl"] - energy["adamantane"]
        )

    dissociation_difference = (
        dissociation["methylene_secondary"] - dissociation["bridgehead_tertiary"]
    )
    abstraction_difference = (
        abstraction["methylene_secondary"] - abstraction["bridgehead_tertiary"]
    )
    radical_difference = _kcal(
        energy["adamantyl_methylene"] - energy["adamantyl_bridgehead"]
    )
    identity_residual = max(
        abs(dissociation_difference - radical_difference),
        abs(abstraction_difference - radical_difference),
    )

    preferred = (
        "bridgehead_tertiary" if dissociation_difference > 0 else "methylene_secondary"
    )
    return {
        "status": "complete",
        "units": "kcal/mol",
        "quantity_type": "adiabatic electronic energies; no zero-point, thermal or entropic terms",
        "stage_1_adiabatic_dissociation_energy": dissociation,
        "stage_2_ethynyl_abstraction_reaction_energy": abstraction,
        "site_difference": {
            "from_dissociation_cycle": dissociation_difference,
            "from_abstraction_cycle": abstraction_difference,
            "from_radical_energies_alone": radical_difference,
            "identity_residual_kcal_per_mol": identity_residual,
            "identity_holds": bool(identity_residual < 1e-6),
            "identity_note": (
                "The hydrogen-atom, ethynyl and acetylene terms cancel exactly in the "
                "site difference, so the two cycles must agree to numerical precision. "
                "Stage 2 therefore adds no new thermodynamic site information."
            ),
            "thermodynamically_preferred_site": preferred,
            "magnitude_kcal_per_mol": abs(dissociation_difference),
        },
        "energies_ev": energy,
        "spin_diagnostics": {
            name: {
                "s2": records[name]["diagnostics"].get("s2"),
                "expected_s2": records[name]["diagnostics"].get("expected_s2"),
                "spin_contamination": records[name]["diagnostics"].get("spin_contamination"),
            }
            for name in records
        },
        "force_residuals_ev_per_angstrom": {
            name: records[name]["max_force_ev_per_angstrom"] for name in records
        },
        "interpretation_limits": [
            "These are thermodynamic reaction energies, not barriers. A thermodynamic site preference is not a kinetic selectivity and does not order the two abstraction rates.",
            "This project's DFT has not reproduced a correct barrier for its calibration reaction, so DFT site energetics are screening values pending a high-level check.",
            "Reaction energies omit zero-point energy, which differs between a broken C-H at a tertiary and a secondary carbon and is not negligible at this scale.",
            "An energy difference expressed in units of k_B T is an energy-scale comparison only. It is not converted here into a success probability or an equilibrium population, because a driven positional operation is not an equilibrium process.",
            "Both radicals are isolated relaxed molecules. The mounted, constrained radical in the 53-atom candidate is a different system.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="New evidence directory (never named 'runs')")
    parser.add_argument("--xc", default="pbe0")
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--density-fit", action="store_true")
    parser.add_argument("--fmax", type=float, default=0.01)
    parser.add_argument("--max-steps", type=int, default=200)
    parser.add_argument("--resume", action="store_true",
                        help="Reuse an existing directory, skipping species already completed")
    parser.add_argument("--only", nargs="*", default=None, help="Restrict to named species")
    parser.add_argument("--no-seed", action="store_true",
                        help="Build radicals from the ideal-lattice cage instead of the relaxed one")
    arguments = parser.parse_args()

    output = Path(arguments.output).resolve()
    if output.name == "runs":
        raise SystemExit("Refusing to write to a directory named 'runs'; .gitignore drops it")
    output.mkdir(parents=True, exist_ok=arguments.resume)

    species = build_species()
    names = list(species) if arguments.only is None else list(arguments.only)
    for name in names:
        if name not in species:
            raise SystemExit(f"Unknown species {name!r}; known: {sorted(species)}")
    # Seeding the radicals from a relaxed cage saves optimizer steps without
    # changing the converged answer, since each radical still relaxes to the
    # same force threshold on its own.
    seed_from_relaxed_cage = "adamantane" in names and not arguments.no_seed

    manifest = {
        "stage": "1_and_2_site_energetics",
        "method": f"{arguments.xc.upper()}-D3(BJ)/{arguments.basis}",
        "density_fit": bool(arguments.density_fit),
        "fmax_ev_per_angstrom": arguments.fmax,
        "species_planned": names,
        "host": {"platform": platform.platform(), "machine": platform.machine()},
        "concurrency_caveat": (
            "Other sessions were computing on this eight-core host during these runs, and this "
            "PySCF build exposes one thread, so wall times are not a controlled speed benchmark."
        ),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    records: dict[str, dict] = {}
    for name in names:
        directory = output / name
        existing = directory / "result.json"
        if arguments.resume and existing.exists():
            record = json.loads(existing.read_text())
            if record.get("completed"):
                records[name] = record
                print(f"[skip] {name}: already complete", flush=True)
                continue
        target = species[name]
        if seed_from_relaxed_cage and name.startswith("adamantyl_"):
            relaxed_cage = output / "adamantane" / "relaxed.xyz"
            if relaxed_cage.exists():
                reseeded = build_species(relaxed_cage=ase_read(relaxed_cage))[name]
                target = reseeded
                print(f"[seed] {name} seeded from relaxed adamantane", flush=True)
        print(f"[run ] {name} ({target.formula}, spin={target.spin})", flush=True)
        record = relax_species(
            target, directory,
            xc=arguments.xc, basis=arguments.basis, density_fit=arguments.density_fit,
            fmax=arguments.fmax, max_steps=arguments.max_steps,
        )
        records[name] = record
        if record["completed"]:
            print(
                f"[done] {name}: E = {record['energy_ev']:.6f} eV, "
                f"{record['optimizer_steps']} steps, "
                f"|F|max = {record['max_force_ev_per_angstrom']:.4f}, "
                f"{record['wall_seconds']:.0f} s, "
                f"nbf = {record['diagnostics'].get('basis_functions')}, "
                f"S2 = {record['diagnostics'].get('s2')}",
                flush=True,
            )
        else:
            print(f"[FAIL] {name}: {record.get('error_type')}: {record.get('error')}", flush=True)

    summary = analyse(records)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary.get("site_difference", summary), indent=2), flush=True)
    print("wrote", output / "summary.json", flush=True)


if __name__ == "__main__":
    main()
