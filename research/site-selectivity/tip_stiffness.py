"""Measure a mounted ethynyl tip's lateral stiffness instead of assuming it.

The headroom claim in this lane -- 4.4 N/m required against 10 to 100 N/m
available -- has an assumption doing all of its work: that a diamondoid mount
supplies 10 to 100 N/m.  That figure is quoted from the positional-assembly
literature and has never been computed for this tool.  If the real tip were 2
N/m the conclusion would invert.  This script computes it.

**Model, and why it is the conservative choice.**  A full Hessian on the 53-atom
candidate is out of reach here, so the mount is reduced to a methyl group:
CH3-C(triple)C-H, propyne, with the methyl hydrogens anchored.  A methyl group is
*far floppier* than an adamantane cage, so the stiffness this returns is a
**lower estimate** of what the real mounted tip provides.  That is the useful
direction: if even a methyl-mounted ethynyl clears the requirement, a
cage-mounted one clears it comfortably.

**Why propyne rather than the propynyl radical.**  The quantity wanted is a
framework mechanical property -- how stiffly the mount resists lateral
displacement of the apex -- not a property of the unpaired electron.  Using the
closed-shell propyne avoids the open-shell SCF hazard entirely, and that hazard is
real for exactly this species: a peer lane measured propynyl at
PBE0-D3(BJ)/def2-SVP with the default ``minao`` guess converging 11.2 kcal/mol
above the solution ``atom``, ``huckel`` and ``1e`` agree on, and
``PySCFCalculator`` exposes no initial-guess setting.  So the radical version
cannot currently be computed reliably through the core calculator at all.  The
substitution is recorded as a limitation rather than hidden: the C(triple)C
framework stiffness is similar between propyne and propynyl, but it is not
identical and this is a proxy.

**Two approximations that push in opposite directions**, stated together so
neither is mistaken for the whole error:

- The methyl mount is floppier than a cage, so the result is *too soft*.
- The anchors are rigid, so the compliance is *too stiff*.

Both are bounded in direction but not in size, so this is an estimate with a
known sign structure, not a converged number.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import write
from ase.optimize import BFGS

from nanodesign.quantum import PySCFCalculator, QuantumSettings
from nanodesign.stationary import characterize_stationary_point

from positional_requirements import assess, mode_wavenumber
from positional_uncertainty import atom_uncertainty, effective_stiffness

# The margin this lane measured, census r5.
MARGIN_ANGSTROM = 2.4949626493635075

# Propyne starting geometry, Angstrom. Standard bond lengths; relaxed before use.
_CC_TRIPLE = 1.21
_C_METHYL = 1.46
_CH_METHYL = 1.09
_CH_ACETYLENIC = 1.06


def build_propyne() -> tuple[Atoms, dict]:
    """CH3-C(triple)C-H along +z, methyl at the origin, apex at the top."""
    methyl_carbon = np.array([0.0, 0.0, 0.0])
    inner_carbon = np.array([0.0, 0.0, _C_METHYL])
    apex_carbon = inner_carbon + np.array([0.0, 0.0, _CC_TRIPLE])
    apex_hydrogen = apex_carbon + np.array([0.0, 0.0, _CH_ACETYLENIC])

    # Three methyl hydrogens, tetrahedral, pointing away from the chain.
    tetrahedral_z = -1.0 / 3.0
    radial = float(np.sqrt(1.0 - tetrahedral_z**2))
    hydrogens = []
    for index in range(3):
        angle = 2.0 * np.pi * index / 3.0
        direction = np.array([radial * np.cos(angle), radial * np.sin(angle), tetrahedral_z])
        hydrogens.append(methyl_carbon + direction * _CH_METHYL)

    positions = np.vstack([methyl_carbon, inner_carbon, apex_carbon, apex_hydrogen, *hydrogens])
    atoms = Atoms("CCCH" + "H" * 3, positions=positions, pbc=False)
    indices = {
        "methyl_carbon": 0,
        "inner_carbon": 1,
        "apex_carbon": 2,
        "apex_hydrogen": 3,
        "mount_hydrogens": [4, 5, 6],
    }
    if atoms.get_chemical_formula() != "C3H4":
        raise ValueError(f"Expected propyne C3H4, built {atoms.get_chemical_formula()}")
    return atoms, indices


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fmax", type=float, default=0.005,
                        help="Tight, because a Hessian is only meaningful at a minimum")
    parser.add_argument("--step", type=float, default=0.005, help="Finite-difference step, Angstrom")
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--xc", default="pbe0")
    arguments = parser.parse_args()

    output = Path(arguments.output).resolve()
    if output.name == "runs":
        raise SystemExit("Refusing to write to a directory named 'runs'")
    output.mkdir(parents=True, exist_ok=False)

    atoms, indices = build_propyne()
    settings = QuantumSettings(
        charge=0, spin=0, xc=arguments.xc, basis=arguments.basis, dispersion="d3bj",
        grid_level=3, conv_tol=1e-10, max_cycle=200, threads=1, memory_mb=3000,
        density_fit=False,
    )
    calculator = PySCFCalculator(settings, event_log=output / "electronic.jsonl")
    atoms.calc = calculator

    write(output / "input.xyz", atoms)
    started = time.monotonic()
    print("relaxing propyne (closed shell, no anchors)...", flush=True)
    optimizer = BFGS(atoms, trajectory=str(output / "relax.traj"), logfile=str(output / "relax.log"))
    relaxed = bool(optimizer.run(fmax=arguments.fmax, steps=300))
    forces = np.asarray(atoms.get_forces(), dtype=float)
    relax_residual = float(np.linalg.norm(forces, axis=1).max())
    write(output / "relaxed.xyz", atoms)
    print(f"  relaxed={relaxed}  |F|max={relax_residual:.5f} eV/A  "
          f"{optimizer.get_number_of_steps()} steps  {time.monotonic() - started:.0f} s", flush=True)

    geometry = {
        "c_c_triple_angstrom": float(atoms.get_distance(indices["inner_carbon"], indices["apex_carbon"])),
        "methyl_c_c_angstrom": float(atoms.get_distance(indices["methyl_carbon"], indices["inner_carbon"])),
        "apex_c_h_angstrom": float(atoms.get_distance(indices["apex_carbon"], indices["apex_hydrogen"])),
    }
    print("  relaxed geometry:", {k: round(v, 4) for k, v in geometry.items()}, flush=True)

    # Anchor the mount. Three fixed atoms remove translation and rotation, so the
    # remaining curvature is internal and the compliance is well defined.
    anchored = atoms.copy()
    anchored.calc = calculator
    anchored.set_constraint(FixAtoms(indices=indices["mount_hydrogens"]))

    free_count = len(anchored) - len(indices["mount_hydrogens"])
    print(f"characterizing Hessian: {3 * free_count} free coordinates, "
          f"{1 + 6 * free_count} force evaluations...", flush=True)
    hessian_started = time.monotonic()
    characterization = characterize_stationary_point(
        anchored,
        step_angstrom=arguments.step,
        force_tolerance_ev_per_angstrom=max(arguments.fmax * 4.0, 0.02),
        max_free_coordinates=60,
    )
    print(f"  done in {time.monotonic() - hessian_started:.0f} s; "
          f"classification {characterization['classification']}; "
          f"negative modes {characterization['negative_mode_count']}", flush=True)
    print("  frequencies cm-1:",
          [round(f, 1) for f in characterization["frequencies_cm1"]], flush=True)

    report: dict = {
        "purpose": "Lateral stiffness of a mounted ethynyl tip, measured rather than assumed",
        "model": "propyne CH3-C#C-H with the three methyl hydrogens anchored",
        "method": f"{arguments.xc.upper()}-D3(BJ)/{arguments.basis}",
        "settings": settings.to_dict(),
        "indices": indices,
        "relaxation": {
            "converged": relaxed,
            "fmax_requested_ev_per_angstrom": arguments.fmax,
            "force_residual_ev_per_angstrom": relax_residual,
            "steps": int(optimizer.get_number_of_steps()),
        },
        "relaxed_geometry_angstrom": geometry,
        "hessian": {
            "classification": characterization["classification"],
            "negative_mode_count": characterization["negative_mode_count"],
            "frequencies_cm1": characterization["frequencies_cm1"],
            "free_atom_indices": characterization["free_atom_indices"],
            "frozen_atom_indices": characterization["frozen_atom_indices"],
            "step_angstrom": arguments.step,
            "asymmetry_max_abs": characterization["hessian_asymmetry_max_abs_ev_per_angstrom2"],
        },
        "margin_angstrom": MARGIN_ANGSTROM,
        "tips": {},
    }

    for label in ("apex_carbon", "apex_hydrogen"):
        index = indices[label]
        stiffness = effective_stiffness(characterization, index)
        quantum = atom_uncertainty(characterization, index, 300.0, quantum=True)
        classical = atom_uncertainty(characterization, index, 300.0, quantum=False)
        softest = stiffness["softest_stiffness_n_per_m"]
        verdict = assess(MARGIN_ANGSTROM, softest, 1e-15, 300.0)
        entry = {
            "principal_stiffness_n_per_m": stiffness["principal_stiffness_n_per_m"],
            "softest_stiffness_n_per_m": softest,
            "softest_direction": stiffness["softest_direction"],
            "clamped_diagonal_block_stiffness_n_per_m": stiffness["clamped_diagonal_block_stiffness_n_per_m"],
            "sigma_quantum_angstrom_300K": quantum["largest_sigma_angstrom"],
            "sigma_classical_angstrom_300K": classical["largest_sigma_angstrom"],
            "rms_displacement_quantum_angstrom_300K": quantum["rms_displacement_angstrom"],
            "quantum_over_classical": quantum["largest_sigma_angstrom"] / classical["largest_sigma_angstrom"],
            "margin_in_quantum_sigma": MARGIN_ANGSTROM / quantum["largest_sigma_angstrom"],
            "required_stiffness_n_per_m": verdict["requirement"]["required_stiffness_classical_n_per_m"],
            "requirement_met": verdict["requirement_met"],
            "stiffness_headroom_factor": verdict["stiffness_headroom_factor"],
            "softest_mode_wavenumber_cm1_if_carbon_mass": mode_wavenumber(softest, 12.011),
        }
        report["tips"][label] = entry
        print(f"\n{label}:", flush=True)
        print(f"  stiffness, compliance-based   {softest:8.2f} N/m  (softest of "
              f"{[round(v, 1) for v in stiffness['principal_stiffness_n_per_m']]})", flush=True)
        print(f"  stiffness, clamped block      "
              f"{min(stiffness['clamped_diagonal_block_stiffness_n_per_m']):8.2f} N/m", flush=True)
        print(f"  sigma at 300 K, quantum       {quantum['largest_sigma_angstrom']:8.4f} A", flush=True)
        print(f"  sigma at 300 K, classical     {classical['largest_sigma_angstrom']:8.4f} A", flush=True)
        print(f"  margin / sigma                {entry['margin_in_quantum_sigma']:8.1f}", flush=True)
        print(f"  requirement {entry['required_stiffness_n_per_m']:.2f} N/m  met={entry['requirement_met']}  "
              f"headroom {entry['stiffness_headroom_factor']:.1f}x", flush=True)

    report["limitations"] = [
        "A methyl mount is far floppier than an adamantane cage, so this stiffness is a lower estimate of the real mounted tip and the headroom is understated.",
        "Rigid anchors make the compliance an upper bound on stiffness, which pushes the opposite way; the two approximations have known signs but unbounded sizes.",
        "Closed-shell propyne substitutes for the ethynyl radical, because PySCFCalculator exposes no initial-guess setting and propynyl has a documented 11.2 kcal/mol wrong-solution hazard at this level. The C#C framework stiffness is similar but not identical.",
        "Harmonic only. The softest direction is a bending mode and bending is anharmonic well before a 2.5 Angstrom displacement, so the tail of the distribution is not described by this curvature.",
        "A stiffness that clears the positional requirement says nothing about whether a barrier exists, whether competing chemistry intervenes, or whether this mount can be built.",
        "Three anchored hydrogens is a boundary condition, not a mechanically validated mount.",
    ]
    report["wall_seconds"] = time.monotonic() - started
    report["threads_honored"] = calculator.diagnostics.get("threads_honored")
    report["basis_functions"] = calculator.diagnostics.get("basis_functions")
    (output / "tip_stiffness.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print("\nwrote", output / "tip_stiffness.json", flush=True)


if __name__ == "__main__":
    main()
