"""Reduced tool models built by truncating the repository's own candidate.

The 53-atom candidate mounts the ethynyl tip on a second adamantane cage. That
handle exists so the tip is a mounted species rather than a free radical, and
it is also most of the cost: 27 of 53 atoms. The question A2 has to answer is
whether a smaller handle is a faithful stand-in.

Two distinct things change when the handle shrinks, and they must not be
conflated:

1. **Electronic substituent effect.** What the handle does to the apex
   radical's appetite for a hydrogen. Measurable as the shift in tip
   H-affinity, ``E(R-CC-H) - E(R-CC*) - E(H)``, across R.
2. **Boundary stiffness and sterics.** The handle also carries the fixed
   anchors, resists the tip being pushed back, and occupies space near the
   substrate. A smaller handle changes those too, and an H-affinity comparison
   says nothing about them.

This module supplies structures for (1) only, plus reduced full candidates for
cost measurement. It deliberately does not claim the reduction is safe; the
report quantifies the part it can measure and names the part it cannot.

Geometry provenance: every atom position here is taken from
``nanodesign.candidates.make_h_abstraction`` at the requested pose, or placed
along the same tool axis with the same nominal bond lengths. No separate
geometry guess is introduced. All structures are unrelaxed starting points.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms

from nanodesign.candidates import make_h_abstraction

# Same nominal lengths as the candidate builder, reused rather than re-guessed.
_HANDLE_CC = 1.46
_METHYL_CH = 1.09
_APEX_CH = 1.06

# Fixed layout of the 53-atom candidate.
SUBSTRATE = range(0, 26)
TIP_APEX = 26
TIP_DISTAL = 27
HANDLE_BRIDGEHEAD = 28
HANDLE = range(28, 53)

HANDLES = ("hydrogen", "methyl", "adamantyl")


def _methyl_hydrogens(carbon: np.ndarray, axis: np.ndarray) -> list[np.ndarray]:
    """Three staggered H positions on a methyl carbon pointing away from ``axis``.

    ``axis`` is the unit vector from the methyl carbon toward the atom it is
    bonded to (the distal ethynyl carbon). The hydrogens sit at the
    tetrahedral angle from that bond.
    """
    reference = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(reference, axis))) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    perpendicular = np.cross(axis, reference)
    perpendicular /= np.linalg.norm(perpendicular)
    second = np.cross(axis, perpendicular)
    tetrahedral = math.radians(109.4712206)
    positions = []
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        radial = math.cos(phi) * perpendicular + math.sin(phi) * second
        direction = math.cos(tetrahedral) * axis + math.sin(tetrahedral) * radial
        positions.append(carbon + _METHYL_CH * direction / np.linalg.norm(direction))
    return positions


def tool_fragment(handle: str, hydrogenated: bool, separation: float = 3.6) -> tuple[Atoms, dict[str, Any]]:
    """Isolated tool ``R-CC*`` or ``R-CC-H``, substrate removed.

    The substrate is absent on purpose: this fragment isolates what the handle
    does to the tip electronically. It is not the candidate and carries no
    anchors, so it says nothing about mechanical boundary conditions.
    """
    if handle not in HANDLES:
        raise ValueError(f"handle must be one of {HANDLES}, not {handle!r}")

    reactant, _, metadata = make_h_abstraction(separation, 0.0)
    positions = reactant.positions
    apex = positions[TIP_APEX]
    distal = positions[TIP_DISTAL]
    axis = distal - apex
    axis = axis / np.linalg.norm(axis)

    symbols = ["C", "C"]
    coordinates = [apex.copy(), distal.copy()]

    if handle == "adamantyl":
        for index in HANDLE:
            symbols.append(reactant.symbols[index])
            coordinates.append(positions[index].copy())
    elif handle == "methyl":
        carbon = distal + _HANDLE_CC * axis
        symbols.append("C")
        coordinates.append(carbon)
        for hydrogen in _methyl_hydrogens(carbon, -axis):
            symbols.append("H")
            coordinates.append(hydrogen)
    else:  # hydrogen
        symbols.append("H")
        coordinates.append(distal + 1.06 * axis)

    info: dict[str, Any] = {
        "handle": handle,
        "hydrogenated": hydrogenated,
        "apex_index": 0,
        "distal_index": 1,
        "geometry_status": "unrelaxed_starting_guess_from_candidate_pose",
        "substrate_removed": True,
        "scope": (
            "Isolated tool fragment for the electronic substituent comparison "
            "only. No anchors, no substrate, no mechanical boundary condition."
        ),
        "separation_angstrom_of_source_pose": separation,
        "source_candidate": metadata["candidate"],
    }

    if hydrogenated:
        # The abstracted H arrives on the apex, anti-parallel to the handle.
        symbols.append("H")
        coordinates.append(apex - _APEX_CH * axis)
        info["added_hydrogen_index"] = len(symbols) - 1
        charge, spin = 0, 0
    else:
        charge, spin = 0, 1

    atoms = Atoms(symbols=symbols, positions=np.array(coordinates), pbc=False)
    moments = np.zeros(len(atoms))
    if spin:
        moments[0] = 1.0  # unpaired electron on the apex carbon
    atoms.set_initial_magnetic_moments(moments)
    info.update(total_charge=charge, spin=spin, spin_multiplicity=spin + 1, formula=atoms.get_chemical_formula())
    atoms.info.update(info)
    return atoms, info


def reduced_candidate(handle: str, separation: float = 3.6, lateral: float = 0.0) -> tuple[Atoms, dict[str, Any]]:
    """The full candidate with its adamantane handle replaced.

    Substrate, pose and tip are untouched; only the mount changes. Anchors:
    the substrate keeps its three fixed carbons. The reduced handles have no
    skeleton to anchor, so the terminal handle atom is fixed instead. That is a
    *different* mechanical boundary condition, not a smaller version of the
    same one, and the report says so rather than treating the models as
    interchangeable mounts.
    """
    if handle not in HANDLES:
        raise ValueError(f"handle must be one of {HANDLES}, not {handle!r}")

    reactant, _, metadata = make_h_abstraction(separation, lateral)
    if handle == "adamantyl":
        # The candidate's own metadata names the electronic state
        # "spin_multiplicity"; the reduced variants below use "spin" (the
        # PySCF convention, N_alpha - N_beta). Normalize so every handle
        # returns one schema, rather than making callers branch on which
        # variant they asked for.
        info = dict(metadata)
        info["handle"] = handle
        info["n_atoms"] = len(reactant)
        info["spin"] = metadata["spin_multiplicity"] - 1
        info["anchor_scheme"] = "three distal skeletal carbons in each cage (original candidate)"
        return reactant, info

    positions = reactant.positions
    apex = positions[TIP_APEX]
    distal = positions[TIP_DISTAL]
    axis = (distal - apex) / np.linalg.norm(distal - apex)

    symbols = [reactant.symbols[i] for i in SUBSTRATE] + ["C", "C"]
    coordinates = [positions[i].copy() for i in SUBSTRATE] + [apex.copy(), distal.copy()]

    if handle == "methyl":
        carbon = distal + _HANDLE_CC * axis
        symbols.append("C")
        coordinates.append(carbon)
        handle_anchor = len(symbols) - 1
        for hydrogen in _methyl_hydrogens(carbon, -axis):
            symbols.append("H")
            coordinates.append(hydrogen)
    else:  # hydrogen
        symbols.append("H")
        coordinates.append(distal + 1.06 * axis)
        handle_anchor = len(symbols) - 1

    atoms = Atoms(symbols=symbols, positions=np.array(coordinates), pbc=False)
    fixed = list(metadata["target_anchor_indices"]) + [handle_anchor]
    atoms.set_constraint(FixAtoms(indices=fixed))
    moments = np.zeros(len(atoms))
    moments[TIP_APEX] = 1.0
    atoms.set_initial_magnetic_moments(moments)

    info = {
        "handle": handle,
        "formula": atoms.get_chemical_formula(),
        "n_atoms": len(atoms),
        "target_carbon": metadata["target_carbon"],
        "transferred_hydrogen": metadata["transferred_hydrogen"],
        "tip_apex": TIP_APEX,
        "tip_distal": TIP_DISTAL,
        "fixed_indices": fixed,
        "total_charge": 0,
        "spin": 1,
        "spin_multiplicity": 2,
        "separation_angstrom": separation,
        "lateral_offset_angstrom": lateral,
        "geometry_status": "unrelaxed_starting_guess",
        "anchor_scheme": (
            "Substrate keeps its three fixed cage carbons; the reduced handle "
            "has a single fixed terminal atom. This is a different mechanical "
            "boundary condition from the two-cage candidate, not a scaled one."
        ),
        "scope": (
            "Cost-measurement and screening model. Reduced handle changes the "
            "mount's stiffness, mass and sterics as well as its electronics."
        ),
    }
    atoms.info.update(info)
    return atoms, info
