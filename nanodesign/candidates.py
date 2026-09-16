"""Atom-conserving, unoptimized starting structures for a small research model.

The cages have the connectivity of adamantane, constructed directly as a
finite fragment of the ideal cubic diamond network.  This does not supply an
optimized tool, a bulk diamond surface, or evidence that H transfer is viable.
"""

from __future__ import annotations

import math
from dataclasses import asdict
import json
from pathlib import Path
from typing import Any

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import write

from .quantum import QuantumSettings


_CAGE_CC = 1.54
_CAGE_CH = 1.09
_ETHYNYL_CC = 1.21
_HANDLE_CC = 1.46
_PRODUCT_CH = 1.06


def create_design(output, separation=3.6, offset=0.0, settings=None):
    """Save an explicit unrelaxed candidate without overwriting existing work."""
    settings = settings or QuantumSettings()
    if settings.charge != 0 or settings.spin != 1:
        raise ValueError("The generated H-abstraction candidate requires a neutral doublet.")
    initial, final, metadata = make_h_abstraction(separation, offset)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    write(out / "initial.xyz", initial)
    write(out / "final.xyz", final)
    design = {
        "schema_version": 1, "length_unit": "angstrom",
        "name": "Adamantane-supported ethynyl H-abstraction candidate",
        "scope": "Finite diamondoid cluster in vacuum with fixed distal carbon anchors. Unrelaxed candidate, not a diamond surface or validated assembly tool.",
        "initial": "initial.xyz", "final": "final.xyz", "fixed_indices": metadata["fixed_indices"],
        "hydrogen_transfer": {"donor": metadata["target_carbon"], "hydrogen": metadata["transferred_hydrogen"], "acceptor": metadata["tip_apex"]},
        "quantum": asdict(settings), "metadata": metadata,
    }
    (out / "design.json").write_text(json.dumps(design, indent=2, allow_nan=False) + "\n")
    return out / "design.json"


def _adamantane() -> tuple[list[str], np.ndarray, list[list[int]]]:
    """Return C10H16 with bridgehead C=0 at the origin and its H on +z.

    Four bridgeheads occupy alternating cube corners; six methylene carbons
    occupy the corresponding external face-axis sites.  The 12 cage edges
    all have length ``_CAGE_CC`` and all cage bond angles are tetrahedral.
    Missing tetrahedral neighbors are passivated with hydrogen.
    """
    corners = np.array(
        [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
        dtype=float,
    )
    axes = np.array(
        [[2, 0, 0], [0, 2, 0], [0, 0, 2], [-2, 0, 0], [0, -2, 0], [0, 0, -2]],
        dtype=float,
    )
    lattice = np.vstack((corners, axes))
    carbon_positions = lattice * (_CAGE_CC / math.sqrt(3.0))
    positions = list(carbon_positions.copy())
    symbols = ["C"] * 10
    bonds: list[list[int]] = []
    neighbors: list[list[int]] = [[] for _ in range(10)]
    for i in range(10):
        for j in range(i + 1, 10):
            if np.isclose(np.linalg.norm(lattice[i] - lattice[j]), math.sqrt(3.0)):
                bonds.append([i, j, 1])
                neighbors[i].append(j)
                neighbors[j].append(i)
    for i in range(10):
        # Opposite sublattices have opposite sets of tetrahedral directions.
        directions = corners if i < 4 else -corners
        occupied = [lattice[j] - lattice[i] for j in neighbors[i]]
        for direction in directions:
            if any(np.array_equal(direction, vector) for vector in occupied):
                continue
            h_index = len(symbols)
            symbols.append("H")
            positions.append(carbon_positions[i] + direction * (_CAGE_CH / math.sqrt(3.0)))
            bonds.append([i, h_index, 1])

    # An orthonormal right-handed basis aligns the chosen bridgehead bond
    # with +z without changing any cage or C-H bond length or angle.
    rotation = np.array(
        [[1 / math.sqrt(2), -1 / math.sqrt(2), 0],
         [1 / math.sqrt(6), 1 / math.sqrt(6), -2 / math.sqrt(6)],
         [1 / math.sqrt(3), 1 / math.sqrt(3), 1 / math.sqrt(3)]]
    )
    transformed = (np.asarray(positions) - carbon_positions[0]) @ rotation.T
    transformed[np.abs(transformed) < 1e-14] = 0.0
    return symbols, transformed, bonds


def make_h_abstraction(
    separation_angstrom: float = 3.6,
    lateral_offset_angstrom: float = 0.0,
) -> tuple[Atoms, Atoms, dict[str, Any]]:
    """Make unoptimized endpoints for one neutral-doublet H-transfer candidate.

    The reactant is adamantane plus an adamantane-bound ethynyl radical,
    ``C10H16 + C10H15-C≡C•``.  The product transfers the *same* bridgehead H
    to the ethynyl apex, leaving a target adamantyl radical.  There are 22 C
    and 31 H atoms in each endpoint, with identical atom ordering.

    ``separation_angstrom`` is the z separation between target carbon and
    apex carbon; at nonzero lateral offset their actual distance is larger.
    ``lateral_offset_angstrom`` translates the complete tool along +x.
    The product C≡C-H guess stays collinear with the tool axis.

    Three distal backbone carbons in each cage are held fixed.  This is a
    finite-cluster boundary condition, not a mechanically validated mount.
    All quoted bond lengths are initial guesses, not computed equilibrium
    values.  Both endpoints require electronic-structure relaxation and
    assessment of competing chemistry before any feasibility claim.
    """
    separation = float(separation_angstrom)
    lateral = float(lateral_offset_angstrom)
    if not math.isfinite(separation) or separation <= 0:
        raise ValueError("separation_angstrom must be finite and greater than zero")
    if not math.isfinite(lateral):
        raise ValueError("lateral_offset_angstrom must be finite")

    cage_symbols, cage_positions, cage_bonds = _adamantane()
    target_carbon = 0
    transferred_hydrogen = 10
    tip_apex = 26
    tip_distal = 27
    handle_bridgehead = 28

    # Invert the second cage so its substituted bridgehead points toward
    # the target.  Omit only its outward H to make C10H15-C≡C•.
    handle_indices = [i for i in range(26) if i != transferred_hydrogen]
    handle_map = {old: 28 + new for new, old in enumerate(handle_indices)}
    handle_origin = np.array([lateral, 0.0, separation + _ETHYNYL_CC + _HANDLE_CC])
    handle_positions = -cage_positions[handle_indices] + handle_origin
    symbols = cage_symbols + ["C", "C"] + [cage_symbols[i] for i in handle_indices]
    positions = np.vstack(
        (cage_positions,
         [lateral, 0.0, separation],
         [lateral, 0.0, separation + _ETHYNYL_CC],
         handle_positions)
    )

    reactant_bonds = [bond.copy() for bond in cage_bonds]
    reactant_bonds += [[tip_apex, tip_distal, 3], [tip_distal, handle_bridgehead, 1]]
    reactant_bonds += [
        [handle_map[i], handle_map[j], order]
        for i, j, order in cage_bonds
        if i in handle_map and j in handle_map
    ]
    product_bonds = [
        bond.copy() for bond in reactant_bonds
        if set(bond[:2]) != {target_carbon, transferred_hydrogen}
    ]
    product_bonds.append([tip_apex, transferred_hydrogen, 1])

    target_anchors = [7, 8, 9]
    handle_anchors = [handle_map[i] for i in target_anchors]
    fixed_indices = target_anchors + handle_anchors
    reactant = Atoms(symbols=symbols, positions=positions, pbc=False)
    reactant.set_constraint(FixAtoms(indices=fixed_indices))
    product = reactant.copy()
    product.positions[transferred_hydrogen] = [lateral, 0.0, separation - _PRODUCT_CH]
    for endpoint, label, radical in (
        (reactant, "reactant", tip_apex), (product, "product", target_carbon)
    ):
        magnetic_moments = np.zeros(len(endpoint))
        magnetic_moments[radical] = 1.0
        endpoint.set_initial_magnetic_moments(magnetic_moments)
        endpoint.info.update(
            candidate="finite_adamantane_ethynyl_h_abstraction",
            endpoint=label,
            total_charge=0,
            spin_multiplicity=2,
            unpaired_electrons=1,
            geometry_status="unoptimized_starting_guess",
        )
        distances = endpoint.get_all_distances()
        np.fill_diagonal(distances, np.inf)
        if float(distances.min()) < 0.7:
            raise ValueError(
                f"{label} starting guess has atoms closer than 0.7 Å; "
                "increase separation or change the lateral offset"
            )

    metadata: dict[str, Any] = {
        "candidate": "finite_adamantane_ethynyl_h_abstraction",
        "description": (
            "Unoptimized, atom-conserving H-abstraction endpoints for two finite "
            "diamondoid clusters: an adamantane target and an adamantane-bound "
            "ethynyl radical tool. Ideal diamond-lattice cages with approximate "
            "linkage lengths; not a bulk diamond surface or a validated assembler."
        ),
        "geometry_status": "unoptimized_starting_guess",
        "formula": "C22H31",
        "total_charge": 0,
        "charge": 0,
        "spin_multiplicity": 2,
        "multiplicity": 2,
        "unpaired_electrons": 1,
        "target_carbon": target_carbon,
        "transferred_hydrogen": transferred_hydrogen,
        "tip_apex": tip_apex,
        "tip_distal": tip_distal,
        "handle_bridgehead": handle_bridgehead,
        "substrate_indices": list(range(26)),
        "tool_indices": list(range(26, 53)),
        "fixed_indices": fixed_indices,
        "target_anchor_indices": target_anchors,
        "tool_anchor_indices": handle_anchors,
        "reactant_radical_index": tip_apex,
        "product_radical_index": target_carbon,
        "reactant_bonds": reactant_bonds,
        "product_bonds": product_bonds,
        "bond_list_format": "[atom_index_1, atom_index_2, nominal_bond_order]; zero-based indices",
        "separation_angstrom": separation,
        "lateral_offset_angstrom": lateral,
        "target_apex_distance_angstrom": math.hypot(separation, lateral),
        "initial_bond_lengths_angstrom": {
            "cage_C_C": _CAGE_CC,
            "cage_C_H": _CAGE_CH,
            "ethynyl_C_C": _ETHYNYL_CC,
            "handle_to_ethynyl_C_C": _HANDLE_CC,
            "product_apex_C_H": _PRODUCT_CH,
        },
        "boundary_condition": "Three distal skeletal carbons fixed in each finite cage",
        "periodic": False,
    }
    return reactant, product, metadata
