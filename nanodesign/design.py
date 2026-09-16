"""Explicit structures, mechanical boundary conditions, and atom identities."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.constraints import FixAtoms
from ase.data import covalent_radii
from ase.io import read

from .quantum import QuantumSettings


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_pair(initial, final, fixed_indices):
    if len(initial) < 2 or not np.array_equal(initial.numbers, final.numbers):
        raise ValueError("Endpoints must preserve every atom and its index/element.")
    if initial.pbc.any() or final.pbc.any():
        raise ValueError("This backend supports isolated finite clusters only.")
    if any(type(i) is not int or i < 0 or i >= len(initial) for i in fixed_indices):
        raise ValueError("fixed_indices must be valid zero-based integer atom indices.")
    if len(set(fixed_indices)) != len(fixed_indices):
        raise ValueError("Duplicate fixed atom indices.")
    if len(fixed_indices) == len(initial):
        raise ValueError("At least one atom must be mobile.")
    for atoms in (initial, final):
        if not np.isfinite(atoms.positions).all():
            raise ValueError("Non-finite atomic coordinates.")
        distance = atoms.get_all_distances()
        np.fill_diagonal(distance, np.inf)
        if distance.min() < 0.55:
            raise ValueError("Atoms closer than 0.55 Å: inspect the input geometry.")
    if fixed_indices and not np.allclose(initial.positions[fixed_indices], final.positions[fixed_indices], atol=1e-8, rtol=0):
        raise ValueError("NEB endpoints must have identical mechanical anchors; tool motion requires separate poses.")


def load_design(path):
    path = Path(path).resolve()
    data = json.loads(path.read_text())
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported design schema_version; expected 1.")
    if data.get("length_unit") != "angstrom":
        raise ValueError("Declare length_unit='angstrom'; convert other coordinate units before import.")
    initial_path = (path.parent / data["initial"]).resolve()
    final_path = (path.parent / data["final"]).resolve()
    initial, final = read(initial_path), read(final_path)
    fixed = data.get("fixed_indices", [])
    validate_pair(initial, final, fixed)
    settings = QuantumSettings(**data["quantum"])
    reaction = data.get("hydrogen_transfer")
    if reaction is not None:
        indices = [reaction[k] for k in ("donor", "hydrogen", "acceptor")]
        if any(type(i) is not int or not 0 <= i < len(initial) for i in indices) or len(set(indices)) != 3:
            raise ValueError("Invalid hydrogen-transfer atom identities.")
        if initial[indices[1]].symbol != "H" or initial[indices[0]].symbol != "C" or initial[indices[2]].symbol != "C":
            raise ValueError("Hydrogen-transfer gate currently supports C–H to C transfer only.")
        if indices[1] in fixed:
            raise ValueError("The transferred hydrogen cannot be a fixed anchor.")
    metadata = data.get("metadata", {})
    if not isinstance(metadata, dict):
        raise ValueError("metadata must be an object.")
    groups = []
    for name in ("substrate_indices", "tool_indices"):
        group = metadata.get(name, [])
        if not isinstance(group, list) or any(type(i) is not int or not 0 <= i < len(initial) for i in group) or len(set(group)) != len(group):
            raise ValueError(f"{name} must contain unique in-range integer atom indices.")
        groups.append(set(group))
    if groups[0] & groups[1]:
        raise ValueError("Substrate and tool groups must not overlap.")
    for name in ("reactant_bonds", "product_bonds"):
        bonds = metadata.get(name, [])
        if not isinstance(bonds, list):
            raise ValueError(f"{name} must be a list of [atom_i, atom_j, nominal_order].")
        seen = set()
        for bond in bonds:
            if not isinstance(bond, list) or len(bond) != 3 or any(type(v) is not int for v in bond):
                raise ValueError(f"Invalid nominal bond in {name}.")
            i, j, order = bond
            pair = tuple(sorted((i, j)))
            if not 0 <= i < len(initial) or not 0 <= j < len(initial) or i == j or order not in (1, 2, 3) or pair in seen:
                raise ValueError(f"Invalid or duplicate nominal bond in {name}.")
            seen.add(pair)
    for atoms in (initial, final):
        atoms.set_constraint(FixAtoms(indices=fixed))
    provenance = {"design_sha256": sha256(path), "initial_sha256": sha256(initial_path), "final_sha256": sha256(final_path)}
    return data, initial, final, settings, provenance


def transfer_distances(atoms, reaction):
    if not reaction:
        return {}
    donor, hydrogen, acceptor = (reaction[k] for k in ("donor", "hydrogen", "acceptor"))
    return {"donor_h_angstrom": float(atoms.get_distance(donor, hydrogen)),
            "acceptor_h_angstrom": float(atoms.get_distance(acceptor, hydrogen))}


def endpoint_identity_ok(atoms, reaction, state):
    if not reaction:
        return None
    distances = transfer_distances(atoms, reaction)
    bound, unbound = (("donor_h_angstrom", "acceptor_h_angstrom") if state == "initial" else ("acceptor_h_angstrom", "donor_h_angstrom"))
    return bool(distances[bound] < 1.30 and distances[unbound] > 1.50)


def topology_screen(atoms, bonds):
    """Conservative geometry flag, NOT an electronic bond-order measurement."""
    if not bonds:
        return {"assessed": False, "preserved": None}
    expected = {tuple(sorted(bond[:2])) for bond in bonds}
    distances = atoms.get_all_distances()
    radii = covalent_radii[atoms.numbers]
    stretched, new_contacts = [], []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            ratio = float(distances[i, j] / (radii[i] + radii[j]))
            if (i, j) in expected and ratio > 1.35:
                stretched.append([i, j, float(distances[i, j])])
            elif (i, j) not in expected and ratio < 1.15:
                new_contacts.append([i, j, float(distances[i, j])])
    return {"assessed": True, "preserved": not (stretched or new_contacts),
            "stretched_expected_bonds": stretched, "unexpected_close_contacts": new_contacts,
            "criterion": "Heuristic only: expected bonds <=1.35*(covalent radii sum), other contacts >=1.15*sum; no bond-order or spin-density inference."}
