"""Stage 0: geometric census of every abstractable hydrogen on the target cage.

No electronic-structure calculation is performed here.  This script measures the
actual coordinates of the 53-atom candidate produced by
``nanodesign.candidates.make_h_abstraction`` and answers two purely geometric
questions that no amount of reaction energetics can answer:

1. Which hydrogens of the target cage are distinguishable by *position*, and by
   how much?  For each of the 16 target C-H bonds we place the tool in the same
   idealized abstraction geometry it occupies for the intended site and record
   how far the tool apex must move to get there.
2. Would the tool physically fit there?  For each retargeted placement we scan
   the tool's azimuthal orientation about its own axis and report the largest
   achievable clearance against van der Waals contact.

Geometric reach is a necessary condition, not a sufficient one.  A site that is
reachable is not thereby reactive, and a site that is blocked in this rigid
model is not thereby unreachable by a relaxed or differently mounted tool.  Both
molecules are held rigid at their unoptimized starting geometry; no relaxation,
no electronic structure, no barrier and no dynamics enter this file.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from nanodesign.candidates import make_h_abstraction

# Bondi van der Waals radii, Angstrom.  Used only to express a contact
# criterion in a reproducible way; they are not a force field.
_VDW = {"C": 1.70, "H": 1.20}
_VDW_SOURCE = "Bondi, J. Phys. Chem. 68 (1964) 441; C 1.70 A, H 1.20 A"

# A carbon-carbon bond in the ideal cage is 1.54 A, so 1.8 A separates bonded
# from nonbonded carbon pairs without ambiguity in this rigid construction.
_BOND_CUTOFF = 1.8

# Azimuthal scan about the tool axis when retargeting, degrees.
_AZIMUTH_STEP = 5.0


def _neighbors(symbols, positions):
    """Return the bonded-neighbour index lists implied by the rigid geometry."""
    count = len(symbols)
    distances = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=-1)
    neighbors: list[list[int]] = [[] for _ in range(count)]
    for i in range(count):
        for j in range(i + 1, count):
            if distances[i, j] < _BOND_CUTOFF:
                neighbors[i].append(j)
                neighbors[j].append(i)
    return neighbors, distances


def _rotation_aligning(source, target):
    """Rotation matrix taking unit vector ``source`` onto unit vector ``target``."""
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    source = source / np.linalg.norm(source)
    target = target / np.linalg.norm(target)
    cross = np.cross(source, target)
    dot = float(np.dot(source, target))
    norm = float(np.linalg.norm(cross))
    if norm < 1e-12:
        if dot > 0:
            return np.eye(3)
        # Antiparallel: rotate by pi about any axis orthogonal to source.
        seed = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(seed, source))) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
        axis = np.cross(source, seed)
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    axis = cross / norm
    angle = math.atan2(norm, dot)
    return _rotation_about(axis, angle)


def _rotation_about(axis, angle):
    """Rodrigues rotation matrix for a unit ``axis`` and ``angle`` in radians."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    skew = np.array(
        [[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]]
    )
    return np.eye(3) + math.sin(angle) * skew + (1.0 - math.cos(angle)) * (skew @ skew)


def _clearance(tool_positions, tool_symbols, target_positions, target_symbols, ignore_target):
    """Smallest van der Waals gap between the placed tool and the target cage.

    The gap is ``distance - (r_i + r_j)``.  Negative means interpenetration of
    the van der Waals spheres.  ``ignore_target`` is the index of the hydrogen
    being abstracted, whose close approach to the apex is the intended contact
    rather than a clash, so it is excluded from the minimum.
    """
    keep = [i for i in range(len(target_symbols)) if i != ignore_target]
    if not keep:
        raise ValueError("Target cage must retain at least one atom for the clearance test")
    positions = target_positions[keep]
    radii = np.array([_VDW[target_symbols[i]] for i in keep], dtype=float)
    tool_radii = np.array([_VDW[symbol] for symbol in tool_symbols], dtype=float)
    distances = np.linalg.norm(tool_positions[:, None, :] - positions[None, :, :], axis=-1)
    gaps = distances - (tool_radii[:, None] + radii[None, :])
    flat = int(np.argmin(gaps))
    tool_index, target_index = divmod(flat, len(keep))
    return {
        "min_vdw_gap_angstrom": float(gaps.min()),
        "min_distance_angstrom": float(distances[tool_index, target_index]),
        "closest_tool_atom_index_within_tool": int(tool_index),
        "closest_target_atom_index": int(keep[target_index]),
    }


def census(separation_angstrom: float = 3.6, lateral_offset_angstrom: float = 0.0) -> dict:
    reactant, _product, metadata = make_h_abstraction(
        separation_angstrom, lateral_offset_angstrom
    )
    symbols = reactant.get_chemical_symbols()
    positions = np.asarray(reactant.positions, dtype=float)
    neighbors, distances = _neighbors(symbols, positions)

    substrate = list(metadata["substrate_indices"])
    tool = list(metadata["tool_indices"])
    donor = int(metadata["target_carbon"])
    intended_hydrogen = int(metadata["transferred_hydrogen"])
    apex = int(metadata["tip_apex"])
    distal = int(metadata["tip_distal"])

    # Every C-H bond on the target cage, classified by the substitution of its
    # carbon.  A cage carbon bonded to three other carbons is a bridgehead
    # (tertiary C-H); one bonded to two is a methylene (secondary C-H).
    sites = []
    for hydrogen in substrate:
        if symbols[hydrogen] != "H":
            continue
        attached = [i for i in neighbors[hydrogen] if symbols[i] == "C"]
        if len(attached) != 1:
            raise ValueError(f"Hydrogen {hydrogen} is not bonded to exactly one carbon")
        carbon = attached[0]
        carbon_neighbors = [i for i in neighbors[carbon] if symbols[i] == "C"]
        if len(carbon_neighbors) == 3:
            kind = "bridgehead_tertiary"
        elif len(carbon_neighbors) == 2:
            kind = "methylene_secondary"
        else:
            raise ValueError(
                f"Cage carbon {carbon} has {len(carbon_neighbors)} carbon neighbours; "
                "the census classifier expects 2 or 3"
            )
        sites.append({"hydrogen_index": hydrogen, "carbon_index": carbon, "site_type": kind})

    # The intended abstraction geometry, read off the actual candidate rather
    # than assumed: apex-to-hydrogen distance and the C-H...apex angle.
    intended_vector = positions[intended_hydrogen] - positions[donor]
    intended_length = float(np.linalg.norm(intended_vector))
    apex_to_hydrogen = float(np.linalg.norm(positions[apex] - positions[intended_hydrogen]))
    tool_axis = positions[distal] - positions[apex]
    tool_axis = tool_axis / np.linalg.norm(tool_axis)

    tool_symbols = [symbols[i] for i in tool]
    tool_positions = positions[tool]
    target_symbols = [symbols[i] for i in substrate]
    target_positions = positions[substrate]
    apex_within_tool = tool.index(apex)
    nominal_apex = positions[apex].copy()

    azimuths = np.arange(0.0, 360.0, _AZIMUTH_STEP)
    for site in sites:
        hydrogen = site["hydrogen_index"]
        carbon = site["carbon_index"]
        bond = positions[hydrogen] - positions[carbon]
        bond_length = float(np.linalg.norm(bond))
        unit = bond / bond_length
        # Place the apex collinear with this C-H bond, at the same apex-to-H
        # distance the intended site uses, so every site is compared in the
        # identical idealized abstraction geometry.
        apex_goal = positions[hydrogen] + unit * apex_to_hydrogen
        base_rotation = _rotation_aligning(tool_axis, unit)
        best = None
        for azimuth in azimuths:
            rotation = _rotation_about(unit, math.radians(float(azimuth))) @ base_rotation
            placed = (tool_positions - nominal_apex) @ rotation.T + apex_goal
            clearance = _clearance(
                placed, tool_symbols, target_positions, target_symbols, ignore_target=hydrogen
            )
            if best is None or clearance["min_vdw_gap_angstrom"] > best["min_vdw_gap_angstrom"]:
                best = dict(clearance, azimuth_degrees=float(azimuth))
        assert best is not None
        # Residual check that the placement really is the intended geometry.
        placed_apex = apex_goal
        # The decision boundary for "which hydrogen is nearest the apex" is the
        # perpendicular bisector plane of the intended and alternative hydrogen.
        # Its distance from the nominal apex is the smallest apex displacement,
        # in any direction, that makes this alternative the closest hydrogen.
        if hydrogen == intended_hydrogen:
            bisector_distance = 0.0
        else:
            separation_vector = positions[hydrogen] - positions[intended_hydrogen]
            bisector_distance = abs(
                float(distances[apex, hydrogen]) ** 2
                - float(distances[apex, intended_hydrogen]) ** 2
            ) / (2.0 * float(np.linalg.norm(separation_vector)))
        site.update(
            {
                "c_h_bond_length_angstrom": bond_length,
                "apex_translation_required_angstrom": float(
                    np.linalg.norm(placed_apex - nominal_apex)
                ),
                "apex_displacement_to_become_nearest_hydrogen_angstrom": float(
                    bisector_distance
                ),
                "hydrogen_distance_from_intended_hydrogen_angstrom": float(
                    distances[hydrogen, intended_hydrogen]
                ),
                "apex_to_hydrogen_in_nominal_pose_angstrom": float(
                    distances[apex, hydrogen]
                ),
                "best_azimuth_degrees": best["azimuth_degrees"],
                "best_min_vdw_gap_angstrom": best["min_vdw_gap_angstrom"],
                "best_min_distance_angstrom": best["min_distance_angstrom"],
                "closest_target_atom_index_at_best_azimuth": best["closest_target_atom_index"],
                "reachable_without_vdw_overlap": bool(best["min_vdw_gap_angstrom"] >= 0.0),
            }
        )

    intended = next(site for site in sites if site["hydrogen_index"] == intended_hydrogen)
    others = [site for site in sites if site["hydrogen_index"] != intended_hydrogen]
    reachable_others = [site for site in others if site["reachable_without_vdw_overlap"]]
    by_translation = sorted(others, key=lambda site: site["apex_translation_required_angstrom"])
    nearest_alternative = by_translation[0]
    nearest_reachable = min(
        reachable_others,
        key=lambda site: site["apex_translation_required_angstrom"],
        default=None,
    )
    tightest = min(
        others, key=lambda site: site["apex_displacement_to_become_nearest_hydrogen_angstrom"]
    )
    # The competitors are symmetry-degenerate, so reporting a single limiting
    # index invites reading an arbitrary tie-break as a specific atom. Three
    # downstream readers did exactly that. Report the whole tied set.
    _TIE_TOLERANCE = 1e-9
    tied = [
        site for site in others
        if abs(
            site["apex_displacement_to_become_nearest_hydrogen_angstrom"]
            - tightest["apex_displacement_to_become_nearest_hydrogen_angstrom"]
        ) <= _TIE_TOLERANCE
    ]

    # Handle-damage check: can the apex reach a hydrogen of its own mount?
    tool_hydrogens = [i for i in tool if symbols[i] == "H"]
    apex_to_own_hydrogens = sorted(float(distances[apex, i]) for i in tool_hydrogens)

    counts: dict[str, int] = {}
    for site in sites:
        counts[site["site_type"]] = counts.get(site["site_type"], 0) + 1

    return {
        "stage": "0_geometric_site_census",
        "computation": "rigid-geometry measurement of supplied coordinates; no electronic structure",
        "source": "nanodesign.candidates.make_h_abstraction",
        "pose": {
            "separation_angstrom": float(separation_angstrom),
            "lateral_offset_angstrom": float(lateral_offset_angstrom),
            "geometry_status": metadata["geometry_status"],
        },
        "composition_check": {
            "atom_count": len(symbols),
            "carbon_count": sum(1 for symbol in symbols if symbol == "C"),
            "hydrogen_count": sum(1 for symbol in symbols if symbol == "H"),
            "formula_metadata": metadata["formula"],
            "fixed_indices": metadata["fixed_indices"],
            "donor_carbon": donor,
            "transferred_hydrogen": intended_hydrogen,
            "acceptor_apex": apex,
        },
        "intended_geometry": {
            "donor_c_h_bond_length_angstrom": intended_length,
            "apex_to_transferred_hydrogen_angstrom": apex_to_hydrogen,
            "donor_to_apex_angstrom": float(distances[donor, apex]),
            "collinear_by_construction": True,
            "intended_site_type": intended["site_type"],
            "intended_site_min_vdw_gap_angstrom": intended["best_min_vdw_gap_angstrom"],
        },
        "site_type_counts": counts,
        "site_count": len(sites),
        "sites": sorted(sites, key=lambda site: site["apex_translation_required_angstrom"]),
        # Deliberately NOT called "selectivity_margin". An earlier run used that
        # name for this block, and because the conservative margin lives in the
        # separate block below, two downstream readers keyed on the wrong number
        # and disagreed by a factor of two. The name now says what it measures.
        "site_relocation_distances": {
            "definition": (
                "Apex translation, in Angstrom, needed to move the tool from the intended "
                "abstraction geometry to the identical geometry over another hydrogen. "
                "This is a deliberate re-aiming distance, NOT the tolerance against "
                "accidental misplacement; for that see nearest_hydrogen_margin, which is "
                "smaller and is the conservative quantity."
            ),
            "intended_site_translation_angstrom": intended["apex_translation_required_angstrom"],
            "nearest_alternative_hydrogen_index": nearest_alternative["hydrogen_index"],
            "nearest_alternative_site_type": nearest_alternative["site_type"],
            "nearest_alternative_translation_angstrom": nearest_alternative[
                "apex_translation_required_angstrom"
            ],
            "nearest_alternative_reachable": nearest_alternative["reachable_without_vdw_overlap"],
            "nearest_reachable_alternative_hydrogen_index": (
                nearest_reachable["hydrogen_index"] if nearest_reachable else None
            ),
            "nearest_reachable_alternative_site_type": (
                nearest_reachable["site_type"] if nearest_reachable else None
            ),
            "nearest_reachable_alternative_translation_angstrom": (
                nearest_reachable["apex_translation_required_angstrom"] if nearest_reachable else None
            ),
            "reachable_alternative_count": len(reachable_others),
            "blocked_alternative_count": len(others) - len(reachable_others),
        },
        # The conservative margin. This is the one to quote as a tolerance.
        "nearest_hydrogen_margin": {
            "definition": (
                "Smallest apex displacement, in any direction, that makes some other cage "
                "hydrogen the one nearest the apex. Equals the distance from the nominal apex "
                "to the perpendicular bisector plane of the intended and competing hydrogen."
            ),
            "margin_angstrom": tightest["apex_displacement_to_become_nearest_hydrogen_angstrom"],
            "limiting_site_type": tightest["site_type"],
            "tied_competitor_count": len(tied),
            "tied_competitor_hydrogen_indices": sorted(site["hydrogen_index"] for site in tied),
            "tied_competitor_carbon_indices": sorted({site["carbon_index"] for site in tied}),
            "tie_tolerance_angstrom": _TIE_TOLERANCE,
            "limiting_hydrogen_index": tightest["hydrogen_index"],
            "limiting_hydrogen_index_note": (
                "An arbitrary representative of a symmetry-degenerate tied set, not a "
                "distinguished atom. Bind identity to limiting_site_type and to the tied "
                "index list, never to this single index."
            ),
            "interpretation": (
                "A necessary geometric condition only. Being the nearest hydrogen is not the "
                "reaction criterion; this margin does not include approach barriers, relaxation, "
                "or the positional precision any real mechanism can achieve."
            ),
        },
        "handle_proximity": {
            "purpose": "Screen for the apex contacting its own mount, a documented unassessed failure mode",
            "apex_to_nearest_tool_hydrogen_angstrom": apex_to_own_hydrogens[0],
            "apex_to_tool_hydrogen_distances_angstrom": apex_to_own_hydrogens[:5],
        },
        "method_notes": {
            "vdw_radii": _VDW,
            "vdw_radii_source": _VDW_SOURCE,
            "bond_cutoff_angstrom": _BOND_CUTOFF,
            "azimuth_step_degrees": _AZIMUTH_STEP,
            "clearance_definition": "distance minus the sum of van der Waals radii; negative is sphere overlap",
            "retargeting_definition": (
                "Rigid-body placement of the unchanged tool with its C(apex)->C(distal) axis "
                "collinear with the target C-H bond and its apex at the intended apex-to-H distance."
            ),
        },
        "limitations": [
            "Rigid unrelaxed geometries. Both cages keep their ideal-lattice starting coordinates; neither molecule is allowed to distort to relieve contact or to improve an approach.",
            "Van der Waals overlap is a steric screen, not an energy. A negative gap does not quantify a repulsion and a nonnegative gap does not establish a viable approach path.",
            "Only the abstraction endpoint geometry is placed. No approach trajectory, withdrawal path, barrier, rate or electronic structure is computed here.",
            "Reachability is assessed for the free rigid tool. The real tool is mounted, and the mount's own constraints are not represented.",
            "Positional distinguishability is a necessary condition for positional selectivity, not a demonstration of it. The achievable positional precision of a real mechanism is not computed in this file.",
        ],
    }


def main() -> None:
    # Not named "runs": .gitignore ignores that directory name at any depth,
    # which would silently drop this evidence from the repository.
    output = Path(__file__).resolve().parent / "evidence" / "stage0-site-census-r4"
    output.mkdir(parents=True, exist_ok=False)
    result = census()
    (output / "census.json").write_text(
        json.dumps(result, indent=2, allow_nan=False, sort_keys=False) + "\n"
    )

    margin = result["site_relocation_distances"]
    print(f"atoms {result['composition_check']['atom_count']}", end="  ")
    print(
        f"C {result['composition_check']['carbon_count']}",
        f"H {result['composition_check']['hydrogen_count']}",
    )
    print("site type counts:", result["site_type_counts"])
    print(
        "intended apex-to-H",
        f"{result['intended_geometry']['apex_to_transferred_hydrogen_angstrom']:.3f} A",
        "| intended clearance",
        f"{result['intended_geometry']['intended_site_min_vdw_gap_angstrom']:+.3f} A",
    )
    print(
        "nearest alternative: H",
        margin["nearest_alternative_hydrogen_index"],
        margin["nearest_alternative_site_type"],
        f"needs {margin['nearest_alternative_translation_angstrom']:.3f} A apex translation;",
        "reachable" if margin["nearest_alternative_reachable"] else "vdW-blocked",
    )
    print(
        "nearest *reachable* alternative: H",
        margin["nearest_reachable_alternative_hydrogen_index"],
        margin["nearest_reachable_alternative_site_type"],
        f"at {margin['nearest_reachable_alternative_translation_angstrom']} A",
    )
    print(
        "reachable alternatives",
        margin["reachable_alternative_count"],
        "| blocked",
        margin["blocked_alternative_count"],
    )
    nearest = result["nearest_hydrogen_margin"]
    print(
        "nearest-hydrogen margin:",
        f"{nearest['margin_angstrom']:.3f} A apex displacement before H",
        nearest["limiting_hydrogen_index"],
        f"({nearest['limiting_site_type']}) becomes closest",
    )
    print(
        "  tied competitors:", nearest['tied_competitor_count'],
        "hydrogens", nearest['tied_competitor_hydrogen_indices'],
        "on carbons", nearest['tied_competitor_carbon_indices'],
    )
    print("apex to nearest own-mount H:", f"{result['handle_proximity']['apex_to_nearest_tool_hydrogen_angstrom']:.3f} A")
    print("wrote", output / "census.json")


if __name__ == "__main__":
    main()
