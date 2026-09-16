"""Screen the post-abstraction state for the tool welding itself to the workpiece.

A failure mode that falls between this project's lanes, because it is not a
selectivity question.  Follow the intended operation to its end: after the
hydrogen transfers, the tool has become a closed-shell alkyne,
H-C(triple)C-adamantyl, and the workpiece has become a 1-adamantyl radical, and
the two are left a few Angstrom apart.

A carbon radical beside an alkyne is not a stable arrangement.  Radical addition
across a C(triple)C forms a C-C sigma bond worth roughly 85 to 90 kcal/mol while
demoting a pi bond worth roughly 60 to 65, so it is exothermic by something like
20 to 25 kcal/mol with only a small barrier for alkyl radicals.  If that applies
here then the intended product is metastable, the deep well is "tool covalently
bonded to workpiece", and every successful abstraction is followed by a race
between withdrawing the tool and destroying it.  Unlike mis-targeting, this
failure yields no further product at all.

**What this file does, and does not do.**  It is a geometric screen only.  It
measures whether the post-abstraction geometry even permits an addition-competent
approach, and how far the radical would have to travel to reach one.  It computes
no energy and no barrier, so it cannot say the addition happens or does not; it
can only say whether geometry obstructs it, and by how much.  The thermodynamic
estimate above is bond-enthalpy reasoning, not a calculation from this
repository.

**The specific geometric question.**  Radical addition to an alkyne requires
attack on the pi system, which lies perpendicular to the C(triple)C axis.  An
approach along the axis points at the in-line position between the sp lobes and
is close to the worst geometry for addition.  So the relevant coordinate is the
*attack angle* between the axis and the vector from the alkyne carbon to the
radical: 90 degrees is addition-competent, 0 or 180 degrees is end-on and
obstructed.  That makes the tolerance here ANGULAR, which is a different
quantity from the lateral distance tolerance the positional analysis computes.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from nanodesign.candidates import make_h_abstraction

# A representative carbon-carbon distance at a radical-addition transition
# state, Angstrom.  Used only to define "addition-competent" reproducibly; it is
# a literature-scale number, not a computed saddle for this system.
_ADDITION_APPROACH = 2.2
_ADDITION_APPROACH_SOURCE = (
    "Typical C...C separation at an alkyl-radical-plus-alkyne addition transition "
    "state, about 2.2 Angstrom. Order-of-magnitude scale only."
)
_VDW = {"C": 1.70, "H": 1.20}


def screen(separation_angstrom: float = 3.6, lateral_offset_angstrom: float = 0.0) -> dict:
    _reactant, product, metadata = make_h_abstraction(separation_angstrom, lateral_offset_angstrom)
    symbols = product.get_chemical_symbols()
    positions = np.asarray(product.positions, dtype=float)

    radical = int(metadata["product_radical_index"])      # cage carbon that lost its H
    apex = int(metadata["tip_apex"])                      # alkyne carbon nearer the workpiece
    distal = int(metadata["tip_distal"])                  # the other alkyne carbon
    transferred = int(metadata["transferred_hydrogen"])   # now bonded to the apex

    axis = positions[distal] - positions[apex]
    axis /= np.linalg.norm(axis)

    results = {}
    for label, carbon in (("apex_carbon", apex), ("distal_carbon", distal)):
        to_radical = positions[radical] - positions[carbon]
        distance = float(np.linalg.norm(to_radical))
        direction = to_radical / distance
        # Angle between the C#C axis and the direction to the radical. Measured
        # from the apex the axis points away from the workpiece, so an end-on
        # approach gives ~180 degrees.
        cosine = float(np.clip(np.dot(direction, axis), -1.0, 1.0))
        attack_angle = math.degrees(math.acos(cosine))
        # How far off the axis the radical sits, and how far along it.
        along = float(np.dot(to_radical, axis))
        perpendicular = float(np.linalg.norm(to_radical - along * axis))
        # Deviation from the addition-competent 90 degrees.
        results[label] = {
            "radical_to_carbon_angstrom": distance,
            "attack_angle_degrees": attack_angle,
            "degrees_from_addition_competent": abs(attack_angle - 90.0),
            "perpendicular_offset_from_axis_angstrom": perpendicular,
            "displacement_along_axis_angstrom": along,
            "addition_competent_target_position": (
                positions[carbon] + _ADDITION_APPROACH * np.array([1.0, 0.0, 0.0])
            ).tolist(),
        }
        # Distance the radical would have to travel to sit perpendicular to the
        # axis at an addition transition-state separation. Minimised over the
        # azimuth about the axis, which is free by symmetry of the screen.
        target_perpendicular = _ADDITION_APPROACH
        travel = math.hypot(along - 0.0, perpendicular - target_perpendicular)
        results[label]["travel_to_addition_geometry_angstrom"] = abs(travel)

    # Is the newly transferred hydrogen physically between them?
    radical_to_transferred = float(np.linalg.norm(positions[radical] - positions[transferred]))
    radical_to_apex = results["apex_carbon"]["radical_to_carbon_angstrom"]
    obstructing = bool(radical_to_transferred < radical_to_apex)
    vdw_gap_to_transferred = radical_to_transferred - (_VDW["C"] + _VDW["H"])

    nearest = min(results, key=lambda key: results[key]["travel_to_addition_geometry_angstrom"])
    return {
        "screen": "post_abstraction_tool_welding_geometry",
        "computation": "rigid-geometry measurement of the product endpoint; no electronic structure",
        "state": "product endpoint from nanodesign.candidates.make_h_abstraction",
        "pose": {
            "separation_angstrom": float(separation_angstrom),
            "lateral_offset_angstrom": float(lateral_offset_angstrom),
        },
        "species_after_transfer": {
            "workpiece": "1-adamantyl radical, unpaired electron on cage carbon "
                         f"{radical}",
            "tool": f"closed-shell alkyne, H on carbon {apex}, C#C between {apex} and {distal}",
        },
        "alkyne_carbons": results,
        "transferred_hydrogen_obstructs": obstructing,
        "radical_to_transferred_hydrogen_angstrom": radical_to_transferred,
        "radical_to_transferred_hydrogen_vdw_gap_angstrom": vdw_gap_to_transferred,
        "closest_addition_target": nearest,
        "minimum_travel_to_addition_geometry_angstrom": results[nearest][
            "travel_to_addition_geometry_angstrom"
        ],
        "geometric_verdict": (
            "obstructed: the radical sits end-on to the C#C axis, which is the wrong "
            "approach for pi addition, and must travel a substantial distance to reach an "
            "addition-competent geometry"
            if results[nearest]["degrees_from_addition_competent"] > 45.0
            else "not obstructed: the radical is already near an addition-competent approach angle"
        ),
        "method_notes": {
            "addition_approach_angstrom": _ADDITION_APPROACH,
            "addition_approach_source": _ADDITION_APPROACH_SOURCE,
            "attack_angle_definition": (
                "Angle between the C#C axis, directed from the named carbon toward the other "
                "alkyne carbon, and the vector from that carbon to the radical centre. "
                "90 degrees is addition-competent; 0 or 180 is end-on."
            ),
            "vdw_radii": _VDW,
        },
        "limitations": [
            "Geometric only. No addition energy, no barrier and no rate is computed, so this cannot establish that welding does or does not occur; it establishes only whether the geometry obstructs it.",
            "The 20 to 25 kcal/mol exothermicity motivating this screen is bond-enthalpy reasoning, not a calculation from this repository. A small-model check is the cheap next step.",
            "Rigid unrelaxed product endpoint. Relaxation of the product, and the real mount's freedom, could change the approach geometry substantially.",
            "The obstruction found here depends on the approach staying collinear, which is exactly what thermal angular wander degrades. The relevant tolerance is angular and has not been computed anywhere in this project.",
            "Only two addition sites are considered, the two alkyne carbons. Addition elsewhere on the tool, and intramolecular chemistry after tool regeneration, are not screened.",
        ],
    }


def main() -> None:
    output = Path(__file__).resolve().parent / "evidence" / "product-state-welding-screen"
    output.mkdir(parents=True, exist_ok=False)
    result = screen()
    (output / "screen.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

    print("post-abstraction welding screen, nominal pose")
    for label, entry in result["alkyne_carbons"].items():
        print(f"  {label}:")
        print(f"    radical-to-carbon        {entry['radical_to_carbon_angstrom']:.3f} A")
        print(f"    attack angle             {entry['attack_angle_degrees']:.1f} deg "
              f"({entry['degrees_from_addition_competent']:.1f} deg from addition-competent)")
        print(f"    perpendicular offset     {entry['perpendicular_offset_from_axis_angstrom']:.3f} A")
        print(f"    travel to addition geom  {entry['travel_to_addition_geometry_angstrom']:.3f} A")
    print(f"  transferred H obstructs: {result['transferred_hydrogen_obstructs']} "
          f"(radical-to-H {result['radical_to_transferred_hydrogen_angstrom']:.3f} A, "
          f"vdW gap {result['radical_to_transferred_hydrogen_vdw_gap_angstrom']:+.3f} A)")
    print(f"  verdict: {result['geometric_verdict']}")
    print("wrote", output / "screen.json")


if __name__ == "__main__":
    main()
