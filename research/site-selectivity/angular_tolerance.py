"""Screen 3: can thermal wander reach the welding geometry, and can retraction?

Screens 1 and 2 established the shape of the problem.  Screen 1 found that the
post-abstraction geometry obstructs the addition that would weld the tool to the
workpiece: the radical sits exactly end-on to the C(triple)C axis at 180 degrees
where 90 is addition-competent, and the transferred hydrogen occupies the
approach vector at a negative van der Waals gap.  Screen 2 measured the
thermodynamics and made the problem worse -- the welding addition is 41.4
kcal/mol downhill at PBE0-D3(BJ)/def2-SVP against 38.4 for the abstraction it
follows, so **the intended product sits above a well deeper than the reaction
that produced it.**

Neither screen says whether the protection survives.  Screen 1 evaluated one
frozen pose, and the two protections both come from collinearity, which thermal
motion degrades.  This file asks the question that decides it: how far is the
addition geometry, measured in units of the tool's own thermal spread?

**It uses both stiffness components, which is the point.**  Positioning the tip
laterally and driving it axially are governed by different components of the same
Hessian -- 7.27 N/m for the doubly degenerate transverse bend and 251.03 N/m for
the axial, stretch-like coordinate, a ratio of 34.5.  Reaching the addition
geometry requires moving *closer* along the stiff axis and rotating through the
soft one, so the two must be treated separately.  Using a single stiffness would
be wrong by a factor of 34 in whichever direction it was taken.

**Angular spread.**  A lateral displacement of the apex at lever arm L tilts the
tool axis by approximately sigma_lateral / L for small angles, which gives the
angular wander directly from the measured bend stiffness and the recorded lever.

**What this cannot do.**  It is geometry and equilibrium spread, with no barrier
and no rate, so it cannot say welding does not occur -- only whether *thermal
wander at the product pose* is a plausible route to it, and whether a retraction
path passes through a competent geometry.  The addition barrier remains unknown,
and a driven retraction is not an equilibrium process.  The 2.2 Angstrom approach
distance is a literature-scale figure carried over from screen 1.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.constants import Boltzmann, electron_volt

from nanodesign.candidates import make_h_abstraction

EV_PER_ANGSTROM2_TO_N_PER_M = electron_volt / 1e-20

# Measured on the propyne model, evidence/tip-stiffness-propyne/.
LATERAL_STIFFNESS_N_PER_M = 7.2715
AXIAL_STIFFNESS_N_PER_M = 251.03
LEVER_ARM_ANGSTROM = 3.0593

# Addition-competence window. The distance is the literature-scale approach from
# screen 1; the angular window is deliberately generous, because a wide window
# makes the "cannot be reached" conclusion harder rather than easier to obtain.
ADDITION_APPROACH_ANGSTROM = 2.2
ADDITION_ANGLE_DEGREES = 90.0
ADDITION_ANGLE_WINDOW_DEGREES = 30.0
ADDITION_DISTANCE_WINDOW_ANGSTROM = 0.5


def _sigma(stiffness_n_per_m: float, temperature_kelvin: float) -> float:
    """Classical rms displacement, Angstrom. Both modes here are far below the
    338 cm^-1 quantum crossover, so the classical form is adequate."""
    stiffness_ev = stiffness_n_per_m / EV_PER_ANGSTROM2_TO_N_PER_M
    return math.sqrt((Boltzmann * temperature_kelvin / electron_volt) / stiffness_ev)


def thermal_spreads(temperature_kelvin: float = 300.0) -> dict:
    lateral = _sigma(LATERAL_STIFFNESS_N_PER_M, temperature_kelvin)
    axial = _sigma(AXIAL_STIFFNESS_N_PER_M, temperature_kelvin)
    angular = math.degrees(lateral / LEVER_ARM_ANGSTROM)
    return {
        "temperature_kelvin": temperature_kelvin,
        "lateral_stiffness_n_per_m": LATERAL_STIFFNESS_N_PER_M,
        "axial_stiffness_n_per_m": AXIAL_STIFFNESS_N_PER_M,
        "stiffness_ratio_axial_over_lateral": AXIAL_STIFFNESS_N_PER_M / LATERAL_STIFFNESS_N_PER_M,
        "lever_arm_angstrom": LEVER_ARM_ANGSTROM,
        "sigma_lateral_angstrom": lateral,
        "sigma_axial_angstrom": axial,
        "sigma_angular_degrees": angular,
        "angular_note": (
            "Small-angle estimate: a lateral apex displacement of sigma at lever arm L tilts the "
            "tool axis by sigma/L radians. Valid because sigma/L is about 0.08."
        ),
    }


def _attack_angle(radical, carbon, axis) -> float:
    direction = radical - carbon
    norm = float(np.linalg.norm(direction))
    cosine = float(np.clip(np.dot(direction / norm, axis), -1.0, 1.0))
    return math.degrees(math.acos(cosine))


def distance_to_welding_geometry(temperature_kelvin: float = 300.0) -> dict:
    """How far the product pose is from addition competence, in units of spread."""
    _reactant, product, metadata = make_h_abstraction()
    positions = np.asarray(product.positions, dtype=float)
    radical = positions[int(metadata["product_radical_index"])]
    apex = positions[int(metadata["tip_apex"])]
    distal = positions[int(metadata["tip_distal"])]
    axis = (distal - apex) / np.linalg.norm(distal - apex)

    separation = float(np.linalg.norm(radical - apex))
    angle = _attack_angle(radical, apex, axis)
    spreads = thermal_spreads(temperature_kelvin)

    # Reaching competence needs BOTH: close the axial gap, and rotate.
    axial_gap = separation - ADDITION_APPROACH_ANGSTROM
    angular_gap = abs(angle - ADDITION_ANGLE_DEGREES) - ADDITION_ANGLE_WINDOW_DEGREES
    return {
        "product_pose": {
            "radical_to_apex_angstrom": separation,
            "attack_angle_degrees": angle,
        },
        "spreads": spreads,
        "axial_approach_required_angstrom": axial_gap,
        "axial_approach_in_sigma": axial_gap / spreads["sigma_axial_angstrom"],
        "rotation_required_degrees": max(angular_gap, 0.0),
        "rotation_in_sigma": max(angular_gap, 0.0) / spreads["sigma_angular_degrees"],
        "both_required": True,
        "note": (
            "Addition competence requires the radical to be both near the approach distance and "
            "near perpendicular. These are independent coordinates governed by different "
            "stiffnesses, so the events must both occur; quoting either alone understates the "
            "protection."
        ),
    }


def retraction_scan(tilt_degrees_list=(0.0, 15.0, 30.0, 45.0, 60.0, 90.0),
                    separations=None) -> dict:
    """Does any retraction path pass through an addition-competent geometry?

    The tool is translated outward along its axis from the product separation,
    optionally tilted about the apex first.  A purely axial retraction should keep
    the radical end-on and only increase the distance, which is the benign case;
    the question is whether a tilted path passes nearer to perpendicular while
    still close enough to react.
    """
    _reactant, product, metadata = make_h_abstraction()
    positions = np.asarray(product.positions, dtype=float)
    radical = positions[int(metadata["product_radical_index"])]
    apex0 = positions[int(metadata["tip_apex"])]
    distal0 = positions[int(metadata["tip_distal"])]
    axis0 = (distal0 - apex0) / np.linalg.norm(distal0 - apex0)

    if separations is None:
        # Includes NEGATIVE offsets, i.e. approach as well as retraction. Scanning
        # only outward would make the screen one-sided: the dangerous direction is
        # approach, and omitting it would obtain a benign result by construction.
        separations = np.arange(-1.8, 4.01, 0.2)

    paths = []
    for tilt in tilt_degrees_list:
        angle = math.radians(float(tilt))
        # Rotate the tool axis about x, which tips it away from the radical.
        rotation = np.array([
            [1.0, 0.0, 0.0],
            [0.0, math.cos(angle), -math.sin(angle)],
            [0.0, math.sin(angle), math.cos(angle)],
        ])
        axis = rotation @ axis0
        competent = []
        closest = None
        for offset in separations:
            apex = apex0 + axis * float(offset)
            separation = float(np.linalg.norm(radical - apex))
            attack = _attack_angle(radical, apex, axis)
            distance_ok = abs(separation - ADDITION_APPROACH_ANGSTROM) <= ADDITION_DISTANCE_WINDOW_ANGSTROM
            angle_ok = abs(attack - ADDITION_ANGLE_DEGREES) <= ADDITION_ANGLE_WINDOW_DEGREES
            record = {
                "retraction_angstrom": float(offset),
                "radical_to_apex_angstrom": separation,
                "attack_angle_degrees": attack,
                "distance_competent": bool(distance_ok),
                "angle_competent": bool(angle_ok),
                "addition_competent": bool(distance_ok and angle_ok),
            }
            if record["addition_competent"]:
                competent.append(record)
            if closest is None or abs(attack - ADDITION_ANGLE_DEGREES) < abs(
                closest["attack_angle_degrees"] - ADDITION_ANGLE_DEGREES
            ):
                closest = record
        paths.append({
            "tilt_degrees": float(tilt),
            "competent_points": len(competent),
            "any_competent": bool(competent),
            "closest_approach_to_perpendicular": closest,
        })

    dangerous = [path for path in paths if path["any_competent"]]
    return {
        "paths": paths,
        "any_path_reaches_competence": bool(dangerous),
        "dangerous_tilts_degrees": [path["tilt_degrees"] for path in dangerous],
        "window": {
            "approach_angstrom": ADDITION_APPROACH_ANGSTROM,
            "distance_window_angstrom": ADDITION_DISTANCE_WINDOW_ANGSTROM,
            "angle_degrees": ADDITION_ANGLE_DEGREES,
            "angle_window_degrees": ADDITION_ANGLE_WINDOW_DEGREES,
        },
        "note": (
            "Rigid-body translation of an unrelaxed tool along its own axis, with the workpiece "
            "fixed. A real retraction relaxes both fragments and need not follow the axis."
        ),
    }


def joint_thermal_cost(temperature_kelvin: float = 300.0,
                       tilt_step_degrees: float = 2.5,
                       offset_step_angstrom: float = 0.05) -> dict:
    """Cheapest thermal route to addition competence, over tilt AND approach jointly.

    The one-sided version of this screen scanned only outward and concluded no
    path reaches competence.  That was an artifact of the scan: including approach
    shows the welding geometry IS reachable in principle, at roughly 40 degrees of
    tilt combined with 1.4 Angstrom of approach.  So the protection is not absolute
    geometric exclusion, and saying so would have been wrong.

    What it is instead: the competent region sits far away in the *joint*
    coordinate.  Tilt and axial approach are independent and governed by different
    stiffnesses, so the combined cost is the quadrature sum of the two in units of
    their own spreads.  That is the honest statement of the protection -- thermal
    and mechanical rather than forbidden.
    """
    spreads = thermal_spreads(temperature_kelvin)
    _reactant, product, metadata = make_h_abstraction()
    positions = np.asarray(product.positions, dtype=float)
    radical = positions[int(metadata["product_radical_index"])]
    apex0 = positions[int(metadata["tip_apex"])]
    distal0 = positions[int(metadata["tip_distal"])]
    axis0 = (distal0 - apex0) / np.linalg.norm(distal0 - apex0)

    competent = []
    for tilt in np.arange(0.0, 90.01, tilt_step_degrees):
        angle = math.radians(float(tilt))
        rotation = np.array([
            [1.0, 0.0, 0.0],
            [0.0, math.cos(angle), -math.sin(angle)],
            [0.0, math.sin(angle), math.cos(angle)],
        ])
        axis = rotation @ axis0
        for offset in np.arange(-1.8, 0.61, offset_step_angstrom):
            apex = apex0 + axis * float(offset)
            separation = float(np.linalg.norm(radical - apex))
            attack = _attack_angle(radical, apex, axis)
            if (abs(separation - ADDITION_APPROACH_ANGSTROM) <= ADDITION_DISTANCE_WINDOW_ANGSTROM
                    and abs(attack - ADDITION_ANGLE_DEGREES) <= ADDITION_ANGLE_WINDOW_DEGREES):
                cost = math.hypot(
                    float(tilt) / spreads["sigma_angular_degrees"],
                    abs(float(offset)) / spreads["sigma_axial_angstrom"],
                )
                competent.append({
                    "tilt_degrees": float(tilt),
                    "axial_offset_angstrom": float(offset),
                    "radical_to_apex_angstrom": separation,
                    "attack_angle_degrees": attack,
                    "tilt_in_sigma": float(tilt) / spreads["sigma_angular_degrees"],
                    "approach_in_sigma": abs(float(offset)) / spreads["sigma_axial_angstrom"],
                    "joint_cost_sigma": cost,
                })
    competent.sort(key=lambda item: item["joint_cost_sigma"])
    return {
        "spreads": spreads,
        "competent_point_count": len(competent),
        "reachable_in_principle": bool(competent),
        "cheapest_route": competent[0] if competent else None,
        "cheapest_five": competent[:5],
        "conclusion": (
            "The welding geometry is reachable in principle, so the protection is not geometric "
            "exclusion. It sits "
            + (f"{competent[0]['joint_cost_sigma']:.1f}" if competent else "an unreachable number of")
            + " combined standard deviations from the product pose, in a joint coordinate of tool "
            "tilt and axial approach, so thermal access at 300 K is negligible. That is a "
            "mechanical protection, contingent on the measured stiffnesses, not an absolute one."
        ),
    }


def main() -> None:
    output = Path(__file__).resolve().parent / "evidence" / "angular-tolerance-screen-r3"
    output.mkdir(parents=True, exist_ok=False)

    distance = distance_to_welding_geometry(300.0)
    scan = retraction_scan()
    joint = joint_thermal_cost(300.0)
    report = {
        "screen": "3_angular_tolerance_and_retraction_path",
        "question": "Can thermal wander or a retraction path reach the welding geometry?",
        "computation": "rigid geometry plus classical equilibrium spreads from a measured Hessian",
        "stiffness_source": "research/site-selectivity/evidence/tip-stiffness-propyne/",
        "thermodynamic_context": {
            "welding_addition_kcal_per_mol": -41.42,
            "intended_abstraction_kcal_per_mol": -38.38,
            "welding_deeper_by_kcal_per_mol": 3.04,
            "level": "PBE0-D3(BJ)/def2-SVP, electronic only, methyl model for the radical",
            "source": "data/validation/si-energy-reproduction/welding-thermo-screen.json",
            "caveat": (
                "Quote the internal comparison, not the absolute. PBE0-D3 is measurably too "
                "exothermic on this chemistry, by 2.01 kcal/mol for primary and 6.15 for tertiary "
                "abstraction in this project's own measurements, and the addition carries a "
                "similar bias. The comparison survives because both sides share it."
            ),
        },
        "thermal_reach": distance,
        "retraction": scan,
        "joint_thermal_cost": joint,
        "limitations": [
            "Geometry and equilibrium spread only. No addition barrier and no rate, so this cannot establish that welding does not occur.",
            "Equilibrium statistics applied to a driven retraction, which is not an equilibrium process.",
            "Rigid unrelaxed fragments throughout; a real retraction relaxes both and need not follow the tool axis.",
            "Stiffnesses come from the propyne proxy: a methyl mount is floppier than the real cage, so the spreads here are overestimates and the protection is understated, while rigid anchors push the other way.",
            "The 2.2 Angstrom approach distance and the 30 degree angular window are literature-scale choices, not computed from an addition saddle. The window was chosen generously so that a negative result is harder to obtain.",
        ],
    }
    (output / "screen.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    spreads = distance["spreads"]
    print("thermal spreads at 300 K, from the measured Hessian")
    print(f"  lateral  sigma {spreads['sigma_lateral_angstrom']:.4f} A   (k = {spreads['lateral_stiffness_n_per_m']:.2f} N/m)")
    print(f"  axial    sigma {spreads['sigma_axial_angstrom']:.4f} A   (k = {spreads['axial_stiffness_n_per_m']:.2f} N/m)")
    print(f"  angular  sigma {spreads['sigma_angular_degrees']:.2f} deg")
    print(f"  stiffness ratio axial/lateral {spreads['stiffness_ratio_axial_over_lateral']:.1f}x")
    print()
    print("distance from the product pose to addition competence")
    print(f"  axial approach required   {distance['axial_approach_required_angstrom']:.3f} A "
          f"= {distance['axial_approach_in_sigma']:.1f} sigma_axial")
    print(f"  rotation required         {distance['rotation_required_degrees']:.1f} deg "
          f"= {distance['rotation_in_sigma']:.1f} sigma_angular")
    print("  and BOTH are required simultaneously")
    print()
    print("axial paths, negative offset is approach and positive is retraction")
    for path in scan["paths"]:
        closest = path["closest_approach_to_perpendicular"]
        print(f"  tilt {path['tilt_degrees']:5.1f} deg: competent points {path['competent_points']}, "
              f"closest angle {closest['attack_angle_degrees']:6.1f} deg at "
              f"{closest['radical_to_apex_angstrom']:.2f} A")
    print(f"  any path reaches competence: {scan['any_path_reaches_competence']}")
    print()
    print("joint cost over tilt AND approach together")
    cheapest = joint["cheapest_route"]
    print(f"  competent grid points: {joint['competent_point_count']}")
    if cheapest:
        print(f"  cheapest route: tilt {cheapest['tilt_degrees']:.1f} deg "
              f"({cheapest['tilt_in_sigma']:.1f} sigma) plus approach "
              f"{abs(cheapest['axial_offset_angstrom']):.2f} A ({cheapest['approach_in_sigma']:.1f} sigma)")
        print(f"  JOINT COST: {cheapest['joint_cost_sigma']:.1f} combined sigma")
    print("wrote", output / "screen.json")


if __name__ == "__main__":
    main()
