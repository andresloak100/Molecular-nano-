"""Does thermal motion let the tip hit the wrong hydrogen?

This is the question the whole repository is built on top of and has not asked.
Everything so far measures *energies*: barriers, reaction energies, site
preferences. But a peer lane established that this tool has no steric
discrimination between adamantane's sixteen C-H sites, and this project's own
calibration puts the thermodynamic preference between site types at a few
kcal/mol at most. If the chemistry does not choose the site, then position
does, and position is not exact: the tip is a mechanical object at finite
temperature and it wanders.

Drexler's criterion for atomically precise manufacturing is an error rate per
operation small enough that a large product is still perfect - conventionally
quoted around 1e-15. That is a statement about a probability distribution of
tip positions, not about a barrier. This module computes it.

THE GEOMETRY. At the nominal 3.6 A pose the target bridgehead hydrogen sits on
the tool axis. The margin is the smallest displacement of the apex, IN ANY
DIRECTION, that puts a different hydrogen nearer than the target: geometrically,
the distance to the closest perpendicular-bisector plane between the target and
a competitor. That is 2.495 A, and the easiest escape is not lateral but tilted
about 120 degrees from the tool axis, down toward the equatorial methylenes.

Taking the purely lateral value instead gives 2.878 A. Both are correct
measures of different things, and only the minimum is safe for a tolerance,
since thermal displacement is three-dimensional and finds the easiest direction.
The lateral figure would overstate the allowance by about 15%.

THE PHYSICS. Displacement is three-dimensional, so the escape region is a
sphere. For an isotropic Gaussian with per-axis sigma,

    P(r > R) = erfc(x) + (2x / sqrt(pi)) exp(-x^2),   x = R / (sigma sqrt 2)

The sphere sits inside the target's Voronoi cell, so leaving the cell requires
at least leaving the sphere: this overestimates the error, which is the
direction an error budget should lean. Setting it below 1e-15 at R = 2.495 A
needs sigma < 0.292 A, which converts through classical equipartition,
sigma^2 = kT/k, into a *required stiffness* of 0.30 eV/A^2 (4.8 N/m) at 298 K.
That number is the design specification this project has never written down.

A CONDITION THAT CARRIES THE WHOLE RESULT, and is not decorative. This is a
criterion about which hydrogen is NEAREST, not which one REACTS. The two
coincide only if the competing barriers are comparable, which is A1's open
question. Below the tunneling crossover temperature they decouple further,
because a more distant site with a narrower barrier can win on width alone. So
the finding is "thermal wander is not the binding risk, GIVEN that nearest
implies reacting" - never "positional control is solved."

TWO CONTRIBUTIONS TO THE STIFFNESS, and only one is measured here.

The *interaction* stiffness is the curvature of the potential energy surface as
the tool slides laterally over the substrate: does the chemistry itself pull the
tip back over the target? This module measures it by scanning the lateral offset
the candidate builder already exposes, and reading the curvature at the minimum.

The *mount* stiffness is what the handle and its anchors contribute, and it is
not measured here. It needs a Hessian, which is expensive and is named as the
refinement rather than guessed at. The two add, so the interaction stiffness
alone is a lower bound on the total, and the sigma computed from it alone is an
upper bound.

WHY COOLING DOES NOT SOLVE IT. Classical equipartition says sigma^2 = kT/k, so
sigma goes to zero as the temperature does. That is wrong at low temperature.
The quantum result,

    sigma^2 = (hbar / (2 mu omega)) * coth(hbar omega / (2 kT))

approaches hbar/(2 mu omega) rather than zero: zero-point motion sets a floor
that no amount of cooling removes. Both are computed, and the floor is reported,
because "run it cold" is the obvious first suggestion and it has a limit.

WHAT THIS IS NOT. A rigid scan, so no relaxation: real curvatures are softer,
which makes the measured interaction stiffness an upper bound on itself. A
harmonic fit near one minimum. A model that asks only whether the tip is *over*
the right hydrogen, not whether being there makes the reaction go - that is a
barrier question and belongs to S1. And an error rate computed this way counts
mis-positioning only; it is a necessary condition for precision, not a
sufficient one.

Run from the repository root:

    python research/candidate-feasibility/positional_control.py --handle hydrogen
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

from nanodesign.candidates import make_h_abstraction  # noqa: E402
from nanodesign.quantum import PySCFCalculator, QuantumSettings  # noqa: E402
from reduced_models import HANDLES, TIP_APEX, reduced_candidate  # noqa: E402
from timing import load_snapshot, timed  # noqa: E402

EVIDENCE = LANE / "evidence"

BOLTZMANN_EV_PER_K = 8.617333262e-5
HBAR_J_S = 1.054571817e-34
EV_TO_J = 1.602176634e-19
ANGSTROM_TO_M = 1e-10
AMU_TO_KG = 1.66053906660e-27
DREXLER_ERROR_TARGET = 1e-15
TEMPERATURES_K = (4.0, 77.0, 195.0, 298.15, 500.0)


def competing_site_geometry(separation: float = 3.6) -> dict:
    """How far the tip may wander before a wrong hydrogen is the nearest one.

    A wrong turn worth recording, because it is the obvious way to do this and
    it is wrong by a factor of four. Ranking competitors by their lateral
    distance from the tool axis picks out three hydrogens 1.452 A off-axis and
    concludes the tip has only a 0.73 A margin. Those three sit at z = -3.66,
    on the far side of the cage, 7.4 A from the apex: the tip cannot reach them
    at any lateral offset. Lateral distance alone ignores whether a site is
    reachable in z.

    The criterion used instead is Voronoi from the apex. Slide the apex
    laterally toward a competitor and find the offset at which that competitor
    becomes equidistant from the apex with the target. The smallest such offset
    over all competitors is the margin. This correctly identifies the six
    equatorial methylene hydrogens, at z = -0.15, as the real rivals.

    Still only a geometric proxy. "Nearest hydrogen" is not "the hydrogen that
    reacts": abstraction also wants a roughly collinear C-C...H-C approach, and
    a tip that is nearest to the target but badly angled may do nothing at all.
    This bounds mis-targeting, it does not predict the product.
    """
    reactant, _, metadata = make_h_abstraction(separation, 0.0)
    positions = reactant.positions
    symbols = reactant.get_chemical_symbols()
    target = metadata["transferred_hydrogen"]
    target_position = positions[target]
    apex = positions[metadata["tip_apex"]]

    competitors = []
    for index in metadata["substrate_indices"]:
        if symbols[index] != "H" or index == target:
            continue
        position = positions[index]
        lateral_vector = position[:2] - apex[:2]
        lateral_distance = float(np.linalg.norm(lateral_vector))
        entry = {
            "index": index,
            "position": position.tolist(),
            "lateral_from_tool_axis_angstrom": float(np.linalg.norm(position[:2])),
            "z_angstrom": float(position[2]),
            "distance_to_apex_angstrom": float(np.linalg.norm(position - apex)),
        }
        # Crossover along an arbitrary unit direction d: with apex(s) = apex + s*d
        # the s^2 terms cancel between the two squared distances, leaving
        #     s = (|apex-C|^2 - |apex-T|^2) / (2 d . (C - T))
        numerator = float(
            np.dot(apex - position, apex - position)
            - np.dot(apex - target_position, apex - target_position)
        )
        separation_vector = position - target_position

        # Lateral-only, kept for comparison with the earlier reading.
        if lateral_distance >= 1e-9:
            lateral_direction = np.array([lateral_vector[0], lateral_vector[1], 0.0]) / lateral_distance
            denominator = 2.0 * float(np.dot(lateral_direction, separation_vector))
            lateral_crossover = numerator / denominator if abs(denominator) > 1e-12 else None
            if lateral_crossover is not None and lateral_crossover <= 0:
                lateral_crossover = None
        else:
            lateral_crossover = None
        entry["lateral_crossover_angstrom"] = lateral_crossover

        # Minimum over ALL displacement directions. s(d) is minimised by aligning
        # d with (C - T) when the numerator is positive, giving
        #     s_min = |apex-C|^2 - |apex-T|^2) / (2 |C - T|)
        # which is the distance from the apex to the perpendicular bisector plane
        # between target and competitor. The apex must cross that plane for the
        # competitor to become nearest, in any direction whatsoever.
        separation_norm = float(np.linalg.norm(separation_vector))
        entry["any_direction_crossover_angstrom"] = (
            numerator / (2.0 * separation_norm) if separation_norm > 1e-12 and numerator > 0 else None
        )
        competitors.append(entry)

    any_direction = [c for c in competitors if c["any_direction_crossover_angstrom"] is not None]
    lateral_only = [c for c in competitors if c["lateral_crossover_angstrom"] is not None]
    if not any_direction:
        raise RuntimeError("No competing hydrogen becomes nearest under any displacement.")
    any_direction.sort(key=lambda entry: entry["any_direction_crossover_angstrom"])
    margin = any_direction[0]["any_direction_crossover_angstrom"]
    lateral_margin = (
        min(c["lateral_crossover_angstrom"] for c in lateral_only) if lateral_only else None
    )

    naive = min(c["lateral_from_tool_axis_angstrom"] for c in competitors) / 2.0
    return {
        "target_hydrogen_index": target,
        "target_is_on_tool_axis": bool(np.linalg.norm(target_position[:2]) < 1e-9),
        "n_competing_substrate_hydrogens": len(competitors),
        "criterion": (
            "Distance from the apex to the nearest perpendicular-bisector plane "
            "between the target hydrogen and any competitor: the smallest "
            "displacement IN ANY DIRECTION that makes a wrong hydrogen nearest."
        ),
        "target_radius_angstrom": margin,
        "nearest_rival_index": any_direction[0]["index"],
        "nearest_rival_z_angstrom": any_direction[0]["z_angstrom"],
        "n_rivals_at_that_offset": sum(
            1 for c in any_direction if abs(c["any_direction_crossover_angstrom"] - margin) < 1e-6
        ),
        "lateral_only_margin_angstrom": lateral_margin,
        "why_the_minimum_and_not_the_lateral_value": (
            "Thermal displacement is three-dimensional, so the tolerance is set "
            "by the easiest escape direction, not by the lateral one. The "
            "lateral figure is larger and would overstate the allowance by "
            "roughly 15% here; the minimum is the conservative quantity and is "
            "what the specification uses."
        ),
        "rejected_naive_radius_angstrom": naive,
        "why_naive_is_wrong": (
            "Half the smallest lateral offset of any competitor. Picks hydrogens "
            "on the far side of the cage that the tip cannot reach at all, and "
            "understates the margin roughly fourfold."
        ),
        "limits": (
            "Nearest-hydrogen is not the same as the hydrogen that reacts. "
            "Abstraction also wants a near-collinear approach, so this bounds "
            "mis-targeting rather than predicting a product."
        ),
        "reachable_rivals": any_direction[:6],
    }


def error_probability(target_radius: float, sigma: float) -> float:
    """Probability that a 3-D isotropic Gaussian displacement leaves the sphere.

    Thermal displacement of the tip is three-dimensional, so the escape region
    is a sphere of radius ``target_radius``, not a disc. For per-axis sigma,

        P(r > R) = erfc(x) + (2x / sqrt(pi)) exp(-x^2),   x = R / (sigma sqrt 2)

    The sphere fits inside the target's Voronoi cell, so leaving the cell
    requires at least leaving the sphere: this is an upper bound on the true
    mis-targeting probability, which is the direction an error budget wants.
    """
    if sigma <= 0:
        return 0.0
    x = target_radius / (sigma * math.sqrt(2.0))
    if x > 30:  # exp(-900) underflows; the probability is zero to any precision
        return 0.0
    return math.erfc(x) + (2.0 * x / math.sqrt(math.pi)) * math.exp(-x * x)


def required_sigma(target_radius: float, error_target: float = DREXLER_ERROR_TARGET) -> float:
    """Largest sigma meeting the error target. Inverted numerically by bisection.

    The 3-D tail has no elementary inverse, so this brackets and bisects rather
    than using a closed form that would only be right in 2-D.
    """
    low, high = 1e-6, target_radius
    # error_probability increases with sigma, so bracket then halve.
    for _ in range(200):
        middle = 0.5 * (low + high)
        if error_probability(target_radius, middle) > error_target:
            high = middle
        else:
            low = middle
    return 0.5 * (low + high)


def classical_sigma(stiffness_ev_per_a2: float, temperature_k: float) -> float:
    """sigma = sqrt(kT/k). One harmonic coordinate, classical equipartition.

    Amplitude scales as 1/sqrt(stiffness), not equally across coordinates.
    """
    if stiffness_ev_per_a2 <= 0:
        return float("inf")
    return math.sqrt(BOLTZMANN_EV_PER_K * temperature_k / stiffness_ev_per_a2)


def quantum_sigma(stiffness_ev_per_a2: float, reduced_mass_amu: float, temperature_k: float) -> dict:
    """Harmonic-oscillator sigma including zero-point motion.

    sigma^2 = (hbar / (2 mu omega)) coth(hbar omega / (2 kT)), which tends to
    hbar/(2 mu omega) as T -> 0 rather than to zero.
    """
    if stiffness_ev_per_a2 <= 0:
        return {"sigma_angstrom": float("inf"), "zero_point_sigma_angstrom": float("inf")}
    k_si = stiffness_ev_per_a2 * EV_TO_J / (ANGSTROM_TO_M ** 2)
    mu_si = reduced_mass_amu * AMU_TO_KG
    omega = math.sqrt(k_si / mu_si)
    zero_point_variance = HBAR_J_S / (2.0 * mu_si * omega)
    x = HBAR_J_S * omega / (2.0 * BOLTZMANN_EV_PER_K * EV_TO_J * temperature_k)
    # coth overflows for large x (cold, stiff); it is 1 to machine precision there.
    coth = 1.0 / math.tanh(x) if x < 350 else 1.0
    variance = zero_point_variance * coth
    return {
        "sigma_angstrom": math.sqrt(variance) / ANGSTROM_TO_M,
        "zero_point_sigma_angstrom": math.sqrt(zero_point_variance) / ANGSTROM_TO_M,
        "omega_rad_per_s": omega,
        "hbar_omega_ev": HBAR_J_S * omega / EV_TO_J,
        "kt_ev": BOLTZMANN_EV_PER_K * temperature_k,
        "quantum_regime": bool(HBAR_J_S * omega / EV_TO_J > BOLTZMANN_EV_PER_K * temperature_k),
    }


# Buildable mount stiffnesses, from the rung-5 lane's bracket. Load geometry,
# not material, spans the range: a strut loaded along its axis is two orders of
# magnitude stiffer than the same material worked in bending, because a
# cantilever softens as the cube of its length.
MOUNT_BRACKET_N_PER_M = (
    ("bending cantilever, soft end", 2.0),
    ("bending cantilever, stiff end", 20.0),
    ("nm-scale axial strut, low", 130.0),
    ("nm-scale axial strut, high", 400.0),
    ("single C-C bond, axial (hard cap)", 450.0),
)
NEWTON_PER_METRE_PER_EV_PER_A2 = 16.02176634


def mount_feasibility(target_radius: float, temperatures=(77.0, 298.15)) -> dict:
    """Does a *buildable* mount meet the positional requirement? Not automatic.

    An earlier version of this lane reported that the positional requirement,
    4.8 N/m at 298 K, is so soft that any mount clears it. That was checked
    against nothing: it compared the requirement to an imagined stiffness. The
    rung-5 lane's bracket of what is actually buildable starts at about 2 N/m
    for a long handle worked in bending, which is BELOW the requirement.

    At 2 N/m and 298 K the mis-targeting rate is about 1e-6 - nine orders of
    magnitude short of the 1e-15 target. The same mount at 77 K gives 3e-25 and
    passes comfortably. So the design constraint is real and has a shape:

        mount stiffly (axial, short) OR operate cold.

    A long cantilever handle at room temperature does not satisfy positional
    control, and nothing else in this repository would have caught that,
    because the requirement and the buildable range were never compared.
    """
    rows = {}
    for name, newtons in MOUNT_BRACKET_N_PER_M:
        stiffness = newtons / NEWTON_PER_METRE_PER_EV_PER_A2
        per_temperature = {}
        for temperature in temperatures:
            sigma = classical_sigma(stiffness, temperature)
            probability = error_probability(target_radius, sigma)
            per_temperature[f"{temperature:g}K"] = {
                "sigma_angstrom": sigma,
                "error_probability": probability,
                "meets_target": probability < DREXLER_ERROR_TARGET,
            }
        rows[name] = {
            "stiffness_n_per_m": newtons,
            "stiffness_ev_per_angstrom_squared": stiffness,
            "temperatures": per_temperature,
        }
    requirements = {
        f"{temperature:g}K": {
            "required_stiffness_n_per_m": (
                BOLTZMANN_EV_PER_K * temperature / required_sigma(target_radius) ** 2
                * NEWTON_PER_METRE_PER_EV_PER_A2
            ),
        }
        for temperature in temperatures
    }
    failures = [
        name for name, row in rows.items()
        if not row["temperatures"]["298.15K"]["meets_target"]
    ]
    return {
        "bracket_source": "rung-5 mechanical-coupling lane; load geometry spans the range, not material",
        "requirements": requirements,
        "mounts": rows,
        "fails_at_298K": failures,
        "verdict": (
            "Positional control is NOT automatically satisfied. Mounts at the "
            "soft end of the buildable range fail at room temperature and pass "
            "cold, so the constraint is: mount stiffly, or operate cold."
        ),
        "condition": (
            "Still conditional on nearest-hydrogen implying reacting-hydrogen, "
            "which is unestablished."
        ),
    }


def lateral_scan(handle: str, offsets, settings: QuantumSettings) -> list[dict]:
    """Energy against lateral tool offset, rigid: no relaxation at any offset."""
    points = []
    for offset in offsets:
        atoms, _ = reduced_candidate(handle, separation=3.6, lateral=float(offset))
        atoms.calc = PySCFCalculator(settings, event_log=EVIDENCE / f"lateral-{handle}.jsonl")
        point: dict = {"lateral_offset_angstrom": float(offset), "n_atoms": len(atoms)}
        try:
            with timed(point, "timing"):
                point["energy_ev"] = float(atoms.get_potential_energy())
            diagnostics = atoms.calc.diagnostics
            point.update(
                scf_converged=diagnostics["scf_converged"],
                s2=diagnostics["s2"],
                basis_functions=diagnostics["basis_functions"],
            )
            print(f"  offset {offset:5.2f} A   E = {point['energy_ev']:.6f} eV"
                  f"   cpu {point['timing']['cpu_seconds']:.1f} s")
        except Exception as exc:  # noqa: BLE001 - a failed point is a result
            point["error"] = f"{type(exc).__name__}: {exc}"
            print(f"  offset {offset:5.2f} A   FAILED: {point['error'][:70]}")
        points.append(point)
    return points


def fit_curvature(points: list[dict], window: float = 0.45) -> dict:
    """Quadratic fit near the minimum; curvature is the lateral stiffness.

    Only points within ``window`` of the minimum are fitted, because a harmonic
    reading of a landscape is only meaningful near its bottom.
    """
    usable = [p for p in points if "energy_ev" in p]
    if len(usable) < 3:
        return {"fitted": False, "reason": f"only {len(usable)} points converged; need 3"}
    minimum = min(usable, key=lambda p: p["energy_ev"])
    centre = minimum["lateral_offset_angstrom"]
    local = [p for p in usable if abs(p["lateral_offset_angstrom"] - centre) <= window]
    if len(local) < 3:
        local = sorted(usable, key=lambda p: abs(p["lateral_offset_angstrom"] - centre))[:3]
    x = np.array([p["lateral_offset_angstrom"] for p in local])
    y = np.array([p["energy_ev"] for p in local])
    coefficients = np.polyfit(x, y, 2)
    curvature = 2.0 * float(coefficients[0])  # E = a x^2 + ... -> k = d2E/dx2 = 2a
    residuals = y - np.polyval(coefficients, x)
    return {
        "fitted": True,
        "stiffness_ev_per_angstrom_squared": curvature,
        "minimum_at_offset_angstrom": centre,
        "points_used": len(local),
        "fit_window_angstrom": window,
        "max_fit_residual_ev": float(np.abs(residuals).max()),
        "energy_range_over_scan_ev": float(max(p["energy_ev"] for p in usable)
                                           - min(p["energy_ev"] for p in usable)),
        "curvature_sign_note": (
            "A negative or near-zero curvature means the interaction does not "
            "localize the tip laterally at all, so every bit of positional "
            "control would have to come from the mount."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handle", default="hydrogen", choices=list(HANDLES))
    parser.add_argument("--max-offset", type=float, default=1.2)
    parser.add_argument("--step", type=float, default=0.2)
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    EVIDENCE.mkdir(parents=True, exist_ok=True)
    out = EVIDENCE / (args.out or f"positional-control-{args.handle}.json")
    if out.exists():
        raise SystemExit(f"{out} exists; evidence is never overwritten. Choose --out.")

    geometry = competing_site_geometry()
    radius = geometry["target_radius_angstrom"]
    needed = required_sigma(radius)

    settings = QuantumSettings(
        xc="pbe0", basis=args.basis, dispersion="d3bj",
        density_fit=True, threads=1, scf_initial_guess="atom",
    )

    report: dict = {
        "schema_version": 1,
        "question": "Does thermal motion of the tip let it reach the wrong hydrogen?",
        "why": (
            "This tool has no steric site discrimination and only a few kcal/mol "
            "of thermodynamic preference, so selectivity rests on position. "
            "Position at finite temperature is a distribution, not a point."
        ),
        "geometry": geometry,
        "specification": {
            "error_target_per_operation": DREXLER_ERROR_TARGET,
            "target_radius_angstrom": radius,
            "required_lateral_sigma_angstrom": needed,
            "derivation": "P(r > R) = exp(-R^2 / 2 sigma^2) for an isotropic 2-D Gaussian",
            "required_stiffness_ev_per_angstrom_squared_at_298K": (
                BOLTZMANN_EV_PER_K * 298.15 / needed ** 2
            ),
        },
        "method": f"PBE0-D3(BJ)/{args.basis}, density fitting, rigid lateral scan, no relaxation",
        "model": f"{args.handle}-handle reduced candidate",
        "measures": "interaction contribution to lateral stiffness only",
        "does_not_measure": (
            "the mount's own stiffness, which adds to it and needs a Hessian; "
            "so sigma computed here is an UPPER bound and the verdict it gives "
            "is pessimistic by an unknown amount"
        ),
        "host_at_start": load_snapshot(),
        "status": "running",
    }
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    print(f"target radius {radius:.3f} A "
          f"({geometry['n_competitors_at_that_distance']} competing H at "
          f"{geometry['nearest_competitor_lateral_angstrom']:.3f} A lateral)")
    print(f"required sigma for {DREXLER_ERROR_TARGET:.0e} error rate: {needed:.4f} A")
    print(f"required stiffness at 298 K: "
          f"{report['specification']['required_stiffness_ev_per_angstrom_squared_at_298K']:.2f} eV/A^2\n")

    offsets = np.arange(0.0, args.max_offset + 1e-9, args.step)
    print(f"lateral scan, {args.handle} handle, {len(offsets)} points:")
    points = lateral_scan(args.handle, offsets, settings)
    report["scan"] = points
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    fit = fit_curvature(points)
    report["fit"] = fit
    if fit["fitted"]:
        stiffness = fit["stiffness_ev_per_angstrom_squared"]
        tool_mass = sum(
            mass for index, mass in enumerate(reduced_candidate(args.handle)[0].get_masses())
            if index >= TIP_APEX
        )
        report["reduced_mass_amu"] = float(tool_mass)
        report["reduced_mass_note"] = (
            "Mass of the displaced tool atoms. The substrate is anchored and "
            "treated as infinitely massive, which understates sigma slightly."
        )
        rows = {}
        for temperature in TEMPERATURES_K:
            classical = classical_sigma(stiffness, temperature)
            quantum = quantum_sigma(stiffness, tool_mass, temperature)
            rows[f"{temperature:g}K"] = {
                "temperature_k": temperature,
                "classical_sigma_angstrom": classical,
                "classical_error_probability": error_probability(radius, classical),
                "quantum_sigma_angstrom": quantum["sigma_angstrom"],
                "quantum_error_probability": error_probability(radius, quantum["sigma_angstrom"]),
                "zero_point_sigma_angstrom": quantum["zero_point_sigma_angstrom"],
                "hbar_omega_ev": quantum.get("hbar_omega_ev"),
                "kt_ev": quantum.get("kt_ev"),
                "quantum_regime": quantum.get("quantum_regime"),
                "meets_target_classical": error_probability(radius, classical) < DREXLER_ERROR_TARGET,
                "meets_target_quantum": error_probability(radius, quantum["sigma_angstrom"]) < DREXLER_ERROR_TARGET,
            }
        report["thermal"] = rows
        zero_point = quantum_sigma(stiffness, tool_mass, 1.0)["zero_point_sigma_angstrom"]
        report["zero_point_floor"] = {
            "sigma_angstrom": zero_point,
            "error_probability_at_floor": error_probability(radius, zero_point),
            "coolable_to_target": error_probability(radius, zero_point) < DREXLER_ERROR_TARGET,
            "meaning": (
                "Cooling reduces sigma only until zero-point motion dominates. "
                "If the error probability at this floor already exceeds the "
                "target, no temperature fixes this stiffness."
            ),
        }

        print(f"\nfitted interaction stiffness: {stiffness:.3f} eV/A^2 "
              f"(minimum at {fit['minimum_at_offset_angstrom']:.2f} A, "
              f"scan spans {fit['energy_range_over_scan_ev']:.4f} eV)")
        print(f"reduced mass {tool_mass:.1f} amu, zero-point sigma {zero_point:.4f} A\n")
        print(f"{'T (K)':>8} {'sigma_cl':>10} {'sigma_qm':>10} {'P(wrong site)':>16} {'meets 1e-15':>12}")
        for name, row in rows.items():
            print(f"{row['temperature_k']:8.2f} {row['classical_sigma_angstrom']:10.4f} "
                  f"{row['quantum_sigma_angstrom']:10.4f} {row['quantum_error_probability']:16.3e} "
                  f"{str(row['meets_target_quantum']):>12}")
    else:
        print(f"\nno curvature fit: {fit['reason']}")

    report["status"] = "completed"
    report["host_at_end"] = load_snapshot()
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
