"""Local curvature checks of a supplied geometry on its attached potential.

This module does not find a saddle or establish reaction connectivity.  It
finite-differences actual calculator forces and keeps the coordinate space
explicit, including all free Cartesian modes.  No toy potential is provided.
"""

from __future__ import annotations

from copy import deepcopy
import math
from numbers import Integral, Real

import numpy as np
from ase import Atoms
from ase.calculators.calculator import compare_atoms
from ase.constraints import FixAtoms
from scipy.constants import atomic_mass, electron_volt, speed_of_light


# sqrt(eV / Angstrom**2 / atomic_mass) is an angular frequency in s**-1.
# Divide by 2*pi and by the speed of light in cm/s for a wavenumber.
_WAVENUMBER_FACTOR = math.sqrt(electron_volt / (1e-20 * atomic_mass)) / (
    2 * math.pi * speed_of_light * 100
)
_ABSENT = object()
# A coordinate-representation guard only, not an electronic-force or chemical
# accuracy target. The actual offset must represent the requested stencil to
# within one part per million before dividing force differences by that step.
_MAX_STEP_REPRESENTATION_RELATIVE_ERROR = 1e-6


def _positive_real(name, value):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _displacement_plan(atoms, coordinate_indices, step):
    """Reject unresolved stencils before touching the attached calculator."""
    probe = atoms.copy()
    reference = atoms.positions
    plan = []
    for flattened_coordinate in coordinate_indices:
        atom_index, axis = divmod(int(flattened_coordinate), 3)
        for sign in (1, -1):
            requested = sign * step
            with np.errstate(over="ignore", invalid="ignore"):
                displaced = reference[atom_index, axis] + requested
                actual = displaced - reference[atom_index, axis]
                relative_error = abs(actual - requested) / step
            if (not np.isfinite(displaced) or not np.isfinite(relative_error)
                    or sign * actual <= 0
                    or relative_error > _MAX_STEP_REPRESENTATION_RELATIVE_ERROR):
                raise ValueError(
                    "step_angstrom is not numerically resolved at "
                    f"atom {atom_index}, axis {axis}; choose a representable step "
                    "or recenter the geometry. This is a coordinate precision guard."
                )
            probe.positions[:] = reference
            probe.positions[atom_index, axis] = displaced
            if "positions" not in compare_atoms(atoms, probe):
                raise ValueError(
                    "step_angstrom is indistinguishable from the reference "
                    "geometry under ASE's position-cache tolerance."
                )
            plan.append({
                "atom_index": atom_index,
                "axis": axis,
                "requested_offset_angstrom": requested,
                "actual_offset_angstrom": float(actual),
                "displaced_coordinate_angstrom": float(displaced),
                "relative_step_representation_error": float(relative_error),
            })
    return plan


def characterize_stationary_point(
    atoms: Atoms,
    *,
    step_angstrom: float = 0.005,
    force_tolerance_ev_per_angstrom: float = 0.03,
    frequency_tolerance_cm1: float = 20.0,
    max_free_coordinates: int = 120,
) -> dict:
    """Characterize the local Cartesian Hessian without changing ``atoms``.

    The caller attaches the desired ASE calculator, including any electronic
    evaluation log destination.  This function requests 1 + 2*d force arrays,
    where d is the number of free Cartesian coordinates; an existing baseline
    cache may avoid one calculation.  Only whole-atom ``FixAtoms`` constraints
    and finite, nonperiodic geometries are supported.  Frozen atoms are omitted
    from differentiation and from the gradient residual.

    All d modes are reported.  Negative signed frequencies denote imaginary
    harmonic frequencies; only frequencies below -frequency_tolerance_cm1
    count as resolved negative modes.  Near-zero modes remain unresolved,
    including numerical translations/rotations in an unanchored molecule.
    The function neither removes external modes nor certifies a minimum or a
    transition state from a mode count.

    Stencils indistinguishable under ASE's position cache, or distorted by
    coordinate rounding, are rejected before force work. Successful results
    include full baseline and displaced force arrays so that the unsymmetrized
    Hessian and its asymmetry can be independently reconstructed.

    Calculations use a geometry copy with the actual attached calculator.  Its
    ASE geometry/result cache and optional ``diagnostics`` are restored in a
    finally block, including on failure.  External logs, files, and backend
    counters intentionally retain evidence of the performed evaluations.  Do
    not use the same calculator concurrently from another thread.
    """
    step = _positive_real("step_angstrom", step_angstrom)
    force_tolerance = _positive_real("force_tolerance_ev_per_angstrom", force_tolerance_ev_per_angstrom)
    frequency_tolerance = _positive_real("frequency_tolerance_cm1", frequency_tolerance_cm1)
    if isinstance(max_free_coordinates, bool) or not isinstance(max_free_coordinates, Integral) or max_free_coordinates < 1:
        raise ValueError("max_free_coordinates must be a positive integer")
    if not isinstance(atoms, Atoms) or len(atoms) == 0:
        raise ValueError("A nonempty ASE Atoms geometry is required")
    if atoms.calc is None:
        raise ValueError("Attach an actual force calculator before characterizing the geometry")
    if np.any(atoms.pbc):
        raise ValueError("Periodic systems are unsupported; supply a finite nonperiodic geometry")
    if not np.all(np.isfinite(atoms.positions)):
        raise ValueError("Atomic coordinates must be finite")

    frozen = set()
    for constraint in atoms.constraints:
        if not isinstance(constraint, FixAtoms):
            raise ValueError("Only FixAtoms constraints are supported by this Cartesian Hessian")
        indices = np.asarray(constraint.get_indices(), dtype=int)
        if np.any(indices < 0) or np.any(indices >= len(atoms)):
            raise ValueError("FixAtoms contains an out-of-range atom index")
        frozen.update(indices.tolist())
    free_indices = np.array([index for index in range(len(atoms)) if index not in frozen], dtype=int)
    dimension = 3 * len(free_indices)
    if dimension == 0:
        raise ValueError("At least one atom must be free")
    if dimension > max_free_coordinates:
        raise ValueError(
            f"Hessian requires {dimension} free coordinates and {1 + 2 * dimension} force requests; "
            f"the configured max_free_coordinates is {max_free_coordinates}. "
            "Use a smaller physical model or explicitly raise the calculation limit."
        )
    free_masses = np.asarray(atoms.get_masses()[free_indices], dtype=float)
    if not np.all(np.isfinite(free_masses)) or np.any(free_masses <= 0):
        raise ValueError("Free-atom masses must be finite and positive")

    coordinate_indices = (3 * free_indices[:, None] + np.arange(3)).ravel()
    plan = _displacement_plan(atoms, coordinate_indices, step)

    calculator = atoms.calc
    # Copy mutable arrays before doing any calculator work.  Saving these ASE
    # fields avoids a displaced cache masquerading as the original geometry.
    saved_atoms = getattr(calculator, "atoms", _ABSENT)
    if saved_atoms is not _ABSENT and saved_atoms is not None:
        saved_atoms = saved_atoms.copy()
    saved_results = deepcopy(calculator.results)
    saved_diagnostics = getattr(calculator, "diagnostics", _ABSENT)
    if saved_diagnostics is not _ABSENT:
        saved_diagnostics = deepcopy(saved_diagnostics)
    working = atoms.copy()
    working.calc = calculator
    reference_positions = atoms.positions.copy()

    def get_full_forces():
        forces = np.asarray(working.get_forces(apply_constraint=False), dtype=float)
        if forces.shape != (len(atoms), 3) or not np.all(np.isfinite(forces)):
            raise ValueError("Calculator returned an invalid or nonfinite force array")
        return forces.copy()

    def calculation_call_id():
        diagnostics = getattr(calculator, "diagnostics", None)
        call_id = diagnostics.get("call_id") if isinstance(diagnostics, dict) else None
        return call_id if isinstance(call_id, str) else None

    displacements = []
    try:
        baseline_full_forces = get_full_forces()
        baseline_call_id = calculation_call_id()
        baseline_forces = baseline_full_forces[free_indices]
        raw_hessian = np.empty((dimension, dimension), dtype=float)
        for column in range(dimension):
            displaced_forces = []
            for displacement in plan[2 * column:2 * column + 2]:
                working.positions[:] = reference_positions
                working.positions[displacement["atom_index"], displacement["axis"]] = displacement["displaced_coordinate_angstrom"]
                # A custom calculator may use a coarser cache rule than ASE's
                # standard comparison. Never accept its unchanged force cache.
                if "forces" in calculator.results and "positions" not in calculator.check_state(working):
                    raise ValueError("step_angstrom is unresolved by the attached calculator's position cache")
                full_forces = get_full_forces()
                displacements.append({
                    **displacement,
                    "forces_ev_per_angstrom": full_forces.tolist(),
                    "calculation_call_id": calculation_call_id(),
                })
                displaced_forces.append(full_forces[free_indices].ravel())
            plus, minus = displaced_forces
            # F_i = -dE/dx_i, so H_ij = -dF_i/dx_j.
            raw_hessian[:, column] = -(plus - minus) / (2 * step)
    finally:
        if saved_atoms is _ABSENT:
            if hasattr(calculator, "atoms"):
                delattr(calculator, "atoms")
        else:
            calculator.atoms = saved_atoms
        calculator.results = saved_results
        if saved_diagnostics is _ABSENT:
            if hasattr(calculator, "diagnostics"):
                delattr(calculator, "diagnostics")
        else:
            calculator.diagnostics = saved_diagnostics

    if not np.all(np.isfinite(raw_hessian)):
        raise ValueError("Finite-difference Hessian contains nonfinite entries")
    asymmetry = raw_hessian - raw_hessian.T
    hessian = (raw_hessian + raw_hessian.T) / 2
    repeated_mass = np.repeat(free_masses, 3)
    inverse_sqrt_mass = 1 / np.sqrt(repeated_mass)
    weighted = hessian * inverse_sqrt_mass[:, None] * inverse_sqrt_mass[None, :]
    eigenvalues, eigenvectors = np.linalg.eigh(weighted)
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * _WAVENUMBER_FACTOR
    negative_count = int(np.count_nonzero(frequencies < -frequency_tolerance))
    unresolved_count = int(np.count_nonzero(np.abs(frequencies) <= frequency_tolerance))
    force_max = float(np.linalg.norm(baseline_forces, axis=1).max())
    stationary = force_max <= force_tolerance
    if not stationary:
        classification = "not_stationary_within_force_tolerance"
    elif negative_count == 0:
        classification = "stationary_point_with_no_resolved_negative_modes"
    elif negative_count == 1:
        classification = "first_order_saddle_candidate"
    else:
        classification = "higher_order_saddle_candidate"
    raw_norm = float(np.linalg.norm(raw_hessian))

    return {
        "method": "central finite differences of attached-calculator forces",
        "settings": {
            "step_angstrom": step,
            "force_tolerance_ev_per_angstrom": force_tolerance,
            "frequency_tolerance_cm1": frequency_tolerance,
            "max_free_coordinates": int(max_free_coordinates),
        },
        "displacement_resolution": {
            "checked": True,
            "cache_visibility_check": "ASE compare_atoms before force work, then attached-calculator check_state for each displacement",
            "max_relative_step_representation_error_allowed": _MAX_STEP_REPRESENTATION_RELATIVE_ERROR,
            "max_relative_step_representation_error_observed": max(d["relative_step_representation_error"] for d in plan),
            "scope": "Numerical coordinate/cache resolution only; not a force-error or chemical-accuracy bound.",
        },
        "finite_difference_evidence": {
            "units": {"length": "angstrom", "force": "eV/angstrom"},
            "reference_positions_angstrom": reference_positions.tolist(),
            "baseline_forces_ev_per_angstrom": baseline_full_forces.tolist(),
            "baseline_calculation_call_id": baseline_call_id,
            "calculation_call_id_scope": "Matches calculator diagnostics and electronic.jsonl where supplied; null means the calculator supplied no string call ID.",
            "force_scope": "Unconstrained forces on all atoms, including fixed anchors",
            "axis_convention": "0=x, 1=y, 2=z; atom_index is zero-based",
            "displacement_order": "Plus then minus for each free Cartesian coordinate in coordinate_order",
            "hessian_reconstruction": "Select free-atom force components; column j is -(F_plus - F_minus)/(2*settings.step_angstrom). Symmetrize only after reconstructing all columns.",
            "displacements": displacements,
        },
        "force_requests": 1 + 2 * dimension,
        "free_atom_indices": free_indices.tolist(),
        "frozen_atom_indices": sorted(frozen),
        "free_coordinate_count": dimension,
        "free_masses_amu": free_masses.tolist(),
        "coordinate_order": "x,y,z for each atom in free_atom_indices",
        "free_gradient_ev_per_angstrom": (-baseline_forces).tolist(),
        "free_force_max_ev_per_angstrom": force_max,
        "free_force_cartesian_rms_ev_per_angstrom": float(np.sqrt(np.mean(baseline_forces**2))),
        "stationary_within_force_tolerance": stationary,
        "classification": classification,
        "hessian_ev_per_angstrom2": hessian.tolist(),
        "hessian_asymmetry_max_abs_ev_per_angstrom2": float(np.abs(asymmetry).max()),
        "hessian_asymmetry_relative_frobenius": float(np.linalg.norm(asymmetry)) / raw_norm if raw_norm else 0.0,
        "hessian_symmetrized": True,
        "mass_weighted_eigenvalues_ev_per_angstrom2_amu": eigenvalues.tolist(),
        "frequencies_cm1": frequencies.tolist(),
        "frequency_sign_convention": "Negative signed wavenumbers denote imaginary harmonic frequencies.",
        "negative_mode_count": negative_count,
        "raw_negative_eigenvalue_count": int(np.count_nonzero(eigenvalues < 0)),
        "unresolved_near_zero_mode_count": unresolved_count,
        "resolved_positive_mode_count": dimension - negative_count - unresolved_count,
        "mode_classifications": [
            "resolved_negative" if frequency < -frequency_tolerance else
            "resolved_positive" if frequency > frequency_tolerance else "unresolved_near_zero"
            for frequency in frequencies
        ],
        "cartesian_modes_per_sqrt_amu": (eigenvectors.T * inverse_sqrt_mass[None, :]).reshape(
            dimension, len(free_indices), 3
        ).tolist(),
        "mode_normalization": "Each Cartesian mode obeys sum(mass_amu * displacement**2) = 1; sign is arbitrary.",
        "external_modes_removed": False,
        "coordinate_space": "free Cartesian coordinates with fixed anchors" if frozen else "all 3N Cartesian coordinates of an unanchored molecule",
        "irc_performed": False,
        "transition_state_verified": False,
        "finite_difference_step_convergence_checked": False,
        "limitations": [
            "One resolved imaginary mode at a stationary geometry is necessary but not sufficient for a first-order reaction transition state; IRC or equivalent endpoint connectivity is absent.",
            "Near-zero modes are unresolved, not verified stable vibrations. Unanchored molecules include numerical translations and rotations; no external-mode projection is performed.",
            "Fixed anchors define the curvature subspace; global translations and rotations must not be removed from an anchored calculation.",
            "Central differences have finite-step and electronic-force errors. Hessian symmetry is a diagnostic, not an accuracy guarantee; repeat with other steps and converged electronic settings.",
            "Local harmonic curvature of the attached potential does not establish free energies, competing reactions, assembly reliability, or model accuracy.",
        ],
    }
