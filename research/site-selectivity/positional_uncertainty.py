"""The other half of a positional-selectivity argument: how far the tip moves.

Stage 0 of this lane measured a *positional margin*: the apex of the tool would
have to be displaced 2.495 Angstrom before a different cage hydrogen became the
nearest one.  A margin is only half an argument.  The other half is the
positional *uncertainty* of a real mounted tip -- the thermal and zero-point
displacement it undergoes while held by a handle of finite stiffness.  This is
the quantity Drexler's feasibility case for positional assembly actually turns
on, and `docs/MODEL.md` states plainly that it is not represented anywhere in
this repository.

This module computes it from a Hessian that has already been calculated, so it
adds analysis rather than new electronic-structure cost.  It consumes the
dictionary returned by ``nanodesign.stationary.characterize_stationary_point``
read-only and does not modify any core module.

Physics, stated explicitly because each choice is a place to be wrong.

For a harmonic potential U = (1/2) x^T H x in Cartesian displacements, the
classical Boltzmann distribution gives a displacement covariance

    Cov_classical = k_B T H^-1                                            (1)

which is mass independent.  The *inverse* is essential and is the whole reason
a diagonal block of H is not good enough: inverting H lets every other atom
relax in response to a tip displacement, whereas the diagonal block holds them
clamped and therefore reports the tip as stiffer than it is.

Quantum mechanically the harmonic ground state has zero-point motion that does
not vanish at low temperature, and for a stiff diamondoid tip at 300 K the
stiff modes have hbar*omega >> k_B T, so zero-point motion *dominates* thermal
motion.  Using equipartition there underestimates the spread.  In mass-weighted
normal coordinates Q_k with angular frequencies omega_k,

    <Q_k^2> = (hbar / 2 omega_k) coth(hbar omega_k / 2 k_B T)             (2)

which tends to k_B T / omega_k^2 in the high-temperature limit, recovering (1).
The Cartesian covariance is then the mode sum

    Cov_ij = sum_k <Q_k^2> d_k,i d_k,j                                    (3)

over the mass-weighted-normalized Cartesian mode vectors d_k.  Both routes are
implemented and cross-checked against each other in the classical limit.

What this does NOT establish.  A harmonic equilibrium displacement distribution
is not an operating error rate for a driven assembly step: the tool is being
pushed along a reaction coordinate, not sampling equilibrium, and a reaction
occurs or does not for reasons that include barriers this repository has not
yet obtained.  Nothing here is converted into a success probability.  The
anchors are held rigid, which makes every stiffness reported here an upper
bound and every displacement a lower bound -- a real mount is softer, so the
true uncertainty is larger than these numbers.
"""

from __future__ import annotations

import math
from numbers import Real

import numpy as np
from scipy.constants import (
    Boltzmann,
    atomic_mass,
    electron_volt,
    hbar,
)

# 1 eV/Angstrom**2 expressed in N/m, the unit tip stiffness is quoted in.
EV_PER_ANGSTROM2_TO_N_PER_M = electron_volt / 1e-20
# Angular frequency in rad/s from a mass-weighted eigenvalue in eV/(A**2 amu).
_OMEGA_FROM_EIGENVALUE = math.sqrt(electron_volt / (1e-20 * atomic_mass))
# <Q**2> in SI (kg m**2) converted to amu*Angstrom**2.
_Q2_SI_TO_AMU_ANGSTROM2 = 1.0 / (atomic_mass * 1e-20)


class PositionalUncertaintyError(ValueError):
    """The supplied curvature data cannot define a displacement distribution."""


def _positive(name: str, value) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value <= 0:
        raise PositionalUncertaintyError(f"{name} must be finite and positive")
    return float(value)


def _kt_ev(temperature_kelvin: float) -> float:
    return Boltzmann * temperature_kelvin / electron_volt


def _validate(characterization: dict, *, require_stationary: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int]]:
    """Pull the Hessian, eigenvalues and mode vectors out, checking each."""
    if not isinstance(characterization, dict):
        raise PositionalUncertaintyError("characterization must be the dict from characterize_stationary_point")
    for key in ("hessian_ev_per_angstrom2", "mass_weighted_eigenvalues_ev_per_angstrom2_amu",
                "cartesian_modes_per_sqrt_amu", "free_atom_indices"):
        if key not in characterization:
            raise PositionalUncertaintyError(f"characterization is missing {key!r}")

    if require_stationary and not characterization.get("stationary_within_force_tolerance", False):
        raise PositionalUncertaintyError(
            "Geometry is not stationary within its force tolerance "
            f"(max free force {characterization.get('free_force_max_ev_per_angstrom')!r} eV/A). "
            "A harmonic displacement distribution is only defined about a minimum; "
            "relax first, or pass require_stationary=False and accept that the result "
            "is not an equilibrium distribution."
        )

    hessian = np.asarray(characterization["hessian_ev_per_angstrom2"], dtype=float)
    eigenvalues = np.asarray(characterization["mass_weighted_eigenvalues_ev_per_angstrom2_amu"], dtype=float)
    modes = np.asarray(characterization["cartesian_modes_per_sqrt_amu"], dtype=float)
    free = [int(index) for index in characterization["free_atom_indices"]]

    dimension = hessian.shape[0]
    if hessian.ndim != 2 or hessian.shape[1] != dimension:
        raise PositionalUncertaintyError("Hessian must be square")
    if not np.all(np.isfinite(hessian)):
        raise PositionalUncertaintyError("Hessian contains nonfinite entries")
    if eigenvalues.shape != (dimension,):
        raise PositionalUncertaintyError("Eigenvalue count does not match the Hessian dimension")
    if modes.shape != (dimension, len(free), 3):
        raise PositionalUncertaintyError("Mode array shape does not match the free-coordinate space")
    return hessian, eigenvalues, modes.reshape(dimension, dimension), free


def _check_positive_definite(eigenvalues: np.ndarray, floor: float) -> None:
    """Fail closed on any mode that cannot carry a bounded displacement."""
    smallest = float(eigenvalues.min())
    if smallest <= 0.0:
        raise PositionalUncertaintyError(
            f"The curvature has a nonpositive mode (smallest eigenvalue {smallest:.6g} "
            "eV/A^2/amu). A saddle or a flat direction has no bounded thermal "
            "displacement distribution; this is a transition structure or an "
            "unanchored molecule, not a mounted tip at a minimum."
        )
    if smallest < floor:
        raise PositionalUncertaintyError(
            f"The softest mode has eigenvalue {smallest:.6g} eV/A^2/amu, below the "
            f"floor of {floor:.6g}. Displacements along it would dominate and are "
            "almost certainly numerical translations or rotations of an unanchored "
            "molecule rather than a real restoring force. Anchor the structure, or "
            "use pair_distance_uncertainty, which is translation and rotation free."
        )


def _mode_amplitudes(eigenvalues: np.ndarray, temperature_kelvin: float, *, quantum: bool) -> np.ndarray:
    """<Q_k^2> per mode, in amu*Angstrom**2."""
    if quantum:
        omega = np.sqrt(eigenvalues) * _OMEGA_FROM_EIGENVALUE
        x = hbar * omega / (2.0 * Boltzmann * temperature_kelvin)
        # coth overflows for stiff modes at low T; 1/tanh is stable and tends to 1.
        amplitude_si = (hbar / (2.0 * omega)) / np.tanh(x)
        return amplitude_si * _Q2_SI_TO_AMU_ANGSTROM2
    return _kt_ev(temperature_kelvin) / eigenvalues


def covariance(characterization: dict, temperature_kelvin: float = 298.15, *,
               quantum: bool = True, require_stationary: bool = True,
               eigenvalue_floor: float = 1e-6) -> dict:
    """Cartesian displacement covariance of every free atom, in Angstrom**2.

    ``quantum=True`` uses equation (2) and therefore includes zero-point motion.
    ``quantum=False`` uses classical equipartition, equation (1), which
    underestimates the spread whenever hbar*omega is not small against k_B T.
    """
    temperature = _positive("temperature_kelvin", temperature_kelvin)
    hessian, eigenvalues, modes, free = _validate(characterization, require_stationary=require_stationary)
    _check_positive_definite(eigenvalues, eigenvalue_floor)

    amplitudes = _mode_amplitudes(eigenvalues, temperature, quantum=quantum)
    matrix = (modes.T * amplitudes) @ modes

    # Independent classical route, for cross-checking the mode sum against a
    # direct inversion of the Hessian.  These must agree when quantum=False.
    classical_direct = _kt_ev(temperature) * np.linalg.inv(hessian)
    if quantum:
        consistency = float(
            np.abs(((modes.T * (_kt_ev(temperature) / eigenvalues)) @ modes) - classical_direct).max()
        )
    else:
        consistency = float(np.abs(matrix - classical_direct).max())

    return {
        "temperature_kelvin": temperature,
        "statistics": "quantum harmonic, includes zero-point motion" if quantum else "classical equipartition",
        "covariance_angstrom2": matrix,
        "free_atom_indices": free,
        "classical_mode_sum_vs_hessian_inverse_max_abs_angstrom2": consistency,
        "softest_mode_eigenvalue_ev_per_angstrom2_amu": float(eigenvalues.min()),
    }


def atom_uncertainty(characterization: dict, atom_index: int, temperature_kelvin: float = 298.15,
                     **options) -> dict:
    """Displacement statistics of one atom, e.g. the tool apex.

    ``atom_index`` is an index into the original geometry, not into the free
    subset.  Returns the principal axes of that atom's 3x3 covariance block,
    the root-mean-square displacement along each, and the total.
    """
    result = covariance(characterization, temperature_kelvin, **options)
    free = result["free_atom_indices"]
    if atom_index not in free:
        raise PositionalUncertaintyError(
            f"Atom {atom_index} is frozen or absent; it has no displacement distribution. "
            f"Free atoms are {free}."
        )
    position = free.index(atom_index)
    block = result["covariance_angstrom2"][3 * position:3 * position + 3, 3 * position:3 * position + 3]
    variances, axes = np.linalg.eigh(block)
    if variances.min() <= 0:
        raise PositionalUncertaintyError("Atom covariance is not positive definite; check the Hessian")
    sigma = np.sqrt(variances)
    return {
        "atom_index": int(atom_index),
        "temperature_kelvin": result["temperature_kelvin"],
        "statistics": result["statistics"],
        "rms_displacement_angstrom": float(math.sqrt(float(np.trace(block)))),
        "principal_sigma_angstrom": sigma.tolist(),
        "largest_sigma_angstrom": float(sigma.max()),
        "principal_axes": axes.T.tolist(),
        "covariance_angstrom2": block.tolist(),
        "softest_direction_stiffness_n_per_m": float(
            EV_PER_ANGSTROM2_TO_N_PER_M / (_kt_ev(result["temperature_kelvin"]) / float(variances.max()))
        ) if not options.get("quantum", True) else None,
    }


def effective_stiffness(characterization: dict, atom_index: int, direction=None,
                        *, require_stationary: bool = True, eigenvalue_floor: float = 1e-6) -> dict:
    """Compliance-based effective stiffness of one atom, in N/m.

    Along a unit direction n the effective stiffness is 1 / (n^T C n) where C is
    the compliance, the inverse Hessian restricted to that atom.  This lets the
    rest of the structure relax, which is what a real handle does; the diagonal
    Hessian block, which clamps every other atom, is reported alongside so the
    difference between the two is visible rather than hidden.

    Quoted against Drexler's own scale: stiffnesses of order 10 to 100 N/m are
    what the positional-assembly argument assumes for a stiff diamondoid mount.
    """
    hessian, eigenvalues, _modes, free = _validate(characterization, require_stationary=require_stationary)
    _check_positive_definite(eigenvalues, eigenvalue_floor)
    if atom_index not in free:
        raise PositionalUncertaintyError(f"Atom {atom_index} is frozen or absent; free atoms are {free}")
    position = free.index(atom_index)
    slice_ = slice(3 * position, 3 * position + 3)

    compliance = np.linalg.inv(hessian)[slice_, slice_]
    clamped = hessian[slice_, slice_]

    eigen_compliance, axes = np.linalg.eigh(compliance)
    if eigen_compliance.min() <= 0:
        raise PositionalUncertaintyError("Atom compliance is not positive definite")
    # Largest compliance is the softest direction, so the smallest stiffness.
    stiffness_principal = np.sort(1.0 / eigen_compliance)[::-1]

    result = {
        "atom_index": int(atom_index),
        "principal_stiffness_n_per_m": (stiffness_principal * EV_PER_ANGSTROM2_TO_N_PER_M).tolist(),
        "softest_stiffness_n_per_m": float(stiffness_principal.min() * EV_PER_ANGSTROM2_TO_N_PER_M),
        "softest_direction": axes[:, int(np.argmax(eigen_compliance))].tolist(),
        "clamped_diagonal_block_stiffness_n_per_m": (
            np.sort(np.linalg.eigvalsh(clamped))[::-1] * EV_PER_ANGSTROM2_TO_N_PER_M
        ).tolist(),
        "note": (
            "Compliance-based values let the rest of the structure relax and are the "
            "physical ones. The clamped diagonal block is always stiffer and is reported "
            "only to show the size of that error."
        ),
    }
    if direction is not None:
        unit = np.asarray(direction, dtype=float)
        norm = float(np.linalg.norm(unit))
        if not math.isfinite(norm) or norm == 0:
            raise PositionalUncertaintyError("direction must be a finite nonzero vector")
        unit = unit / norm
        result["requested_direction"] = unit.tolist()
        result["requested_direction_stiffness_n_per_m"] = float(
            EV_PER_ANGSTROM2_TO_N_PER_M / float(unit @ compliance @ unit)
        )
    return result


def pair_distance_uncertainty(characterization: dict, first: int, second: int, positions,
                              temperature_kelvin: float = 298.15, *, quantum: bool = True,
                              require_stationary: bool = True,
                              eigenvalue_floor: float = 1e-8) -> dict:
    """RMS fluctuation of the distance between two atoms, in Angstrom.

    A difference coordinate is exactly invariant to rigid translation and, to
    first order, to rigid rotation.  So unlike an absolute displacement this is
    well defined for an unanchored molecule, and it is the natural measure for a
    transfer coordinate such as donor-to-hydrogen or hydrogen-to-apex.
    """
    temperature = _positive("temperature_kelvin", temperature_kelvin)
    _hessian, eigenvalues, modes, free = _validate(characterization, require_stationary=require_stationary)
    for index in (first, second):
        if index not in free:
            raise PositionalUncertaintyError(f"Atom {index} is frozen or absent; free atoms are {free}")
    if first == second:
        raise PositionalUncertaintyError("first and second must be different atoms")

    coordinates = np.asarray(positions, dtype=float)
    separation = coordinates[second] - coordinates[first]
    length = float(np.linalg.norm(separation))
    if not math.isfinite(length) or length == 0:
        raise PositionalUncertaintyError("The two atoms coincide; the distance has no direction")
    unit = separation / length

    # d|r2 - r1| = u . (dx2 - dx1), so the projection vector is +u on atom 2
    # and -u on atom 1.  Modes that translate the whole molecule give zero here.
    projection = np.zeros(len(free) * 3)
    projection[3 * free.index(second):3 * free.index(second) + 3] = unit
    projection[3 * free.index(first):3 * free.index(first) + 3] = -unit

    # A near-zero mode carries a near-divergent amplitude, so even a tiny
    # numerical overlap with it can contaminate the sum.  Rigid translations and
    # rotations have exactly zero overlap with a difference coordinate in the
    # ideal case, but that is a symmetry accident of the ideal case and is not
    # something to depend on.  Modes below the floor are excluded and the
    # variance they would have contributed is reported, so the exclusion is
    # auditable rather than invisible.
    usable = eigenvalues > max(eigenvalue_floor, 0.0)
    amplitudes = np.zeros_like(eigenvalues)
    amplitudes[usable] = _mode_amplitudes(eigenvalues[usable], temperature, quantum=quantum)
    overlaps = modes @ projection
    variance = float(np.sum(amplitudes[usable] * overlaps[usable] ** 2))

    excluded = ~usable
    excluded_variance = 0.0
    if np.any(excluded & (eigenvalues > 0)):
        soft = eigenvalues[excluded & (eigenvalues > 0)]
        soft_overlaps = overlaps[excluded & (eigenvalues > 0)]
        excluded_variance = float(
            np.sum(_mode_amplitudes(soft, temperature, quantum=quantum) * soft_overlaps**2)
        )
    if variance < 0:
        raise PositionalUncertaintyError("Computed a negative variance; check the curvature data")

    return {
        "atom_indices": [int(first), int(second)],
        "distance_angstrom": length,
        "temperature_kelvin": temperature,
        "statistics": "quantum harmonic, includes zero-point motion" if quantum else "classical equipartition",
        "sigma_angstrom": math.sqrt(variance),
        "variance_angstrom2": variance,
        "excluded_nonpositive_modes": int(np.count_nonzero(eigenvalues <= 0)),
        "excluded_soft_modes": int(np.count_nonzero(usable == False) - np.count_nonzero(eigenvalues <= 0)),
        "excluded_mode_variance_angstrom2": excluded_variance,
        "excluded_mode_variance_fraction": (
            excluded_variance / variance if variance > 0 else None
        ),
        "eigenvalue_floor_ev_per_angstrom2_amu": float(max(eigenvalue_floor, 0.0)),
        "note": (
            "A difference coordinate removes rigid translation exactly and rigid rotation "
            "to first order, so this is defined without anchors. Nonpositive and soft modes "
            "are excluded, counted, and the variance they would have contributed is reported, "
            "rather than silently dropped or silently included with a divergent amplitude."
        ),
    }


def compare_to_margin(sigma_angstrom: float, margin_angstrom: float) -> dict:
    """Express a measured positional margin in units of the positional spread.

    Deliberately reports a ratio and nothing else. It does **not** return an
    error rate or a success probability: the tip in an assembly step is driven
    along a reaction coordinate rather than sampling an equilibrium
    distribution, the harmonic model has no barrier in it, and the anchors here
    are rigid so the true spread is larger than the computed one. Converting
    this ratio into an operating reliability would be exactly the unsupported
    step the project's coordination notes warn against.
    """
    sigma = _positive("sigma_angstrom", sigma_angstrom)
    margin = _positive("margin_angstrom", margin_angstrom)
    return {
        "sigma_angstrom": sigma,
        "margin_angstrom": margin,
        "margin_in_sigma": margin / sigma,
        "interpretation": (
            "The margin is this many root-mean-square positional spreads wide, in a rigid-anchor "
            "harmonic equilibrium model. This is a length comparison only. It is not an error "
            "rate, not a success probability, and not a statement that the operation works."
        ),
        "known_omissions": [
            "Anchors are rigid, so the real mount is softer and the real spread larger.",
            "No barrier, reaction coordinate or driving force is represented.",
            "Harmonic only; anharmonicity grows exactly where displacements are large.",
            "Equilibrium statistics applied to a driven operation.",
        ],
    }
