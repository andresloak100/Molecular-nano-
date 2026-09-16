"""Validation of the positional-uncertainty layer against closed-form physics.

Every expected value here is an independent analytic result, not a previously
recorded output of this code, so these are checks rather than regressions.  The
canonical case is Drexler's own: a 10 N/m restraint at 300 K gives a 0.2035 A
root-mean-square displacement, which is the number the positional-assembly
argument rests on.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.constants import Boltzmann, atomic_mass, electron_volt, hbar

from positional_uncertainty import (
    EV_PER_ANGSTROM2_TO_N_PER_M,
    PositionalUncertaintyError,
    atom_uncertainty,
    compare_to_margin,
    covariance,
    effective_stiffness,
    pair_distance_uncertainty,
)


def characterization(hessian, masses_amu, *, stationary=True):
    """Build the same dict ``characterize_stationary_point`` returns.

    Constructed here the same way the core module does it, so the tests exercise
    the real data contract rather than a convenience shape.
    """
    hessian = np.asarray(hessian, dtype=float)
    masses = np.asarray(masses_amu, dtype=float)
    repeated = np.repeat(masses, 3)
    inverse_sqrt = 1.0 / np.sqrt(repeated)
    weighted = hessian * inverse_sqrt[:, None] * inverse_sqrt[None, :]
    eigenvalues, eigenvectors = np.linalg.eigh(weighted)
    dimension = hessian.shape[0]
    modes = (eigenvectors.T * inverse_sqrt[None, :]).reshape(dimension, len(masses), 3)
    return {
        "hessian_ev_per_angstrom2": hessian.tolist(),
        "mass_weighted_eigenvalues_ev_per_angstrom2_amu": eigenvalues.tolist(),
        "cartesian_modes_per_sqrt_amu": modes.tolist(),
        "free_atom_indices": list(range(len(masses))),
        "free_masses_amu": masses.tolist(),
        "stationary_within_force_tolerance": stationary,
        "free_force_max_ev_per_angstrom": 0.0 if stationary else 2.38,
    }


def isotropic(stiffness_n_per_m, mass_amu=12.011):
    """One free atom held isotropically at a given stiffness in N/m."""
    k = stiffness_n_per_m / EV_PER_ANGSTROM2_TO_N_PER_M
    return characterization(np.eye(3) * k, [mass_amu])


def test_unit_conversion_matches_si():
    # 1 eV / Angstrom**2 = 1.602176634e-19 J / 1e-20 m**2.
    assert EV_PER_ANGSTROM2_TO_N_PER_M == pytest.approx(16.02176634, rel=1e-12)


def test_drexler_canonical_amplitude():
    """10 N/m at 300 K must give sigma = sqrt(kT/k) = 0.2035 Angstrom."""
    expected = math.sqrt(Boltzmann * 300.0 / 10.0) * 1e10
    assert expected == pytest.approx(0.20352, abs=1e-5)

    result = atom_uncertainty(isotropic(10.0), 0, 300.0, quantum=False)
    for sigma in result["principal_sigma_angstrom"]:
        assert sigma == pytest.approx(expected, rel=1e-10)
    # Three independent directions, so the total RMS is sqrt(3) times one axis.
    assert result["rms_displacement_angstrom"] == pytest.approx(expected * math.sqrt(3.0), rel=1e-10)


def test_classical_amplitude_is_mass_independent():
    """Equipartition on the potential energy cannot depend on the mass."""
    light = atom_uncertainty(isotropic(25.0, mass_amu=1.008), 0, 298.15, quantum=False)
    heavy = atom_uncertainty(isotropic(25.0, mass_amu=200.0), 0, 298.15, quantum=False)
    assert light["largest_sigma_angstrom"] == pytest.approx(heavy["largest_sigma_angstrom"], rel=1e-12)


def test_classical_mode_sum_equals_hessian_inverse():
    """The two independent classical routes must agree; (3) must reduce to (1)."""
    generator = np.random.default_rng(20260916)
    root = generator.normal(size=(9, 9))
    hessian = root @ root.T + 6.0 * np.eye(9)
    data = characterization(hessian, [12.011, 1.008, 15.999])
    result = covariance(data, 298.15, quantum=False)
    assert result["classical_mode_sum_vs_hessian_inverse_max_abs_angstrom2"] < 1e-12


def test_quantum_reduces_to_classical_at_high_temperature():
    """coth(x) -> 1/x, so (2) -> kT/omega**2 when kT dominates hbar*omega."""
    data = isotropic(1.0, mass_amu=200.0)
    classical = atom_uncertainty(data, 0, 40000.0, quantum=False)["largest_sigma_angstrom"]
    quantum = atom_uncertainty(data, 0, 40000.0, quantum=True)["largest_sigma_angstrom"]
    assert quantum == pytest.approx(classical, rel=1e-3)
    assert quantum > classical  # zero-point motion always adds


def test_zero_point_motion_survives_zero_temperature():
    """A stiff mode is zero-point dominated; classical equipartition vanishes."""
    data = isotropic(300.0, mass_amu=1.008)
    cold_quantum = atom_uncertainty(data, 0, 1.0, quantum=True)["largest_sigma_angstrom"]
    cold_classical = atom_uncertainty(data, 0, 1.0, quantum=False)["largest_sigma_angstrom"]
    # At 1 K the ratio is sqrt(hbar*omega / 2kT) = sqrt(1617) ~ 40, not arbitrary.
    assert cold_quantum / cold_classical == pytest.approx(
        math.sqrt(hbar * math.sqrt(300.0 / (1.008 * atomic_mass)) / (2.0 * Boltzmann * 1.0)),
        rel=1e-6,
    )

    # Analytic ground-state amplitude of a one-dimensional oscillator,
    # sigma = sqrt(hbar / (2 m omega)), computed entirely in SI.
    mass = 1.008 * atomic_mass
    omega = math.sqrt(300.0 / mass)
    expected = math.sqrt(hbar / (2.0 * mass * omega)) * 1e10
    assert cold_quantum == pytest.approx(expected, rel=1e-6)


def test_stiff_mode_is_zero_point_dominated_at_room_temperature():
    """The claim in the module docstring, checked rather than asserted."""
    data = isotropic(300.0, mass_amu=1.008)
    quantum = atom_uncertainty(data, 0, 298.15, quantum=True)["largest_sigma_angstrom"]
    classical = atom_uncertainty(data, 0, 298.15, quantum=False)["largest_sigma_angstrom"]
    assert quantum > 2.0 * classical


def test_compliance_is_softer_than_the_clamped_block():
    """Letting the rest of the structure relax must never look stiffer."""
    # Two coupled atoms: a strong spring between them and a weak tether on each.
    coupling, tether = 40.0, 3.0
    block = np.zeros((6, 6))
    for axis in range(3):
        block[axis, axis] = coupling + tether
        block[3 + axis, 3 + axis] = coupling + tether
        block[axis, 3 + axis] = -coupling
        block[3 + axis, axis] = -coupling
    data = characterization(block, [12.011, 12.011])
    result = effective_stiffness(data, 0)
    compliance_based = min(result["principal_stiffness_n_per_m"])
    clamped = min(result["clamped_diagonal_block_stiffness_n_per_m"])
    assert compliance_based < clamped
    # Closed form for this topology. Per axis the 2x2 block is
    # [[c+t, -c], [-c, c+t]], whose determinant is t(t+2c), so the compliance
    # seen by one atom is (c+t)/(t(t+2c)) and the relaxed stiffness is its
    # reciprocal. Limits check out: c -> infinity gives 2t (rigidly bound pair
    # pulling on two parallel tethers) and c -> 0 gives t (isolated tether).
    relaxed = tether * (tether + 2.0 * coupling) / (coupling + tether)
    expected = relaxed * EV_PER_ANGSTROM2_TO_N_PER_M
    assert compliance_based == pytest.approx(expected, rel=1e-9)
    assert clamped == pytest.approx((coupling + tether) * EV_PER_ANGSTROM2_TO_N_PER_M, rel=1e-9)


def test_requested_direction_stiffness():
    data = isotropic(50.0)
    result = effective_stiffness(data, 0, direction=[1.0, 0.0, 0.0])
    assert result["requested_direction_stiffness_n_per_m"] == pytest.approx(50.0, rel=1e-9)


def test_fails_closed_on_a_saddle():
    """A negative curvature has no bounded thermal distribution."""
    hessian = np.diag([5.0, 5.0, -5.0])
    data = characterization(hessian, [12.011])
    with pytest.raises(PositionalUncertaintyError, match="nonpositive mode"):
        covariance(data, 298.15)


def test_fails_closed_on_a_soft_translation():
    """An unanchored molecule's near-zero modes must not be silently used."""
    hessian = np.diag([20.0, 20.0, 1e-12])
    data = characterization(hessian, [12.011])
    with pytest.raises(PositionalUncertaintyError, match="softest mode"):
        covariance(data, 298.15)


def test_fails_closed_on_a_nonstationary_geometry():
    data = isotropic(20.0)
    data["stationary_within_force_tolerance"] = False
    with pytest.raises(PositionalUncertaintyError, match="not stationary"):
        covariance(data, 298.15)
    # ...but an explicit override is available and honoured.
    assert covariance(data, 298.15, require_stationary=False)["temperature_kelvin"] == 298.15


def test_frozen_atom_has_no_distribution():
    data = isotropic(20.0)
    with pytest.raises(PositionalUncertaintyError, match="frozen or absent"):
        atom_uncertainty(data, 7, 298.15)


def test_diatomic_bond_length_fluctuation_matches_analytic_result():
    """A diatomic stretch is exactly solvable: sigma = sqrt(hbar / 2 mu omega).

    This also exercises the translation invariance of the difference
    coordinate: the free diatomic's Hessian is singular along the translations,
    yet the bond-length fluctuation is finite and correct.
    """
    force_constant = 36.0  # eV/Angstrom**2, roughly an H-H stretch
    mass = 1.008
    axis = np.array([0.0, 0.0, 1.0])
    projector = np.outer(axis, axis)
    hessian = np.zeros((6, 6))
    hessian[:3, :3] = force_constant * projector
    hessian[3:, 3:] = force_constant * projector
    hessian[:3, 3:] = -force_constant * projector
    hessian[3:, :3] = -force_constant * projector
    data = characterization(hessian, [mass, mass])

    result = pair_distance_uncertainty(
        data, 0, 1, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]], 1.0, quantum=True
    )
    reduced = (mass / 2.0) * atomic_mass
    omega = math.sqrt(force_constant * electron_volt / 1e-20 / reduced)
    expected = math.sqrt(hbar / (2.0 * reduced * omega)) * 1e10
    assert result["sigma_angstrom"] == pytest.approx(expected, rel=1e-8)
    # Four of the six modes are rigid translations/rotations and contribute zero.
    assert result["excluded_nonpositive_modes"] >= 1


def test_pair_distance_is_translation_invariant():
    force_constant = 36.0
    axis = np.array([0.0, 0.0, 1.0])
    projector = np.outer(axis, axis)
    hessian = np.zeros((6, 6))
    hessian[:3, :3] = force_constant * projector
    hessian[3:, 3:] = force_constant * projector
    hessian[:3, 3:] = -force_constant * projector
    hessian[3:, :3] = -force_constant * projector
    data = characterization(hessian, [1.008, 1.008])
    here = pair_distance_uncertainty(data, 0, 1, [[0, 0, 0], [0, 0, 0.74]], 298.15)
    moved = pair_distance_uncertainty(
        data, 0, 1, [[5.0, -2.0, 3.0], [5.0, -2.0, 3.74]], 298.15
    )
    assert here["sigma_angstrom"] == pytest.approx(moved["sigma_angstrom"], rel=1e-12)


def test_margin_comparison_reports_a_ratio_and_no_probability():
    result = compare_to_margin(0.1, 2.495)
    assert result["margin_in_sigma"] == pytest.approx(24.95, rel=1e-12)
    # No *field* may carry a probability or rate. The prose deliberately uses
    # those words to deny them, so the check is on the returned keys and on the
    # numeric values, not on the disclaimer text.
    for key in result:
        assert not any(
            token in key for token in ("probability", "success", "error_rate", "reliab", "yield")
        )
    numeric = [value for value in result.values() if isinstance(value, float)]
    assert numeric and all(value > 1.0 or value in (0.1, 2.495) for value in numeric)
    assert result["known_omissions"]


def test_rejects_nonsense_inputs():
    with pytest.raises(PositionalUncertaintyError):
        compare_to_margin(0.0, 1.0)
    with pytest.raises(PositionalUncertaintyError):
        compare_to_margin(0.1, float("nan"))
    with pytest.raises(PositionalUncertaintyError):
        covariance({"hessian_ev_per_angstrom2": [[1.0]]}, 298.15)
    with pytest.raises(PositionalUncertaintyError):
        atom_uncertainty(isotropic(20.0), 0, -5.0)


def test_end_to_end_against_the_real_core_hessian():
    """Drive the real ``characterize_stationary_point`` contract, not a mock.

    The other tests build the characterization dict themselves, which cannot
    catch a mismatch with the core module's actual output. This one attaches a
    genuine ASE calculator to a relaxed dimer and pushes the real returned
    dictionary through this module. Lennard-Jones is used deliberately: it costs
    microseconds, so the interface is exercised without an electronic-structure
    calculation and without competing for a loaded host. It is a contract test,
    not a chemical result.
    """
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    from ase.optimize import BFGS

    from nanodesign.stationary import characterize_stationary_point

    dimer = Atoms("Ar2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 3.8]], pbc=False)
    dimer.calc = LennardJones(sigma=3.4, epsilon=0.01)
    BFGS(dimer, logfile=None).run(fmax=1e-6, steps=200)

    data = characterize_stationary_point(dimer, step_angstrom=1e-3, force_tolerance_ev_per_angstrom=1e-4)
    assert data["free_coordinate_count"] == 6

    # An unanchored dimer has five near-zero modes, so an absolute displacement
    # is undefined and the module must refuse rather than return a huge number.
    # Which branch fires depends on whether finite differences put those modes
    # a hair above or below zero; here they land at about -8e-19, so the
    # nonpositive branch catches it. Either refusal is correct.
    with pytest.raises(PositionalUncertaintyError, match="nonpositive mode|softest mode"):
        covariance(data, 298.15)

    # The bond length is still well defined, because it is a difference.
    result = pair_distance_uncertainty(data, 0, 1, dimer.positions, 298.15, quantum=True)
    assert result["sigma_angstrom"] > 0.0
    assert math.isfinite(result["sigma_angstrom"])

    # Cross-check against the exact diatomic amplitude at this temperature,
    # using the stiffest mode the core module itself reports so the two paths
    # are independent. The full expression is required, not the zero-point
    # limit: this dimer's stretch is only about 26 cm^-1, so at 298 K
    # coth(hbar*omega / 2kT) is about 16 and the motion is nearly classical.
    # Using sqrt(hbar / 2 mu omega) here would be wrong by a factor of 4.
    eigenvalue = max(data["mass_weighted_eigenvalues_ev_per_angstrom2_amu"])
    omega = math.sqrt(eigenvalue * electron_volt / (1e-20 * atomic_mass))
    reduced = (39.948 / 2.0) * atomic_mass
    thermal = 1.0 / math.tanh(hbar * omega / (2.0 * Boltzmann * 298.15))
    expected = math.sqrt(hbar / (2.0 * reduced * omega) * thermal) * 1e10
    assert result["sigma_angstrom"] == pytest.approx(expected, rel=1e-6)
    assert thermal > 10.0  # this mode really is in the classical regime

    # And it must sit just above the purely classical amplitude sqrt(kT/k),
    # with k = lambda * mu recovered from the mass-weighted eigenvalue.
    force_constant = eigenvalue * (39.948 / 2.0)  # eV/Angstrom**2
    classical = math.sqrt(Boltzmann * 298.15 / electron_volt / force_constant)
    assert result["sigma_angstrom"] == pytest.approx(classical, rel=1e-3)
    assert result["sigma_angstrom"] > classical


def test_stiffness_and_spread_are_mutually_consistent():
    """Cross-check the two paths against the classical relation k = kT / sigma^2.

    An earlier version returned a stiffness from atom_uncertainty with k_B T and
    the variance inverted, giving 25.7 N/m for a system built at 10 N/m. Nothing
    caught it because the two paths were never compared. They are now.
    """
    for built in (5.0, 10.0, 47.5, 300.0):
        data = isotropic(built)
        stiffness = effective_stiffness(data, 0)["softest_stiffness_n_per_m"]
        assert stiffness == pytest.approx(built, rel=1e-9)

        sigma = atom_uncertainty(data, 0, 300.0, quantum=False)["largest_sigma_angstrom"]
        # Classically k = k_B T / sigma^2, which must reproduce the same value.
        from positional_uncertainty import EV_PER_ANGSTROM2_TO_N_PER_M as CONVERT

        implied = (Boltzmann * 300.0 / electron_volt) / sigma**2 * CONVERT
        assert implied == pytest.approx(built, rel=1e-9)


def test_atom_uncertainty_reports_no_numeric_stiffness():
    """Stiffness has exactly one source in this module, and it is not here."""
    result = atom_uncertainty(isotropic(10.0), 0, 300.0, quantum=False)
    assert isinstance(result["stiffness"], str)
    assert not any(
        isinstance(value, float) and "stiff" in key
        for key, value in result.items()
    )
