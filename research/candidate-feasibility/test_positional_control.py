"""Checks on the positional-control physics.

These functions produce the only number in this lane that is a design
specification rather than a measurement, so they are checked against closed
forms and against a brute-force geometric search rather than trusted.

    python -m pytest research/candidate-feasibility/test_positional_control.py -q
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

from nanodesign.candidates import make_h_abstraction  # noqa: E402
from positional_control import (  # noqa: E402
    BOLTZMANN_EV_PER_K,
    DREXLER_ERROR_TARGET,
    classical_sigma,
    competing_site_geometry,
    error_probability,
    fit_curvature,
    quantum_sigma,
    required_sigma,
)


def test_required_sigma_round_trips_through_the_error_probability():
    """The sigma the specification demands must land exactly on the target."""
    for radius in (0.5, 1.0, 2.878, 5.0):
        sigma = required_sigma(radius)
        assert error_probability(radius, sigma) == pytest.approx(DREXLER_ERROR_TARGET, rel=1e-9)


def test_error_probability_matches_the_two_dimensional_gaussian_integral():
    """exp(-R^2/2s^2) must equal the numerically integrated 2-D tail."""
    sigma, radius = 0.3, 0.75
    # P(r > R) for an isotropic 2-D Gaussian, by direct radial integration.
    r = np.linspace(radius, radius + 40 * sigma, 400_000)
    density = (r / sigma ** 2) * np.exp(-(r ** 2) / (2 * sigma ** 2))
    assert np.trapezoid(density, r) == pytest.approx(error_probability(radius, sigma), rel=1e-6)


def test_classical_sigma_is_equipartition():
    """Amplitude goes as 1/sqrt(stiffness), not equally across coordinates."""
    temperature = 298.15
    assert classical_sigma(4.0, temperature) == pytest.approx(
        classical_sigma(1.0, temperature) / 2.0
    )
    assert classical_sigma(2.0, temperature) ** 2 == pytest.approx(
        BOLTZMANN_EV_PER_K * temperature / 2.0
    )


def test_quantum_sigma_approaches_the_classical_result_when_hot_and_soft():
    """coth(x) -> 1/x for small x, so the quantum result must recover kT/k."""
    stiffness, mass, temperature = 0.05, 200.0, 1000.0
    quantum = quantum_sigma(stiffness, mass, temperature)
    assert quantum["hbar_omega_ev"] < 0.05 * BOLTZMANN_EV_PER_K * temperature
    assert quantum["sigma_angstrom"] == pytest.approx(
        classical_sigma(stiffness, temperature), rel=2e-3
    )


def test_quantum_sigma_has_a_zero_point_floor_that_cooling_cannot_beat():
    """The whole point: cooling does not take sigma to zero."""
    stiffness, mass = 5.0, 27.0
    cold = quantum_sigma(stiffness, mass, 1.0)
    colder = quantum_sigma(stiffness, mass, 0.01)
    floor = cold["zero_point_sigma_angstrom"]
    assert cold["sigma_angstrom"] == pytest.approx(floor, rel=1e-9)
    assert colder["sigma_angstrom"] == pytest.approx(floor, rel=1e-9)
    # And the classical result wrongly goes to zero over the same range.
    assert classical_sigma(stiffness, 0.01) < floor / 10.0


def test_quantum_sigma_never_falls_below_the_classical_one():
    """Zero-point motion adds to thermal motion; it cannot subtract."""
    for temperature in (4.0, 77.0, 298.15, 1000.0):
        quantum = quantum_sigma(3.0, 27.0, temperature)
        assert quantum["sigma_angstrom"] >= classical_sigma(3.0, temperature) - 1e-12


def test_cold_stiff_case_does_not_overflow():
    """coth of a large argument must be handled, not raised."""
    result = quantum_sigma(50.0, 12.0, 0.1)
    assert math.isfinite(result["sigma_angstrom"])
    assert result["quantum_regime"] is True


def test_geometry_margin_matches_a_brute_force_azimuthal_search():
    """The analytic crossover must agree with sliding the apex numerically."""
    geometry = competing_site_geometry()
    reactant, _, metadata = make_h_abstraction(3.6, 0.0)
    positions = reactant.positions
    symbols = reactant.get_chemical_symbols()
    apex = positions[metadata["tip_apex"]]
    target = positions[metadata["transferred_hydrogen"]]
    rivals = np.array([
        positions[i] for i in metadata["substrate_indices"]
        if symbols[i] == "H" and i != metadata["transferred_hydrogen"]
    ])

    azimuths = np.linspace(0.0, 2 * np.pi, 181)
    offsets = np.arange(0.0, 6.0, 0.005)
    directions = np.stack([np.cos(azimuths), np.sin(azimuths), np.zeros_like(azimuths)], axis=1)
    swept = apex[None, None, :] + offsets[None, :, None] * directions[:, None, :]
    to_target = np.linalg.norm(swept - target, axis=2)
    to_rivals = np.linalg.norm(swept[:, :, None, :] - rivals[None, None, :, :], axis=3).min(axis=2)
    crossed = to_rivals <= to_target
    first = np.where(crossed.any(axis=1), np.argmax(crossed, axis=1), len(offsets) - 1)
    assert offsets[first].min() == pytest.approx(geometry["target_radius_angstrom"], abs=0.01)


def test_the_naive_lateral_criterion_is_recorded_and_rejected():
    """The wrong-turn must stay visible, since it is the obvious thing to do."""
    geometry = competing_site_geometry()
    assert geometry["rejected_naive_radius_angstrom"] < geometry["target_radius_angstrom"] / 3
    # The rivals it would have picked are behind the cage, not reachable.
    assert geometry["nearest_rival_z_angstrom"] > -1.0


def test_the_real_rival_is_an_equatorial_methylene_hydrogen():
    geometry = competing_site_geometry()
    assert geometry["nearest_rival_index"] in range(14, 20)
    assert geometry["n_rivals_at_that_offset"] >= 1


def test_curvature_fit_recovers_a_known_parabola():
    stiffness = 0.8
    points = [
        {"lateral_offset_angstrom": x, "energy_ev": 0.5 * stiffness * x ** 2 - 3.0}
        for x in (-0.4, -0.2, 0.0, 0.2, 0.4)
    ]
    fit = fit_curvature(points)
    assert fit["fitted"] is True
    assert fit["stiffness_ev_per_angstrom_squared"] == pytest.approx(stiffness, rel=1e-8)
    assert fit["minimum_at_offset_angstrom"] == pytest.approx(0.0)


def test_curvature_fit_reports_failure_rather_than_guessing():
    fit = fit_curvature([{"lateral_offset_angstrom": 0.0, "energy_ev": 1.0}])
    assert fit["fitted"] is False
    assert "need 3" in fit["reason"]


def test_curvature_fit_survives_failed_scan_points():
    """A failed calculation must be skipped, not crash the fit."""
    points = [
        {"lateral_offset_angstrom": -0.2, "energy_ev": 0.02},
        {"lateral_offset_angstrom": 0.0, "energy_ev": 0.0},
        {"lateral_offset_angstrom": 0.2, "error": "QuantumCalculationError: no convergence"},
        {"lateral_offset_angstrom": 0.4, "energy_ev": 0.08},
    ]
    fit = fit_curvature(points)
    assert fit["fitted"] is True
    assert fit["points_used"] == 3
