"""Checks on the requirement inverter, including against real committed evidence.

The last test is the important one: it reads a DFT Hessian another lane already
computed and committed, runs both the forward and inverse paths over it, and
confirms they agree with closed-form diatomic physics.  It costs no
electronic-structure time because the calculation was already paid for.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from scipy.constants import Boltzmann, atomic_mass, electron_volt, hbar, speed_of_light

from positional_requirements import (
    EV_PER_ANGSTROM2_TO_N_PER_M,
    PositionalRequirementError,
    assess,
    crossing_probability,
    crossover_wavenumber,
    quantum_correction,
    required_sigma,
    required_stiffness,
)

# The margin measured by this lane's stage-0 census, r5.
MARGIN = 2.4949626493635075
REPOSITORY = Path(__file__).resolve().parents[3]


def test_probability_and_sigma_invert_each_other():
    for target in (1e-2, 1e-6, 1e-12, 1e-15):
        sigma = required_sigma(MARGIN, target)
        assert crossing_probability(sigma, MARGIN) == pytest.approx(target, rel=1e-9)


def test_one_sigma_crossing_probability_is_the_standard_normal_tail():
    # At sigma = margin the one-sided tail is the familiar 15.87 percent.
    assert crossing_probability(1.0, 1.0) == pytest.approx(0.158655254, rel=1e-8)
    assert crossing_probability(1.0, 2.0) == pytest.approx(0.022750132, rel=1e-8)


def test_classical_requirement_is_exactly_kt_over_sigma_squared():
    result = required_stiffness(MARGIN, 1e-15, 300.0)
    sigma = result["required_sigma_angstrom"]
    expected = (Boltzmann * 300.0 / electron_volt) / sigma**2 * EV_PER_ANGSTROM2_TO_N_PER_M
    assert result["required_stiffness_classical_n_per_m"] == pytest.approx(expected, rel=1e-12)


def test_the_measured_margin_needs_only_a_few_newtons_per_metre():
    """The headline inversion, pinned so a regression would be visible."""
    result = required_stiffness(MARGIN, 1e-15, 300.0)
    assert result["required_sigma_angstrom"] == pytest.approx(0.3142, abs=1e-4)
    assert result["required_stiffness_classical_n_per_m"] == pytest.approx(4.196, abs=1e-3)
    # Far below what a stiff diamondoid mount is expected to supply.
    assert result["required_stiffness_classical_n_per_m"] < 10.0


def test_requirement_tightens_with_a_more_demanding_target():
    loose = required_stiffness(MARGIN, 1e-3, 300.0)["required_stiffness_classical_n_per_m"]
    tight = required_stiffness(MARGIN, 1e-15, 300.0)["required_stiffness_classical_n_per_m"]
    assert tight > loose


def test_requirement_tightens_with_temperature():
    cold = required_stiffness(MARGIN, 1e-15, 77.0)["required_stiffness_classical_n_per_m"]
    warm = required_stiffness(MARGIN, 1e-15, 300.0)["required_stiffness_classical_n_per_m"]
    assert warm == pytest.approx(cold * 300.0 / 77.0, rel=1e-9)


def test_quantum_correction_limits():
    soft = quantum_correction(1.0, 300.0)
    assert soft["sigma_quantum_over_sigma_classical"] == pytest.approx(1.0, abs=1e-6)
    assert "near-classical" in soft["regime"]

    stiff = quantum_correction(20000.0, 300.0)
    # For large x the ratio grows as sqrt(x).
    x = stiff["hbar_omega_over_2kT"]
    assert stiff["sigma_quantum_over_sigma_classical"] == pytest.approx(math.sqrt(x), rel=1e-6)
    assert "zero-point dominated" in stiff["regime"]


def test_crossover_is_a_real_root():
    result = crossover_wavenumber(300.0, 0.10)
    wavenumber = result["crossover_wavenumber_cm1"]
    assert wavenumber == pytest.approx(338.0, abs=1.0)
    assert quantum_correction(wavenumber, 300.0)["sigma_quantum_over_sigma_classical"] == pytest.approx(1.10, rel=1e-6)
    # A colder machine pushes the crossover down, so more modes need the quantum form.
    assert crossover_wavenumber(77.0, 0.10)["crossover_wavenumber_cm1"] < wavenumber


def test_assess_reports_headroom_and_meets_the_requirement_at_30_n_per_m():
    result = assess(MARGIN, 30.0, 1e-15, 300.0)
    assert result["requirement_met"]
    assert result["stiffness_headroom_factor"] == pytest.approx(30.0 / 4.196, rel=1e-3)
    assert result["achievable_sigma_angstrom_classical"] == pytest.approx(0.1175, abs=1e-4)
    assert result["margin_in_achievable_sigma"] == pytest.approx(21.2, abs=0.1)


def test_assess_fails_a_mount_that_is_too_soft():
    result = assess(MARGIN, 1.0, 1e-15, 300.0)
    assert not result["requirement_met"]
    assert result["stiffness_headroom_factor"] < 1.0


def test_no_field_claims_to_be_an_error_rate():
    result = required_stiffness(MARGIN, 1e-15, 300.0)
    # The disclaimer key deliberately names the thing it denies, so exclude it
    # and check that nothing else offers a rate or reliability figure.
    assert "not_an_error_rate" in result
    for key, value in result.items():
        if key == "not_an_error_rate":
            continue
        assert "error_rate" not in key
        assert "reliab" not in key
        # And no numeric field is a bare probability masquerading as a result.
        if isinstance(value, float):
            assert key != "probability"


def test_rejects_ill_posed_requests():
    for bad in (0.0, 0.5, 1.0, -1e-3, float("nan")):
        with pytest.raises(PositionalRequirementError):
            required_sigma(MARGIN, bad)
    with pytest.raises(PositionalRequirementError):
        required_stiffness(-1.0, 1e-9)
    with pytest.raises(PositionalRequirementError):
        quantum_correction(0.0)


def test_quantum_requirement_needs_a_mass_and_matches_classical_when_soft():
    """With a heavy effective mass the mode is soft, so both routes agree."""
    result = required_stiffness(MARGIN, 1e-15, 300.0, effective_mass_amu=500.0)
    assert result["required_stiffness_quantum_n_per_m"] == pytest.approx(
        result["required_stiffness_classical_n_per_m"], rel=0.02
    )
    assert result["required_mode_wavenumber_cm1"] < crossover_wavenumber(300.0)["crossover_wavenumber_cm1"]
    assert "near-classical" in result["mode_regime"]


def test_against_the_committed_h2_dft_hessian():
    """Reuse a Hessian another lane already computed; spend nothing, verify hard.

    Confirms three independent things agree on real DFT data: the forward
    quantum amplitude from the stored Hessian, the closed-form diatomic
    zero-point amplitude, and this module's quantum/classical ratio.
    """
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from ase.io import read

    from positional_uncertainty import pair_distance_uncertainty

    stored = REPOSITORY / "data/validation/h2-integration/modes/result.json"
    geometry = REPOSITORY / "data/validation/h2-integration/modes/input.extxyz"
    if not stored.exists() or not geometry.exists():
        pytest.skip("The reused H2 evidence is not present in this checkout")

    record = json.loads(stored.read_text())
    stationary = record["stationary"]
    assert record["quantum_settings"]["xc"] == "pbe0"
    atoms = read(geometry)

    temperature = 298.15
    quantum = pair_distance_uncertainty(
        stationary, 0, 1, atoms.positions, temperature, quantum=True
    )["sigma_angstrom"]
    classical = pair_distance_uncertainty(
        stationary, 0, 1, atoms.positions, temperature, quantum=False
    )["sigma_angstrom"]

    # 1. The stretch frequency the core module reported, used independently.
    wavenumber = max(stationary["frequencies_cm1"])
    assert wavenumber == pytest.approx(4383.9, abs=0.1)

    # 2. Closed-form diatomic zero-point amplitude at that frequency.
    omega = wavenumber * 2.0 * math.pi * speed_of_light * 100.0
    reduced = (1.008 / 2.0) * atomic_mass
    analytic = math.sqrt(hbar / (2.0 * reduced * omega)) * 1e10
    assert quantum == pytest.approx(analytic, rel=1e-4)
    assert quantum == pytest.approx(0.087348, abs=1e-5)

    # 3. This module's ratio must reproduce the measured ratio of the two paths.
    ratio = quantum_correction(wavenumber, temperature)["sigma_quantum_over_sigma_classical"]
    assert quantum / classical == pytest.approx(ratio, rel=1e-6)
    assert ratio == pytest.approx(3.25, abs=0.02)

    # The whole point: at room temperature this bond is zero-point dominated, so
    # the classical formula is not slightly off, it misses most of the spread.
    assert quantum > 3.0 * classical
    cold = pair_distance_uncertainty(
        stationary, 0, 1, atoms.positions, 4.0, quantum=True
    )["sigma_angstrom"]
    assert cold == pytest.approx(quantum, rel=1e-9)


def test_cooling_widens_the_margin_for_a_soft_mount():
    """Independently reproduces a peer lane's cryogenic extension.

    A 30 N/m mount moving 50 amu is a ~101 cm^-1 mode. Its quantum spread
    saturates at a small zero-point floor, so cooling keeps widening the margin
    instead of hitting a wall, even though the classical formula becomes badly
    optimistic down there.
    """
    from positional_requirements import cooling_assessment, mode_wavenumber

    assert mode_wavenumber(30.0, 50.0) == pytest.approx(101.0, abs=0.5)

    result = cooling_assessment(MARGIN, 30.0, 50.0, (298.15, 77.0, 4.0))
    rows = {row["temperature_kelvin"]: row for row in result["temperatures"]}

    warm = rows[298.15]
    assert warm["sigma_classical_angstrom"] == pytest.approx(0.1171, abs=1e-4)
    assert warm["sigma_quantum_angstrom"] == pytest.approx(0.1183, abs=1e-4)
    assert warm["margin_in_quantum_sigma"] == pytest.approx(21.0, abs=0.5)

    cold = rows[4.0]
    assert cold["sigma_classical_angstrom"] == pytest.approx(0.0136, abs=1e-4)
    assert cold["sigma_quantum_angstrom"] == pytest.approx(0.0578, abs=1e-4)
    assert cold["margin_in_quantum_sigma"] == pytest.approx(43.0, abs=1.0)
    # The classical formula is optimistic by over 4x at 4 K.
    assert cold["classical_understates_by_factor"] == pytest.approx(4.26, abs=0.05)

    # The headline: cooling helps rather than hurts, for this mount.
    assert result["cooling_materially_helps"]
    assert result["margin_improvement_factor_on_cooling"] == pytest.approx(43.2/21.1, rel=0.05)
    assert result["margin_in_sigma_floor"] == pytest.approx(43.2, abs=1.0)


def test_a_stiff_mode_gains_nothing_from_cooling():
    """The contrasting case the conditional rule exists to separate."""
    from positional_requirements import cooling_assessment

    # A bond-stiffness mode on a light mass is already at its zero-point floor.
    result = cooling_assessment(MARGIN, 570.0, 0.504, (298.15, 4.0))
    rows = {row["temperature_kelvin"]: row for row in result["temperatures"]}
    warm, cold = rows[298.15], rows[4.0]
    assert warm["sigma_quantum_angstrom"] == pytest.approx(cold["sigma_quantum_angstrom"], rel=1e-6)
    # Already at its zero-point floor at room temperature: the residual gain is
    # about 1e-8, so cooling buys nothing a designer could use.
    assert not result["cooling_materially_helps"]
    assert result["margin_improvement_factor_on_cooling"] == pytest.approx(1.0, abs=1e-6)


def test_crossover_scales_linearly_with_temperature():
    """Peer-reported table: 336 cm^-1 at 298 K down to 4.5 at 4 K."""
    for temperature, expected in ((298.15, 336.0), (150.0, 169.0), (77.0, 87.0), (4.0, 4.5)):
        got = crossover_wavenumber(temperature, 0.10)["crossover_wavenumber_cm1"]
        assert got == pytest.approx(expected, rel=0.02)


def test_six_degenerate_competitors_barely_change_the_requirement():
    """The census finds a six-fold tied shell, so the honest target is a union bound.

    The point of computing it is to show it does not matter: the requirement
    grows only logarithmically in the number of competitors, so accounting for
    all six moves it by a few percent rather than by a factor of six.
    """
    single = required_stiffness(MARGIN, 1e-15, 300.0, competitor_count=1)
    six = required_stiffness(MARGIN, 1e-15, 300.0, competitor_count=6)

    assert six["per_competitor_target_probability"] == pytest.approx(1e-15 / 6)
    assert six["required_sigma_angstrom"] < single["required_sigma_angstrom"]
    # sigma = margin / (sqrt(2) * erfcinv(2e-15/6)); erfcinv(3.333e-16) = 5.772,
    # so margin/sigma = 8.161 and sigma = 0.3057 A. Checked by hand.
    assert six["required_sigma_angstrom"] == pytest.approx(0.3057, abs=1e-4)
    assert six["required_stiffness_classical_n_per_m"] == pytest.approx(4.431, abs=0.01)
    # A sixfold tighter probability budget costs under ten percent in stiffness.
    ratio = (
        six["required_stiffness_classical_n_per_m"]
        / single["required_stiffness_classical_n_per_m"]
    )
    assert 1.0 < ratio < 1.10
    # And it is still far below any real mount, which is the conclusion that matters.
    assert six["required_stiffness_classical_n_per_m"] < 10.0


def test_competitor_count_is_validated():
    for bad in (0, -1, 2.5, True):
        with pytest.raises(PositionalRequirementError):
            required_stiffness(MARGIN, 1e-15, 300.0, competitor_count=bad)


def test_maximum_operating_temperature_brackets_room_temperature():
    """The measured tip, expressed as an operating specification.

    Both criteria are quoted because they bracket rather than disagree: the union
    over bisector half-spaces is the exact failure event, and the sphere-exit
    criterion is a strict upper bound on it, so it yields the lower T_max.
    """
    from positional_requirements import maximum_operating_temperature

    measured = 7.2715  # apex carbon, propyne model, compliance-based

    exact = maximum_operating_temperature(measured, 4.431, 300.0)
    conservative = maximum_operating_temperature(measured, 4.825, 298.15)

    assert exact["maximum_operating_temperature_kelvin"] == pytest.approx(492.3, abs=1.0)
    assert conservative["maximum_operating_temperature_kelvin"] == pytest.approx(449.3, abs=1.0)
    # The conservative criterion must give the lower limit, never the higher.
    assert (
        conservative["maximum_operating_temperature_kelvin"]
        < exact["maximum_operating_temperature_kelvin"]
    )
    for result in (exact, conservative):
        assert result["meets_requirement_at_reference"]
        assert result["margin_above_room_temperature_kelvin"] > 100.0


def test_operating_temperature_is_linear_in_stiffness():
    from positional_requirements import maximum_operating_temperature

    single = maximum_operating_temperature(10.0, 5.0, 300.0)
    double = maximum_operating_temperature(20.0, 5.0, 300.0)
    assert single["maximum_operating_temperature_kelvin"] == pytest.approx(600.0)
    assert double["maximum_operating_temperature_kelvin"] == pytest.approx(1200.0)
    with pytest.raises(PositionalRequirementError):
        maximum_operating_temperature(10.0, 0.0, 300.0)


def test_the_union_bound_is_contained_in_the_sphere_exit_event():
    """Why the two criteria bracket rather than compete.

    Crossing a bisector plane at perpendicular distance d puts the apex at radius
    at least d, since d is the closest approach of that plane to the origin. So
    the union of half-spaces is a subset of the sphere-exit event and its
    probability can never be larger -- which is why the sphere-exit criterion
    demands the higher stiffness.
    """
    margin = MARGIN
    for competitors in (1, 6):
        union = required_stiffness(margin, 1e-15, 300.0, competitor_count=competitors)
        sigma_union = union["required_sigma_angstrom"]
        # 3D sphere-exit tail at the same sigma must be at least the union's target.
        # P(|r| > d) for an isotropic Gaussian, via the chi-3 survival function.
        from scipy.stats import chi

        sphere_tail = float(chi.sf(margin / sigma_union, df=3))
        assert sphere_tail >= 1e-15
