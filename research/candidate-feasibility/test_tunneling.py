"""Checks on the tunneling corrections.

Every function here is a closed form, so each is checked against an
independently derived expression rather than against itself. The constant in
the crossover temperature is the load-bearing one: the claim that this
reaction is in the tunneling regime at room temperature stands or falls on it.

    python -m pytest research/candidate-feasibility/test_tunneling.py -q
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import pytest

LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(LANE))

from tunneling import (  # noqa: E402
    angular_frequency,
    bell_kappa,
    crossover_temperature,
    effective_barrier_reduction_kcal,
    isotope_effect,
    survey,
    wigner_kappa,
)

# hc/k_B = 1.43877 K cm, so T_c = 1.43877 * wavenumber / (2 pi).
SECOND_RADIATION_CONSTANT_K_CM = 1.438776877


def test_crossover_temperature_matches_the_spectroscopic_closed_form():
    """T_c = hc nu / (2 pi k_B), derived independently of the module's constants."""
    for wavenumber in (500.0, 1000.0, 1500.0, 2000.0):
        expected = SECOND_RADIATION_CONSTANT_K_CM * wavenumber / (2.0 * math.pi)
        assert crossover_temperature(wavenumber) == pytest.approx(expected, rel=1e-5)


def test_the_headline_number():
    """A 1500i cm^-1 saddle puts the crossover above room temperature."""
    t_c = crossover_temperature(1500.0)
    assert t_c == pytest.approx(343.5, abs=0.5)
    assert t_c > 298.15


def test_crossover_scales_linearly_with_the_imaginary_frequency():
    assert crossover_temperature(2000.0) == pytest.approx(2.0 * crossover_temperature(1000.0))


def test_angular_frequency_conversion():
    """omega = 2 pi c nu, checked against a hand value."""
    assert angular_frequency(1000.0) == pytest.approx(2.0 * math.pi * 2.99792458e13, rel=1e-9)


def test_wigner_is_the_small_u_expansion_of_bell():
    """(u/2)/sin(u/2) = 1 + u^2/24 + O(u^4), so they agree when u is small."""
    wavenumber, temperature = 200.0, 1000.0  # deliberately deep in the classical regime
    bell = bell_kappa(wavenumber, temperature)
    assert bell["valid"]
    assert bell["u"] < 0.5
    assert wigner_kappa(wavenumber, temperature) == pytest.approx(bell["kappa"], rel=1e-4)


def test_wigner_diverges_from_bell_once_the_correction_is_large():
    """Wigner must not be trusted when it stops being small; show it breaks."""
    wavenumber, temperature = 1200.0, 298.15
    bell = bell_kappa(wavenumber, temperature)
    assert bell["valid"]
    assert bell["kappa"] > 5.0
    assert wigner_kappa(wavenumber, temperature) < bell["kappa"] * 0.75


def test_bell_refuses_below_the_crossover_temperature():
    """The divergence must be reported as invalidity, not returned as a number."""
    wavenumber = 1500.0
    t_c = crossover_temperature(wavenumber)
    below = bell_kappa(wavenumber, t_c * 0.99)
    assert below["valid"] is False
    assert below["kappa"] is None
    assert "dominant" in below["reason"]
    above = bell_kappa(wavenumber, t_c * 1.05)
    assert above["valid"] is True
    assert above["kappa"] > 1.0


def test_bell_blows_up_as_the_crossover_is_approached_from_above():
    wavenumber = 1000.0
    t_c = crossover_temperature(wavenumber)
    near = bell_kappa(wavenumber, t_c * 1.001)["kappa"]
    far = bell_kappa(wavenumber, t_c * 2.0)["kappa"]
    assert near > far > 1.0
    assert near > 100.0


def test_u_equals_two_pi_exactly_at_the_crossover():
    wavenumber = 1337.0
    at = bell_kappa(wavenumber, crossover_temperature(wavenumber))
    assert at["u"] == pytest.approx(2.0 * math.pi, rel=1e-9)
    assert at["valid"] is False


def test_tunneling_is_always_an_enhancement():
    for wavenumber in (200.0, 500.0, 800.0):
        for temperature in (298.15, 500.0, 1000.0):
            assert wigner_kappa(wavenumber, temperature) >= 1.0
            bell = bell_kappa(wavenumber, temperature)
            if bell["valid"]:
                assert bell["kappa"] >= 1.0


def test_equivalent_barrier_lowering_inverts_the_arrhenius_factor():
    """The reported barrier error must reproduce kappa through exp(dE/kT)."""
    temperature, kappa = 298.15, 10.0
    lowering = effective_barrier_reduction_kcal(kappa, temperature)
    boltzmann_kcal = 1.380649e-23 * temperature / (4184.0 / 6.02214076e23)
    assert math.exp(lowering / boltzmann_kcal) == pytest.approx(kappa, rel=1e-9)


def test_deuterium_tunnels_less_so_the_isotope_effect_exceeds_one():
    result = isotope_effect(1200.0, 298.15)
    assert result["deuterium_wavenumber_cm"] == pytest.approx(1200.0 / math.sqrt(2.0))
    assert result["tunneling_contribution_to_kie"] > 1.0
    assert result["hydrogen_kappa"] > result["deuterium_kappa"]


def test_isotope_effect_reports_rather_than_fabricates_below_crossover():
    result = isotope_effect(2000.0, 77.0)
    assert result["tunneling_contribution_to_kie"] is None


def test_survey_flags_room_temperature_against_the_crossover():
    rows = {row["imaginary_wavenumber_cm"]: row for row in survey()["rows"]}
    assert rows[1000.0]["room_temperature_is_below_crossover"] is False
    assert rows[1500.0]["room_temperature_is_below_crossover"] is True
    # Cryogenic operation is below the crossover even for a soft barrier.
    assert rows[500.0]["temperatures"]["77K"]["bell_valid"] is False
