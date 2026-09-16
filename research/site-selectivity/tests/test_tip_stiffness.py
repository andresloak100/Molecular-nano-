"""The frequency-based corroboration of the measured mount stiffness.

A stiffness and an observed frequency together imply a reduced mass, and a real
vibration's reduced mass cannot exceed the molecule's total mass. That makes this
a falsification test rather than a plausibility check: it rules the clamped
Hessian block out on physical grounds rather than on a preference for the
compliance treatment.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from scipy.constants import atomic_mass, speed_of_light

from tip_stiffness import build_propyne, reduced_mass_check

EVIDENCE = Path(__file__).resolve().parents[1] / "evidence" / "tip-stiffness-propyne"
PROPYNE_MASS = 40.064


def test_propyne_model_is_built_correctly():
    atoms, indices = build_propyne()
    assert atoms.get_chemical_formula() == "C3H4"
    assert len(indices["mount_hydrogens"]) == 3
    # Three anchored atoms is the minimum that removes translation and rotation.
    assert atoms.get_chemical_symbols()[indices["apex_carbon"]] == "C"
    assert atoms.get_chemical_symbols()[indices["apex_hydrogen"]] == "H"


def test_reduced_mass_formula_against_a_known_case():
    """H2 at 4383.9 cm^-1 with k = 570.7 N/m must give mu = 0.504 amu."""
    result = reduced_mass_check(570.7, 4383.9, 2.016)
    assert result["implied_reduced_mass_amu"] == pytest.approx(0.504, rel=2e-3)
    assert result["physically_possible"]


def test_the_measured_stiffness_is_corroborated_by_its_frequency():
    if not (EVIDENCE / "tip_stiffness.json").exists():
        pytest.skip("propyne tip-stiffness evidence is not present in this checkout")
    record = json.loads((EVIDENCE / "tip_stiffness.json").read_text())
    apex = record["tips"]["apex_carbon"]
    softest = min(f for f in record["hessian"]["frequencies_cm1"] if f > 1.0)
    assert softest == pytest.approx(88.3, abs=0.5)

    compliance = reduced_mass_check(apex["softest_stiffness_n_per_m"], softest, PROPYNE_MASS)
    assert apex["softest_stiffness_n_per_m"] == pytest.approx(7.27, abs=0.02)
    # A transverse bend of a linear fragment should be a sizeable fraction of the
    # molecular mass, and this lands at 0.39 of it.
    assert compliance["implied_reduced_mass_amu"] == pytest.approx(15.8, abs=0.2)
    assert 0.2 < compliance["reduced_mass_over_molecular_mass"] < 0.6
    assert compliance["physically_possible"]


def test_the_clamped_stiffness_is_falsified_by_the_same_frequency():
    """The strong form: clamping is not merely wrong, it is impossible here."""
    if not (EVIDENCE / "tip_stiffness.json").exists():
        pytest.skip("propyne tip-stiffness evidence is not present in this checkout")
    record = json.loads((EVIDENCE / "tip_stiffness.json").read_text())
    apex = record["tips"]["apex_carbon"]
    softest = min(f for f in record["hessian"]["frequencies_cm1"] if f > 1.0)
    clamped = min(apex["clamped_diagonal_block_stiffness_n_per_m"])

    assert clamped == pytest.approx(79.8, abs=0.5)
    check = reduced_mass_check(clamped, softest, PROPYNE_MASS)
    # 173.6 amu for a 40.06 amu molecule: 4.3 times the whole thing.
    assert check["implied_reduced_mass_amu"] == pytest.approx(173.6, abs=1.0)
    assert check["reduced_mass_over_molecular_mass"] == pytest.approx(4.33, abs=0.05)
    assert not check["physically_possible"]
    # And the overestimate factor that a naive treatment would have introduced.
    assert clamped / apex["softest_stiffness_n_per_m"] == pytest.approx(11.0, abs=0.2)


def test_the_measured_stiffness_still_meets_the_requirement_but_barely():
    if not (EVIDENCE / "tip_stiffness.json").exists():
        pytest.skip("propyne tip-stiffness evidence is not present in this checkout")
    record = json.loads((EVIDENCE / "tip_stiffness.json").read_text())
    apex = record["tips"]["apex_carbon"]
    assert apex["requirement_met"]
    # Met, but under a factor of two: the earlier sevenfold claim is withdrawn.
    assert 1.5 < apex["stiffness_headroom_factor"] < 2.2


def test_the_soft_mode_is_below_the_quantum_crossover():
    """Consistency with the crossover table: a soft mode needs no quantum term."""
    if not (EVIDENCE / "tip_stiffness.json").exists():
        pytest.skip("propyne tip-stiffness evidence is not present in this checkout")
    from positional_requirements import crossover_wavenumber

    record = json.loads((EVIDENCE / "tip_stiffness.json").read_text())
    softest = min(f for f in record["hessian"]["frequencies_cm1"] if f > 1.0)
    assert softest < crossover_wavenumber(300.0)["crossover_wavenumber_cm1"]
    apex = record["tips"]["apex_carbon"]
    assert apex["quantum_over_classical"] == pytest.approx(1.009, abs=0.005)
