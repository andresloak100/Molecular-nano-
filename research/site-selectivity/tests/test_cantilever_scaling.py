"""The tip behaves as a cantilever, so stiffness scales as the inverse cube of length."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

from tip_stiffness import cantilever_scaling, maximum_tip_length

EVIDENCE = Path(__file__).resolve().parents[1] / "evidence" / "tip-stiffness-propyne"


def measured():
    record = json.loads((EVIDENCE / "tip_stiffness.json").read_text())
    atoms = read(EVIDENCE / "relaxed.xyz")
    positions = np.asarray(atoms.positions)
    indices = record["indices"]
    anchor_z = positions[indices["mount_hydrogens"]][:, 2].mean()
    return {
        "k_carbon": record["tips"]["apex_carbon"]["softest_stiffness_n_per_m"],
        "k_hydrogen": record["tips"]["apex_hydrogen"]["softest_stiffness_n_per_m"],
        "lever_carbon": float(positions[indices["apex_carbon"]][2] - anchor_z),
        "lever_hydrogen": float(positions[indices["apex_hydrogen"]][2] - anchor_z),
    }


def test_the_two_measured_reference_points_obey_cubic_scaling():
    """The load-bearing check: one Hessian gives two levers, and they agree."""
    if not (EVIDENCE / "tip_stiffness.json").exists():
        pytest.skip("propyne tip-stiffness evidence is not present")
    data = measured()
    assert data["lever_carbon"] == pytest.approx(3.0593, abs=1e-3)
    assert data["lever_hydrogen"] == pytest.approx(4.1332, abs=1e-3)

    predicted = cantilever_scaling(data["k_carbon"], data["lever_carbon"], data["lever_hydrogen"])
    assert predicted == pytest.approx(data["k_hydrogen"], rel=0.02)
    # And the cubic exponent beats the neighbouring integers decisively.
    ratio = data["k_carbon"] / data["k_hydrogen"]
    lengths = data["lever_hydrogen"] / data["lever_carbon"]
    assert abs(ratio - lengths**3) < abs(ratio - lengths**2)
    assert abs(ratio - lengths**3) < abs(ratio - lengths**4)


def test_a_stiffness_headroom_is_only_its_cube_root_as_a_length_budget():
    if not (EVIDENCE / "tip_stiffness.json").exists():
        pytest.skip("propyne tip-stiffness evidence is not present")
    data = measured()
    result = maximum_tip_length(data["k_carbon"], data["lever_carbon"], 4.431)
    assert result["stiffness_headroom_factor"] == pytest.approx(1.641, abs=0.01)
    # 1.64x in stiffness is only 1.18x in length.
    assert result["maximum_lever_angstrom"] == pytest.approx(3.611, abs=0.01)
    assert result["spare_length_fraction"] == pytest.approx(0.18, abs=0.01)
    assert result["spare_length_angstrom"] == pytest.approx(0.55, abs=0.02)
    assert result["spare_length_fraction"] == pytest.approx(
        result["stiffness_headroom_factor"] ** (1 / 3) - 1, rel=1e-9
    )


def test_scaling_is_monotonic_and_validated():
    assert cantilever_scaling(10.0, 3.0, 6.0) == pytest.approx(10.0 / 8.0)
    assert cantilever_scaling(10.0, 3.0, 1.5) == pytest.approx(80.0)
    for bad in ((0.0, 1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 1.0, 0.0)):
        with pytest.raises(ValueError):
            cantilever_scaling(*bad)
    with pytest.raises(ValueError):
        maximum_tip_length(1.0, 1.0, 0.0)
