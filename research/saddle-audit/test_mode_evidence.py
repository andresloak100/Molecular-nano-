"""Synthetic algebraic fixtures only; no electronic structure or fitted accuracy."""
from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from mode_evidence import transfer_projection

INDICES = {"acceptor_carbon": 0, "transferring_hydrogen": 1, "donor_carbon": 2}
SYMBOLS = ["C", "H", "C"]


def fixture():
    return {"geometry_angstrom": [[1., 0, 0], [0., 0, 0], [-1., 0, 0]],
            "free_atom_indices": [0, 1, 2], "free_masses_amu": [12., 1., 12.],
            "frequencies_cm1": [-100.],
            "cartesian_modes_per_sqrt_amu": [[[0., 0, 0], [1., 0, 0], [0., 0, 0]]],
            "settings": {"frequency_tolerance_cm1": 20.}}


def score(data):
    return transfer_projection(data, INDICES, SYMBOLS)["negative_modes"][0]["squared_mass_metric_overlap"]


def test_hydrogen_axial_overlap_includes_carbon_mass():
    out = transfer_projection(fixture(), INDICES, SYMBOLS)
    assert score(fixture()) == pytest.approx(.96)
    assert out["negative_modes"][0]["hydrogen_mass_weighted_displacement_share"] == 1
    assert not out["establishes_connectivity"]
    assert not out["establishes_electronic_state"]


@pytest.mark.parametrize("axial", [0., 1e-14, -1e-14])
def test_transverse_noise_does_not_pass(axial):
    data = fixture()
    data["cartesian_modes_per_sqrt_amu"][0][1] = [axial, 1., 0]
    out = transfer_projection(data, INDICES, SYMBOLS)
    assert out["negative_modes"][0]["squared_mass_metric_overlap"] <= 1e-25
    assert not out["single_negative_mode_passes_overlap_heuristic"]


@pytest.mark.parametrize("scale", [-7., -.1, .2, 9.])
def test_sign_and_amplitude_invariance(scale):
    data = fixture()
    data["cartesian_modes_per_sqrt_amu"] = (np.asarray(data["cartesian_modes_per_sqrt_amu"])*scale).tolist()
    assert score(data) == pytest.approx(score(fixture()))


def test_rigid_coordinate_rotation_and_translation():
    data = fixture()
    axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
    cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    angle = .77
    rot = np.eye(3)*np.cos(angle) + (1-np.cos(angle))*np.outer(axis,axis) + np.sin(angle)*cross
    data["geometry_angstrom"] = (np.asarray(data["geometry_angstrom"]) @ rot.T + [3., -9., 5.]).tolist()
    data["cartesian_modes_per_sqrt_amu"] = (np.asarray(data["cartesian_modes_per_sqrt_amu"]) @ rot.T).tolist()
    assert score(data) == pytest.approx(score(fixture()))


def test_free_atom_row_reordering():
    data = fixture(); order = [2, 0, 1]
    data["free_atom_indices"] = order
    data["free_masses_amu"] = np.array(data["free_masses_amu"])[order].tolist()
    data["cartesian_modes_per_sqrt_amu"] = np.array(data["cartesian_modes_per_sqrt_amu"])[:,order].tolist()
    assert score(data) == pytest.approx(score(fixture()))


def test_optimal_gradient_direction_has_unit_overlap():
    data = fixture()
    data["cartesian_modes_per_sqrt_amu"] = [[[1/12,0,0],[-2,0,0],[1/12,0,0]]]
    assert score(data) == pytest.approx(1.)


@pytest.mark.parametrize("direction", [[[1,2,3]]*3, [[0,1,0],[0,0,0],[0,-1,0]]])
def test_rigid_body_modes_have_zero_projection(direction):
    data = fixture(); data["cartesian_modes_per_sqrt_amu"] = [direction]
    assert score(data) == pytest.approx(0., abs=1e-25)


@pytest.mark.parametrize("mass", [0., -1., float("nan"), float("inf")])
def test_invalid_mass_rejected(mass):
    data = fixture(); data["free_masses_amu"][1] = mass
    with pytest.raises(ValueError): score(data)


def test_zero_negative_mode_rejected():
    data = fixture(); data["cartesian_modes_per_sqrt_amu"] = np.zeros((1,3,3)).tolist()
    with pytest.raises(ValueError): score(data)


def test_coincident_reaction_atoms_rejected():
    data = fixture(); data["geometry_angstrom"][0] = data["geometry_angstrom"][1]
    with pytest.raises(ValueError): score(data)


def test_frozen_reaction_atom_is_explicitly_unassessed():
    data = fixture(); data["free_atom_indices"] = [0,2]
    data["free_masses_amu"] = [12,12]
    data["cartesian_modes_per_sqrt_amu"] = [[[1,0,0],[0,0,0]]]
    out = transfer_projection(data, INDICES, SYMBOLS)
    assert not out["assessed"]


def test_overlap_is_not_hydrogen_participation():
    data = fixture()
    data["geometry_angstrom"][2] = [2.,0,0]  # Both carbons on one side.
    data["cartesian_modes_per_sqrt_amu"] = [[[1/12,0,0],[0,0,0],[-1/12,0,0]]]
    out = transfer_projection(data, INDICES, SYMBOLS)["negative_modes"][0]
    assert out["squared_mass_metric_overlap"] == pytest.approx(1.)
    assert out["hydrogen_mass_weighted_displacement_share"] == 0


def test_input_record_is_not_changed():
    data = fixture(); before = deepcopy(data)
    transfer_projection(data, INDICES, SYMBOLS)
    assert data == before
