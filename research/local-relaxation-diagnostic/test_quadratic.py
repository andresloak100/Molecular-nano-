"""Closed-form Cartesian quadratic diagnostics; no production or quantum imports."""

import json
import math

import numpy as np
import pytest

from diagnose import DiagnosticError, diagnose_quadratic


DEFAULT_CONTROLS = {
    "curvature_floor_ev_per_angstrom2": 1e-8,
    "max_condition_number": 1e8,
    "max_atom_shift_angstrom": 2.0,
    "expected_index": 0,
}


def diagnose(hessian, gradient, **controls):
    report = diagnose_quadratic(
        hessian, gradient, **(DEFAULT_CONTROLS | controls),
    )
    assert report["energy_error_bound_established"] is False
    assert report["actual_relaxation_performed"] is False
    # Refused calculations must also be portable as strict JSON evidence.
    json.dumps(report, allow_nan=False)
    return report


def codes(report):
    return {item["code"] for item in report["findings"]}


def assert_refused(report, code):
    assert report["status"] == "refused"
    assert code in codes(report)
    assert report["estimate"] is None


def assert_estimate(report, displacement, *, signed, absolute, positive, negative):
    estimate = report["estimate"]
    np.testing.assert_allclose(
        estimate["displacement_free_atoms_angstrom"], displacement, atol=1e-13,
    )
    assert estimate["signed_quadratic_energy_change_ev"] == pytest.approx(signed, abs=1e-13)
    assert estimate["absolute_mode_contribution_sum_ev"] == pytest.approx(absolute, abs=1e-13)
    assert estimate["positive_curvature_energy_change_ev"] == pytest.approx(positive, abs=1e-13)
    assert estimate["negative_curvature_energy_change_ev"] == pytest.approx(negative, abs=1e-13)
    assert estimate["linearized_gradient_residual_norm_ev_per_angstrom"] < 1e-12


def test_minimum_has_downhill_signed_correction_and_known_displacement():
    report = diagnose(np.diag([2.0, 4, 6]).tolist(), [[0.2, -0.4, 0]])

    assert report["status"] == "estimated"
    assert_estimate(report, [[-0.1, 0.1, 0]], signed=-0.03, absolute=0.03, positive=-0.03, negative=0)
    spectrum = report["spectrum"]
    np.testing.assert_allclose(spectrum["eigenvalues_ev_per_angstrom2"], [2, 4, 6])
    assert spectrum["negative_count"] == spectrum["soft_count"] == 0
    assert spectrum["absolute_condition_number"] == pytest.approx(3)


def test_saddle_energy_cancellation_does_not_erase_displacement_or_absolute_contributions():
    report = diagnose(
        np.diag([-2.0, 2, 4]).tolist(), [[0.2, 0.2, 0]], expected_index=1,
    )

    assert report["status"] == "estimated"
    assert_estimate(report, [[0.1, -0.1, 0]], signed=0, absolute=0.02, positive=-0.01, negative=0.01)
    assert report["estimate"]["max_atom_shift_angstrom"] == pytest.approx(math.sqrt(0.02))
    assert report["spectrum"]["negative_count"] == 1


def test_negative_curvature_correction_can_raise_energy_toward_saddle():
    report = diagnose(
        np.diag([-2.0, 2, 4]).tolist(), [[0.2, 0, 0]], expected_index=1,
    )

    assert report["status"] == "estimated"
    assert_estimate(report, [[0.1, 0, 0]], signed=0.01, absolute=0.01, positive=0, negative=0.01)


def test_equal_gradient_does_not_imply_equal_displacement_on_stiff_and_soft_surfaces():
    stiff = diagnose(np.diag([1.0, 2, 3]).tolist(), [[0.1, 0, 0]])
    soft = diagnose(np.diag([0.1, 2, 3]).tolist(), [[0.1, 0, 0]])

    assert stiff["status"] == soft["status"] == "estimated"
    assert_estimate(stiff, [[-0.1, 0, 0]], signed=-0.005, absolute=0.005, positive=-0.005, negative=0)
    assert_estimate(soft, [[-1.0, 0, 0]], signed=-0.05, absolute=0.05, positive=-0.05, negative=0)


def test_rotated_coupled_hessian_preserves_the_closed_form_quadratic_solution():
    # Orthogonal rotation with exact 3/5 and 4/5 entries; expected modal
    # displacement is (-0.1, 0.1, -0.1), rotated into Cartesian coordinates.
    rotation = np.array([[0.6, -0.8, 0], [0.8, 0.6, 0], [0, 0, 1]])
    # Supply the exactly symmetric Cartesian matrix explicitly, avoiding a
    # rounding asymmetry from computing two nominally equal off-diagonals.
    hessian = [[3.92, -1.44, 0], [-1.44, 3.08, 0], [0, 0, 7]]
    gradient = rotation @ np.array([0.2, -0.5, 0.7])
    report = diagnose(hessian, [gradient.tolist()])

    assert report["status"] == "estimated"
    assert_estimate(report, [[-0.14, -0.02, -0.1]], signed=-0.07, absolute=0.07, positive=-0.07, negative=0)
    np.testing.assert_allclose(report["spectrum"]["eigenvalues_ev_per_angstrom2"], [2, 5, 7], atol=1e-13)


def test_rms_atom_shift_uses_atom_vector_lengths_and_preserves_two_atom_mapping():
    report = diagnose(
        np.diag([2.0, 2, 2, 4, 4, 4]).tolist(),
        [[-0.2, 0, 0], [-1.2, 0, 0]],
    )

    assert report["status"] == "estimated"
    assert_estimate(report, [[0.1, 0, 0], [0.3, 0, 0]], signed=-0.19, absolute=0.19, positive=-0.19, negative=0)
    assert report["estimate"]["max_atom_shift_angstrom"] == pytest.approx(0.3)
    assert report["estimate"]["rms_atom_shift_angstrom"] == pytest.approx(math.sqrt(0.05))


def test_positive_definite_zero_gradient_is_a_zero_estimate():
    report = diagnose(np.diag([1.0, 2, 3]).tolist(), [[0, 0, 0]])

    assert report["status"] == "estimated"
    assert_estimate(report, [[0, 0, 0]], signed=0, absolute=0, positive=0, negative=0)
    assert report["estimate"]["max_atom_shift_angstrom"] == 0
    assert report["estimate"]["rms_atom_shift_angstrom"] == 0


@pytest.mark.parametrize("soft_eigenvalue", [0.0, 1e-10, -1e-10, 1e-6, -1e-6])
def test_curvature_at_or_below_floor_refuses_inversion_without_partial_estimate(soft_eigenvalue):
    report = diagnose(
        np.diag([soft_eigenvalue, 1, 2]).tolist(), [[0.1, 0, 0]],
        curvature_floor_ev_per_angstrom2=1e-6, max_condition_number=1e14,
    )

    assert_refused(report, "soft_curvature")
    assert report["spectrum"]["soft_count"] == 1


def test_singular_hessian_is_refused_even_when_gradient_is_zero():
    report = diagnose(np.diag([0.0, 1, 2]).tolist(), [[0, 0, 0]])

    assert_refused(report, "soft_curvature")


def test_excessive_absolute_condition_number_refuses_inversion():
    report = diagnose(
        np.diag([-1e-4, 1, 2]).tolist(), [[0.1, 0, 0]],
        expected_index=1, max_condition_number=1e4,
    )

    assert_refused(report, "ill_conditioned")
    assert report["spectrum"]["soft_count"] == 0
    assert report["spectrum"]["absolute_condition_number"] == pytest.approx(2e4)


def test_condition_number_exactly_at_limit_is_allowed():
    report = diagnose(
        np.diag([1.0, 2, 2]).tolist(), [[0.1, 0, 0]], max_condition_number=2,
    )

    assert report["status"] == "estimated"
    assert report["spectrum"]["absolute_condition_number"] == pytest.approx(2)


@pytest.mark.parametrize("hessian,expected_index,actual_index", [
    (np.diag([-2.0, 2, 4]).tolist(), 0, 1),
    (np.diag([2.0, 2, 4]).tolist(), 1, 0),
    (np.diag([-2.0, -3, 4]).tolist(), 1, 2),
])
def test_wrong_stationary_index_retains_estimate_but_is_inconclusive(hessian, expected_index, actual_index):
    report = diagnose(hessian, [[0.2, 0, 0]], expected_index=expected_index)

    assert report["status"] == "inconclusive"
    assert "index_mismatch" in codes(report)
    assert report["estimate"]["displacement_free_atoms_angstrom"] is not None
    assert report["spectrum"]["negative_count"] == actual_index


def test_displacement_limit_retains_unclipped_estimate_but_is_inconclusive():
    report = diagnose(
        np.diag([1.0, 2, 3]).tolist(), [[0.2, 0, 0]], max_atom_shift_angstrom=0.1,
    )

    assert report["status"] == "inconclusive"
    assert "displacement_limit_exceeded" in codes(report)
    assert_estimate(report, [[-0.2, 0, 0]], signed=-0.02, absolute=0.02, positive=-0.02, negative=0)
    assert report["estimate"]["max_atom_shift_angstrom"] == pytest.approx(0.2)


def test_index_and_displacement_findings_are_both_retained():
    report = diagnose(
        np.diag([-2.0, 2, 3]).tolist(), [[0.4, 0, 0]],
        expected_index=0, max_atom_shift_angstrom=0.1,
    )

    assert report["status"] == "inconclusive"
    assert {"index_mismatch", "displacement_limit_exceeded"} <= codes(report)


@pytest.mark.parametrize("hessian,gradient", [
    ([], []),
    ([[1, 0], [0, 1]], [[0, 0, 0]]),
    ([[1, 0, 0], [0, 1], [0, 0, 1]], [[0, 0, 0]]),
    ([[1, 0.1, 0], [0, 1, 0], [0, 0, 1]], [[0, 0, 0]]),
    (np.eye(3).tolist(), [[0, 0]]),
    (np.eye(3).tolist(), [[0, 0, 0], [0, 0, 0]]),
    ([[True, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 0, 0]]),
    (np.eye(3).tolist(), [[False, 0, 0]]),
    ([[float("nan"), 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 0, 0]]),
    ([[float("inf"), 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 0, 0]]),
    (np.eye(3).tolist(), [[0, float("nan"), 0]]),
    (np.eye(3).tolist(), [[0, 0, float("-inf")]]),
])
def test_malformed_hessian_or_gradient_raises_diagnostic_error(hessian, gradient):
    with pytest.raises(DiagnosticError):
        diagnose_quadratic(hessian, gradient, **DEFAULT_CONTROLS)


@pytest.mark.parametrize("controls", [
    {"curvature_floor_ev_per_angstrom2": True},
    {"curvature_floor_ev_per_angstrom2": 0},
    {"curvature_floor_ev_per_angstrom2": float("nan")},
    {"max_condition_number": True},
    {"max_condition_number": 0.5},
    {"max_condition_number": float("inf")},
    {"max_atom_shift_angstrom": True},
    {"max_atom_shift_angstrom": 0},
    {"max_atom_shift_angstrom": -1},
    {"expected_index": True},
    {"expected_index": 1.0},
    {"expected_index": -1},
    {"expected_index": 2},
])
def test_malformed_controls_raise_diagnostic_error(controls):
    with pytest.raises(DiagnosticError):
        diagnose_quadratic(
            np.eye(3).tolist(), [[0, 0, 0]], **(DEFAULT_CONTROLS | controls),
        )
