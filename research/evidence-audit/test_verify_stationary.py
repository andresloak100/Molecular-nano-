"""Adversarial saved-record tests; all force values below are declared synthetic.

No production, ASE, or electronic-structure calculation creates these fixtures.
Synthetic JSON is written only into pytest temporary directories, never evidence.
"""
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


VERIFIER = Path(__file__).with_name("verify_stationary.py")
SPEC = importlib.util.spec_from_file_location("independent_stationary_verifier", VERIFIER)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
# Independent SI expression; deliberately do not import the verifier's constant.
CM1_FACTOR = math.sqrt(1.602176634e-19 / 1.66053906892e-27) / (2 * math.pi * 299792458 * 1e-8)


def install_curvature(record, raw):
    """Derive a fixture's saved summaries from a declared mathematical matrix."""
    symmetric = (raw + raw.T) / 2
    mass = np.repeat(record["free_masses_amu"], 3)
    weighted = symmetric / np.sqrt(np.outer(mass, mass))
    eigenvalues, eigenvectors = np.linalg.eigh(weighted)
    modes = eigenvectors.T / np.sqrt(mass)
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * CM1_FACTOR
    threshold = record["settings"]["frequency_tolerance_cm1"]
    labels = ["resolved_negative" if f < -threshold else "resolved_positive" if f > threshold
              else "unresolved_near_zero" for f in frequencies]
    asymmetry = raw - raw.T
    norm = np.linalg.norm(raw)
    negative = labels.count("resolved_negative")
    stationary = record["stationary_within_force_tolerance"]
    classification = "not_stationary_within_force_tolerance" if not stationary else (
        "stationary_point_with_no_resolved_negative_modes" if negative == 0 else
        "first_order_saddle_candidate" if negative == 1 else "higher_order_saddle_candidate")
    record.update({
        "hessian_ev_per_angstrom2": symmetric.tolist(),
        "hessian_asymmetry_max_abs_ev_per_angstrom2": float(np.abs(asymmetry).max()),
        "hessian_asymmetry_relative_frobenius": float(np.linalg.norm(asymmetry) / norm) if norm else 0.,
        "hessian_symmetrized": True,
        "mass_weighted_eigenvalues_ev_per_angstrom2_amu": eigenvalues.tolist(),
        "cartesian_modes_per_sqrt_amu": modes.reshape(6, 2, 3).tolist(),
        "frequencies_cm1": frequencies.tolist(),
        "negative_mode_count": negative,
        "resolved_positive_mode_count": labels.count("resolved_positive"),
        "unresolved_near_zero_mode_count": labels.count("unresolved_near_zero"),
        "raw_negative_eigenvalue_count": int((eigenvalues < 0).sum()),
        "mode_classifications": labels,
        "classification": classification,
    })


def synthetic_record(*, rounded=False, degenerate=False, zero=False):
    """F(x)=B-K(x-x0): a test field, explicitly not a molecular potential.

    Middle atom is fixed; free atoms have masses 2 and 5 amu. K includes
    nonconservative coupling, so an audit must reconstruct before symmetrizing.
    Expected raw columns use the analytic linear map and realized step ratio,
    rather than reproducing the verifier's subtraction of saved force arrays.
    """
    reference = np.array([[.125, -.75, 1.5], [2., -1., 3.], [3.75, .625, -.25]])
    if rounded:
        reference[0, 0] = 1e7
    if degenerate or zero:
        reference[[0, 2]] = 0.
    free, frozen = [0, 2], [1]
    coordinates = [0, 1, 2, 6, 7, 8]
    mass = np.repeat([2., 5.], 3)
    matrix = np.diag([-4., 5., 6., 11., 12., 13., 7., 8., 9.])
    matrix[0, 6] = matrix[6, 0] = .6
    matrix[1, 7] = matrix[7, 1] = -.4
    matrix[0, 7], matrix[7, 0] = .9, -.3
    matrix[0, 3] = matrix[3, 0] = 2.
    if degenerate:
        matrix[np.ix_(coordinates, coordinates)] = np.diag(mass * [1., 1., 2., 3., 4., 5.])
    if zero:
        matrix[:] = 0.
    baseline = np.array([[.004, -.005, .002], [20., -11., 7.], [.007, .002, -.003]])
    step = .005
    displacements, realized_widths = [], []
    for column, coordinate in enumerate(coordinates):
        atom, axis = divmod(coordinate, 3)
        offsets = []
        for sign in (1, -1):
            requested = sign * step
            displaced = reference[atom, axis] + requested
            actual = displaced - reference[atom, axis]
            delta = np.zeros(9)
            delta[coordinate] = actual
            forces = baseline - (matrix @ delta).reshape(3, 3)
            displacements.append({
                "atom_index": atom, "axis": axis,
                "requested_offset_angstrom": requested,
                "actual_offset_angstrom": float(actual),
                "displaced_coordinate_angstrom": float(displaced),
                "relative_step_representation_error": float(abs(actual - requested) / step),
                "forces_ev_per_angstrom": forces.tolist(),
                "calculation_call_id": f"synthetic-only-{2 * column + (1 if sign > 0 else 2)}",
            })
            offsets.append(actual)
        realized_widths.append(offsets[0] - offsets[1])
    free_forces = baseline[free]
    force_max = float(np.linalg.norm(free_forces, axis=1).max())
    record = {
        "fixture_notice": "Synthetic linear force field for software tests; no scientific evidence",
        "settings": {"step_angstrom": step, "force_tolerance_ev_per_angstrom": .03,
                     "frequency_tolerance_cm1": 20., "max_free_coordinates": 120},
        "free_atom_indices": free, "frozen_atom_indices": frozen,
        "free_coordinate_count": 6, "force_requests": 13,
        "free_masses_amu": [2., 5.],
        "coordinate_order": "x,y,z for each atom in free_atom_indices",
        "free_gradient_ev_per_angstrom": (-free_forces).tolist(),
        "free_force_max_ev_per_angstrom": force_max,
        "free_force_cartesian_rms_ev_per_angstrom": float(np.sqrt(np.mean(free_forces ** 2))),
        "stationary_within_force_tolerance": True,
        "transition_state_verified": False,
        "external_modes_removed": False,
        "coordinate_space": "free Cartesian coordinates with fixed anchors",
        "geometry_angstrom": reference.tolist(),
        "finite_difference_evidence": {
            "units": {"length": "angstrom", "force": "eV/angstrom"},
            "reference_positions_angstrom": reference.tolist(),
            "baseline_forces_ev_per_angstrom": baseline.tolist(),
            "baseline_calculation_call_id": "synthetic-only-baseline",
            "force_scope": "Unconstrained forces on all atoms, including fixed anchors",
            "axis_convention": "0=x, 1=y, 2=z; atom_index is zero-based",
            "displacement_order": "Plus then minus for each free Cartesian coordinate in coordinate_order",
            "displacements": displacements,
        },
        "displacement_resolution": {
            "checked": True,
            "max_relative_step_representation_error_allowed": 1e-6,
            "max_relative_step_representation_error_observed": max(
                row["relative_step_representation_error"] for row in displacements),
        },
    }
    raw = matrix[np.ix_(coordinates, coordinates)] * (np.array(realized_widths) / (2 * step))
    install_curvature(record, raw)
    return json.loads(json.dumps(record, allow_nan=False))


def assert_failed(record, *, code=None):
    report = MODULE.verify_stationary(record)
    assert report["status"] == "failed", report
    assert report["numerical_reconstruction_verified"] is False
    assert report["scientific_model_validated"] is False
    if code:
        assert report["findings"][0]["code"] == code, report
    return report


def test_roundtrip_anchored_unequal_mass_asymmetric_record_and_input_immutability():
    record = synthetic_record()
    before = deepcopy(record)
    report = MODULE.verify_stationary(record)
    assert report["status"] == "passed", report
    assert record == before
    assert report["details"]["free_coordinates"] == 6
    assert report["details"]["displacements"] == 12
    assert report["details"]["hessian_asymmetry_max_abs_ev_per_angstrom2"] == pytest.approx(1.2)
    assert report["details"]["stationary_within_recorded_force_tolerance"] is True
    assert report["details"]["call_ids"]["external_event_linkage"] == "not_checked"
    assert all(report[key] is False for key in ("scientific_model_validated", "electronic_state_verified",
                                               "transition_state_verified", "connectivity_verified"))


def test_arbitrary_mode_sign_and_degenerate_basis_rotation_are_valid():
    record = synthetic_record(degenerate=True)
    modes = np.array(record["cartesian_modes_per_sqrt_amu"])
    angle = .37
    first, second = modes[:2].copy()
    modes[0] = math.cos(angle) * first + math.sin(angle) * second
    modes[1] = -math.sin(angle) * first + math.cos(angle) * second
    modes[4] *= -1
    record["cartesian_modes_per_sqrt_amu"] = modes.tolist()
    report = MODULE.verify_stationary(record)
    assert report["status"] == "passed", report


def test_zero_raw_hessian_has_zero_relative_asymmetry():
    report = MODULE.verify_stationary(synthetic_record(zero=True))
    assert report["status"] == "passed", report
    assert report["details"]["hessian_asymmetry_relative_frobenius"] == 0
    assert report["details"]["unresolved_mode_count"] == 6


def test_valid_rounded_stencil_uses_requested_denominator():
    record = synthetic_record(rounded=True)
    row = record["finite_difference_evidence"]["displacements"][0]
    assert row["actual_offset_angstrom"] != row["requested_offset_angstrom"]
    assert 1e-8 < row["relative_step_representation_error"] < 1e-6
    assert MODULE.verify_stationary(record)["status"] == "passed"


def test_actual_offset_denominator_is_rejected_even_with_consistent_eigenpairs():
    record = synthetic_record(rounded=True)
    rows = record["finite_difference_evidence"]["displacements"]
    wrong_raw = np.empty((6, 6))
    for column in range(6):
        plus, minus = rows[2 * column:2 * column + 2]
        fp = np.array(plus["forces_ev_per_angstrom"])[[0, 2]].ravel()
        fm = np.array(minus["forces_ev_per_angstrom"])[[0, 2]].ravel()
        wrong_raw[:, column] = -(fp - fm) / (plus["actual_offset_angstrom"] - minus["actual_offset_angstrom"])
    install_curvature(record, wrong_raw)
    assert_failed(record, code="arithmetic_mismatch")


def test_antisymmetric_force_corruption_cannot_hide_behind_unchanged_symmetric_hessian():
    record = synthetic_record()
    rows = record["finite_difference_evidence"]["displacements"]
    h = record["settings"]["step_angstrom"]
    # Add +2 to raw[0,4] and -2 to raw[4,0]: symmetric Hessian is unchanged.
    rows[8]["forces_ev_per_angstrom"][0][0] -= 4 * h
    rows[0]["forces_ev_per_angstrom"][2][1] += 4 * h
    report = assert_failed(record, code="arithmetic_mismatch")
    assert "asymmetry" in report["findings"][0]["message"]


@pytest.mark.parametrize("corruption", ["include_anchor", "duplicate_free", "overlap", "uncovered", "boolean", "out_of_range"])
def test_invalid_or_wrong_free_atom_partition_is_rejected(corruption):
    record = synthetic_record()
    if corruption == "include_anchor":
        record["free_atom_indices"], record["frozen_atom_indices"] = [0, 1], [2]
    elif corruption == "duplicate_free":
        record["free_atom_indices"] = [0, 0]
    elif corruption == "overlap":
        record["frozen_atom_indices"] = [0, 1]
    elif corruption == "uncovered":
        record["frozen_atom_indices"] = []
    elif corruption == "boolean":
        record["free_atom_indices"] = [False, 2]
    else:
        record["free_atom_indices"] = [0, 3]
    assert_failed(record)


@pytest.mark.parametrize("target", ["baseline", "displaced"])
def test_full_force_archive_cannot_drop_anchor_rows(target):
    record = synthetic_record()
    evidence = record["finite_difference_evidence"]
    if target == "baseline":
        evidence["baseline_forces_ev_per_angstrom"].pop(1)
    else:
        for row in evidence["displacements"]:
            row["forces_ev_per_angstrom"].pop(1)
    assert_failed(record, code="invalid_shape")


def test_swapped_unequal_masses_rejected():
    record = synthetic_record()
    record["free_masses_amu"].reverse()
    assert_failed(record, code="arithmetic_mismatch")


@pytest.mark.parametrize("corruption", ["missing", "duplicate_plus", "axis", "atom", "requested_sign"])
def test_stencil_order_and_completeness(corruption):
    record = synthetic_record()
    rows = record["finite_difference_evidence"]["displacements"]
    if corruption == "missing":
        rows.pop()
    elif corruption == "duplicate_plus":
        rows[1] = deepcopy(rows[0])
    elif corruption == "axis":
        rows[0]["axis"] = 2
    elif corruption == "atom":
        rows[0]["atom_index"] = 1
    else:
        rows[0]["requested_offset_angstrom"] *= -1
    assert_failed(record)


@pytest.mark.parametrize("corruption", ["zero_actual", "forged_error", "forged_maximum", "relaxed_limit", "unresolved_coordinate"])
def test_resolution_metadata_cannot_hide_invalid_offsets(corruption):
    record = synthetic_record(rounded=True)
    evidence = record["finite_difference_evidence"]
    row = evidence["displacements"][0]
    if corruption == "zero_actual":
        row["actual_offset_angstrom"] = 0.
    elif corruption == "forged_error":
        row["relative_step_representation_error"] = 0.
    elif corruption == "forged_maximum":
        record["displacement_resolution"]["max_relative_step_representation_error_observed"] = 0.
    elif corruption == "relaxed_limit":
        record["displacement_resolution"]["max_relative_step_representation_error_allowed"] = .01
    else:
        # Algebraically consistent metadata for a lost binary64 displacement.
        evidence["reference_positions_angstrom"][0][0] = 1e16
        record["geometry_angstrom"][0][0] = 1e16
        row["displaced_coordinate_angstrom"] = 1e16
        row["actual_offset_angstrom"] = 0.
        row["relative_step_representation_error"] = 1.
    assert_failed(record)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, "0.0"])
def test_invalid_force_values_are_rejected(value):
    record = synthetic_record()
    record["finite_difference_evidence"]["displacements"][0]["forces_ev_per_angstrom"][1][1] = value
    assert_failed(record, code="nonfinite_json" if type(value) is float and not math.isfinite(value) else "invalid_array")


@pytest.mark.parametrize("field,value", [("step_angstrom", True), ("step_angstrom", 0.),
                                         ("force_tolerance_ev_per_angstrom", -1.),
                                         ("frequency_tolerance_cm1", float("inf"))])
def test_invalid_scalar_settings_are_rejected(field, value):
    record = synthetic_record()
    record["settings"][field] = value
    assert_failed(record, code="nonfinite_json" if type(value) is float and not math.isfinite(value) else "invalid_number")


@pytest.mark.parametrize("corruption", ["wrong_gradient", "wrong_rms", "promoted_ts", "classification", "mode_count", "units", "geometry", "mode_scale"])
def test_saved_summary_corruptions(corruption):
    record = synthetic_record()
    if corruption == "wrong_gradient":
        record["free_gradient_ev_per_angstrom"][0][0] *= -1
    elif corruption == "wrong_rms":
        record["free_force_cartesian_rms_ev_per_angstrom"] = 20.
    elif corruption == "promoted_ts":
        record["transition_state_verified"] = True
    elif corruption == "classification":
        record["classification"] = "stationary_point_with_no_resolved_negative_modes"
    elif corruption == "mode_count":
        record["negative_mode_count"] = True
    elif corruption == "units":
        record["finite_difference_evidence"]["units"]["length"] = "bohr"
    elif corruption == "geometry":
        record["geometry_angstrom"][1][0] += 1.
    else:
        record["cartesian_modes_per_sqrt_amu"][0][0][0] *= 2
    assert_failed(record)


def test_missing_legacy_evidence_is_unavailable_not_reconstructed():
    record = synthetic_record()
    del record["finite_difference_evidence"]
    report = MODULE.verify_stationary(record)
    assert report["status"] == "unavailable"
    assert report["numerical_reconstruction_verified"] is False
    assert report["findings"][0]["code"] == "raw_force_evidence_unavailable"


def test_null_call_ids_allow_math_but_never_claim_external_linkage():
    record = synthetic_record()
    evidence = record["finite_difference_evidence"]
    evidence["baseline_calculation_call_id"] = None
    for row in evidence["displacements"]:
        row["calculation_call_id"] = None
    report = MODULE.verify_stationary(record)
    assert report["status"] == "passed", report
    assert report["details"]["call_ids"] == {
        "present": 0, "missing_or_null": 13, "duplicated_values": [], "external_event_linkage": "not_checked"}


def test_duplicate_call_ids_are_disclosed_without_rejecting_valid_math():
    record = synthetic_record()
    record["finite_difference_evidence"]["displacements"][1]["calculation_call_id"] = "synthetic-only-1"
    report = MODULE.verify_stationary(record)
    assert report["status"] == "passed", report
    assert report["details"]["call_ids"]["duplicated_values"] == ["synthetic-only-1"]
    assert report["details"]["call_ids"]["external_event_linkage"] == "not_checked"


def test_nonstring_call_id_is_invalid():
    record = synthetic_record()
    record["finite_difference_evidence"]["baseline_calculation_call_id"] = 123
    assert_failed(record, code="invalid_call_id")


def write_record(tmp_path, data):
    path = tmp_path / "synthetic-test-record.json"
    path.write_text(json.dumps(data, allow_nan=False))
    return path


def run_cli(path, *args):
    return subprocess.run([sys.executable, "-B", str(VERIFIER), str(path), *map(str, args)],
                          capture_output=True, text=True, check=False)


def test_cli_source_hash_read_only_and_no_synthetic_science_claim(tmp_path):
    path = write_record(tmp_path, synthetic_record())
    original = path.read_bytes()
    output = tmp_path / "verification.json"
    result = run_cli(path, "--output", output)
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert json.loads(output.read_text()) == report
    assert path.read_bytes() == original
    assert report["source"]["sha256"] == hashlib.sha256(original).hexdigest()
    assert report["audit_implementation_sha256"] == hashlib.sha256(VERIFIER.read_bytes()).hexdigest()
    assert report["scientific_model_validated"] is False


@pytest.mark.parametrize("existing_is_input", [False, True])
def test_cli_never_overwrites_output_or_input(tmp_path, existing_is_input):
    path = write_record(tmp_path, synthetic_record())
    output = path if existing_is_input else tmp_path / "existing.json"
    if not existing_is_input:
        output.write_text("existing report must survive")
    original_input, original_output = path.read_bytes(), output.read_bytes()
    result = run_cli(path, "--output", output)
    assert result.returncode != 0
    assert path.read_bytes() == original_input
    assert output.read_bytes() == original_output


def test_cli_missing_legacy_evidence_exits_unavailable(tmp_path):
    record = synthetic_record()
    del record["finite_difference_evidence"]
    result = run_cli(write_record(tmp_path, record))
    assert result.returncode == 2, result.stderr
    assert json.loads(result.stdout)["status"] == "unavailable"


def test_cli_corrupt_record_exits_failed(tmp_path):
    record = synthetic_record()
    record["free_coordinate_count"] = True
    result = run_cli(write_record(tmp_path, record))
    assert result.returncode == 1, result.stderr
    assert json.loads(result.stdout)["status"] == "failed"


def test_cli_invalid_json_reports_failure_without_touching_file(tmp_path):
    path = tmp_path / "invalid.json"
    path.write_text('{"unfinished":')
    original = path.read_bytes()
    result = run_cli(path)
    assert result.returncode == 1
    assert json.loads(result.stdout)["status"] == "failed"
    assert path.read_bytes() == original


@pytest.mark.parametrize("field", ["stationary", "transition_state_characterization"])
def test_wrapper_block_reconstruction_does_not_certify_parent_completion(tmp_path, field):
    path = write_record(tmp_path, {"status": "failed", field: synthetic_record()})
    report = MODULE.inspect_file(path)
    assert report["status"] == "passed", report
    assert report["source"]["selected_field"] == field
    assert report["source"]["parent_status"] == "failed"
    assert report["source"]["parent_completion_assessed"] is False
    assert report["scientific_model_validated"] is False


@pytest.mark.parametrize("ambiguity", ["two_wrappers", "root_and_wrapper"])
def test_ambiguous_characterization_sources_are_rejected(tmp_path, ambiguity):
    record = synthetic_record()
    if ambiguity == "two_wrappers":
        data = {"stationary": record, "transition_state_characterization": deepcopy(record)}
    else:
        data = dict(record, stationary=deepcopy(record))
    report = MODULE.inspect_file(write_record(tmp_path, data))
    assert report["status"] == "failed"
    assert report["findings"][0]["code"] == "ambiguous_record"


def test_null_characterization_is_unavailable(tmp_path):
    report = MODULE.inspect_file(write_record(tmp_path, {"status": "failed", "stationary": None}))
    assert report["status"] == "unavailable"
    assert report["numerical_reconstruction_verified"] is False
    assert report["source"]["parent_completion_assessed"] is False


@pytest.mark.parametrize("token", ["NaN", "Infinity", "1e999"])
def test_cli_rejects_nonfinite_metadata_including_overflowed_json_exponent(tmp_path, token):
    path = write_record(tmp_path, synthetic_record())
    # The corruption lives in metadata, outside every force/mode array.
    text = path.read_text()
    path.write_text(text[:-1] + ', "synthetic_metadata_corruption": ' + token + '}')
    original = path.read_bytes()
    result = run_cli(path)
    assert result.returncode == 1, result.stderr
    report = json.loads(result.stdout)
    assert report["status"] == "failed"
    assert report["findings"][0]["code"] == "nonfinite_json"
    assert path.read_bytes() == original
