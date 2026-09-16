"""Independent evidence-contract checks; all records are synthetic."""
from copy import deepcopy

import numpy as np
import pytest

from nanodesign import electronic_state as state
from e4_fixtures import build_snapshot as _build_snapshot


def build_snapshot(**kwargs):
    """Prove the unmodified control is valid before testing one changed field."""
    snapshot = _build_snapshot(**kwargs)
    state.validate_snapshot(snapshot)
    return snapshot


def assert_no_metrics(report):
    assert report["comparison_status"] != "compared"
    assert report["metrics"] is None
    assert report["electronic_state_identity_verified"] is False
    assert report["ground_state_verified"] is False


def test_distinct_guess_and_requested_resources_are_recorded_but_comparable():
    left = build_snapshot(call_id="e4-left")
    right = deepcopy(left)
    right["call_id"] = "e4-right"
    right["settings"].update(scf_initial_guess="atom", threads=2, memory_mb=512)
    right["scf_energy_hartree"] = -0.7
    report = state.compare_snapshots(left, right)
    assert report["comparison_status"] == "compared"
    assert report["input_digests"]["left"] != report["input_digests"]["right"]
    assert report["electronic_state_identity_verified"] is False
    assert report["ground_state_verified"] is False
    assert left["settings"]["scf_initial_guess"] == "minao"
    assert right["settings"]["scf_initial_guess"] == "atom"


@pytest.mark.parametrize("field", ["grid_level", "density_fit", "conv_tol", "dispersion"])
def test_two_equally_missing_physical_or_numerical_settings_do_not_compare(field):
    snapshot = build_snapshot()
    del snapshot["settings"][field]
    assert_no_metrics(state.compare_snapshots(snapshot, deepcopy(snapshot)))


@pytest.mark.parametrize("field,value", [("grid_level", True), ("density_fit", "false"),
                                         ("conv_tol", -1.0), ("dispersion", "unidentified")])
def test_two_equally_invalid_physical_or_numerical_settings_do_not_compare(field, value):
    snapshot = build_snapshot()
    snapshot["settings"][field] = value
    assert_no_metrics(state.compare_snapshots(snapshot, deepcopy(snapshot)))


@pytest.mark.parametrize("change", ["geometry", "basis", "ao_order", "ao_convention", "method", "version"])
def test_matching_shapes_do_not_bypass_context_binding(change):
    left, right = build_snapshot(), build_snapshot(call_id="e4-other")
    if change == "geometry":
        right["geometry"]["positions_bohr"][0][0] = 1e-5
    elif change == "basis":
        right["basis"]["shells"][0]["exponents"][0] = 1.01
    elif change == "ao_order":
        right["basis"]["ao_labels"][0], right["basis"]["ao_labels"][1] = right["basis"]["ao_labels"][1], right["basis"]["ao_labels"][0]
    elif change == "ao_convention":
        right["basis"]["representation"] = "cartesian"
    elif change == "method":
        right["settings"]["xc"] = "b3lyp"
    elif change == "version":
        right["producer"]["version"] = "different-synthetic-version"
    assert_no_metrics(state.compare_snapshots(left, right))


@pytest.mark.parametrize("bad", ["fractional_occupation", "wrong_electron_count", "nonorthogonal_orbitals",
                                  "indefinite_overlap", "nonfinite_coefficient", "identity_claim"])
def test_invalid_mathematical_evidence_yields_no_comparison(bad):
    snapshot = build_snapshot()
    if bad == "fractional_occupation":
        snapshot["source_occupations"]["alpha"][0] = 0.5
    elif bad == "wrong_electron_count":
        snapshot["channels"]["alpha"]["n_electrons"] = 1
    elif bad == "nonorthogonal_orbitals":
        snapshot["channels"]["alpha"]["occupied_coefficients"][0][0] *= 2
    elif bad == "indefinite_overlap":
        snapshot["ao_overlap"][0][0] = -2
    elif bad == "nonfinite_coefficient":
        snapshot["channels"]["alpha"]["occupied_coefficients"][0][0] = float("nan")
    elif bad == "identity_claim":
        snapshot["electronic_state_identity_verified"] = True
    assert_no_metrics(state.compare_snapshots(snapshot, deepcopy(snapshot)))


def test_save_read_is_lossless_and_existing_evidence_cannot_be_replaced(tmp_path):
    snapshot = build_snapshot()
    destination = tmp_path / "state.json"
    info = state.save_snapshot(destination, snapshot)
    original = destination.read_bytes()
    assert state.read_snapshot(destination) == snapshot
    assert info["size_bytes"] == len(original)
    with pytest.raises(state.ElectronicStateError):
        state.save_snapshot(destination, snapshot)
    assert destination.read_bytes() == original


@pytest.mark.parametrize("reverse", [False, True])
def test_absolute_overlap_tolerance_cannot_hide_large_metric_relative_change(reverse):
    # Both metrics meet the declared conditioning guard. Their entries differ
    # by only 1e-9, yet that doubles the metric of one occupied direction.
    left, right = build_snapshot(), build_snapshot(call_id="e4-metric-change")
    for snapshot, small_eigenvalue in ((left, 1e-9), (right, 2e-9)):
        snapshot["ao_overlap"] = np.diag([1.0, small_eigenvalue, 1.0]).tolist()
        snapshot["channels"]["alpha"]["occupied_coefficients"] = [
            [1.0, 0.0], [0.0, float(1.0 / np.sqrt(small_eigenvalue))], [0.0, 0.0]
        ]
        snapshot["channels"]["beta"]["occupied_coefficients"] = [[0.0], [0.0], [1.0]]
        state.validate_snapshot(snapshot)
    if reverse:
        left, right = right, left
    assert_no_metrics(state.compare_snapshots(left, right))
