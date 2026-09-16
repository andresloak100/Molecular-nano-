"""Synthetic finite-difference evidence checks; no quantum solver is invoked."""

from copy import deepcopy
import hashlib
import importlib.util
from pathlib import Path
import sys

import pytest

from nanodesign.quantum import QuantumSettings


_SPEC = importlib.util.spec_from_file_location(
    "force_energy_consistency_under_test", Path(__file__).with_name("consistency.py")
)
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)
analyze_stencils = _MODULE.analyze_stencils


def quadratic(x):
    return x * x


def quadratic_force(x):
    return -2.0 * x


def cubic(x):
    return x * x * x


def cubic_force(x):
    return -3.0 * x * x


def evaluation(x, identity, *, energy=quadratic, force=quadratic_force):
    return {
        "context": {
            "atomic_numbers": [1, 6], "fixed_indices": [],
            "pbc": [False, False, False],
            "quantum_settings": QuantumSettings(spin=1, threads=1).to_dict(),
            "software_versions": {"synthetic": "1"},
            "resolved_numerics": None,
            "energy_scope": "total_potential_energy",
            "units": {"length": "angstrom", "energy": "eV", "force": "eV/angstrom"},
            "force_scope": "raw_unconstrained",
        },
        "positions_angstrom": [[x, 0.0, 0.0], [2.0, 0.0, 0.0]],
        "source_sha256": hashlib.sha256(identity.encode()).hexdigest(),
        "source_record_id": identity,
        "calculation_call_id": "call-" + identity,
        "status": "completed", "scf_converged": True, "gradient_completed": True,
        "energy_ev": energy(x),
        "forces_ev_per_angstrom": [[force(x), 0.0, 0.0], [0.0, 0.0, 0.0]],
        "state_evidence_id": None,
    }


def stencil(hminus=0.125, hplus=None, *, energy=quadratic, force=quadratic_force):
    hplus = hminus if hplus is None else hplus
    return {
        "atom_index": 0, "axis": 0,
        "requested_minus_offset_angstrom": -hminus,
        "requested_plus_offset_angstrom": hplus,
        "reference": evaluation(0.5, "reference", energy=energy, force=force),
        "minus": evaluation(0.5 - hminus, f"minus-{hminus}", energy=energy, force=force),
        "plus": evaluation(0.5 + hplus, f"plus-{hplus}", energy=energy, force=force),
    }


def document(*stencils):
    return {"schema_version": 1, "stencils": list(stencils) if stencils else [stencil()]}


def points(item):
    return [item[role] for role in ("reference", "minus", "plus")]


def assert_invalid(value):
    report = analyze_stencils(value)
    assert report["status"] == "invalid"
    assert report["findings"]
    return report


def assert_no_certification(report):
    for name in (
        "electronic_state_identity_verified", "electronic_branch_continuity_verified",
        "scientific_model_validated", "numerical_convergence_certified",
    ):
        assert report[name] is False


def test_quadratic_force_has_correct_negative_gradient_sign():
    report = analyze_stencils(document(), force_tolerance_ev_per_angstrom=1e-10)
    assert report["status"] == "comparison_produced"
    row = report["comparisons"][0]
    assert row["analytic_force_ev_per_angstrom"] == -1.0
    assert row["finite_difference_force_ev_per_angstrom"] == pytest.approx(-1.0)
    assert row["residual_ev_per_angstrom"] == pytest.approx(0.0, abs=1e-12)
    assert row["absolute_residual_ev_per_angstrom"] == pytest.approx(0.0, abs=1e-12)
    assert row["within_declared_tolerance"] is True
    assert row["estimator"]
    assert_no_certification(report)


def test_nonuniform_quadratic_stencil_uses_reference_energy_and_actual_distances():
    report = analyze_stencils(document(stencil(0.125, 0.25)))
    assert report["status"] == "comparison_produced"
    row = report["comparisons"][0]
    # A naive secant between endpoints would incorrectly give -1.125.
    assert row["finite_difference_force_ev_per_angstrom"] == pytest.approx(-1.0)
    assert row["actual_minus_distance_angstrom"] == 0.125
    assert row["actual_plus_distance_angstrom"] == 0.25
    assert row["within_declared_tolerance"] is None


def test_cubic_step_sensitivity_keeps_one_identical_shared_baseline():
    first = stencil(0.125, energy=cubic, force=cubic_force)
    second = stencil(0.25, energy=cubic, force=cubic_force)
    assert first["reference"] == second["reference"]
    report = analyze_stencils(document(first, second))
    assert report["status"] == "comparison_produced"
    rows = report["comparisons"]
    assert [row["finite_difference_force_ev_per_angstrom"] for row in rows] == pytest.approx(
        [-0.75 - 0.125**2, -0.75 - 0.25**2]
    )
    group = report["step_sensitivity"][0]
    assert group["atom_index"] == group["axis"] == 0
    assert group["distinct_stencils"] == 2
    assert group["finite_difference_force_spread_ev_per_angstrom"] == pytest.approx(0.25**2 - 0.125**2)
    assert group["analytic_force_spread_ev_per_angstrom"] == pytest.approx(0.0)
    assert_no_certification(report)


def test_inconsistent_analytic_force_reports_signed_residual_and_failed_tolerance():
    item = stencil()
    item["reference"]["forces_ev_per_angstrom"][0][0] = -0.8
    report = analyze_stencils(document(item), force_tolerance_ev_per_angstrom=0.05)
    assert report["status"] == "comparison_produced"
    row = report["comparisons"][0]
    assert row["residual_ev_per_angstrom"] == pytest.approx(-0.2)
    assert row["absolute_residual_ev_per_angstrom"] == pytest.approx(0.2)
    assert row["within_declared_tolerance"] is False


def test_common_energy_offset_does_not_change_force_comparison():
    original = document()
    shifted = deepcopy(original)
    for point in points(shifted["stencils"][0]):
        point["energy_ev"] += 1024.0
    before = analyze_stencils(original)["comparisons"][0]
    after = analyze_stencils(shifted)["comparisons"][0]
    assert after["finite_difference_force_ev_per_angstrom"] == before["finite_difference_force_ev_per_angstrom"]
    assert after["residual_ev_per_angstrom"] == before["residual_ev_per_angstrom"]


@pytest.mark.parametrize("missing", ["point", "energy"])
def test_missing_point_or_energy_is_unavailable(missing):
    item = stencil()
    if missing == "point":
        item["plus"] = None
    else:
        item["plus"]["energy_ev"] = None
    report = analyze_stencils(document(item))
    assert report["status"] == "unavailable"
    assert report["comparisons"] == []
    assert_no_certification(report)


@pytest.mark.parametrize("status", ["failed", "not_run"])
def test_failed_or_unrun_point_cannot_supply_energy_even_when_number_is_present(status):
    item = stencil()
    item["plus"]["status"] = status
    report = analyze_stencils(document(item))
    assert report["status"] == "unavailable"
    assert report["comparisons"] == []


def test_partial_report_preserves_available_comparison_without_imputing_missing_one():
    second = stencil(0.25)
    second["plus"] = None
    report = analyze_stencils(document(stencil(), second))
    assert report["status"] == "partial"
    assert len(report["comparisons"]) == 1
    assert_no_certification(report)


@pytest.mark.parametrize("mismatch", ["settings", "atom_order", "units"])
def test_changed_method_atom_order_or_units_is_invalid(mismatch):
    item = stencil()
    context = item["plus"]["context"]
    if mismatch == "settings":
        context["quantum_settings"]["basis"] = "sto-3g"
    elif mismatch == "atom_order":
        context["atomic_numbers"] = [6, 1]
    else:
        context["units"]["energy"] = "Hartree"
    assert_invalid(document(item))


@pytest.mark.parametrize("invalid", ["boolean_energy", "nan_energy", "boolean_coordinate"])
def test_invalid_numeric_evidence_returns_invalid_without_raising(invalid):
    item = stencil()
    if invalid == "boolean_energy":
        item["plus"]["energy_ev"] = True
    elif invalid == "nan_energy":
        item["plus"]["energy_ev"] = float("nan")
    else:
        item["reference"]["positions_angstrom"][0][1] = False
    assert_invalid(document(item))


def test_duplicate_calculation_call_cannot_describe_different_displacements():
    item = stencil()
    item["plus"]["calculation_call_id"] = item["minus"]["calculation_call_id"]
    assert_invalid(document(item))


def test_identical_baseline_call_can_be_preserved_in_another_source_snapshot():
    first, second = stencil(), stencil(0.25)
    second["reference"]["source_record_id"] = "reference-in-resumed-checkpoint"
    second["reference"]["source_sha256"] = hashlib.sha256(b"resumed-checkpoint").hexdigest()
    report = analyze_stencils(document(first, second))
    assert report["status"] == "comparison_produced"
    assert len(report["comparisons"]) == 2


def test_reusing_baseline_source_identity_with_different_content_is_invalid():
    first, second = stencil(), stencil(0.25)
    second["reference"]["energy_ev"] += 0.1
    assert_invalid(document(first, second))


def test_masked_force_scope_is_not_a_raw_energy_gradient_comparison():
    item = stencil()
    for point in points(item):
        point["context"]["force_scope"] = "constrained_masked"
    assert_invalid(document(item))


def test_raw_anchor_load_is_preserved_while_only_free_coordinate_is_displaced():
    item = stencil()
    for point in points(item):
        point["context"]["fixed_indices"] = [1]
        point["forces_ev_per_angstrom"][1] = [123.0, -8.0, 3.0]
    value = document(item)
    original = deepcopy(value)
    report = analyze_stencils(value)
    assert report["status"] == "comparison_produced"
    assert report["comparisons"][0]["analytic_force_ev_per_angstrom"] == -1.0
    assert value == original
    assert value["stencils"][0]["reference"]["forces_ev_per_angstrom"][1] == [123.0, -8.0, 3.0]


def test_equal_opaque_state_ids_do_not_verify_branch_continuity():
    item = stencil()
    for point in points(item):
        point["state_evidence_id"] = "same-opaque-state-label"
    report = analyze_stencils(document(item))
    assert report["status"] == "comparison_produced"
    assert_no_certification(report)


@pytest.mark.parametrize("movement", ["other_coordinate", "fixed_coordinate"])
def test_unrequested_or_frozen_coordinate_motion_is_invalid(movement):
    item = stencil()
    if movement == "other_coordinate":
        item["plus"]["positions_angstrom"][1][2] += 0.01
    else:
        for point in points(item):
            point["context"]["fixed_indices"] = [0]
    assert_invalid(document(item))


def test_unrepresentable_requested_displacement_is_invalid():
    item = stencil(1e-30)
    assert item["minus"]["positions_angstrom"] == item["reference"]["positions_angstrom"]
    assert_invalid(document(item))


def test_boolean_quantum_integer_cannot_alias_valid_setting():
    item = stencil()
    item["plus"]["context"]["quantum_settings"]["charge"] = False
    assert_invalid(document(item))


def test_finite_offsets_with_overflowing_combined_span_are_invalid():
    linear_energy = lambda x: x
    linear_force = lambda x: -1.0
    item = {
        "atom_index": 0, "axis": 0,
        "requested_minus_offset_angstrom": -1e308,
        "requested_plus_offset_angstrom": 1e308,
        "reference": evaluation(0.0, "overflow-reference", energy=linear_energy, force=linear_force),
        "minus": evaluation(-1e308, "overflow-minus", energy=linear_energy, force=linear_force),
        "plus": evaluation(1e308, "overflow-plus", energy=linear_energy, force=linear_force),
    }
    # Every supplied scalar is finite, but a+b is not representable; zeroed
    # weights from division by an infinite span must not manufacture a force.
    assert_invalid(document(item))


def test_rounding_resolution_warning_is_separate_from_literal_residual_tolerance():
    item = stencil(0.125)
    for point in points(item):
        point["energy_ev"] = 1e16
        point["forces_ev_per_angstrom"][0][0] = 0.0
    report = analyze_stencils(document(item), force_tolerance_ev_per_angstrom=1e-3)
    assert report["status"] == "comparison_produced"
    row = report["comparisons"][0]
    assert row["finite_difference_force_ev_per_angstrom"] == 0.0
    assert row["residual_ev_per_angstrom"] == 0.0
    assert row["within_declared_tolerance"] is True
    assert row["binary64_input_rounding_sensitivity_ev_per_angstrom"] == pytest.approx(8.0)
    assert row["input_rounding_exceeds_declared_tolerance"] is True
    assert row["numerical_resolution_findings"]
    assert_no_certification(report)


@pytest.mark.parametrize("field", ["source_sha256", "source_record_id"])
def test_nullable_source_identity_is_unavailable_instead_of_invented(field):
    item = stencil()
    item["plus"][field] = None
    report = analyze_stencils(document(item))
    assert report["status"] == "unavailable"
    assert report["comparisons"] == []
    assert_no_certification(report)


@pytest.mark.parametrize("settings", [None, {"basis": "sto-3g"}])
def test_null_or_incomplete_quantum_settings_are_unavailable(settings):
    item = stencil()
    item["plus"]["context"]["quantum_settings"] = settings
    report = analyze_stencils(document(item))
    assert report["status"] == "unavailable"
    assert report["comparisons"] == []
    assert_no_certification(report)
