"""Independent F2 polynomial/provenance cases; no quantum calculations."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys

import pytest


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


oracle = load("f2_polynomial_oracle", HERE / "analytic_fixtures.py")
f1 = load("f2_reviewed_f1", HERE.parent / "force-energy-consistency/consistency.py")


def no_scientific_promotion(report):
    for key in ("electronic_state_identity_verified", "electronic_branch_continuity_verified",
                "scientific_model_validated", "numerical_convergence_certified"):
        assert report[key] is False


def compare(document, tolerance=None):
    before = deepcopy(document)
    report = f1.analyze_stencils(document, force_tolerance_ev_per_angstrom=tolerance)
    assert document == before
    no_scientific_promotion(report)
    return report


@pytest.mark.parametrize("a,b", [(0.125, 0.125), (0.125, 0.0625), (0.5, 0.125), (0.0625, 0.5)])
def test_exact_quadratic_force_at_reference_even_for_unequal_steps(a, b):
    document = oracle.stencil_document(((a, b),))
    report = compare(document, 1e-12)
    assert report["status"] == "comparison_produced"
    row = report["comparisons"][0]
    assert row["analytic_force_ev_per_angstrom"] == -1.25
    assert row["finite_difference_force_ev_per_angstrom"] == pytest.approx(-1.25, abs=1e-12)
    assert row["within_declared_tolerance"] is True
    assert row["actual_minus_distance_angstrom"] == a
    assert row["actual_plus_distance_angstrom"] == b
    assert document["stencils"][0]["reference"]["forces_ev_per_angstrom"][1] == [-3.0, 0.0, 0.0]
    if a != b:
        point = document["stencils"][0]
        wrong_secant = -(point["plus"]["energy_ev"] - point["minus"]["energy_ev"]) / (a + b)
        assert abs(wrong_secant + 1.25) >= 0.125


@pytest.mark.parametrize("a,b", [(0.125, 0.0625), (0.0625, 0.125), (0.125, 0.125)])
def test_asymmetric_cubic_and_quartic_truncation_matches_symbolic_oracle(a, b):
    report = compare(oracle.stencil_document(((a, b),), cubic=8, quartic=16))
    actual = report["comparisons"][0]["finite_difference_force_ev_per_angstrom"]
    assert actual == pytest.approx(oracle.exact_three_point_force(a=a, b=b, cubic=8, quartic=16), abs=1e-12)
    assert report["comparisons"][0]["within_declared_tolerance"] is None


def test_three_step_cubic_trend_is_observation_not_a_certified_bound():
    steps = [0.125, 0.0625, 0.03125]
    report = compare(oracle.stencil_document(tuple((h, h) for h in steps), cubic=8))
    errors = [row["absolute_residual_ev_per_angstrom"] for row in report["comparisons"]]
    assert errors == pytest.approx([8*h*h for h in steps], abs=1e-12)
    assert errors[0] / errors[1] == errors[1] / errors[2] == 4
    assert report["step_sensitivity"][0]["distinct_stencils"] == 3
    assert report["step_sensitivity"][0]["finite_difference_force_spread_ev_per_angstrom"] == pytest.approx(max(errors)-min(errors))


def test_wrong_force_sign_fails_numerical_tolerance_without_changing_evidence():
    document = oracle.stencil_document()
    document["stencils"][0]["reference"]["forces_ev_per_angstrom"][0][0] *= -1
    row = compare(document, .01)["comparisons"][0]
    assert row["absolute_residual_ev_per_angstrom"] == pytest.approx(2.5)
    assert row["within_declared_tolerance"] is False


def test_large_saved_energy_offset_cannot_recover_lost_difference_digits():
    document = oracle.stencil_document(offset=1e16)
    point = document["stencils"][0]
    assert point["reference"]["energy_ev"] == point["plus"]["energy_ev"] == point["minus"]["energy_ev"]
    row = compare(document, 1e-3)["comparisons"][0]
    assert row["finite_difference_force_ev_per_angstrom"] == 0
    assert row["absolute_residual_ev_per_angstrom"] == 1.25
    assert row["within_declared_tolerance"] is False


def test_unresolved_energy_equality_gets_explicit_resolution_warning_even_if_residual_is_zero():
    document = oracle.stencil_document(offset=1e16, slope=0, cubic=8)
    row = compare(document, .001)["comparisons"][0]
    assert oracle.exact_three_point_force(a=.125, b=.125, slope=0, cubic=8) == -.125
    assert row["finite_difference_force_ev_per_angstrom"] == row["analytic_force_ev_per_angstrom"] == 0
    assert row["within_declared_tolerance"] is True  # Literal residual statement only.
    assert row["binary64_input_rounding_sensitivity_ev_per_angstrom"] == 8
    assert row["input_rounding_exceeds_declared_tolerance"] is True
    assert row["numerical_resolution_findings"]


def test_finite_extreme_offset_sum_cannot_silently_collapse_derivative_weights():
    document = oracle.stencil_document()
    stencil = document["stencils"][0]
    stencil["requested_minus_offset_angstrom"] = -1e308
    stencil["requested_plus_offset_angstrom"] = 1e308
    for role, x in (("reference", 0.0), ("minus", -1e308), ("plus", 1e308)):
        point = stencil[role]
        point["positions_angstrom"][0][0] = x
        point["energy_ev"] = x
        point["forces_ev_per_angstrom"][0][0] = -1
    report = compare(document)
    # The inputs are finite, but a+b overflows. A guarded rejection or a
    # correctly scaled evaluation is acceptable; silently zero weights are not.
    if report["status"] == "invalid":
        assert report["comparisons"] == []
    else:
        assert report["status"] == "comparison_produced"
        assert report["comparisons"][0]["finite_difference_force_ev_per_angstrom"] == pytest.approx(-1)


def test_representable_constant_energy_and_coordinate_origin_do_not_change_result():
    document = oracle.stencil_document(((.125, .0625),), cubic=8)
    original = compare(document)["comparisons"][0]
    for role in ("reference", "minus", "plus"):
        point = document["stencils"][0][role]
        point["energy_ev"] += 4096
        for row in point["positions_angstrom"]:
            for axis in range(3):
                row[axis] += 16
    shifted = compare(document)["comparisons"][0]
    assert shifted["finite_difference_force_ev_per_angstrom"] == original["finite_difference_force_ev_per_angstrom"]
    assert shifted["absolute_residual_ev_per_angstrom"] == original["absolute_residual_ev_per_angstrom"]


@pytest.mark.parametrize("claimed_id", [None, "same-converged-guess", "verified-electronic-state"])
def test_state_label_never_becomes_branch_verification(claimed_id):
    report = compare(oracle.stencil_document(state_evidence_id=claimed_id), 1e-12)
    assert report["status"] == "comparison_produced"
    assert report["comparisons"][0]["within_declared_tolerance"] is True


@pytest.mark.parametrize("missing", ["plus_record", "minus_energy", "baseline_energy", "baseline_force", "failed_plus"])
def test_missing_or_failed_input_never_becomes_an_accepted_difference(missing):
    document = oracle.stencil_document()
    point = document["stencils"][0]
    if missing == "plus_record":
        point["plus"] = None
    elif missing == "minus_energy":
        point["minus"]["energy_ev"] = None
    elif missing == "baseline_energy":
        point["reference"]["energy_ev"] = None
    elif missing == "baseline_force":
        point["reference"]["forces_ev_per_angstrom"] = None
    else:
        point["plus"]["status"] = "failed"
    report = compare(document)
    assert report["status"] == "unavailable"
    assert report["comparisons"] == []


@pytest.mark.parametrize("fault", ["unit", "energy_scope", "masked_force", "settings", "guess", "off_axis",
                                   "collapsed", "wrong_sign", "frozen", "call_identity", "boolean_energy", "force_shape"])
def test_incompatible_evidence_is_rejected_before_numerical_comparison(fault):
    document = oracle.stencil_document()
    point = document["stencils"][0]
    if fault == "unit":
        for role in ("reference", "minus", "plus"):
            point[role]["context"]["units"]["energy"] = "Hartree"
    elif fault == "energy_scope":
        point["plus"]["context"]["energy_scope"] = "SCF-only; excludes added dispersion"
    elif fault == "masked_force":
        for role in ("reference", "minus", "plus"):
            point[role]["context"]["force_scope"] = "constrained_masked"
    elif fault == "settings":
        point["plus"]["context"]["quantum_settings"]["basis"] = "def2-svp"
    elif fault == "guess":
        point["plus"]["context"]["quantum_settings"]["scf_initial_guess"] = "atom"
    elif fault == "off_axis":
        point["plus"]["positions_angstrom"][0][1] += .125
    elif fault == "collapsed":
        point["plus"]["positions_angstrom"] = deepcopy(point["reference"]["positions_angstrom"])
    elif fault == "wrong_sign":
        point["minus"]["positions_angstrom"] = deepcopy(point["plus"]["positions_angstrom"])
    elif fault == "frozen":
        for role in ("reference", "minus", "plus"):
            point[role]["context"]["fixed_indices"] = [0, 1]
    elif fault == "call_identity":
        point["plus"]["calculation_call_id"] = point["minus"]["calculation_call_id"]
    elif fault == "boolean_energy":
        point["plus"]["energy_ev"] = True
    elif fault == "force_shape":
        point["reference"]["forces_ev_per_angstrom"][0].pop()
    report = compare(document)
    assert report["status"] == "invalid"
    assert report["comparisons"] == []
    assert report["findings"]


def test_missing_second_step_keeps_first_result_without_imputing_convergence():
    document = oracle.stencil_document(((.125, .125), (.0625, .0625)), cubic=8)
    document["stencils"][1]["minus"] = None
    report = compare(document)
    assert report["status"] == "partial"
    assert len(report["comparisons"]) == 1
    assert report["comparisons"][0]["finite_difference_force_ev_per_angstrom"] == -1.375


def test_file_reader_rejects_duplicate_keys_instead_of_silently_overwriting(tmp_path):
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema_version":1,"schema_version":1,"stencils":[]}')
    with pytest.raises(ValueError):
        f1.read_stencil_document(path)


def test_file_reader_preserves_valid_fixture_and_source_bytes(tmp_path):
    document = oracle.stencil_document()
    path = tmp_path / "valid.json"
    path.write_text(json.dumps(document))
    before = path.read_bytes()
    readback = f1.read_stencil_document(path)
    assert readback == document
    assert path.read_bytes() == before
