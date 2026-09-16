"""Synthetic grid/metadata fixtures only; no quantum calculations or real claims."""

import copy
import hashlib
import importlib.util
import math
from pathlib import Path
import struct

import pytest


spec = importlib.util.spec_from_file_location("resolved_numerics_review", Path(__file__).with_name("resolved_numerics.py"))
resolved = importlib.util.module_from_spec(spec)
spec.loader.exec_module(resolved)


@pytest.fixture
def evidence():
    case = {"id": "synthetic-only", "atom_count": 2, "symbols": ["H", "C"],
            "settings": {"grid_level": 3, "density_fit": True}}
    grids = [[0.0, 0.0, 1.0, 0.125], [1.0, 0.0, 0.0, 0.25]]
    cpu = {
        "synthetic_fixture": True,
        "resolved_numerics": {
            "precision": {"energy": "float64", "forces": "float64", "density": "float64"},
            "orbital_basis": {"sha256": hashlib.sha256(b"synthetic orbital basis").hexdigest()},
            "auxiliary_basis": {"status": "resolved", "sha256": hashlib.sha256(b"synthetic auxiliary basis").hexdigest()},
            "quadrature": {
                "level": 3, "atom_grid": {"H": {"radial": 50, "angular": 194}, "C": {"radial": 75, "angular": 302}},
                "pruning": "synthetic-pruning-v1", "radial_method": "synthetic-radial-v1",
                "becke_scheme": "synthetic-becke-v1", "point_count": 2,
                "canonicalization": resolved.GRID_CANONICALIZATION,
                "signature_sha256": resolved.canonical_grid_signature(grids),
            },
        },
    }
    return cpu, copy.deepcopy(cpu), case


def test_matching_complete_metadata_passes_without_mutation(evidence):
    before = copy.deepcopy(evidence)
    assert resolved.validate_resolved_numerics(*evidence) == []
    assert evidence == before


def test_canonical_signature_byte_contract_and_order_invariance():
    rows = [[1.0, 0.0, 0.0, 0.25], [0.0, 0.0, 1.0, 0.125]]
    expected = hashlib.sha256(struct.pack("<Q", 2)
                              + struct.pack("<dddd", 0.0, 0.0, 1.0, 0.125)
                              + struct.pack("<dddd", 1.0, 0.0, 0.0, 0.25)).hexdigest()
    assert resolved.canonical_grid_signature(rows) == expected
    assert resolved.canonical_grid_signature(list(reversed(rows))) == expected


def test_signed_zero_is_normalized_before_sorting_and_serialization():
    positive = [[0.0, 1.0, 0.0, 0.0]]
    signed = [[-0.0, 1.0, -0.0, -0.0]]
    assert resolved.canonical_grid_signature(positive) == resolved.canonical_grid_signature(signed)


def test_grid_signature_preserves_coordinate_weight_pairing():
    original = [[0, 0, 1, 0.125], [1, 0, 0, 0.25]]
    reassigned = [[0, 0, 1, 0.25], [1, 0, 0, 0.125]]
    assert resolved.canonical_grid_signature(original) != resolved.canonical_grid_signature(reassigned)


def test_single_ulp_weight_change_is_not_rounded_away():
    original = [[0, 0, 1, 0.125]]
    changed = [[0, 0, 1, math.nextafter(0.125, 1.0)]]
    assert resolved.canonical_grid_signature(original) != resolved.canonical_grid_signature(changed)


def test_duplicate_rows_are_preserved_and_count_is_bound():
    row = [0, 0, 1, 0.125]
    assert resolved.canonical_grid_signature([row]) != resolved.canonical_grid_signature([row, row])


@pytest.mark.parametrize("rows", [None, [], [[0, 0, 0]], [[0, 0, 0, True]],
                                 [[0, 0, 0, "1"]], [[0, 0, 0, float("nan")]],
                                 [[0, 0, 0, float("inf")]], [[0, 0, 0, 10**400]]])
def test_grid_helper_rejects_malformed_or_nonfinite_rows(rows):
    with pytest.raises(ValueError):
        resolved.canonical_grid_signature(rows)


@pytest.mark.parametrize("field", ["precision", "orbital_basis", "auxiliary_basis", "quadrature"])
@pytest.mark.parametrize("target", [0, 1], ids=["cpu", "gpu"])
def test_required_resolved_objects_cannot_be_omitted(evidence, field, target):
    del evidence[target]["resolved_numerics"][field]
    errors = resolved.validate_resolved_numerics(*evidence)
    assert any(item.startswith("unavailable:") and field in item for item in errors)


def test_capture_template_with_nulls_is_unavailable(evidence):
    evidence[1]["resolved_numerics"] = {"precision": None, "orbital_basis": None, "auxiliary_basis": None, "quadrature": None}
    errors = resolved.validate_resolved_numerics(*evidence)
    assert len(errors) == 4
    assert all(item.startswith("unavailable:") for item in errors)


@pytest.mark.parametrize("placeholder", ["", "default", "AUTO", "unknown", "unresolved", "not_recorded", "not recorded"])
@pytest.mark.parametrize("field", ["pruning", "radial_method", "becke_scheme"])
def test_placeholder_recipes_are_not_resolved_choices(evidence, placeholder, field):
    evidence[1]["resolved_numerics"]["quadrature"][field] = placeholder
    assert any(item.startswith("unavailable:") and field in item for item in resolved.validate_resolved_numerics(*evidence))


@pytest.mark.parametrize("field", ["energy", "forces", "density"])
def test_precision_must_be_explicit_float64_for_every_quantity(evidence, field):
    evidence[1]["resolved_numerics"]["precision"][field] = "float32"
    assert any(item.startswith("mismatch:") and field in item for item in resolved.validate_resolved_numerics(*evidence))


@pytest.mark.parametrize("field", ["energy", "forces", "density"])
def test_missing_precision_is_unavailable(evidence, field):
    del evidence[1]["resolved_numerics"]["precision"][field]
    assert any(item.startswith("unavailable:") and field in item for item in resolved.validate_resolved_numerics(*evidence))


@pytest.mark.parametrize("field", ["level", "point_count"])
def test_boolean_quadrature_counts_are_malformed(evidence, field):
    evidence[1]["resolved_numerics"]["quadrature"][field] = True
    assert any(item.startswith("malformed:") and field in item for item in resolved.validate_resolved_numerics(*evidence))


@pytest.mark.parametrize("field,value", [("radial", True), ("angular", False), ("radial", 0), ("angular", 50.0)])
def test_element_grid_counts_are_positive_integers(evidence, field, value):
    evidence[1]["resolved_numerics"]["quadrature"]["atom_grid"]["H"][field] = value
    assert any(item.startswith("malformed:") and field in item for item in resolved.validate_resolved_numerics(*evidence))


def test_every_case_element_requires_resolved_grid_counts(evidence):
    del evidence[1]["resolved_numerics"]["quadrature"]["atom_grid"]["C"]
    assert any("missing case elements ['C']" in item for item in resolved.validate_resolved_numerics(*evidence))


def test_atom_grid_requires_final_object_schema_not_count_list(evidence):
    evidence[1]["resolved_numerics"]["quadrature"]["atom_grid"]["H"] = [50, 194]
    assert any(item.startswith("malformed:") and "atom_grid.H" in item for item in resolved.validate_resolved_numerics(*evidence))


@pytest.mark.parametrize("field", ["pruning", "radial_method", "becke_scheme", "point_count", "atom_grid", "signature_sha256"])
def test_resolved_cpu_gpu_disagreements_cannot_pass(evidence, field):
    quadrature = evidence[1]["resolved_numerics"]["quadrature"]
    if field == "point_count":
        quadrature[field] += 1
    elif field == "atom_grid":
        quadrature[field]["H"]["radial"] += 1
    elif field == "signature_sha256":
        quadrature[field] = hashlib.sha256(b"different synthetic grid").hexdigest()
    else:
        quadrature[field] = "different-resolved-id"
    assert any(item.startswith("mismatch:") and field in item for item in resolved.validate_resolved_numerics(*evidence))


def test_matching_wrong_requested_grid_level_is_rejected(evidence):
    for record in evidence[:2]:
        record["resolved_numerics"]["quadrature"]["level"] = 4
    assert any("requested grid_level" in item for item in resolved.validate_resolved_numerics(*evidence))


def test_canonicalization_identity_is_required(evidence):
    evidence[1]["resolved_numerics"]["quadrature"]["canonicalization"] = "rounded-grid-v1"
    assert any("canonicalization" in item and item.startswith("mismatch:") for item in resolved.validate_resolved_numerics(*evidence))


@pytest.mark.parametrize("field", ["orbital_basis", "auxiliary_basis"])
def test_expanded_basis_digest_mismatch_is_explicit(evidence, field):
    evidence[1]["resolved_numerics"][field]["sha256"] = hashlib.sha256(b"different synthetic basis").hexdigest()
    assert any(item.startswith("mismatch:") and field in item for item in resolved.validate_resolved_numerics(*evidence))


def test_direct_calculation_requires_no_auxiliary_basis(evidence):
    evidence[2]["settings"]["density_fit"] = False
    for record in evidence[:2]:
        record["resolved_numerics"]["auxiliary_basis"] = {"status": "not_applicable"}
    assert resolved.validate_resolved_numerics(*evidence) == []
    evidence[1]["resolved_numerics"]["auxiliary_basis"]["sha256"] = "a" * 64
    assert any(item.startswith("malformed:") and "auxiliary_basis.sha256" in item for item in resolved.validate_resolved_numerics(*evidence))


def test_fitted_calculation_cannot_claim_auxiliary_basis_not_applicable(evidence):
    evidence[1]["resolved_numerics"]["auxiliary_basis"] = {"status": "not_applicable"}
    errors = resolved.validate_resolved_numerics(*evidence)
    assert any("auxiliary_basis.status" in item for item in errors)
    assert any("auxiliary_basis.sha256" in item for item in errors)


def test_raw_artifact_hashes_can_differ_without_overriding_canonical_signature(evidence):
    for index, record in enumerate(evidence[:2]):
        quadrature = record["resolved_numerics"]["quadrature"]
        quadrature["grid_points_sha256"] = hashlib.sha256(f"synthetic points {index}".encode()).hexdigest()
        quadrature["grid_weights_sha256"] = hashlib.sha256(f"synthetic weights {index}".encode()).hexdigest()
    assert resolved.validate_resolved_numerics(*evidence) == []
    del evidence[1]["resolved_numerics"]["quadrature"]["signature_sha256"]
    assert any(item.startswith("unavailable:") and "signature_sha256" in item for item in resolved.validate_resolved_numerics(*evidence))


@pytest.mark.parametrize("bad_digest", ["default", "a" * 63, "z" * 64, True, None])
def test_invalid_basis_hashes_never_default_to_requested_basis(evidence, bad_digest):
    evidence[1]["resolved_numerics"]["orbital_basis"]["sha256"] = bad_digest
    assert any("orbital_basis.sha256" in item for item in resolved.validate_resolved_numerics(*evidence))


@pytest.mark.parametrize("field,value", [("symbols", None), ("symbols", ["H", "Qq"]), ("atom_count", True),
                                         ("atom_count", 3), ("settings", None)])
def test_case_identity_must_be_complete_and_well_typed(evidence, field, value):
    evidence[2][field] = value
    errors = resolved.validate_resolved_numerics(*evidence)
    assert errors
    assert all("case." in item for item in errors)


@pytest.mark.parametrize("field,value", [("density_fit", 1), ("density_fit", None), ("grid_level", True), ("grid_level", 10)])
def test_requested_settings_cannot_be_inferred(evidence, field, value):
    evidence[2]["settings"][field] = value
    assert any(f"case.settings.{field}" in item for item in resolved.validate_resolved_numerics(*evidence))
