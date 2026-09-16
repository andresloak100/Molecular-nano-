"""Synthetic evidence fixtures only: these are not quantum-calculation results."""

import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


MODULE_PATH = Path(__file__).with_name("equivalence.py")
spec = importlib.util.spec_from_file_location("gpu_evidence_equivalence", MODULE_PATH)
equivalence = importlib.util.module_from_spec(spec)
spec.loader.exec_module(equivalence)


@pytest.fixture
def evidence():
    settings = {
        "charge": 0, "spin": 1, "xc": "pbe0", "basis": "def2-svp",
        "dispersion": "d3bj", "grid_level": 3, "conv_tol": 1e-9,
        "max_cycle": 150, "threads": 1, "memory_mb": 2000,
        "density_fit": True, "scf_initial_guess": "atom",
    }
    plan = {
        "schema_version": 1, "kind": "cpu_gpu_equivalence_protocol",
        "synthetic_fixture": True,
        "cases": [{"id": "synthetic-only", "geometry_sha256": "a" * 64,
                   "atom_count": 2, "settings": settings}],
        "tolerances": {"energy_hartree": 1e-6,
                       "max_force_component_ev_per_angstrom": 1e-4, "s2": 1e-5},
    }
    cpu = {
        "schema_version": 1, "synthetic_fixture": True,
        "case_id": "synthetic-only", "geometry_sha256": "a" * 64,
        "settings": copy.deepcopy(settings), "backend": "pyscf_cpu",
        "execution_device": "cpu", "scf_initial_guess": "atom",
        "scf_converged": True, "gradient_completed": True,
        "energy_hartree": -1.0, "forces_ev_per_angstrom": [[0.01, -0.02, 0.03], [-0.01, 0.02, -0.03]],
        "s2": 0.75, "grid_response": True, "auxiliary_basis_response": True,
        "dispersion_evaluations": 1, "versions": {"pyscf": "synthetic"},
    }
    gpu = copy.deepcopy(cpu)
    gpu.update(backend="gpu4pyscf", execution_device="nvidia_cuda",
               gpu_synchronized=True,
               versions={"pyscf": "synthetic", "gpu4pyscf": "synthetic", "cupy": "synthetic"})
    return plan, cpu, gpu


def test_complete_records_pass_only_numerical_gate_without_mutation(evidence):
    original = copy.deepcopy(evidence)
    result = equivalence.compare_records(*evidence)
    assert result["record_validation_passed"] is True
    assert result["numerical_parity_passed"] is True
    assert result["scientific_validated"] is False
    assert result["electronic_state_identity_verified"] is False
    assert result["ground_state_verified"] is False
    assert set(result["differences"].values()) == {0}
    assert result["mismatches"] == []
    assert evidence == original


@pytest.mark.parametrize("metric", ["energy_hartree", "max_force_component_ev_per_angstrom", "s2"])
def test_over_threshold_is_explicit_with_measured_difference(evidence, metric):
    plan, cpu, gpu = evidence
    delta = plan["tolerances"][metric] * 2
    if metric == "max_force_component_ev_per_angstrom":
        gpu["forces_ev_per_angstrom"][1][2] += delta
    else:
        gpu[metric] += delta
    result = equivalence.compare_records(*evidence)
    assert result["record_validation_passed"] is True
    assert result["numerical_parity_passed"] is False
    assert result["differences"][metric] == pytest.approx(delta)
    assert any(metric in message and "exceeds" in message for message in result["mismatches"])


def test_declared_tolerance_boundary_is_inclusive(evidence):
    plan, cpu, gpu = evidence
    cpu["energy_hartree"] = 0.0
    gpu["energy_hartree"] = plan["tolerances"]["energy_hartree"]
    assert equivalence.compare_records(*evidence)["numerical_parity_passed"] is True


@pytest.mark.parametrize("target", [1, 2], ids=["cpu", "gpu"])
@pytest.mark.parametrize("field", [
    "schema_version", "case_id", "geometry_sha256", "settings", "backend", "execution_device",
    "scf_initial_guess", "scf_converged", "gradient_completed", "energy_hartree",
    "forces_ev_per_angstrom", "s2", "grid_response", "auxiliary_basis_response",
    "dispersion_evaluations", "versions",
])
def test_missing_required_evidence_never_defaults_to_success(evidence, target, field):
    del evidence[target][field]
    result = equivalence.compare_records(*evidence)
    assert result["numerical_parity_passed"] is False
    assert result["record_validation_passed"] is False
    assert result["mismatches"]


@pytest.mark.parametrize("field,value", [
    ("backend", "pyscf_cpu"), ("execution_device", "cpu"),
    ("gpu_synchronized", False), ("gpu_synchronized", 1),
    ("scf_initial_guess", "minao"), ("geometry_sha256", "b" * 64),
    ("case_id", "other-case"), ("gradient_completed", False),
    ("scf_converged", 1), ("grid_response", False),
    ("auxiliary_basis_response", False), ("dispersion_evaluations", 2),
    ("dispersion_evaluations", True), ("schema_version", True),
    ("energy_hartree", True), ("s2", None),
    ("versions", {"pyscf": "synthetic", "gpu4pyscf": "synthetic"}),
    ("versions", {"pyscf": "synthetic", "gpu4pyscf": "synthetic", "cupy": ""}),
])
def test_conflicting_execution_or_missing_support_is_rejected(evidence, field, value):
    evidence[2][field] = value
    result = equivalence.compare_records(*evidence)
    assert result["numerical_parity_passed"] is False
    assert any(field in message for message in result["mismatches"])


def test_missing_gpu_synchronization_is_rejected(evidence):
    del evidence[2]["gpu_synchronized"]
    assert any("gpu_synchronized" in item for item in equivalence.compare_records(*evidence)["mismatches"])


@pytest.mark.parametrize("forces", [[], [[1, 2, 3]], [[0, 0], [0, 0]], [[0, 0, True], [0, 0, 0]], "array"])
def test_invalid_force_array_shape_or_types_are_rejected(evidence, forces):
    evidence[2]["forces_ev_per_angstrom"] = forces
    assert any("forces_ev_per_angstrom" in item for item in equivalence.compare_records(*evidence)["mismatches"])


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_values_anywhere_cannot_appear_in_output(evidence, value):
    evidence[2]["extra_diagnostic"] = {"bad": value}
    result = equivalence.compare_records(*evidence)
    assert result["numerical_parity_passed"] is False
    assert any("nonfinite" in message for message in result["mismatches"])
    json.dumps(result, allow_nan=False)


def test_finite_input_difference_overflow_is_rejected(evidence):
    evidence[1]["energy_hartree"] = 1e308
    evidence[2]["energy_hartree"] = -1e308
    result = equivalence.compare_records(*evidence)
    assert result["numerical_parity_passed"] is False
    assert any("overflowed" in message for message in result["mismatches"])
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("change", ["omit_guess", "change_guess", "unknown_setting", "int_to_bool"])
def test_settings_must_match_exactly(evidence, change):
    settings = evidence[2]["settings"]
    if change == "omit_guess":
        del settings["scf_initial_guess"]
    elif change == "change_guess":
        settings["scf_initial_guess"] = "minao"
    elif change == "unknown_setting":
        settings["unknown"] = True
    else:
        settings["spin"] = True
    result = equivalence.compare_records(*evidence)
    assert any("settings differ" in item for item in result["mismatches"])


@pytest.mark.parametrize("change", ["missing_tolerance", "bool_tolerance", "negative_tolerance", "empty_cases",
                                   "duplicate_case", "bad_hash", "bool_atom_count", "missing_setting",
                                   "invalid_guess", "bool_schema", "bad_kind"])
def test_invalid_plans_are_rejected_before_comparison(evidence, change):
    plan = evidence[0]
    if change == "missing_tolerance":
        del plan["tolerances"]["s2"]
    elif change == "bool_tolerance":
        plan["tolerances"]["s2"] = True
    elif change == "negative_tolerance":
        plan["tolerances"]["s2"] = -1
    elif change == "empty_cases":
        plan["cases"] = []
    elif change == "duplicate_case":
        plan["cases"].append(copy.deepcopy(plan["cases"][0]))
    elif change == "bad_hash":
        plan["cases"][0]["geometry_sha256"] = "unknown"
    elif change == "bool_atom_count":
        plan["cases"][0]["atom_count"] = True
    elif change == "missing_setting":
        del plan["cases"][0]["settings"]["scf_initial_guess"]
    elif change == "invalid_guess":
        plan["cases"][0]["settings"]["scf_initial_guess"] = "invented"
    elif change == "bool_schema":
        plan["schema_version"] = True
    else:
        plan["kind"] = "benchmark-results"
    result = equivalence.compare_records(*evidence)
    assert result["record_validation_passed"] is False
    assert result["numerical_parity_passed"] is False
    assert all(item.startswith("plan") for item in result["mismatches"])


@pytest.mark.parametrize("density_fit,dispersion", [(False, None), (True, None), (False, "d3zero")])
def test_response_and_dispersion_counts_follow_each_explicit_case(evidence, density_fit, dispersion):
    plan, cpu, gpu = evidence
    for settings in (plan["cases"][0]["settings"], cpu["settings"], gpu["settings"]):
        settings.update(density_fit=density_fit, dispersion=dispersion)
    for record in (cpu, gpu):
        record.update(auxiliary_basis_response=density_fit, dispersion_evaluations=0 if dispersion is None else 1)
    assert equivalence.compare_records(*evidence)["numerical_parity_passed"] is True


def test_single_matching_case_is_selected_from_plan(evidence):
    unrelated = copy.deepcopy(evidence[0]["cases"][0])
    unrelated["id"] = "other-case"
    evidence[0]["cases"].insert(0, unrelated)
    assert equivalence.compare_records(*evidence)["case_id"] == "synthetic-only"


@pytest.mark.parametrize("bad_json", ['{"a": 1, "a": 2}', '{"a": NaN}', '{"a": Infinity}', '{'])
def test_json_reader_rejects_ambiguous_or_malformed_json(tmp_path, bad_json):
    path = tmp_path / "invalid.json"
    path.write_text(bad_json)
    with pytest.raises(ValueError):
        equivalence.read_json(path)


def test_cli_reads_explicit_files_and_reports_failures(evidence, tmp_path):
    paths = [tmp_path / name for name in ("plan.json", "cpu.json", "gpu.json")]
    for path, payload in zip(paths, evidence):
        path.write_text(json.dumps(payload))
    command = [sys.executable, str(MODULE_PATH), "compare", "--plan", str(paths[0]),
               "--cpu", str(paths[1]), "--gpu", str(paths[2])]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    assert completed.returncode == 0
    assert json.loads(completed.stdout)["numerical_parity_passed"] is True
    paths[2].write_text('{"energy_hartree": NaN}')
    failed = subprocess.run(command, text=True, capture_output=True, check=False)
    assert failed.returncode == 1
    result = json.loads(failed.stdout)
    assert result["numerical_parity_passed"] is result["scientific_validated"] is False
    assert "Unable to read evidence" in result["mismatches"][0]
