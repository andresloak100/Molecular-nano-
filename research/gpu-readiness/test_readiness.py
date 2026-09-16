"""Offline inventory and exact-input protocol checks; no electronic calculations."""
import json
from pathlib import Path

import pytest

import gpu_readiness as readiness


def inventory(monkeypatch, packages):
    monkeypatch.setattr(readiness, "package_versions", lambda: packages)
    monkeypatch.setattr(readiness.shutil, "which", lambda name: None)
    monkeypatch.setattr(readiness, "memory_bytes", lambda: 24 * 1024**3)


def test_inventory_never_silently_imports_or_runs_solver(monkeypatch):
    inventory(monkeypatch, {"pyscf": "2.14.0"})
    monkeypatch.setattr(readiness, "run_command", lambda *a, **k: pytest.fail("Unexpected runtime execution"))
    result = readiness.environment_report()
    assert result["runtime_probe"]["status"] == "not_requested"
    assert result["nvidia_smi"]["status"] == "command_not_found"
    assert result["prerequisites_observed"] is False
    assert result["quantum_evaluations"] == 0
    assert result["production_gpu_adapter_implemented"] is False
    assert result["gpu_numerical_parity_measured"] is False


@pytest.mark.parametrize("capability,expected", [([8, 0], True), ([6, 1], False)])
def test_importable_gpu_is_only_a_prerequisite(monkeypatch, capability, expected):
    inventory(monkeypatch, {name: "1" for name in ["pyscf", "gpu4pyscf-cuda12x", "cupy-cuda12x", "numpy", "ase", "dftd3"]})
    payload = {"cuda": {"status": "available", "devices": [{"compute_capability": capability}]}}
    monkeypatch.setattr(readiness, "run_command", lambda *a, **k: {"returncode": 0, "stdout": "import message\nG1_JSON=" + json.dumps(payload)})
    result = readiness.environment_report(probe_runtime=True)
    assert result["prerequisites_observed"] is expected
    assert result["gpu_speedup_measured"] is False
    assert result["production_gpu_adapter_implemented"] is False


def test_broken_runtime_is_not_ready(monkeypatch):
    inventory(monkeypatch, {})
    monkeypatch.setattr(readiness, "run_command", lambda *a, **k: {"returncode": None, "error": "timeout"})
    result = readiness.environment_report(probe_runtime=True)
    assert result["runtime_probe"]["status"] == "failed"
    assert result["prerequisites_observed"] is False


def test_plan_preserves_exact_geometry_and_sources(tmp_path):
    plan, payloads = readiness.build_plan()
    assert plan["planned_energy_force_evaluations"] == 10
    assert plan["quantum_evaluations_completed"] == 0
    assert [case["atom_count"] for case in plan["cases"]] == [2, 5, 3, 3, 8]
    assert {case["settings"]["scf_initial_guess"] for case in plan["cases"]} == {"minao", "atom"}
    for case in plan["cases"]:
        raw = (readiness.ROOT / case["coordinate_source"]).read_bytes()
        assert payloads[case["geometry"]] == raw
        assert case["geometry_sha256"] == readiness.digest(raw)
        assert "scf_initial_guess" in case["planned_changes_from_record"]
    destination = tmp_path / "protocol"
    readiness.save_plan(destination, plan, payloads)
    assert json.loads((destination / "plan.json").read_text()) == plan
    for path, expected_digest in plan["bundled_files_sha256"].items():
        assert readiness.digest((destination / path).read_bytes()) == expected_digest
    with pytest.raises(FileExistsError):
        readiness.save_plan(destination, plan, payloads)


def test_large_candidate_requires_explicit_plan_option():
    plan, _ = readiness.build_plan(include_candidate=True)
    assert plan["planned_energy_force_evaluations"] == 12
    candidate = plan["cases"][-1]
    assert candidate["atom_count"] == 53
    assert candidate["symbols"].count("C") == 22
    assert candidate["symbols"].count("H") == 31
    assert candidate["fixed_indices_context"] == [7, 8, 9, 35, 36, 37]
    assert candidate["hydrogen_transfer_context"] == {"donor": 0, "hydrogen": 10, "acceptor": 26}


def test_existing_evidence_cannot_be_overwritten(tmp_path):
    output = tmp_path / "report.json"
    readiness.write_new_json(output, {"original": True})
    before = output.read_bytes()
    with pytest.raises(FileExistsError):
        readiness.write_new_json(output, {"original": False})
    assert output.read_bytes() == before


def test_input_symlinks_cannot_escape_repo(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    external = tmp_path / "outside.xyz"
    external.write_text("external")
    (root / "geometry.xyz").symlink_to(external)
    with pytest.raises(ValueError, match="escapes"):
        readiness.read_inside(root, "geometry.xyz")


def test_generated_protocol_matches_checker_schema_without_computation():
    """Synthetic output values test plumbing, not the H2 energy or GPU execution."""
    from copy import deepcopy
    from equivalence import compare_records
    plan, _ = readiness.build_plan()
    case = plan["cases"][0]
    cpu = {"schema_version": 1, "synthetic_fixture": True,
           "case_id": case["id"], "geometry_sha256": case["geometry_sha256"],
           "settings": case["settings"], "backend": "pyscf_cpu", "execution_device": "cpu",
           "scf_initial_guess": case["settings"]["scf_initial_guess"],
           "scf_converged": True, "gradient_completed": True, "energy_hartree": -1.0,
           "forces_ev_per_angstrom": [[0.0, 0.0, 0.0] for _ in range(case["atom_count"])],
           "s2": 0.0, "grid_response": True, "auxiliary_basis_response": False,
           "dispersion_evaluations": 0, "versions": {"pyscf": "synthetic"}}
    gpu = deepcopy(cpu)
    gpu.update(backend="gpu4pyscf", execution_device="nvidia_cuda", gpu_synchronized=True,
               versions={"pyscf": "synthetic", "gpu4pyscf": "synthetic", "cupy": "synthetic"})
    result = compare_records(plan, cpu, gpu)
    assert result["record_validation_passed"] is True
    assert result["numerical_parity_passed"] is True
    assert result["scientific_validated"] is False
