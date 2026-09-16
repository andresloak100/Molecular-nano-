"""Corrupt temporary archive copies to verify audit failures, never raw evidence."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

SPEC = importlib.util.spec_from_file_location("e1_audit", Path(__file__).with_name("audit_evidence.py"))
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@pytest.fixture
def archive(tmp_path):
    root = tmp_path / "validation"
    for directory in ("pbe0-svp", "pbe0-tzvp", "h2-integration", "h-abstraction-direct-initial", "h-abstraction-df-initial"):
        shutil.copytree(MODULE.DEFAULT_ROOT / directory, root / directory)
    return root


def alter(root, path, mutate):
    file = root / path
    value = json.loads(file.read_text())
    mutate(value)
    file.write_text(json.dumps(value))


def check_failed(root, name, contains):
    report = MODULE.audit(root)
    assert report["status"] == "failed"
    result = next(c for c in report["checks"] if c["check"] == name)
    assert result["status"] == "failed"
    assert contains in result["error"]
    assert report["scientific_model_validated"] is False


def test_actual_archives_and_recorded_limits():
    report = MODULE.audit()
    assert report["status"] == "passed", report["checks"]
    assert len(report["checks"]) == 4
    h2 = next(c["details"] for c in report["checks"] if c["check"] == "h2_characterization")
    assert h2["raw_displaced_forces_available"] is False
    assert h2["unresolved_modes_per_step"] == [5, 5]
    assert h2["max_mode_equation_residual"] < 1e-8


def test_detects_input_byte_change(archive):
    file = archive / "pbe0-svp/inputs/methane.xyz"
    file.write_text(file.read_text() + "\n")
    check_failed(archive, "pbe0_svp", "SHA-256 mismatch")


def test_detects_aggregate_arithmetic_change(archive):
    alter(archive, "pbe0-tzvp/benchmark.json", lambda d: d["computed"].update(reference_geometry_energy_ev=2))
    check_failed(archive, "pbe0_tzvp", "reference_geometry_energy")


def test_rejects_contradictory_reference_diagnostics(archive):
    alter(archive, "pbe0-svp/species/methane.json", lambda d: d["quantum_diagnostics"].update(scf_converged=False))
    alter(archive, "pbe0-svp/benchmark.json", lambda d: d["species"]["methane"]["quantum_diagnostics"].update(scf_converged=False))
    check_failed(archive, "pbe0_svp", "diagnostic SCF convergence")


def test_rejects_contradictory_relaxation_force(archive):
    alter(archive, "h2-integration/relax/result.json", lambda d: d["structure"].update(free_force_max_ev_per_angstrom=999))
    check_failed(archive, "h2_characterization", "relaxed force maximum")


def test_detects_hessian_change(archive):
    alter(archive, "h2-integration/modes/result.json", lambda d: d["stationary"]["hessian_ev_per_angstrom2"][0].__setitem__(0, 3))
    check_failed(archive, "h2_characterization", "Hessian eigenvalues")


def test_detects_mode_normalization_change(archive):
    alter(archive, "h2-integration/modes/result.json", lambda d: d["stationary"]["cartesian_modes_per_sqrt_amu"][0][0].__setitem__(0, 2))
    check_failed(archive, "h2_characterization", "orthonormal")


def test_detects_wrong_frequency_scale(archive):
    alter(archive, "h2-integration/modes/result.json", lambda d: d["stationary"]["frequencies_cm1"].__setitem__(5, 4000))
    check_failed(archive, "h2_characterization", "frequency unit conversion")


def test_detects_wrong_step_summary(archive):
    alter(archive, "h2-integration/step-comparison.json", lambda d: d.update(stretch_frequency_change_cm1=0))
    check_failed(archive, "h2_characterization", "H2 stretch change")


def test_rejects_nonfinite_evidence(archive):
    alter(archive, "h2-integration/modes/result.json", lambda d: d["stationary"].update(free_force_max_ev_per_angstrom=float("nan")))
    check_failed(archive, "h2_characterization", "non-finite JSON")


def test_rejects_missing_completed_displacement(archive):
    file = archive / "h2-integration/modes/electronic.jsonl"
    rows = [json.loads(line) for line in file.read_text().splitlines()]
    last = next(row for row in reversed(rows) if row["event"] == "calculation_completed")
    rows.remove(last)
    file.write_text("\n".join(json.dumps(row) for row in rows) + "\n")
    check_failed(archive, "h2_characterization", "completion events missing")


@pytest.mark.parametrize("mutation,expected", [
    ("failed_completion", "completion flags indicate failure"),
    ("missing_stages", "per-call SCF/gradient evidence missing"),
    ("duplicate_completion", "completion events missing or duplicated"),
])
def test_rejects_incomplete_event_evidence(archive, mutation, expected):
    file = archive / "h2-integration/modes/electronic.jsonl"
    rows = [json.loads(line) for line in file.read_text().splitlines()]
    if mutation == "failed_completion":
        for row in rows:
            if row["event"] == "calculation_completed":
                row["scf_converged"] = row["gradient_completed"] = False
    elif mutation == "missing_stages":
        rows = [row for row in rows if row["event"] == "calculation_completed"]
    else:
        rows.append(next(row for row in rows if row["event"] == "calculation_completed"))
    file.write_text("\n".join(json.dumps(row) for row in rows) + "\n")
    check_failed(archive, "h2_characterization", expected)


def test_rejects_incomplete_base_diagnostics(archive):
    alter(archive, "h2-integration/modes/result.json", lambda d: d["quantum_diagnostics"].update(scf_converged=False))
    check_failed(archive, "h2_characterization", "base electronic evidence incomplete")


def test_rejects_inconsistent_base_event_energy(archive):
    file = archive / "h2-integration/modes/electronic.jsonl"
    rows = [json.loads(line) for line in file.read_text().splitlines()]
    for row in rows:
        if row["event"] == "calculation_completed":
            row["energy_eV"] = 999
    file.write_text("\n".join(json.dumps(row) for row in rows) + "\n")
    check_failed(archive, "h2_characterization", "base-event energy")


def test_rejects_unsupported_saddle_classification(archive):
    alter(archive, "h2-integration/modes/result.json", lambda d: d["stationary"].update(classification="first_order_saddle_candidate"))
    check_failed(archive, "h2_characterization", "stationary classification")


def test_rejects_validation_flag_promotion(archive):
    alter(archive, "h2-integration/modes/result.json", lambda d: d["validation"].update(transition_state_validated=True))
    check_failed(archive, "h2_characterization", "promoted to validated")


def test_refuses_scientific_promotion(archive):
    alter(archive, "h2-integration/modes/result.json", lambda d: d["stationary"].update(transition_state_verified=True))
    check_failed(archive, "h2_characterization", "promoted to validated")


def test_missing_archive_is_a_failed_check(archive):
    (archive / "h-abstraction-direct-initial/result.json").unlink()
    check_failed(archive, "direct_df", "FileNotFoundError")


def test_refuses_path_traversal_and_symlink_escape(archive, tmp_path):
    access = MODULE.Archive(archive)
    with pytest.raises(ValueError, match="escapes"):
        access.raw("../outside.txt")
    outside = tmp_path / "outside.txt"
    outside.write_text("not archive evidence")
    (archive / "escape").symlink_to(outside)
    with pytest.raises(ValueError, match="escapes"):
        access.raw("escape")


def test_cli_failure_status_and_report_no_overwrite(archive, tmp_path):
    output = tmp_path / "report.json"
    alter(archive, "h2-integration/step-comparison.json", lambda d: d.update(stretch_frequency_change_cm1=0))
    command = [sys.executable, str(Path(MODULE.__file__)), "--validation-root", str(archive), "--output", str(output)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 1
    original = output.read_bytes()
    assert json.loads(original)["status"] == "failed"
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert output.read_bytes() == original
