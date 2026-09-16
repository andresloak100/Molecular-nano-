"""Independent public-API checks using a synthetic ASE calculator only."""

from dataclasses import asdict
import hashlib
import json
from pathlib import Path

from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import write
import numpy as np
import pytest

from nanodesign import state_scan
from nanodesign.quantum import QuantumSettings
from mock_backend import install_probe


def file_bytes(directory):
    return {str(path.relative_to(directory)): path.read_bytes()
            for path in Path(directory).rglob("*") if path.is_file()}


@pytest.fixture
def scan_fixture(tmp_path, monkeypatch):
    first = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    selected = Atoms("H2", positions=[[0.3, 0.1, 0], [0.3, 0.1, 1.1]])
    selected.set_constraint(FixAtoms(indices=[0]))
    # Coordinate metadata must not override the caller's quantum state.
    selected.info.update(charge=7, spin=7)
    source = tmp_path / "input.traj"
    write(source, [first, selected])
    settings = QuantumSettings(charge=0, spin=0, basis="sto-3g", dispersion=None,
                               grid_level=1, threads=1, memory_mb=256)
    probe = install_probe(monkeypatch, state_scan)
    return source, selected, settings, probe


def test_selected_frame_and_settings_survive_deleted_source(tmp_path, scan_fixture):
    source, selected, settings, probe = scan_fixture
    output = tmp_path / "scan"
    guesses = ("huckel", "atom", "minao", "1e")
    state_scan.create_state_scan(source, output, settings, image=-1, guesses=guesses)
    assert probe["calls"] == []
    source.unlink()
    state_scan.run_state_scan(output, max_jobs=4)
    assert len(probe["calls"]) == len(guesses)
    assert len({id(call["calculator"]) for call in probe["calls"]}) == len(guesses)
    for call, guess in zip(probe["calls"], guesses):
        np.testing.assert_array_equal(call["positions"], selected.positions)
        np.testing.assert_array_equal(call["numbers"], selected.numbers)
        assert call["settings"] == {**asdict(settings), "scf_initial_guess": guess}


def test_failed_attempt_consumes_budget_and_is_not_implicitly_retried(tmp_path, scan_fixture):
    source, _, settings, probe = scan_fixture
    output = tmp_path / "scan"
    probe["faults"]["minao"] = "failure"
    state_scan.create_state_scan(source, output, settings, guesses=("minao", "atom"))
    state_scan.run_state_scan(output, max_jobs=1)
    assert [call["settings"]["scf_initial_guess"] for call in probe["calls"]] == ["minao"]
    failed_bytes = file_bytes(output)
    assert any(b"S3 deliberately failed" in raw for raw in failed_bytes.values())
    state_scan.run_state_scan(output, max_jobs=1)
    assert [call["settings"]["scf_initial_guess"] for call in probe["calls"]] == ["minao", "atom"]
    state_scan.run_state_scan(output, max_jobs=4)
    assert len(probe["calls"]) == 2
    assert any(b"S3 deliberately failed" in raw for raw in file_bytes(output).values())


def test_reporting_pending_and_finished_evidence_is_read_only(tmp_path, scan_fixture):
    source, _, settings, probe = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, guesses=("minao", "atom"))
    before = file_bytes(output)
    state_scan.state_scan_report(output)
    assert file_bytes(output) == before
    assert probe["calls"] == []
    state_scan.run_state_scan(output, max_jobs=2)
    before = file_bytes(output)
    state_scan.state_scan_report(output)
    assert file_bytes(output) == before
    assert len(probe["calls"]) == 2


def test_interrupted_attempt_is_preserved_and_not_retried(tmp_path, scan_fixture):
    source, _, settings, probe = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, guesses=("minao", "atom", "1e"))
    probe["faults"]["atom"] = "interrupt"
    with pytest.raises(KeyboardInterrupt):
        state_scan.run_state_scan(output, max_jobs=3)
    assert [call["settings"]["scf_initial_guess"] for call in probe["calls"]] == ["minao", "atom"]
    saved = file_bytes(output)
    assert any(b"S3 deliberately interrupted" in raw for raw in saved.values())
    state_scan.state_scan_report(output)
    assert file_bytes(output) == saved
    state_scan.run_state_scan(output, max_jobs=3)
    assert [call["settings"]["scf_initial_guess"] for call in probe["calls"]] == ["minao", "atom", "1e"]
    assert any(b"S3 deliberately interrupted" in raw for raw in file_bytes(output).values())


@pytest.mark.parametrize("fault", ["nonfinite_energy", "malformed_forces", "missing_gradient", "wrong_settings"])
def test_convergence_flag_does_not_accept_contradictory_outputs(tmp_path, scan_fixture, fault):
    source, _, settings, probe = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, guesses=("minao", "atom"))
    probe["faults"]["minao"] = fault
    report = state_scan.run_state_scan(output, max_jobs=2)
    assert report["failed_guesses"] == ["minao"]
    assert report["converged_guesses"] == ["atom"]
    assert report["all_requested_guesses_converged"] is False
    assert report["energy_spread"]["member_guesses"] == ["atom"]
    assert report["energy_spread"]["energy_spread_ev"] is None
    assert report["ground_state_verified"] is False
    assert report["electronic_state_identity_verified"] is False


@pytest.mark.parametrize("first_fails", [False, True])
def test_older_attempt_tampering_is_rejected_after_later_success(tmp_path, scan_fixture, first_fails):
    source, _, settings, probe = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, guesses=("minao", "atom", "1e"))
    if first_fails:
        probe["faults"]["minao"] = "failure"
    report = state_scan.run_state_scan(output, max_jobs=3)
    oldest = output / report["rows"][0]["output"] / "result.json"
    record = json.loads(oldest.read_text())
    record["S3_tamper_marker"] = "altered after later guesses finished"
    oldest.write_text(json.dumps(record))
    before = file_bytes(output)
    with pytest.raises(ValueError):
        state_scan.state_scan_report(output)
    with pytest.raises(ValueError):
        state_scan.run_state_scan(output, max_jobs=4)
    assert file_bytes(output) == before
    assert len(probe["calls"]) == 3


@pytest.mark.parametrize("completed_count", [0, 1, 2])
def test_equal_subset_energies_do_not_claim_complete_state_scan(tmp_path, scan_fixture, completed_count):
    source, _, settings, probe = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings)
    if completed_count:
        state_scan.run_state_scan(output, max_jobs=completed_count)
    report = state_scan.state_scan_report(output)
    assert len(report["converged_guesses"]) == completed_count
    assert len(report["pending_guesses"]) == 4 - completed_count
    assert report["energy_spread"]["scope"] == "converged_subset"
    spread = report["energy_spread"]["energy_spread_ev"]
    assert spread == (0.0 if completed_count == 2 else None)
    assert report["all_requested_attempts_finished"] is False
    assert report["all_requested_guesses_converged"] is False
    assert report["complete_four_guess_check"] is False
    assert report["ground_state_verified"] is False
    assert report["electronic_state_identity_verified"] is False


def test_fixed_atom_forces_are_saved_without_masking(tmp_path, scan_fixture):
    source, selected, settings, _ = scan_fixture
    assert selected.constraints
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, guesses=("minao",))
    report = state_scan.run_state_scan(output)
    assert report["converged_guesses"] == ["minao"]
    result = json.loads((output / report["rows"][0]["output"] / "result.json").read_text())
    np.testing.assert_allclose(result["forces_ev_per_angstrom"],
                               [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])


def test_backend_geometry_mutation_is_rejected_and_next_guess_is_isolated(tmp_path, scan_fixture):
    source, selected, settings, probe = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, guesses=("minao", "atom"))
    probe["faults"]["minao"] = "mutate_atoms"
    report = state_scan.run_state_scan(output, max_jobs=2)
    assert report["failed_guesses"] == ["minao"]
    assert report["converged_guesses"] == ["atom"]
    second = next(call for call in probe["calls"] if call["settings"]["scf_initial_guess"] == "atom")
    np.testing.assert_array_equal(second["positions"], selected.positions)


@pytest.mark.parametrize("result_state", ["completed", "partial_json", "missing"])
def test_recovery_preserves_crash_evidence_and_never_repeats_guess(tmp_path, scan_fixture, result_state):
    source, _, settings, probe = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, guesses=("minao", "atom"))
    first = state_scan.run_state_scan(output, max_jobs=1)
    path = output / first["rows"][0]["output"] / "result.json"
    ledger_path = output / "scan.json"
    ledger = json.loads(ledger_path.read_text())
    # Simulate a process exiting before the terminal ledger update. These are
    # isolated mock artifacts, never real calculation records.
    ledger["attempts"]["minao"]["status"] = "running"
    ledger["attempts"]["minao"].pop("result_sha256")
    if result_state == "partial_json":
        path.write_bytes(b'{"status": "run')
    elif result_state == "missing":
        path.unlink()
    ledger_path.write_text(json.dumps(ledger))
    before = file_bytes(output)
    report = state_scan.state_scan_report(output)
    assert report["rows"][0]["status"] == "running"
    assert report["converged_guesses"] == []
    assert file_bytes(output) == before
    result_bytes = path.read_bytes() if path.exists() else None
    final = state_scan.run_state_scan(output, max_jobs=2)
    assert [call["settings"]["scf_initial_guess"] for call in probe["calls"]] == ["minao", "atom"]
    assert (path.read_bytes() if path.exists() else None) == result_bytes
    if result_state == "completed":
        assert final["converged_guesses"] == ["minao", "atom"]
    else:
        assert final["abandoned_guesses"] == ["minao"]
        assert final["converged_guesses"] == ["atom"]
    assert final["ground_state_verified"] is False


@pytest.mark.parametrize("digest", ["absent", None, ""])
def test_completed_attempt_requires_recorded_digest(tmp_path, scan_fixture, digest):
    source, _, settings, _ = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, guesses=("minao",))
    state_scan.run_state_scan(output)
    path = output / "scan.json"
    ledger = json.loads(path.read_text())
    if digest == "absent":
        del ledger["attempts"]["minao"]["result_sha256"]
    else:
        ledger["attempts"]["minao"]["result_sha256"] = digest
    path.write_text(json.dumps(ledger))
    with pytest.raises(ValueError):
        state_scan.state_scan_report(output)


@pytest.mark.parametrize("field", ["image", "quantum_settings", "geometry"])
def test_binding_types_cannot_alias_frozen_values(tmp_path, scan_fixture, field):
    source, _, settings, _ = scan_fixture
    output = tmp_path / "scan"
    state_scan.create_state_scan(source, output, settings, image=1, guesses=("minao",))
    report = state_scan.run_state_scan(output)
    path = output / report["rows"][0]["output"] / "result.json"
    record = json.loads(path.read_text())
    if field == "image":
        record["image"] = True  # Python equality alone treats True as 1.
    elif field == "quantum_settings":
        record["quantum_settings"]["charge"] = False
    else:
        record["geometry"]["positions_angstrom"][0][2] = False
    raw = json.dumps(record).encode()
    path.write_bytes(raw)
    ledger_path = output / "scan.json"
    ledger = json.loads(ledger_path.read_text())
    # Refresh only the outer file checksum to isolate semantic binding checks.
    # This models malformed writer output, not authentic scientific evidence.
    ledger["attempts"]["minao"]["result_sha256"] = hashlib.sha256(raw).hexdigest()
    ledger_path.write_text(json.dumps(ledger))
    with pytest.raises(ValueError):
        state_scan.state_scan_report(output)
