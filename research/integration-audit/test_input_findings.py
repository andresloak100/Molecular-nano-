"""Bounded Q1 input/campaign integrity probes; no quantum calculations.

These desired-behavior checks reproduced three defects before owner repairs.
All associated strict xfail markers were removed after verified repairs.
Run with .venv/bin/python -m pytest -q research/integration-audit/test_input_findings.py
"""
from dataclasses import asdict
import json
from pathlib import Path

import pytest
from ase import Atoms
from ase.io import write

from nanodesign import campaign, design
from nanodesign.quantum import QuantumSettings


def make_design(directory):
    directory.mkdir()
    initial = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    final = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.76]])
    write(directory / "initial.xyz", initial)
    write(directory / "final.xyz", final)
    document = {
        "schema_version": 1, "length_unit": "angstrom",
        "initial": "initial.xyz", "final": "final.xyz", "fixed_indices": [],
        "quantum": asdict(QuantumSettings(spin=0, threads=1)),
    }
    path = directory / "design.json"
    path.write_text(json.dumps(document))
    return path


def fake_runner(design_path, output, **controls):
    """A self-consistent synthetic record; never calls a physical calculator."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    _, _, _, settings, hashes = design.load_design(design_path)
    result = {
        "schema_version": 1, "status": "completed", "input_hashes": hashes,
        "stage": controls["stage"], "state": controls["state"],
        "quantum_settings": asdict(settings),
        "optimization": {
            "fmax_ev_per_angstrom": controls["fmax"],
            "max_steps_per_stage": controls["steps"], "images": controls["images"],
        },
        "structure": {
            "energy_ev": -1.0, "free_force_max_ev_per_angstrom": 0.0,
            "forces_ev_per_angstrom": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "quantum_diagnostics": {
                "synthetic_test_only": True,
                "scf_converged": True, "gradient_completed": True,
            },
        },
        "validation": {"design_validated": False},
    }
    (output / "result.json").write_text(json.dumps(result))
    return result


def test_campaign_snapshots_valid_symlink_with_different_extension(tmp_path):
    path = make_design(tmp_path / "source")
    target = path.parent / "initial.xyz"
    renamed = path.parent / "actual-initial.extxyz"
    target.rename(renamed)
    target.symlink_to(renamed.name)
    # This is an accepted, readable production design before snapshotting.
    _, initial, _, _, _ = design.load_design(path)
    assert initial.get_distance(0, 1) == pytest.approx(0.74)
    report = campaign.create_campaign([path], tmp_path / "campaign")
    assert report["counts"] == {"pending": 1}
    copied = design.load_design(tmp_path / "campaign/designs/pose-0001/design.json")[1]
    assert copied.get_distance(0, 1) == pytest.approx(0.74)


@pytest.mark.parametrize("changed_source", ["coordinates", "design"])
def test_load_design_rejects_or_preserves_consistency_during_source_edit(tmp_path, monkeypatch, changed_source):
    path = make_design(tmp_path / "source")
    original_design_hash = design.sha256(path)
    initial_path = path.parent / "initial.xyz"
    original_coordinate_hash = design.sha256(initial_path)
    actual_read = design.read
    changed = False

    def read_and_edit(source, *args, **kwargs):
        nonlocal changed
        atoms = actual_read(source, *args, **kwargs)
        if not changed:
            changed = True
            if changed_source == "coordinates":
                replacement = atoms.copy()
                replacement.positions[1, 2] = 0.90
                write(initial_path, replacement)
            else:
                document = json.loads(path.read_text())
                document["quantum"]["basis"] = "sto-3g"
                path.write_text(json.dumps(document))
        return atoms

    monkeypatch.setattr(design, "read", read_and_edit)
    try:
        _, loaded, _, settings, hashes = design.load_design(path)
    except ValueError:
        # Rejecting the concurrent edit is also a valid outcome.
        return
    assert changed
    assert loaded.get_distance(0, 1) == pytest.approx(0.74)
    assert settings.basis == "def2-svp"
    # A returned record must hash the exact bytes actually parsed, or reject.
    assert hashes["initial_sha256"] == original_coordinate_hash
    assert hashes["design_sha256"] == original_design_hash


@pytest.mark.parametrize("operation", ["changed", "missing"])
def test_campaign_still_checks_historical_recorded_attempts(tmp_path, monkeypatch, operation):
    path = make_design(tmp_path / "source")
    root = tmp_path / "campaign"
    campaign.create_campaign([path], root)

    def failure(*args, **kwargs):
        result = fake_runner(*args, **kwargs)
        result["status"] = "failed"
        result["error"] = {"type": "SyntheticError", "message": "Original failure evidence"}
        (Path(args[1]) / "result.json").write_text(json.dumps(result))
        raise RuntimeError("Synthetic failure, no quantum work")

    monkeypatch.setattr(campaign, "run", failure)
    assert campaign.run_campaign(root)["counts"] == {"failed": 1}
    historical = root / "runs/pose-0001/attempt-0001/result.json"
    state = json.loads((root / "campaign.json").read_text())
    assert state["attempts"]["pose-0001"][0]["result_sha256"] == design.sha256(historical)
    monkeypatch.setattr(campaign, "run", fake_runner)
    assert campaign.run_campaign(root, retry_incomplete=True)["counts"] == {"completed": 1}
    if operation == "changed":
        historical.write_text('{"replaced": "original failure record lost"}')
    else:
        historical.unlink()
    with pytest.raises(ValueError, match="changed|missing|integrity"):
        campaign.campaign_report(root)


def test_positive_control_changed_current_coordinates_are_rejected(tmp_path):
    path = make_design(tmp_path / "source")
    root = tmp_path / "campaign"
    campaign.create_campaign([path], root)
    copied = root / "designs/pose-0001/initial.xyz"
    copied.write_text(copied.read_text() + "\n")
    with pytest.raises(ValueError, match="integrity"):
        campaign.campaign_report(root)


def test_positive_control_changed_latest_result_is_rejected(tmp_path, monkeypatch):
    path = make_design(tmp_path / "source")
    root = tmp_path / "campaign"
    campaign.create_campaign([path], root)
    monkeypatch.setattr(campaign, "run", fake_runner)
    assert campaign.run_campaign(root)["counts"] == {"completed": 1}
    current = root / "runs/pose-0001/attempt-0001/result.json"
    current.write_text(current.read_text() + "\n")
    with pytest.raises(ValueError, match="changed"):
        campaign.campaign_report(root)
