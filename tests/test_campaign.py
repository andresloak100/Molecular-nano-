from dataclasses import asdict
from copy import deepcopy
import json

import numpy as np
import pytest
from ase import Atoms
from ase.io import write

import nanodesign.campaign as campaign
from nanodesign.cli import main
from nanodesign.design import load_design
from nanodesign.quantum import QuantumSettings


def designs(tmp_path, count=2):
    paths = []
    for i in range(count):
        directory = tmp_path / f"source-{i}"
        directory.mkdir()
        write(directory / "initial.xyz", Atoms("H2", positions=[[0, 0, 0], [0, 0, .74 + .01*i]]))
        write(directory / "final.xyz", Atoms("H2", positions=[[0, 0, 0], [0, 0, .76 + .01*i]]))
        data = {"schema_version": 1, "length_unit": "angstrom", "initial": "initial.xyz",
                "final": "final.xyz", "quantum": asdict(QuantumSettings(spin=0, threads=1)), "fixed_indices": []}
        path = directory / "design.json"
        path.write_text(json.dumps(data))
        paths.append(path)
    return paths


def fake_run(design, output, **controls):
    output.mkdir(parents=True, exist_ok=False)
    settings = load_design(design)[3]
    result = {"status": "completed", "input_hashes": load_design(design)[4],
              "stage": controls["stage"], "state": controls["state"],
              "quantum_settings": asdict(settings),
              "optimization": {"fmax_ev_per_angstrom": controls.get("fmax", .03),
                    "max_steps_per_stage": controls.get("steps", 200), "images": controls.get("images", 7)},
              "structure": {"energy_ev": -1.0, "endpoint_identity_ok": False,
                  "free_force_max_ev_per_angstrom": 0.0, "forces_ev_per_angstrom": [[0, 0, 0], [0, 0, 0]],
                  "quantum_diagnostics": {"scf_converged": True, "gradient_completed": True}},
              "validation": {"design_validated": False}}
    (output / "result.json").write_text(json.dumps(result))
    return result


def test_snapshot_survives_source_changes_and_copies_exact_coordinates(tmp_path):
    paths = designs(tmp_path)
    original = (paths[0].parent / "initial.xyz").read_bytes()
    root = tmp_path / "campaign"
    report = campaign.create_campaign(paths, root)
    assert report["counts"] == {"pending": 2}
    assert (root / "designs/pose-0001/initial.xyz").read_bytes() == original
    paths[0].write_text("destroyed outside snapshot")
    assert campaign.campaign_report(root)["counts"] == {"pending": 2}


def test_bounded_resumption_does_not_recompute_completed(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path), root)
    calls = []
    def runner(*args, **kwargs):
        calls.append(args)
        return fake_run(*args, **kwargs)
    monkeypatch.setattr(campaign, "run", runner)
    assert campaign.run_campaign(root)["counts"] == {"completed": 1, "pending": 1}
    assert campaign.run_campaign(root)["counts"] == {"completed": 2}
    report = campaign.run_campaign(root)
    assert len(calls) == 2
    assert report["automatic_ranking"] is False
    assert report["design_validated"] is False
    assert report["rows"][0]["hydrogen_basin_preserved"] is False


def test_retry_is_explicit_and_preserves_failed_attempt(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    def failure(*args, **kwargs):
        result = fake_run(*args, **kwargs)
        result["status"] = "failed"
        (args[1] / "result.json").write_text(json.dumps(result))
        raise RuntimeError("SCF failed")
    monkeypatch.setattr(campaign, "run", failure)
    assert campaign.run_campaign(root)["counts"] == {"failed": 1}
    first = root / "runs/pose-0001/attempt-0001/result.json"
    evidence = first.read_bytes()
    monkeypatch.setattr(campaign, "run", fake_run)
    assert campaign.run_campaign(root)["counts"] == {"failed": 1}
    report = campaign.run_campaign(root, retry_incomplete=True)
    assert report["counts"] == {"completed": 1}
    assert report["rows"][0]["attempts"] == 2
    assert first.read_bytes() == evidence


def test_keyboard_interrupt_persists_attempt(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt
    monkeypatch.setattr(campaign, "run", interrupted)
    with pytest.raises(KeyboardInterrupt):
        campaign.run_campaign(root)
    assert campaign.campaign_report(root)["counts"] == {"interrupted": 1}


@pytest.mark.parametrize("completed", [True, False])
def test_recovers_worker_death_without_repeating_finished_work(tmp_path, monkeypatch, completed):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    state = json.loads((root / "campaign.json").read_text())
    output = "runs/pose-0001/attempt-0001"
    state["attempts"]["pose-0001"] = [{"status": "running", "output": output}]
    (root / "campaign.json").write_text(json.dumps(state))
    if completed:
        fake_run(root / "designs/pose-0001/design.json", root / output, stage="singlepoint", state="initial")
    monkeypatch.setattr(campaign, "run", lambda *a, **kw: pytest.fail("Unexpected recomputation"))
    report = campaign.run_campaign(root)
    assert report["counts"] == {"completed" if completed else "abandoned": 1}


@pytest.mark.parametrize("artifact", ["plan.json", "designs/pose-0001/initial.xyz", "designs/pose-0001/design.json"])
def test_changed_inputs_rejected_before_work(tmp_path, monkeypatch, artifact):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    with (root / artifact).open("a") as stream:
        stream.write("\n")
    monkeypatch.setattr(campaign, "run", lambda *a, **kw: pytest.fail("Unexpected calculation"))
    with pytest.raises(ValueError, match="integrity"):
        campaign.run_campaign(root)


def test_changed_completed_result_rejected(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    monkeypatch.setattr(campaign, "run", fake_run)
    campaign.run_campaign(root)
    with (root / "runs/pose-0001/attempt-0001/result.json").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="changed"):
        campaign.campaign_report(root)


def test_missing_completed_evidence_is_an_error(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    monkeypatch.setattr(campaign, "run", fake_run)
    campaign.run_campaign(root)
    (root / "runs/pose-0001/attempt-0001/result.json").unlink()
    with pytest.raises(ValueError, match="missing"):
        campaign.campaign_report(root)
    with pytest.raises(ValueError, match="missing"):
        campaign.run_campaign(root)


def test_truncated_unfinished_evidence_preserved_and_explicitly_retried(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    state = json.loads((root / "campaign.json").read_text())
    output = "runs/pose-0001/attempt-0001"
    state["attempts"]["pose-0001"] = [{"status": "running", "output": output}]
    (root / "campaign.json").write_text(json.dumps(state))
    (root / output).mkdir(parents=True)
    broken = root / output / "result.json"
    broken.write_text('{"status":')
    monkeypatch.setattr(campaign, "run", fake_run)
    result = campaign.run_campaign(root)
    assert result["counts"] == {"abandoned": 1}
    assert result["rows"][0]["evidence_error"]
    assert "energy_ev" not in result["rows"][0]
    assert campaign.run_campaign(root, retry_incomplete=True)["counts"] == {"completed": 1}
    assert broken.read_text() == '{"status":'


def test_false_completion_without_artifact_fails(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    monkeypatch.setattr(campaign, "run", lambda *a, **kw: {"status": "completed"})
    report = campaign.run_campaign(root)
    assert report["counts"] == {"failed": 1}
    assert "energy_ev" not in report["rows"][0]


@pytest.mark.parametrize("change", ["energy_missing", "nonfinite", "gradient_missing", "wrong_method", "wrong_controls", "false_relaxation"])
def test_completed_evidence_is_verified_before_recovery(tmp_path, monkeypatch, change):
    root = tmp_path / "campaign"
    stage = "relax" if change == "false_relaxation" else "singlepoint"
    campaign.create_campaign(designs(tmp_path, 1), root, stage=stage)
    state = json.loads((root / "campaign.json").read_text())
    output = "runs/pose-0001/attempt-0001"
    state["attempts"]["pose-0001"] = [{"status": "running", "output": output}]
    (root / "campaign.json").write_text(json.dumps(state))
    result = fake_run(root / "designs/pose-0001/design.json", root / output, stage=stage, state="initial")
    if change == "energy_missing":
        result["structure"].pop("energy_ev")
    elif change == "nonfinite":
        result["structure"]["energy_ev"] = float("nan")
    elif change == "gradient_missing":
        result["structure"]["quantum_diagnostics"]["gradient_completed"] = False
    elif change == "wrong_method":
        result["quantum_settings"]["basis"] = "sto-3g"
    elif change == "wrong_controls":
        result["optimization"]["fmax_ev_per_angstrom"] = 3.0
    else:
        result["geometry_converged"] = False
    (root / output / "result.json").write_text(json.dumps(result))
    monkeypatch.setattr(campaign, "run", lambda *a, **kw: pytest.fail("Invalid evidence triggered work"))
    with pytest.raises(ValueError):
        campaign.run_campaign(root)


@pytest.mark.parametrize("wrong_barrier", [False, True])
def test_path_report_keeps_barrier_context_and_rejects_wrong_value(tmp_path, monkeypatch, wrong_barrier):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root, stage="path", images=5)
    def runner(*args, **kwargs):
        result = fake_run(*args, **kwargs)
        structure = result.pop("structure")
        structure["endpoint_identity_ok"] = True
        structure["topology_screen"] = {"preserved": True}
        result.update(endpoints_converged=True, neb_converged=True,
            endpoints=[dict(deepcopy(structure), state=state, geometry_converged=True) for state in ("initial", "final")],
            images=[dict(deepcopy(structure), energy_ev=-1.0 + energy) for energy in (0, .1, .2, .1, -.1)],
            relative_energies_ev=[0, .1, .2, .1, -.1], peak_image=2,
            candidate_electronic_barrier_ev=99.0 if wrong_barrier else .2,
            barrier_interpretation="Candidate only; no verified transition state or free-energy barrier.")
        (args[1] / "result.json").write_text(json.dumps(result))
        return result
    monkeypatch.setattr(campaign, "run", runner)
    if wrong_barrier:
        with pytest.raises(ValueError, match="barrier"):
            campaign.run_campaign(root)
    else:
        result = campaign.run_campaign(root)
        row = result["rows"][0]
        assert row["candidate_electronic_barrier_ev"] == .2
        assert "no verified" in row["barrier_interpretation"]
        assert len(row["endpoints"]) == 2
        assert result["quantum_settings"]["basis"] == "def2-svp"
        assert row["electronic_state_identity_verified"] is False


def test_lock_prevents_second_worker(tmp_path, monkeypatch):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    monkeypatch.setattr(campaign, "run", lambda *a, **kw: pytest.fail("Concurrent run"))
    with campaign._lock(root):
        with pytest.raises(RuntimeError, match="Another worker"):
            campaign.run_campaign(root)


def test_mixed_method_campaign_rejected_before_output(tmp_path):
    paths = designs(tmp_path)
    data = json.loads(paths[1].read_text())
    data["quantum"]["basis"] = "sto-3g"
    paths[1].write_text(json.dumps(data))
    with pytest.raises(ValueError, match="share"):
        campaign.create_campaign(paths, tmp_path / "campaign")
    assert not (tmp_path / "campaign").exists()


@pytest.mark.parametrize("kwargs", [{"steps": True}, {"images": 3}, {"fmax": True}, {"fmax": float("nan")}, {"stage": "invented"}])
def test_invalid_controls_rejected(tmp_path, kwargs):
    with pytest.raises(ValueError):
        campaign.create_campaign(designs(tmp_path, 1), tmp_path / "campaign", **kwargs)
    assert not (tmp_path / "campaign").exists()


@pytest.mark.parametrize("max_jobs", [0, True, 101, 1.5])
def test_job_budget_checked_before_execution(tmp_path, max_jobs):
    with pytest.raises(ValueError):
        campaign.run_campaign(tmp_path, max_jobs=max_jobs)


def test_generated_grid_has_distinct_poses_and_no_calculations(tmp_path, monkeypatch):
    monkeypatch.setattr(campaign, "run", lambda *a, **kw: pytest.fail("Grid creation launched a calculation"))
    root = tmp_path / "grid"
    result = campaign.create_pose_campaign(root, [3.4, 3.6], [-.2, .2])
    assert result["counts"] == {"pending": 4}
    assert len({tuple(row["pose"].values()) for row in result["rows"]}) == 4
    for entry in json.loads((root / "plan.json").read_text())["designs"]:
        assert len(load_design(root / entry["design"])[1]) == 53


def test_cli_creates_reads_and_runs_campaign(tmp_path, monkeypatch, capsys):
    root = tmp_path / "campaign"
    paths = designs(tmp_path, 1)
    assert main(["campaign-create", str(paths[0]), "--out", str(root)]) == 0
    assert json.loads(capsys.readouterr().out)["counts"] == {"pending": 1}
    monkeypatch.setattr(campaign, "run", fake_run)
    assert main(["campaign-run", str(root)]) == 0
    assert json.loads(capsys.readouterr().out)["counts"] == {"completed": 1}
    assert main(["campaign-report", str(root)]) == 0
    assert json.loads(capsys.readouterr().out)["automatic_ranking"] is False


@pytest.mark.quantum
def test_real_h2_campaign_and_resumption(tmp_path):
    root = tmp_path / "campaign"
    campaign.create_campaign(designs(tmp_path, 1), root)
    result = campaign.run_campaign(root)
    assert result["counts"] == {"completed": 1}
    assert np.isfinite(result["rows"][0]["energy_ev"])
    assert result["rows"][0]["free_force_max_ev_per_angstrom"] > 0
    evidence = (root / result["rows"][0]["result"]).read_bytes()
    assert campaign.run_campaign(root)["rows"][0]["attempts"] == 1
    assert (root / result["rows"][0]["result"]).read_bytes() == evidence
