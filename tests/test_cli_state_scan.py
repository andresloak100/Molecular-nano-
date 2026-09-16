"""State-scan CLI contracts against real saved evidence and synthetic forces."""

from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write
from ase.units import Hartree

from nanodesign import state_scan
from nanodesign.cli import main
from nanodesign.quantum import QuantumSettings


class SyntheticCalculator(Calculator):
    """Records requested state/coordinates without importing a quantum solver."""

    implemented_properties = ["energy", "forces"]
    instances = []
    calls = []
    failures = {}

    def __init__(self, settings, *, event_log=None):
        super().__init__()
        self.settings = settings
        self.diagnostics = {"synthetic_test_only": True}
        type(self).instances.append(self)

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        guess = self.settings.scf_initial_guess
        type(self).calls.append({"guess": guess, "settings": asdict(self.settings),
                                 "positions": atoms.positions.copy()})
        if guess in self.failures:
            raise self.failures[guess]("Synthetic evaluation stopped")
        energy = {"minao": -1.0, "atom": -1.2, "1e": -1.2, "huckel": -1.0}[guess]
        expected_s2 = self.settings.spin / 2 * (self.settings.spin / 2 + 1)
        self.results = {"energy": energy, "forces": np.zeros((len(atoms), 3))}
        self.diagnostics.update(
            settings=asdict(self.settings), scf_initial_guess=guess,
            scf_converged=True, gradient_completed=True,
            total_energy_hartree=energy / Hartree, hartree_eV=float(Hartree),
            s2=expected_s2, expected_s2=expected_s2, s2_deviation=0.0,
        )


@pytest.fixture(autouse=True)
def synthetic_backend(monkeypatch):
    SyntheticCalculator.instances = []
    SyntheticCalculator.calls = []
    SyntheticCalculator.failures = {}
    monkeypatch.setattr(state_scan, "PySCFCalculator", SyntheticCalculator)


@pytest.fixture
def inputs(tmp_path):
    source = tmp_path / "source coordinates" / "two frames.extxyz"
    source.parent.mkdir()
    frames = [Atoms("H2", positions=[[0, 0, 0], [0, 0, distance]],
                    info={"charge": 0, "spin": 0}) for distance in (0.74, 1.2)]
    write(source, frames)
    settings = {"charge": 1, "spin": 1, "basis": "sto-3g", "dispersion": None,
                "grid_level": 1, "threads": 1}
    settings_path = tmp_path / "separate settings.json"
    settings_path.write_text(json.dumps(settings))
    return source, settings_path, settings


def cli_json(capsys, *arguments):
    code = main([str(argument) for argument in arguments])
    captured = capsys.readouterr()
    assert captured.err == ""
    return code, json.loads(captured.out)


def create_scan(tmp_path, inputs, capsys, *, guesses=None):
    source, settings_path, _ = inputs
    output = tmp_path / "saved scan"
    guess_arguments = [] if guesses is None else ["--guesses", *guesses]
    code, report = cli_json(capsys, "state-scan-create", source, "--out", output,
                            "--settings", settings_path, *guess_arguments)
    assert code == 0
    assert report["attempted_guesses"] == []
    assert report["requested_guesses"] == (list(guesses) if guesses is not None
                                            else ["minao", "atom", "1e", "huckel"])
    assert SyntheticCalculator.instances == []
    return output


def saved_bytes(directory):
    return {str(path.relative_to(directory)): path.read_bytes()
            for path in directory.rglob("*") if path.is_file()}


def assert_uncertified(report):
    for key in ("ground_state_verified", "electronic_state_identity_verified",
                "method_accuracy_validated", "geometry_optimized"):
        assert report[key] is False


@pytest.mark.parametrize("capture", [False, True])
def test_create_freezes_optional_electronic_capture_without_running_solver(tmp_path, inputs, capsys, capture):
    source, settings_path, _ = inputs
    output = tmp_path / "capture choice"
    arguments = ["--capture-electronic-state"] if capture else []
    code, report = cli_json(capsys, "state-scan-create", source, "--settings", settings_path,
                            "--out", output, *arguments)
    assert code == 0 and report["capture_electronic_state"] is capture
    assert json.loads((output / "plan.json").read_text())["capture_electronic_state"] is capture
    assert not (output / "attempts").exists()
    assert SyntheticCalculator.instances == []


@pytest.mark.parametrize("wrapped", [False, True], ids=["raw-settings", "design-quantum"])
def test_create_freezes_explicit_state_selected_frame_and_guess_order_without_work(
        tmp_path, inputs, capsys, wrapped):
    source, settings_path, settings = inputs
    if wrapped:
        # The structure is the CLI argument, independent of design-relative paths.
        settings_path.write_text(json.dumps({"quantum": settings,
                                             "initial": "missing-design-structure.xyz"}))
    output = tmp_path / "explicit scan"
    guesses = ["huckel", "1e", "atom"]
    code, report = cli_json(capsys, "state-scan-create", source, "--settings", settings_path,
                            "--out", output, "--image", 0, "--guesses", *guesses)
    assert code == 0
    assert report["status"] == "pending"
    assert report["requested_guesses"] == report["pending_guesses"] == guesses
    assert report["attempted_guesses"] == []
    assert report["quantum_settings"] == asdict(QuantumSettings(**settings))
    assert report["input"]["image"] == 0
    assert_uncertified(report)
    plan = json.loads((output / "plan.json").read_text())
    assert plan["geometry"]["positions_angstrom"] == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]
    assert plan["input"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert (output / plan["input"]["path"]).read_bytes() == source.read_bytes()
    assert not (output / "attempts").exists()
    assert SyntheticCalculator.instances == []

    source.unlink()
    settings_path.unlink()
    before = saved_bytes(output)
    code, inspected = cli_json(capsys, "state-scan-report", output)
    assert code == 0 and inspected == report
    assert saved_bytes(output) == before
    assert SyntheticCalculator.instances == []


def test_create_requires_settings_argument_before_creating_output(tmp_path, inputs, capsys):
    source, _, _ = inputs
    output = tmp_path / "missing settings"
    with pytest.raises(SystemExit) as error:
        main(["state-scan-create", str(source), "--out", str(output)])
    assert error.value.code == 2
    assert "--settings" in capsys.readouterr().err
    assert not output.exists()
    assert SyntheticCalculator.instances == []


@pytest.mark.parametrize("image", ["99", "-99"])
def test_nonexistent_frame_is_a_clean_cli_error_without_output(tmp_path, inputs, capsys, image):
    source, settings_path, _ = inputs
    output = tmp_path / "nonexistent frame"
    code = main(["state-scan-create", str(source), "--out", str(output),
                 "--settings", str(settings_path), "--image", image])
    captured = capsys.readouterr()
    assert code == 2 and captured.out == ""
    assert captured.err and "Traceback" not in captured.err
    assert not output.exists()
    assert SyntheticCalculator.instances == []


@pytest.mark.parametrize("raw", [
    "{}", '{"spin": 1}', '{"charge": 1}',
    '{"quantum": {}}', '{"quantum": {"spin": 1}}', '{"quantum": {"charge": 1}}',
    "[]", '{"quantum": null}', '{"charge":',
], ids=["no-state", "no-charge", "no-spin", "nested-no-state", "nested-no-charge",
        "nested-no-spin", "non-object", "null-quantum", "malformed-json"])
def test_invalid_or_implicit_settings_fail_before_output(tmp_path, inputs, capsys, raw):
    source, settings_path, _ = inputs
    settings_path.write_text(raw)
    output = tmp_path / "rejected scan"
    code = main(["state-scan-create", str(source), "--out", str(output),
                 "--settings", str(settings_path)])
    captured = capsys.readouterr()
    assert code == 2
    assert captured.out == "" and captured.err
    assert "Traceback" not in captured.err
    assert not output.exists()
    assert SyntheticCalculator.instances == []


def test_bounded_run_uses_frozen_frame_and_settings_and_never_repeats_attempts(tmp_path, inputs, capsys):
    output = create_scan(tmp_path, inputs, capsys)
    source, settings_path, settings = inputs
    source.unlink()
    settings_path.unlink()
    code, partial = cli_json(capsys, "state-scan-run", output)
    assert code == 0
    assert partial["status"] == "partial"
    assert partial["attempted_guesses"] == partial["converged_guesses"] == ["minao"]
    assert partial["pending_guesses"] == ["atom", "1e", "huckel"]
    assert partial["all_requested_attempts_finished"] is False
    assert partial["complete_four_guess_check"] is False
    assert partial["energy_spread"]["energy_spread_ev"] is None
    original_result = (output / "attempts/minao/result.json").read_bytes()

    code, resumed = cli_json(capsys, "state-scan-run", output, "--max-jobs", 2)
    assert code == 0
    assert resumed["converged_guesses"] == ["minao", "atom", "1e"]
    assert resumed["pending_guesses"] == ["huckel"]
    assert resumed["energy_spread"]["energy_spread_ev"] == pytest.approx(0.2)
    code, finished = cli_json(capsys, "state-scan-run", output)
    assert code == 0 and finished["status"] == "finished"
    assert finished["complete_four_guess_check"] is True
    assert_uncertified(finished)
    assert cli_json(capsys, "state-scan-run", output)[0] == 0
    assert (output / "attempts/minao/result.json").read_bytes() == original_result
    assert [call["guess"] for call in SyntheticCalculator.calls] == ["minao", "atom", "1e", "huckel"]
    assert len(SyntheticCalculator.instances) == 4
    for call in SyntheticCalculator.calls:
        assert call["settings"] == {**asdict(QuantumSettings(**settings)), "scf_initial_guess": call["guess"]}
        np.testing.assert_array_equal(call["positions"], [[0, 0, 0], [0, 0, 1.2]])


def test_failed_attempt_consumes_budget_and_stays_visible_during_continuation(tmp_path, inputs, capsys):
    output = create_scan(tmp_path, inputs, capsys, guesses=("minao", "atom"))
    SyntheticCalculator.failures = {"minao": RuntimeError}
    code, failed = cli_json(capsys, "state-scan-run", output)
    assert code == 2
    assert failed["failed_guesses"] == ["minao"] and failed["pending_guesses"] == ["atom"]
    assert [call["guess"] for call in SyntheticCalculator.calls] == ["minao"]
    failure_path = output / "attempts/minao/result.json"
    failure_bytes = failure_path.read_bytes()
    assert json.loads(failure_bytes)["error"]["type"] == "RuntimeError"

    SyntheticCalculator.failures = {}
    code, resumed = cli_json(capsys, "state-scan-run", output)
    assert code == 2
    assert resumed["all_requested_attempts_finished"] is True
    assert resumed["all_requested_guesses_converged"] is False
    assert resumed["failed_guesses"] == ["minao"] and resumed["converged_guesses"] == ["atom"]
    assert cli_json(capsys, "state-scan-run", output)[0] == 2
    before = saved_bytes(output)
    code, report = cli_json(capsys, "state-scan-report", output)
    assert code == 0 and report == resumed
    assert saved_bytes(output) == before
    assert failure_path.read_bytes() == failure_bytes
    assert [call["guess"] for call in SyntheticCalculator.calls] == ["minao", "atom"]
    assert_uncertified(report)


def test_interrupt_is_saved_returns_130_and_continuation_skips_it(tmp_path, inputs, capsys):
    output = create_scan(tmp_path, inputs, capsys, guesses=("minao", "atom", "1e"))
    SyntheticCalculator.failures = {"atom": KeyboardInterrupt}
    assert main(["state-scan-run", str(output), "--max-jobs", "3"]) == 130
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "interrupted" in captured.err.lower() and "Traceback" not in captured.err
    interrupted_path = output / "attempts/atom/result.json"
    interrupted_bytes = interrupted_path.read_bytes()
    assert json.loads(interrupted_bytes)["status"] == "interrupted"
    before = saved_bytes(output)
    code, report = cli_json(capsys, "state-scan-report", output)
    assert code == 0
    assert report["converged_guesses"] == ["minao"]
    assert report["interrupted_guesses"] == ["atom"] and report["pending_guesses"] == ["1e"]
    assert saved_bytes(output) == before

    SyntheticCalculator.failures = {}
    code, resumed = cli_json(capsys, "state-scan-run", output)
    assert code == 2
    assert resumed["converged_guesses"] == ["minao", "1e"]
    assert resumed["interrupted_guesses"] == ["atom"] and resumed["pending_guesses"] == []
    assert interrupted_path.read_bytes() == interrupted_bytes
    assert [call["guess"] for call in SyntheticCalculator.calls] == ["minao", "atom", "1e"]


def test_read_only_report_does_not_recover_crash_and_run_preserves_abandoned_bytes(tmp_path, inputs, capsys):
    output = create_scan(tmp_path, inputs, capsys, guesses=("minao", "atom"))
    ledger_path = output / "scan.json"
    ledger = json.loads(ledger_path.read_text())
    ledger["attempts"]["minao"] = {"status": "running", "output": "attempts/minao"}
    ledger_path.write_text(json.dumps(ledger))
    attempt_path = output / "attempts/minao/result.json"
    attempt_path.parent.mkdir(parents=True)
    incomplete = b'{"status": "running", "unfinished":'
    attempt_path.write_bytes(incomplete)
    before = saved_bytes(output)
    code, report = cli_json(capsys, "state-scan-report", output)
    assert code == 0 and report["rows"][0]["status"] == "running"
    assert saved_bytes(output) == before
    assert SyntheticCalculator.instances == []

    code, recovered = cli_json(capsys, "state-scan-run", output)
    assert code == 2
    assert recovered["abandoned_guesses"] == ["minao"]
    assert recovered["converged_guesses"] == ["atom"]
    assert attempt_path.read_bytes() == incomplete
    assert [call["guess"] for call in SyntheticCalculator.calls] == ["atom"]
    assert cli_json(capsys, "state-scan-report", output)[0] == 0


@pytest.mark.parametrize("command", ["state-scan-report", "state-scan-run"])
def test_changed_saved_result_is_rejected_without_new_work_or_repair(tmp_path, inputs, capsys, command):
    output = create_scan(tmp_path, inputs, capsys, guesses=("minao", "atom"))
    assert cli_json(capsys, "state-scan-run", output)[0] == 0
    result_path = output / "attempts/minao/result.json"
    result_path.write_bytes(result_path.read_bytes() + b"\n")
    before = saved_bytes(output)
    code = main([command, str(output)])
    captured = capsys.readouterr()
    assert code == 2 and captured.out == ""
    assert "changed" in captured.err and "Traceback" not in captured.err
    assert saved_bytes(output) == before
    assert [call["guess"] for call in SyntheticCalculator.calls] == ["minao"]


@pytest.mark.parametrize("max_jobs", ["0", "5", "1.5"])
def test_invalid_job_bound_is_rejected_without_attempts(tmp_path, inputs, capsys, max_jobs):
    output = create_scan(tmp_path, inputs, capsys)
    before = saved_bytes(output)
    arguments = ["state-scan-run", str(output), "--max-jobs", max_jobs]
    if max_jobs == "1.5":
        with pytest.raises(SystemExit) as error:
            main(arguments)
        assert error.value.code == 2
    else:
        assert main(arguments) == 2
    captured = capsys.readouterr()
    assert captured.out == "" and captured.err
    assert saved_bytes(output) == before
    assert SyntheticCalculator.instances == []


def test_module_entrypoint_reports_saved_failure_from_another_directory(tmp_path, inputs, capsys):
    output = create_scan(tmp_path, inputs, capsys, guesses=("minao",))
    SyntheticCalculator.failures = {"minao": RuntimeError}
    assert cli_json(capsys, "state-scan-run", output)[0] == 2
    source, settings_path, _ = inputs
    source.unlink()
    settings_path.unlink()
    before = saved_bytes(output)
    environment = os.environ.copy()
    project = str(Path(__file__).resolve().parents[1])
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, (project, environment.get("PYTHONPATH"))))
    process = subprocess.run(
        [sys.executable, "-m", "nanodesign", "state-scan-report", str(output)],
        cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=30,
    )
    assert process.returncode == 0, process.stderr
    assert process.stderr == ""
    report = json.loads(process.stdout)
    assert report["failed_guesses"] == ["minao"] and report["pending_guesses"] == []
    assert report["scan"] == str(output.resolve())
    assert_uncertified(report)
    assert saved_bytes(output) == before
