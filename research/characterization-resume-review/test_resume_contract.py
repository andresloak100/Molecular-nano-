"""Independent H1 interruption/provenance checks; no electronic calculations."""
from copy import deepcopy
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path

from ase.constraints import FixAtoms
import numpy as np
import pytest

from nanodesign.quantum import QuantumSettings


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fixtures = load_module("h2_review_fixtures", HERE / "fixtures.py")
h1 = load_module("h2_reviewed_h1", ROOT / "research/characterization-resume/prototype.py")
SETTINGS = QuantumSettings(spin=1, basis="sto-3g", dispersion=None, grid_level=0, threads=1)
CONTEXT = {"fixture": "H2 review only; not chemistry", "input_hashes": {"source": "synthetic-source-A"}}


@pytest.fixture(autouse=True)
def mock_only_backend(monkeypatch):
    # H1 deliberately accepts one exact supported backend class. Replace it
    # with our declared analytic fixture, never with a real quantum calculator.
    monkeypatch.setattr(h1, "PySCFCalculator", fixtures.CountingForceFixture)


def run_fixture(output, *, source=None, fail_at=None, failure=KeyboardInterrupt,
                settings=SETTINGS, context=None, atoms=None, **controls):
    atoms = fixtures.anchored_structure() if atoms is None else atoms
    atoms.calc = calculator = fixtures.CountingForceFixture(
        atoms.positions, fail_at=fail_at, failure=failure, label=Path(output).name,
        settings=settings)
    result = h1.run_characterization_checkpoint(
        atoms, output, settings=settings, input_context=CONTEXT if context is None else context,
        resume_from=source, **controls)
    return result, calculator


def interrupted(output, fail_at, failure=KeyboardInterrupt, source=None):
    atoms = fixtures.anchored_structure()
    atoms.calc = calculator = fixtures.CountingForceFixture(
        atoms.positions, fail_at=fail_at, failure=failure, label=Path(output).name,
        settings=SETTINGS)
    before = atoms.positions.copy()
    with pytest.raises(failure, match="synthetic interruption"):
        h1.run_characterization_checkpoint(atoms, output, settings=SETTINGS,
                                          input_context=CONTEXT, resume_from=source)
    np.testing.assert_array_equal(atoms.positions, before)
    assert atoms.calc is calculator and calculator.atoms is None and calculator.results == {}
    assert calculator.diagnostics == {"synthetic_fixture": True}
    result = json.loads((output / "result.json").read_text())
    checkpoint, _ = h1.load_checkpoint(output)
    return result, checkpoint, calculator


def expected_geometry_sequence(atoms, step=.005):
    """Independent ordering oracle, built without H1's manifest helper."""
    reference = atoms.positions.copy()
    sequence = [reference.copy()]
    for atom in (0, 2):
        for axis in range(3):
            for sign in (1, -1):
                point = reference.copy()
                point[atom, axis] += sign * step
                sequence.append(point)
    return sequence


@pytest.mark.parametrize("fail_at", [1, 2, 3], ids=["baseline", "first-plus", "first-minus"])
@pytest.mark.parametrize("failure,status", [(KeyboardInterrupt, "interrupted"), (RuntimeError, "failed")])
def test_each_interruption_reuses_exact_successful_prefix_only(tmp_path, fail_at, failure, status):
    source = tmp_path / "interrupted"
    saved, payload, initial_calculator = interrupted(source, fail_at, failure)
    kept = fail_at - 1
    assert saved["status"] == payload["status"] == status
    assert len(payload["records"]) == kept
    assert saved["checkpoint"]["new_force_evaluations"] == fail_at
    failed = payload["failed_calculation"]
    assert failed["reusable"] is False
    assert failed["request"] == payload["manifest"]["requests"][kept]
    assert failed["quantum_diagnostics"]["deliberate_failure"] is True
    assert failed["quantum_diagnostics"]["scf_converged"] is False
    assert failed["quantum_diagnostics"]["gradient_completed"] is False
    failed_call_id = failed["quantum_diagnostics"]["call_id"]
    assert failed_call_id == f"interrupted-{fail_at}"
    assert failed_call_id not in {record["calculation_call_id"] for record in payload["records"]}
    assert failed["error"]["type"] == failure.__name__
    before = fixtures.snapshot_tree(source)
    result, resumed = run_fixture(tmp_path / "resumed", source=source)
    assert fixtures.snapshot_tree(source) == before
    assert result["checkpoint"]["reused_force_records"] == kept
    assert result["checkpoint"]["new_force_evaluations"] == len(resumed.calls) == 13 - kept
    assert result["checkpoint"]["completed_force_records"] == result["checkpoint"]["planned_force_records"] == 13
    assert len(initial_calculator.calls) + len(resumed.calls) == 14
    oracle = expected_geometry_sequence(fixtures.anchored_structure())
    np.testing.assert_array_equal(resumed.calls, oracle[kept:])
    final_payload, _ = h1.load_checkpoint(tmp_path / "resumed")
    assert final_payload["records"][:kept] == payload["records"]
    for record in final_payload["records"]:
        assert len(record["forces_ev_per_angstrom"]) == 3
        np.testing.assert_array_equal(record["forces_ev_per_angstrom"][1], [3, 0, -2])
    assert result["validation"]["transition_state_validated"] is False
    assert result["validation"]["electronic_state_identity_verified"] is False


def test_successive_failed_attempts_preserve_lineage_and_finish_without_repeating_saved_calls(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    _, first_payload, first_calc = interrupted(first, 2)
    first_bytes = fixtures.snapshot_tree(first)
    _, second_payload, second_calc = interrupted(second, 2, RuntimeError, source=first)
    second_bytes = fixtures.snapshot_tree(second)
    assert len(second_payload["records"]) == 2
    assert second_payload["records"][0] == first_payload["records"][0]
    result, final_calc = run_fixture(tmp_path / "third", source=second)
    assert len(first_calc.calls) + len(second_calc.calls) + len(final_calc.calls) == 15
    assert len(final_calc.calls) == result["checkpoint"]["new_force_evaluations"] == 11
    assert fixtures.snapshot_tree(first) == first_bytes
    assert fixtures.snapshot_tree(second) == second_bytes
    final_payload, _ = h1.load_checkpoint(tmp_path / "third")
    assert final_payload["records"][:2] == second_payload["records"]
    assert final_payload["resume_source"]["attempt_id"] == second_payload["attempt_id"]


def test_completed_checkpoint_replays_with_zero_new_solver_calls(tmp_path):
    source = tmp_path / "complete"
    _, original = run_fixture(source)
    assert len(original.calls) == 13
    before = fixtures.snapshot_tree(source)
    result, calculator = run_fixture(tmp_path / "replay", source=source)
    assert calculator.calls == []
    assert result["checkpoint"]["reused_force_records"] == 13
    assert result["checkpoint"]["new_force_evaluations"] == 0
    assert fixtures.snapshot_tree(source) == before


def test_nonstationary_baseline_is_saved_and_resume_stops_without_more_calls(tmp_path):
    source = tmp_path / "nonstationary"
    for output, prior, expected_calls in [(source, None, 1), (tmp_path / "rejected-again", source, 0)]:
        atoms = fixtures.anchored_structure()
        shifted_equilibrium = atoms.positions.copy()
        shifted_equilibrium[0, 0] += .1
        atoms.calc = calculator = fixtures.CountingForceFixture(shifted_equilibrium, settings=SETTINGS)
        before = fixtures.snapshot_tree(source) if prior else None
        with pytest.raises(h1.CheckpointError, match="Baseline is not stationary"):
            h1.run_characterization_checkpoint(atoms, output, settings=SETTINGS,
                                              input_context=CONTEXT, resume_from=prior)
        assert len(calculator.calls) == expected_calls
        payload, _ = h1.load_checkpoint(output)
        result = json.loads((output / "result.json").read_text())
        assert len(payload["records"]) == 1
        assert payload["status"] == "failed"
        assert result["initial_force_guard_passed"] is False
        assert result["checkpoint"]["new_force_evaluations"] == expected_calls
        # Twelve slots are absent, but the force guard permits zero new work
        # with these unchanged inputs. A missing-slot count is not a job plan.
        assert len(payload["manifest"]["requests"]) - len(payload["records"]) == 12
        if prior:
            assert fixtures.snapshot_tree(source) == before


def test_failed_atomic_replacement_keeps_last_durable_prefix_and_retries_unpersisted_work(tmp_path, monkeypatch):
    source = tmp_path / "write-failure"
    atoms = fixtures.anchored_structure()
    atoms.calc = calculator = fixtures.CountingForceFixture(atoms.positions, settings=SETTINGS)
    real_replace = h1.os.replace

    def fail_after_first_displacement(temporary, destination):
        if Path(destination).name == "force_checkpoint.json":
            pending = json.loads(Path(temporary).read_text())["payload"]
            if len(pending["records"]) >= 2:
                raise OSError("synthetic persistent checkpoint replacement failure")
        return real_replace(temporary, destination)

    with monkeypatch.context() as patch:
        patch.setattr(h1.os, "replace", fail_after_first_displacement)
        with pytest.raises(OSError, match="synthetic persistent"):
            h1.run_characterization_checkpoint(atoms, source, settings=SETTINGS, input_context=CONTEXT)
    durable, _ = h1.load_checkpoint(source)
    assert durable["status"] == "running"
    assert len(durable["records"]) == 1
    assert len(calculator.calls) == 2
    assert not list(source.glob(".checkpoint-*"))
    before = fixtures.snapshot_tree(source)
    result, resumed = run_fixture(tmp_path / "after-write-failure", source=source)
    assert result["checkpoint"]["reused_force_records"] == 1
    assert len(resumed.calls) == result["checkpoint"]["new_force_evaluations"] == 12
    np.testing.assert_array_equal(resumed.calls[0], calculator.calls[-1])
    assert fixtures.snapshot_tree(source) == before


@pytest.mark.parametrize("change", ["positions", "elements", "masses", "anchors", "cell", "charges", "magmoms",
                                    "quantum", "guess", "context", "step", "fmax", "frequency", "limit", "versions"])
def test_changed_resume_identity_is_rejected_before_force_work(tmp_path, monkeypatch, change):
    source = tmp_path / "source"
    interrupted(source, 3)
    before = fixtures.snapshot_tree(source)
    atoms = fixtures.anchored_structure()
    settings = SETTINGS
    context = deepcopy(CONTEXT)
    controls = {}
    if change == "positions":
        atoms.positions[0, 0] += .001
    elif change == "elements":
        atoms.numbers[0] = 2
    elif change == "masses":
        masses = atoms.get_masses(); masses[0] += .1; atoms.set_masses(masses)
    elif change == "anchors":
        atoms.set_constraint(FixAtoms(indices=[0]))
    elif change == "cell":
        atoms.set_cell([5, 5, 5])
    elif change == "charges":
        atoms.set_initial_charges([.1, 0, 0])
    elif change == "magmoms":
        atoms.set_initial_magnetic_moments([1, 0, 0])
    elif change == "quantum":
        settings = replace(SETTINGS, basis="def2-svp")
    elif change == "guess":
        settings = replace(SETTINGS, scf_initial_guess="atom")
    elif change == "context":
        context["input_hashes"]["source"] = "synthetic-source-B"
    elif change == "versions":
        versions = h1._versions()
        monkeypatch.setattr(h1, "_versions", lambda: {**versions, "fixture_changed": True})
    else:
        controls = {"step": {"step": .004}, "fmax": {"fmax": .02},
                    "frequency": {"frequency_tolerance": 25}, "limit": {"max_free_coordinates": 100}}[change]
    atoms.calc = calculator = fixtures.CountingForceFixture(atoms.positions, settings=settings)
    with pytest.raises(h1.CheckpointError, match="manifest mismatch"):
        h1.run_characterization_checkpoint(atoms, tmp_path / "rejected", settings=settings,
                                          input_context=context, resume_from=source, **controls)
    assert calculator.calls == []
    assert not (tmp_path / "rejected").exists()
    assert fixtures.snapshot_tree(source) == before


def rewrite_envelope(path, payload):
    # Independently construct a checksum-consistent invalid record, so tests
    # exercise semantic guards instead of only the outer byte-integrity guard.
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    path.write_text(json.dumps({"payload": payload, "sha256": hashlib.sha256(canonical).hexdigest()}))


@pytest.mark.parametrize("fault", ["truncated", "digest", "wrong-sign", "duplicate-slot", "failed-slot",
                                  "missing-force", "boolean-force", "settings", "call-id", "false-completion"])
def test_bad_checkpoint_cannot_trigger_new_force_evaluation(tmp_path, fault):
    source = tmp_path / "source"
    interrupted(source, 4)
    path = source / "force_checkpoint.json"
    payload = json.loads(path.read_text())["payload"]
    if fault == "truncated":
        path.write_text('{"payload":')
    elif fault == "digest":
        document = json.loads(path.read_text()); document["sha256"] = "0" * 64
        path.write_text(json.dumps(document))
    else:
        if fault == "wrong-sign":
            payload["records"][1]["request"]["requested_offset_angstrom"] *= -1
        elif fault == "duplicate-slot":
            payload["records"][2] = deepcopy(payload["records"][1])
        elif fault == "failed-slot":
            payload["records"][0]["status"] = "failed"
        elif fault == "missing-force":
            payload["records"][0]["forces_ev_per_angstrom"].pop()
        elif fault == "boolean-force":
            payload["records"][0]["forces_ev_per_angstrom"][0][0] = True
        elif fault == "settings":
            payload["records"][0]["quantum_diagnostics"]["settings"]["scf_initial_guess"] = "atom"
        elif fault == "call-id":
            payload["records"][0]["calculation_call_id"] = "inconsistent"
        elif fault == "false-completion":
            payload["status"] = "completed"
        rewrite_envelope(path, payload)
    before = fixtures.snapshot_tree(source)
    atoms = fixtures.anchored_structure()
    atoms.calc = calculator = fixtures.CountingForceFixture(atoms.positions, settings=SETTINGS)
    with pytest.raises(h1.CheckpointError):
        h1.run_characterization_checkpoint(atoms, tmp_path / "rejected", settings=SETTINGS,
                                          input_context=CONTEXT, resume_from=source)
    assert calculator.calls == []
    assert fixtures.snapshot_tree(source) == before


def test_active_source_worker_is_rejected_without_mutating_it(tmp_path):
    import fcntl
    source = tmp_path / "source"
    interrupted(source, 3)
    before = fixtures.snapshot_tree(source)
    atoms = fixtures.anchored_structure()
    atoms.calc = calculator = fixtures.CountingForceFixture(atoms.positions, settings=SETTINGS)
    with (source / ".checkpoint.lock").open("rb") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            with pytest.raises(h1.CheckpointError, match="active worker"):
                h1.run_characterization_checkpoint(atoms, tmp_path / "rejected", settings=SETTINGS,
                                                  input_context=CONTEXT, resume_from=source)
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
    assert calculator.calls == []
    assert fixtures.snapshot_tree(source) == before
