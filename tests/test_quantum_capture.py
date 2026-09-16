"""Electronic-state output lifecycle, using synthetic solvers and capture only.

These tests execute the actual ASE calculator orchestration but replace every
PySCF entry point and snapshot producer. No electronic-structure work is run.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
from uuid import UUID

from ase import Atoms
from ase.constraints import FixAtoms
import numpy as np
import pytest

import nanodesign
from nanodesign.quantum import PySCFCalculator, QuantumCalculationError, QuantumSettings
from nanodesign.stationary import characterize_stationary_point


@pytest.fixture
def synthetic_backend(monkeypatch):
    state = SimpleNamespace(
        threads=7, initial_threads=7, restore_error=None,
        converged=True, energy=-1.0, spin_error=False,
        scf_error=None, capture_error=None, save_error=None,
        gradient_error=None, gradient_error_on=None,
        solvers=[], solves=[], captures=[], saves=[], gradient_calls=0, trace=[],
    )

    def num_threads(value=None):
        if value is None:
            return state.threads
        state.trace.append(("threads", value))
        if value == state.initial_threads and state.restore_error is not None:
            raise state.restore_error
        state.threads = value
        return value

    def molecule(**kwargs):
        # Fixtures use hydrogen only, so atom count equals nuclear charge.
        assert all(symbol == "H" for symbol, _ in kwargs["atom"])
        electrons = len(kwargs["atom"]) - kwargs["charge"]
        return SimpleNamespace(
            natm=len(kwargs["atom"]), nelectron=electrons,
            nelec=((electrons + kwargs["spin"]) // 2, (electrons - kwargs["spin"]) // 2),
            spin=kwargs["spin"], charge=kwargs["charge"], atom=deepcopy(kwargs["atom"]),
            nao_nr=lambda: len(kwargs["atom"]),
        )

    class Gradient:
        def __init__(self, mean_field):
            self.mean_field = mean_field

        def kernel(self):
            state.gradient_calls += 1
            state.trace.append(("gradient", state.gradient_calls))
            if state.gradient_error is not None and (
                state.gradient_error_on is None or state.gradient_calls == state.gradient_error_on
            ):
                raise state.gradient_error
            return np.zeros((self.mean_field.mol.natm, 3))

    class MeanField:
        def __init__(self, mol, kind, fitted=False):
            self.mol, self.kind, self.fitted = mol, kind, fitted
            self.converged, self.cycles = state.converged, 1
            self.grids = SimpleNamespace(level=None)
            state.solvers.append(self)

        def density_fit(self):
            # A distinct final solver proves capture sees the post-DF object.
            fitted = MeanField(self.mol, self.kind, fitted=True)
            fitted.__dict__.update({key: value for key, value in self.__dict__.items()
                                   if key not in {"fitted"}})
            return fitted

        def kernel(self):
            state.solves.append(self)
            state.trace.append(("scf", len(state.solves)))
            if state.scf_error is not None:
                raise state.scf_error
            if getattr(self, "callback", None):
                self.callback({"cycle": 0, "e_tot": state.energy})
            return state.energy

        def spin_square(self):
            if state.spin_error:
                return float("nan"), 1.0
            spin = self.mol.spin / 2
            return spin * (spin + 1), self.mol.spin + 1

        def nuc_grad_method(self):
            return Gradient(self)

    pyscf = ModuleType("pyscf")
    pyscf.__path__ = []
    pyscf.__version__ = "synthetic-test-only"
    pyscf.gto = SimpleNamespace(M=molecule)
    pyscf.lib = SimpleNamespace(num_threads=num_threads, param=SimpleNamespace(BOHR=0.52917721092))
    pyscf.dft = SimpleNamespace(
        RKS=lambda mol: MeanField(mol, "RKS"), UKS=lambda mol: MeanField(mol, "UKS"),
        libxc=SimpleNamespace(libxc_version=lambda: "synthetic", test_deriv_order=lambda *a, **kw: None,
                              is_nlc=lambda code: False),
    )
    scf = ModuleType("pyscf.scf")
    scf.__path__ = []
    dispersion = ModuleType("pyscf.scf.dispersion")
    dispersion.parse_dft = lambda code: (code, None, None)
    scf.dispersion = dispersion
    pyscf.scf = scf
    for name, module in (("pyscf", pyscf), ("pyscf.scf", scf), ("pyscf.scf.dispersion", dispersion)):
        monkeypatch.setitem(sys.modules, name, module)

    electronic_state = ModuleType("nanodesign.electronic_state")

    def capture_snapshot(mean_field, *, settings, call_id):
        state.trace.append(("capture", call_id))
        state.captures.append((mean_field, deepcopy(settings), call_id))
        if state.capture_error is not None:
            raise state.capture_error
        return {"test_fixture": True, "call_id": call_id, "settings": deepcopy(settings),
                "atoms": deepcopy(mean_field.mol.atom), "scope": "converged SCF only"}

    def save_snapshot(path, snapshot):
        path = Path(path)
        state.trace.append(("save", snapshot["call_id"]))
        if state.save_error is not None:
            raise state.save_error
        raw = (json.dumps(snapshot, allow_nan=False) + "\n").encode()
        with path.open("xb") as stream:
            stream.write(raw)
        reference = {"path": str(path.resolve()), "sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}
        state.saves.append((path, deepcopy(snapshot), reference))
        return reference

    electronic_state.capture_snapshot = capture_snapshot
    electronic_state.save_snapshot = save_snapshot
    monkeypatch.setitem(sys.modules, "nanodesign.electronic_state", electronic_state)
    monkeypatch.setattr(nanodesign, "electronic_state", electronic_state, raising=False)
    return state


def attached(tmp_path, *, capture=True, spin=0, density_fit=False):
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    settings = QuantumSettings(spin=spin, charge=spin, basis="sto-3g", dispersion=None,
                               threads=2, density_fit=density_fit, scf_initial_guess="huckel")
    kwargs = {"electronic_state_directory": tmp_path / "states"} if capture else {}
    atoms.calc = calculator = PySCFCalculator(settings, event_log=tmp_path / "events.jsonl", **kwargs)
    return atoms, calculator


def event_records(tmp_path):
    path = tmp_path / "events.jsonl"
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def assert_snapshot(calculator, directory, *, accepted):
    reference = calculator.diagnostics["electronic_state_snapshot"]
    call_id = calculator.diagnostics["call_id"]
    UUID(call_id)
    assert reference["call_id"] == call_id
    assert reference["phase"] == "converged_scf_only"
    assert reference["whole_evaluation_accepted"] is accepted
    path = Path(reference["path"])
    assert path == directory / f"{call_id}.json"
    raw = path.read_bytes()
    assert reference["sha256"] == hashlib.sha256(raw).hexdigest()
    assert reference["size_bytes"] == len(raw)
    snapshot = json.loads(raw)
    assert snapshot["call_id"] == call_id
    # Later call acceptance is diagnostic metadata; the saved SCF bytes remain
    # unchanged and cannot retrospectively certify the gradient or whole call.
    assert "whole_evaluation_accepted" not in snapshot
    return path, snapshot


def test_default_calculation_does_not_capture_or_create_state_directory(tmp_path, synthetic_backend):
    atoms, calculator = attached(tmp_path, capture=False)
    atoms.get_forces()
    assert calculator.diagnostics["calculation_status"] == "completed"
    assert "electronic_state_snapshot" not in calculator.diagnostics
    assert synthetic_backend.captures == synthetic_backend.saves == []
    assert not (tmp_path / "states").exists()
    assert synthetic_backend.threads == synthetic_backend.initial_threads


@pytest.mark.parametrize("spin", [0, 1], ids=["rks", "uks"])
@pytest.mark.parametrize("density_fit", [False, True], ids=["direct", "density-fit"])
def test_capture_binds_final_solver_settings_geometry_and_call_id_and_respects_cache(
    tmp_path, synthetic_backend, spin, density_fit,
):
    atoms, calculator = attached(tmp_path, spin=spin, density_fit=density_fit)
    atoms.get_forces()
    first_id = calculator.diagnostics["call_id"]
    first_path, first_snapshot = assert_snapshot(calculator, tmp_path / "states", accepted=True)
    first_raw = first_path.read_bytes()
    events_before = (tmp_path / "events.jsonl").read_bytes()
    captured_solver, captured_settings, captured_id = synthetic_backend.captures[0]
    assert captured_solver is synthetic_backend.solves[0]
    assert captured_solver.kind == ("RKS" if spin == 0 else "UKS")
    assert captured_solver.fitted is density_fit
    assert captured_settings == calculator.settings.to_dict() == first_snapshot["settings"]
    assert captured_id == first_id
    assert synthetic_backend.trace.index(("capture", first_id)) < synthetic_backend.trace.index(("gradient", 1))
    atoms.get_potential_energy()
    atoms.get_forces()
    assert len(synthetic_backend.solves) == len(synthetic_backend.captures) == 1
    assert first_path.read_bytes() == first_raw
    assert (tmp_path / "events.jsonl").read_bytes() == events_before

    atoms.positions[1, 0] += 0.1
    atoms.get_forces()
    second_path, second_snapshot = assert_snapshot(calculator, tmp_path / "states", accepted=True)
    assert calculator.diagnostics["call_id"] != first_id
    assert second_snapshot["atoms"] != first_snapshot["atoms"]
    assert second_path != first_path and first_path.read_bytes() == first_raw
    assert len(synthetic_backend.solves) == len(synthetic_backend.saves) == 2
    assert {path.name for path in (tmp_path / "states").iterdir()} == {first_path.name, second_path.name}
    assert calculator.diagnostics["calculation_status"] == "completed"
    assert {record["call_id"] for record in event_records(tmp_path)} == {first_id, calculator.diagnostics["call_id"]}


@pytest.mark.parametrize("failure", ["invalid_geometry", "nonconverged"])
def test_later_geometry_failure_clears_old_reference_without_erasing_old_artifact(tmp_path, synthetic_backend, failure):
    atoms, calculator = attached(tmp_path)
    atoms.get_forces()
    old_path, _ = assert_snapshot(calculator, tmp_path / "states", accepted=True)
    old_id, old_bytes = calculator.diagnostics["call_id"], old_path.read_bytes()
    if failure == "invalid_geometry":
        atoms.positions[1] = atoms.positions[0]
        expected_error = ValueError
    else:
        atoms.positions[1, 0] += 0.1
        synthetic_backend.converged = False
        expected_error = QuantumCalculationError
    with pytest.raises(expected_error):
        atoms.get_forces()
    assert calculator.results == {}
    assert calculator.diagnostics["calculation_status"] == "failed"
    assert calculator.diagnostics["call_id"] != old_id
    assert "electronic_state_snapshot" not in calculator.diagnostics
    assert len(synthetic_backend.captures) == 1
    assert list((tmp_path / "states").iterdir()) == [old_path]
    assert old_path.read_bytes() == old_bytes


@pytest.mark.parametrize("failure", ["unconverged", "nonfinite_energy", "nonfinite_spin"])
def test_invalid_scf_or_spin_never_creates_snapshot(tmp_path, synthetic_backend, failure):
    if failure == "unconverged":
        synthetic_backend.converged = False
    elif failure == "nonfinite_energy":
        synthetic_backend.energy = float("nan")
    else:
        synthetic_backend.spin_error = True
    atoms, calculator = attached(tmp_path)
    with pytest.raises(QuantumCalculationError):
        atoms.get_forces()
    assert synthetic_backend.captures == [] and synthetic_backend.gradient_calls == 0
    assert calculator.results == {}
    assert calculator.diagnostics["calculation_status"] == "failed"
    assert "electronic_state_snapshot" not in calculator.diagnostics
    assert not (tmp_path / "states").exists()
    assert synthetic_backend.threads == synthetic_backend.initial_threads


def test_gradient_failure_preserves_current_scf_evidence_but_no_accepted_results(tmp_path, synthetic_backend):
    synthetic_backend.gradient_error = RuntimeError("synthetic gradient failure")
    atoms, calculator = attached(tmp_path)
    with pytest.raises(QuantumCalculationError, match="synthetic gradient failure"):
        atoms.get_forces()
    assert_snapshot(calculator, tmp_path / "states", accepted=False)
    assert calculator.diagnostics["scf_converged"] is True
    assert calculator.diagnostics["gradient_completed"] is False
    assert calculator.diagnostics["calculation_status"] == "failed"
    assert calculator.results == {}
    assert synthetic_backend.threads == synthetic_backend.initial_threads
    events = event_records(tmp_path)
    assert events[-1]["event"] == "calculation_failed"
    assert all(record["event"] != "calculation_completed" for record in events)


@pytest.mark.parametrize("stage", ["capture", "save"])
def test_requested_capture_io_failure_stops_before_gradient(tmp_path, synthetic_backend, stage):
    setattr(synthetic_backend, f"{stage}_error", OSError("synthetic capture destination unavailable"))
    atoms, calculator = attached(tmp_path)
    with pytest.raises(QuantumCalculationError, match="capture destination unavailable"):
        atoms.get_forces()
    assert calculator.results == {}
    assert calculator.diagnostics["calculation_status"] == "failed"
    assert "capture destination unavailable" in calculator.diagnostics["capture_error"]
    assert "electronic_state_snapshot" not in calculator.diagnostics
    assert synthetic_backend.gradient_calls == 0 and synthetic_backend.saves == []
    assert synthetic_backend.threads == synthetic_backend.initial_threads


@pytest.mark.parametrize("interruption", [KeyboardInterrupt, SystemExit])
@pytest.mark.parametrize("stage", ["scf", "save", "gradient"])
def test_interruption_cleans_results_and_preserves_only_written_current_scf_evidence(
    tmp_path, synthetic_backend, interruption, stage,
):
    setattr(synthetic_backend, f"{stage}_error", interruption("synthetic interruption"))
    atoms, calculator = attached(tmp_path)
    with pytest.raises(interruption, match="synthetic interruption"):
        atoms.get_forces()
    assert calculator.results == {}
    assert calculator.diagnostics["calculation_status"] == "interrupted"
    if stage == "gradient":
        assert_snapshot(calculator, tmp_path / "states", accepted=False)
    else:
        assert "electronic_state_snapshot" not in calculator.diagnostics
        assert synthetic_backend.saves == []
    assert synthetic_backend.threads == synthetic_backend.initial_threads
    assert event_records(tmp_path)[-1]["event"] == "calculation_interrupted"


def test_completion_event_follows_thread_restoration_and_log_failure_rejects_call(tmp_path, synthetic_backend, monkeypatch):
    atoms, calculator = attached(tmp_path)
    original_event = calculator._event
    observed_completion = []

    def event(name, *args, **kwargs):
        if name == "calculation_completed":
            observed_completion.append(synthetic_backend.threads)
            raise OSError("synthetic completion log failure")
        return original_event(name, *args, **kwargs)

    monkeypatch.setattr(calculator, "_event", event)
    with pytest.raises(QuantumCalculationError, match="completion log failure"):
        atoms.get_forces()
    assert observed_completion == [synthetic_backend.initial_threads]
    assert_snapshot(calculator, tmp_path / "states", accepted=False)
    assert calculator.diagnostics["gradient_completed"] is True
    assert calculator.diagnostics["calculation_status"] == "failed"
    assert calculator.results == {}
    assert event_records(tmp_path)[-1]["event"] == "calculation_failed"


def test_thread_restoration_failure_rejects_call_without_completion_event(tmp_path, synthetic_backend):
    synthetic_backend.restore_error = RuntimeError("synthetic thread restoration failure")
    atoms, calculator = attached(tmp_path)
    with pytest.raises(QuantumCalculationError, match="thread restoration failure"):
        atoms.get_forces()
    assert calculator.results == {}
    assert calculator.diagnostics["calculation_status"] == "failed"
    assert_snapshot(calculator, tmp_path / "states", accepted=False)
    assert calculator.diagnostics["thread_restoration_error"]["type"] == "RuntimeError"
    assert all(record["event"] != "calculation_completed" for record in event_records(tmp_path))


@pytest.mark.parametrize("displaced_failure", [False, True])
def test_stationary_restores_baseline_cache_and_reference_while_retaining_new_artifacts(
    tmp_path, synthetic_backend, displaced_failure,
):
    atoms, calculator = attached(tmp_path)
    atoms.set_constraint(FixAtoms(indices=[0]))
    baseline_forces = atoms.get_forces()
    baseline_energy = atoms.get_potential_energy()
    baseline_diagnostics = deepcopy(calculator.diagnostics)
    baseline_path, _ = assert_snapshot(calculator, tmp_path / "states", accepted=True)
    baseline_raw = baseline_path.read_bytes()
    if displaced_failure:
        synthetic_backend.gradient_error = RuntimeError("synthetic displaced gradient failure")
        synthetic_backend.gradient_error_on = 2
        with pytest.raises(QuantumCalculationError, match="displaced gradient failure"):
            characterize_stationary_point(atoms, max_free_coordinates=3)
        assert len(synthetic_backend.solves) == 2
        expected_files = 2
    else:
        characterization = characterize_stationary_point(atoms, max_free_coordinates=3)
        assert characterization["finite_difference_evidence"]["baseline_calculation_call_id"] == baseline_diagnostics["call_id"]
        displacements = characterization["finite_difference_evidence"]["displacements"]
        assert len(displacements) == 6
        assert len({item["calculation_call_id"] for item in displacements}) == 6
        assert len(synthetic_backend.solves) == 7
        expected_files = 7
    assert calculator.diagnostics == baseline_diagnostics
    assert baseline_path.read_bytes() == baseline_raw
    assert len(list((tmp_path / "states").glob("*.json"))) == expected_files
    solves_before_cache_reads = len(synthetic_backend.solves)
    np.testing.assert_array_equal(atoms.get_forces(), baseline_forces)
    assert atoms.get_potential_energy() == baseline_energy
    assert len(synthetic_backend.solves) == solves_before_cache_reads
