"""Characterization integration contracts using analytical forces only in tests."""

from dataclasses import asdict
import json

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms
from ase.io import read, write
import numpy as np
import pytest

from nanodesign import workflow
from nanodesign.design import sha256
from nanodesign.quantum import QuantumSettings


class TrackedHarmonic(Calculator):
    """A known Hessian with deliberately large forces on the frozen anchor."""

    implemented_properties = ["energy", "forces"]
    instances = []
    fail_at = None
    failure_type = RuntimeError

    def __init__(self, settings=None, event_log=None):
        super().__init__()
        self.settings, self.event_log = settings, event_log
        self.positions_evaluated = []
        self.diagnostics = {"test_surface": True, "evaluations": 0}
        self.instances.append(self)

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.positions_evaluated.append(atoms.positions.copy())
        number = len(self.positions_evaluated)
        self.diagnostics = {"test_surface": True, "evaluations": number}
        if number == self.fail_at:
            self.diagnostics["failed"] = True
            raise self.failure_type("test force evaluation failed")
        delta = atoms.positions[1] - np.array([1, 0, 0])
        self.results = {
            "energy": float(delta @ delta / 2),
            "forces": np.array([[0, 0, 8], -delta]),
        }


@pytest.fixture
def setup(tmp_path, monkeypatch):
    TrackedHarmonic.instances = []
    monkeypatch.setattr(workflow, "PySCFCalculator", TrackedHarmonic)
    atoms = Atoms("HeH", positions=[[0, 0, 2], [1, 0, 0]])
    write(tmp_path / "initial.xyz", atoms)
    write(tmp_path / "final.xyz", atoms)
    design = tmp_path / "design.json"
    design.write_text(json.dumps({
        "schema_version": 1, "length_unit": "angstrom",
        "initial": "initial.xyz", "final": "final.xyz", "fixed_indices": [0],
        "quantum": asdict(QuantumSettings()),
    }))
    return design, atoms, tmp_path / "characterization"


def test_constraint_replacement_cost_and_baseline_diagnostics_are_accurate(setup, tmp_path):
    design, atoms, out = setup
    # The design owns boundary conditions even when a supplied file stores a
    # conflicting constraint and cached outputs from another calculation.
    atoms.set_constraint(FixAtoms(indices=[0, 1]))
    atoms.calc = SinglePointCalculator(atoms, energy=-999, forces=np.full((2, 3), 99.0))
    source = tmp_path / "source.extxyz"
    write(source, atoms)
    result = workflow.run_characterization(design, source, out, step=0.002, fmax=0.01)
    calculator = TrackedHarmonic.instances[0]
    assert result["status"] == "completed"
    assert result["initial_free_force_max_ev_per_angstrom"] == 0
    assert result["initial_force_guard_passed"] is True
    assert result["expected_force_evaluations"] == len(calculator.positions_evaluated) == 7
    assert result["stationary"]["force_requests"] == 7
    assert result["stationary"]["free_atom_indices"] == [1]
    assert result["stationary"]["frozen_atom_indices"] == [0]
    assert result["quantum_diagnostics"]["evaluations"] == 1
    assert result["characterization_settings"] == result["stationary"]["settings"]
    assert calculator.event_log == out / "electronic.jsonl"
    assert all(np.array_equal(positions[0], atoms.positions[0]) for positions in calculator.positions_evaluated)
    archived = read(out / "input.extxyz")
    assert archived.calc is None
    assert archived.constraints[0].get_indices().tolist() == [0]
    np.testing.assert_allclose(result["stationary"]["hessian_ev_per_angstrom2"], np.eye(3), atol=1e-12)
    saved = json.loads((out / "result.json").read_text())
    assert saved == result
    assert saved["input_hashes"]["structure_sha256"] == sha256(source)
    assert saved["structure_source"] == str(source.resolve())
    assert saved["validation"]["transition_state_validated"] is False


@pytest.mark.parametrize("suffix", [".extxyz", ".traj"])
def test_selected_frame_and_hash_share_snapshot_despite_source_edit(setup, tmp_path, monkeypatch, suffix):
    design, atoms, out = setup
    other = atoms.copy()
    other.positions[1, 0] = 1.2
    source = tmp_path / ("frames" + suffix)
    write(source, [other, atoms])
    captured_hash = sha256(source)
    real_reader = workflow.read_coordinate_snapshot

    def parse_then_edit(path, raw, **kwargs):
        selected = real_reader(path, raw, **kwargs)
        write(path, other)
        return selected

    monkeypatch.setattr(workflow, "read_coordinate_snapshot", parse_then_edit)
    result = workflow.run_characterization(design, source, out, image=1, step=.002)
    assert result["status"] == "completed"
    assert result["structure_image"] == 1
    assert result["input_hashes"]["structure_sha256"] == captured_hash
    assert sha256(source) != captured_hash
    np.testing.assert_array_equal(read(out / "input.extxyz").positions, atoms.positions)


def test_force_guard_archives_residual_and_requested_settings_without_displacements(setup, tmp_path):
    design, atoms, out = setup
    atoms.positions[1, 0] = 1.2
    atoms.set_constraint(FixAtoms(indices=[0, 1]))
    source = tmp_path / "nonstationary.extxyz"
    write(source, atoms)
    with pytest.raises(ValueError, match="not stationary"):
        workflow.run_characterization(design, source, out, step=0.002, fmax=0.01, frequency_tolerance=15)
    saved = json.loads((out / "result.json").read_text())
    assert saved["status"] == "failed"
    assert saved["initial_free_force_max_ev_per_angstrom"] == pytest.approx(0.2)
    assert saved["initial_force_guard_passed"] is False
    assert saved["characterization_settings"]["force_tolerance_ev_per_angstrom"] == 0.01
    assert saved["characterization_settings"]["step_angstrom"] == 0.002
    assert saved["characterization_settings"]["frequency_tolerance_cm1"] == 15
    assert saved["quantum_diagnostics"]["evaluations"] == 1
    assert len(TrackedHarmonic.instances[0].positions_evaluated) == 1
    assert "stationary" not in saved


@pytest.mark.parametrize("options", [
    {"step": True}, {"step": "0.005"}, {"step": np.nan},
    {"fmax": True}, {"fmax": -1}, {"frequency_tolerance": True},
    {"max_free_coordinates": True}, {"image": True}, {"image": ":"}, {"image": 0.5},
])
def test_invalid_controls_are_rejected_before_calculator_or_output_creation(setup, tmp_path, options):
    design, _, out = setup
    with pytest.raises(ValueError):
        workflow.run_characterization(design, tmp_path / "initial.xyz", out, **options)
    assert not TrackedHarmonic.instances
    assert not out.exists()


def test_cost_limit_rejects_before_calculator_creation(setup, tmp_path):
    design, _, out = setup
    with pytest.raises(ValueError, match="3 free coordinates"):
        workflow.run_characterization(design, tmp_path / "initial.xyz", out, max_free_coordinates=2)
    assert not TrackedHarmonic.instances
    assert not out.exists()


@pytest.mark.parametrize("failure_type,status", [(RuntimeError, "failed"), (KeyboardInterrupt, "interrupted")])
@pytest.mark.parametrize("fail_at", [1, 3])
def test_initial_and_displaced_failures_preserve_reproducible_record(setup, tmp_path, monkeypatch, failure_type, status, fail_at):
    design, _, out = setup
    monkeypatch.setattr(TrackedHarmonic, "failure_type", failure_type)
    monkeypatch.setattr(TrackedHarmonic, "fail_at", fail_at)
    with pytest.raises(failure_type):
        workflow.run_characterization(design, tmp_path / "initial.xyz", out, step=0.003, fmax=0.012)
    saved = json.loads((out / "result.json").read_text())
    assert saved["status"] == status
    assert saved["characterization_settings"]["step_angstrom"] == 0.003
    assert saved["characterization_settings"]["force_tolerance_ev_per_angstrom"] == 0.012
    assert saved["quantum_diagnostics"]["evaluations"] == 1
    assert saved["validation"]["transition_state_validated"] is False
    if fail_at == 1:
        assert saved["quantum_diagnostics"]["failed"] is True
        assert "initial_force_guard_passed" not in saved
    else:
        assert saved["initial_force_guard_passed"] is True
    assert "stationary" not in saved


def test_calculator_initialization_failure_is_not_left_running(setup, tmp_path, monkeypatch):
    design, _, out = setup

    def broken_calculator(*args, **kwargs):
        raise RuntimeError("initialization failed")

    monkeypatch.setattr(workflow, "PySCFCalculator", broken_calculator)
    with pytest.raises(RuntimeError, match="initialization"):
        workflow.run_characterization(design, tmp_path / "initial.xyz", out)
    saved = json.loads((out / "result.json").read_text())
    assert saved["status"] == "failed"
    assert "quantum_diagnostics" not in saved


def test_selected_trajectory_frame_is_both_used_and_identified(setup, tmp_path):
    design, atoms, out = setup
    displaced = atoms.copy()
    displaced.positions[1, 0] = 1.2
    source = tmp_path / "frames.extxyz"
    write(source, [displaced, atoms])
    result = workflow.run_characterization(design, source, out, image=1)
    assert result["structure_image"] == 1
    assert result["input_hashes"]["structure_sha256"] == sha256(source)
    np.testing.assert_array_equal(TrackedHarmonic.instances[0].positions_evaluated[0], atoms.positions)


def test_changed_anchor_is_rejected_before_any_calculation(setup, tmp_path):
    design, atoms, out = setup
    atoms.positions[0, 0] = 0.01
    source = tmp_path / "wrong-anchor.extxyz"
    write(source, atoms)
    with pytest.raises(ValueError, match="anchors"):
        workflow.run_characterization(design, source, out)
    assert not TrackedHarmonic.instances
    assert not out.exists()
