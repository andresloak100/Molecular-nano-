"""Synthetic workflow checks. These are software tests, not molecular evidence.

Run from the repository root:
    .venv/bin/python -m pytest -q research/path-contract-audit/test_path_contracts.py
"""
from dataclasses import asdict
import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write

from nanodesign import workflow
from nanodesign.quantum import QuantumSettings


class LoadedDoubleWell(Calculator):
    """A known barrier plus opposing forces on two fixed anchors."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, settings=None, event_log=None):
        super().__init__()
        self.diagnostics = {"test_surface": True, "physical_model": False}

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x, y, z = atoms.positions[2]
        anchor_shift = atoms.positions[0, 0] - atoms.positions[1, 0]
        forces = np.zeros((3, 3))
        forces[0, 0], forces[1, 0] = 2.0, -2.0
        forces[2] = [-0.8 * x * (x * x - 1), -y, -z]
        self.results = {
            "energy": 0.2 * (x * x - 1) ** 2 + 0.5 * (y*y + z*z) - 2.0 * anchor_shift,
            "forces": forces,
        }


def design_fixture(root, x_initial=-1.0, x_final=1.0):
    initial = Atoms("HeHeH", positions=[[0, 0, 3], [0, 0, -3], [x_initial, 0, 0]])
    final = initial.copy()
    final.positions[2, 0] = x_final
    write(root / "initial.xyz", initial)
    write(root / "final.xyz", final)
    design = root / "design.json"
    design.write_text(json.dumps({
        "schema_version": 1, "length_unit": "angstrom",
        "initial": "initial.xyz", "final": "final.xyz",
        "quantum": asdict(QuantumSettings()), "fixed_indices": [0, 1],
        "metadata": {"substrate_indices": [0, 2], "tool_indices": [1]},
    }))
    return design


def test_path_preserves_anchors_and_opposing_holding_loads(tmp_path, monkeypatch):
    design = design_fixture(tmp_path)
    monkeypatch.setattr(workflow, "PySCFCalculator", LoadedDoubleWell)
    result = workflow.run(design, tmp_path / "run", stage="path", fmax=0.01, steps=100)
    assert result["candidate_electronic_barrier_ev"] == pytest.approx(0.2)
    assert result["validation"]["numerical_stage_converged"] is True
    assert result["validation"]["design_validated"] is False
    structures = read(tmp_path / "run" / "path.extxyz", index=":")
    assert len(structures) == 7
    for atoms, summary in zip(structures, result["images"]):
        np.testing.assert_array_equal(atoms.positions[:2], [[0, 0, 3], [0, 0, -3]])
        # Equal/opposite support forces must remain visible even though their sum is zero.
        np.testing.assert_allclose(summary["anchor_force_sum_ev_per_angstrom"], [0, 0, 0])
        groups = summary["anchor_force_groups"]
        np.testing.assert_allclose(groups["substrate"]["external_holding_force_ev_per_angstrom"], [-2, 0, 0])
        np.testing.assert_allclose(groups["tool"]["external_holding_force_ev_per_angstrom"], [2, 0, 0])
        assert summary["free_force_max_ev_per_angstrom"] < 2.0
        assert summary["forces_ev_per_angstrom"][0] == [2.0, 0.0, 0.0]


def test_unrelaxed_endpoint_cannot_produce_path_barrier(tmp_path, monkeypatch):
    design = design_fixture(tmp_path, x_initial=-0.5)
    monkeypatch.setattr(workflow, "PySCFCalculator", LoadedDoubleWell)
    with pytest.raises(ValueError, match="endpoint did not converge"):
        workflow.run(design, tmp_path / "run", stage="path", fmax=1e-6, steps=1)
    saved = json.loads((tmp_path / "run" / "result.json").read_text())
    assert saved["status"] == "failed"
    assert saved["endpoints"][0]["geometry_converged"] is False
    assert "candidate_electronic_barrier_ev" not in saved
    assert saved["validation"]["numerical_stage_converged"] is False
    assert not (tmp_path / "run" / "band-pre.traj").exists()


@pytest.mark.parametrize("failure_type", [OSError, KeyboardInterrupt])
def test_input_snapshot_failure_gets_terminal_status(tmp_path, monkeypatch, failure_type):
    """A snapshot failure must not leave an apparently running orphan result."""
    design = design_fixture(tmp_path)
    monkeypatch.setattr(workflow, "PySCFCalculator", LoadedDoubleWell)
    original_write = workflow.write

    def interrupted_snapshot(path, *args, **kwargs):
        if str(path).endswith("input-final.extxyz"):
            raise failure_type("injected input snapshot failure")
        return original_write(path, *args, **kwargs)

    monkeypatch.setattr(workflow, "write", interrupted_snapshot)
    with pytest.raises(failure_type):
        workflow.run(design, tmp_path / "run", stage="path")
    saved = json.loads((tmp_path / "run" / "result.json").read_text())
    assert saved["status"] == ("interrupted" if failure_type is KeyboardInterrupt else "failed")
    assert saved["error"]["type"] == failure_type.__name__
    assert saved["validation"]["numerical_stage_converged"] is False
    assert "candidate_electronic_barrier_ev" not in saved


def test_calculator_setup_failure_gets_terminal_status(tmp_path, monkeypatch):
    design = design_fixture(tmp_path)

    def broken_setup(*args, **kwargs):
        raise RuntimeError("injected calculator initialization failure")

    monkeypatch.setattr(workflow, "PySCFCalculator", broken_setup)
    with pytest.raises(RuntimeError, match="initialization failure"):
        workflow.run(design, tmp_path / "run", stage="path")
    saved = json.loads((tmp_path / "run" / "result.json").read_text())
    assert saved["status"] == "failed"
    assert saved["error"]["type"] == "RuntimeError"
    assert saved["validation"]["numerical_stage_converged"] is False
    assert "candidate_electronic_barrier_ev" not in saved
