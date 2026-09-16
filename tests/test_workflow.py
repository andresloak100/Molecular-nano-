"""Workflow safety properties, using a known analytic surface only in tests."""
from dataclasses import asdict
import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write

from nanodesign import workflow
from nanodesign.design import validate_pair, load_design, endpoint_identity_ok, topology_screen
from nanodesign.quantum import QuantumSettings


class DoubleWell(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, settings=None):
        super().__init__()
        self.diagnostics = {"test_surface": True}

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x = atoms.positions[1, 0]
        forces = np.zeros((2, 3))
        forces[1, 0] = -0.8 * x * (x * x - 1)
        self.results = {"energy": 0.2 * (x * x - 1)**2, "forces": forces}


def write_test_design(tmp_path):
    initial = Atoms("HeH", positions=[[0, 0, 3], [-1, 0, 0]])
    final = initial.copy()
    final.positions[1, 0] = 1
    write(tmp_path / "initial.xyz", initial)
    write(tmp_path / "final.xyz", final)
    path = tmp_path / "design.json"
    path.write_text(json.dumps({"schema_version": 1, "length_unit": "angstrom", "initial": "initial.xyz", "final": "final.xyz",
                                "quantum": asdict(QuantumSettings()), "fixed_indices": [0]}))
    return path


def test_rejects_moving_anchors_and_overlap():
    a = Atoms("CH", positions=[[0, 0, 0], [1.1, 0, 0]])
    b = a.copy()
    b.positions[0, 1] = 0.1
    with pytest.raises(ValueError, match="anchors"):
        validate_pair(a, b, [0])
    b.positions[1] = b.positions[0]
    with pytest.raises(ValueError, match="closer"):
        validate_pair(a, b, [])


def test_neb_known_barrier_and_no_accuracy_claim(tmp_path, monkeypatch):
    path = write_test_design(tmp_path)
    monkeypatch.setattr(workflow, "PySCFCalculator", DoubleWell)
    result = workflow.run(path, tmp_path / "run", stage="path", fmax=0.01, steps=100)
    assert result["status"] == "completed"
    assert result["candidate_electronic_barrier_ev"] == pytest.approx(0.2, abs=1e-8)
    assert result["reaction_energy_ev"] == pytest.approx(0.0, abs=1e-8)
    assert result["validation"]["numerical_stage_converged"]
    assert result["validation"]["design_validated"] is False
    saved = json.loads((tmp_path / "run" / "result.json").read_text())
    assert len(saved["input_hashes"]["design_sha256"]) == 64
    assert (tmp_path / "run" / "path.extxyz").is_file()


def test_failure_saved_without_barrier(tmp_path, monkeypatch):
    path = write_test_design(tmp_path)

    class Broken(DoubleWell):
        def calculate(self, *args, **kwargs):
            raise RuntimeError("SCF did not converge")

    monkeypatch.setattr(workflow, "PySCFCalculator", Broken)
    with pytest.raises(RuntimeError, match="SCF"):
        workflow.run(path, tmp_path / "run")
    result = json.loads((tmp_path / "run" / "result.json").read_text())
    assert result["status"] == "failed"
    assert "candidate_electronic_barrier_ev" not in result
    assert not result["validation"]["numerical_stage_converged"]
    with pytest.raises(FileExistsError):
        workflow.run(path, tmp_path / "run")


def test_endpoint_collapse_cannot_pass(tmp_path, monkeypatch):
    path = write_test_design(tmp_path)
    write(tmp_path / "final.xyz", load_design(path)[1])
    monkeypatch.setattr(workflow, "PySCFCalculator", DoubleWell)
    with pytest.raises(ValueError, match="same structure"):
        workflow.run(path, tmp_path / "run", stage="path")


def test_transfer_identity_checks_both_bonds():
    a = Atoms("CHC", positions=[[0, 0, 0], [1.09, 0, 0], [3.6, 0, 0]])
    r = {"donor": 0, "hydrogen": 1, "acceptor": 2}
    assert endpoint_identity_ok(a, r, "initial")
    assert not endpoint_identity_ok(a, r, "final")
    a.positions[1, 0] = 2.54
    assert endpoint_identity_ok(a, r, "final")
    assert not endpoint_identity_ok(a, r, "initial")


def test_rejects_ambiguous_coordinate_units(tmp_path):
    path = write_test_design(tmp_path)
    data = json.loads(path.read_text())
    data["length_unit"] = "bohr"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="length_unit"):
        load_design(path)


def test_interruption_persisted(tmp_path, monkeypatch):
    path = write_test_design(tmp_path)

    class Interrupted(DoubleWell):
        def calculate(self, *args, **kwargs):
            raise KeyboardInterrupt()

    monkeypatch.setattr(workflow, "PySCFCalculator", Interrupted)
    with pytest.raises(KeyboardInterrupt):
        workflow.run(path, tmp_path / "run")
    result = json.loads((tmp_path / "run" / "result.json").read_text())
    assert result["status"] == "interrupted"
    assert not result["validation"]["numerical_stage_converged"]


@pytest.mark.parametrize("metadata", [{"tool_indices": [0.0]}, {"tool_indices": [True]}, {"tool_indices": [99]}, {"tool_indices": [0], "substrate_indices": [0]}, {"reactant_bonds": [[0, 0, 1]]}])
def test_bad_metadata_rejected_before_computation(tmp_path, metadata):
    path = write_test_design(tmp_path)
    data = json.loads(path.read_text())
    data["metadata"] = metadata
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        load_design(path)


def test_geometric_screen_detects_dissociated_bond():
    atoms = Atoms("CH", positions=[[0, 0, 0], [1.09, 0, 0]])
    assert topology_screen(atoms, [[0, 1, 1]])["preserved"]
    atoms.positions[1, 0] = 2.0
    assert not topology_screen(atoms, [[0, 1, 1]])["preserved"]
