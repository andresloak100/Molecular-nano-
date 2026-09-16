"""F2 actual H1 producer-to-F1 adapter checks using an explicit analytic mock."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
import pytest

from nanodesign.quantum import QuantumSettings


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


oracle = load("f2_bridge_oracle", HERE / "analytic_fixtures.py")
h1 = load("f2_bridge_h1", HERE.parent / "characterization-resume/prototype.py")
f1 = load("f2_bridge_f1", HERE.parent / "force-energy-consistency/consistency.py")
SETTINGS = QuantumSettings(spin=0, basis="sto-3g", dispersion=None, threads=1)
CONVERSION = 27.211386245988  # Explicit conversion recorded by this synthetic producer.


class AnalyticBackend(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, label, *, missing_total=False, fail_at=None):
        super().__init__()
        self.settings = SETTINGS
        self.label = label
        self.missing_total = missing_total
        self.fail_at = fail_at
        self.calls = []
        self.diagnostics = {"synthetic_fixture": True}

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {}
        self.calls.append(atoms.positions.tolist())
        self.diagnostics = {
            "call_id": f"{self.label}-{len(self.calls)}", "settings": SETTINGS.to_dict(),
            "synthetic_fixture": True, "scf_converged": False, "gradient_completed": False,
            "hartree_eV": CONVERSION,
        }
        if len(self.calls) == self.fail_at:
            raise RuntimeError("synthetic F2 interrupted acquisition")
        values = oracle.polynomial_record(atoms.positions.tolist(), cubic=8)
        # SCF-only energy is deliberately unsuitable as the total-energy record.
        self.diagnostics["dft_energy_hartree"] = values["energy_ev"] / CONVERSION + 1.5
        if not self.missing_total:
            self.diagnostics["total_energy_hartree"] = values["energy_ev"] / CONVERSION
        self.diagnostics.update(scf_converged=True, gradient_completed=True)
        import numpy as np
        self.results = {"energy": values["energy_ev"], "forces": np.array(values["forces_ev_per_angstrom"])}


def tree_bytes(directory):
    return {path.relative_to(directory).as_posix(): path.read_bytes()
            for path in directory.rglob("*") if path.is_file()}


def produce(output, monkeypatch, *, step=.125, missing_total=False, fail_at=None):
    monkeypatch.setattr(h1, "PySCFCalculator", AnalyticBackend)
    atoms = Atoms("H2", positions=deepcopy(oracle.REFERENCE))
    atoms.set_constraint(FixAtoms(indices=[1]))
    atoms.calc = backend = AnalyticBackend(output.name, missing_total=missing_total, fail_at=fail_at)
    def run():
        h1.run_characterization_checkpoint(atoms, output, settings=SETTINGS,
            input_context={"synthetic_fixture": "F2 polynomial producer; not chemistry"}, step=step, fmax=2)
    if fail_at:
        with pytest.raises(RuntimeError, match="synthetic F2 interrupted"):
            run()
    else:
        run()
    return backend


def test_real_checkpoint_schema_converts_total_energy_and_multiple_steps_without_solver_work(tmp_path, monkeypatch):
    paths = [tmp_path / "coarse", tmp_path / "fine"]
    backends = [produce(path, monkeypatch, step=step) for path, step in zip(paths, [.125, .0625])]
    source_bytes = [tree_bytes(path) for path in paths]
    before_calls = [len(backend.calls) for backend in backends]
    document = f1.from_h1_checkpoints(paths)
    report = f1.analyze_stencils(document)
    assert report["status"] == "comparison_produced"
    assert len(report["comparisons"]) == 6
    for index, row in enumerate(report["comparisons"]):
        step = .125 if index < 3 else .0625
        expected = oracle.exact_three_point_force(a=step, b=step, cubic=8) if row["axis"] == 0 else 0
        assert row["finite_difference_force_ev_per_angstrom"] == pytest.approx(expected, abs=1e-10)
    assert document["stencils"][0]["reference"]["energy_ev"] == pytest.approx(-100)
    assert document["stencils"][0]["reference"]["forces_ev_per_angstrom"][1] == [-3.0, 0.0, 0.0]
    assert [len(backend.calls) for backend in backends] == before_calls == [7, 7]
    assert [tree_bytes(path) for path in paths] == source_bytes
    assert report["electronic_branch_continuity_verified"] is False
    assert report["scientific_model_validated"] is False


def test_scf_only_record_does_not_substitute_for_missing_total_energy(tmp_path, monkeypatch):
    path = tmp_path / "no-total"
    backend = produce(path, monkeypatch, missing_total=True)
    before = tree_bytes(path)
    document = f1.from_h1_checkpoints([path])
    assert document["stencils"][0]["reference"]["energy_ev"] is None
    report = f1.analyze_stencils(document)
    assert report["status"] == "unavailable"
    assert report["comparisons"] == []
    assert len(backend.calls) == 7
    assert tree_bytes(path) == before


def test_interrupted_checkpoint_preserves_failure_reason_and_missing_slots(tmp_path, monkeypatch):
    path = tmp_path / "failed"
    backend = produce(path, monkeypatch, fail_at=3)
    before = tree_bytes(path)
    document = f1.from_h1_checkpoints([path])
    report = f1.analyze_stencils(document)
    assert "synthetic F2 interrupted acquisition" in json.dumps(document["source_checkpoints"])
    assert report["status"] == "unavailable"
    assert len(report["unavailable_stencils"]) == 3
    assert report["comparisons"] == []
    assert len(backend.calls) == 3
    assert tree_bytes(path) == before
