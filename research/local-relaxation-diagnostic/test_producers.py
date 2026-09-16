"""Exercise actual producers with analytical forces; no quantum calculation."""
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
from ase.io import write
import numpy as np
import pytest

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from diagnose import diagnose_saved_result, main
from nanodesign import workflow
from nanodesign.quantum import QuantumSettings


OPTIONS = dict(curvature_floor_ev_per_angstrom2=.001, max_condition_number=1e5,
               max_atom_shift_angstrom=.1, expected_index=0,
               max_relative_hessian_asymmetry=.001)
REFERENCE = np.array([[.123456789123, -.023456789123, .034567891234],
                      [.864197532468, .063456789123, -.074567891234]])
HESSIAN = np.array([[2., .5, 0.], [.5, 3., .25], [0., .25, 4.]])
GRADIENT = np.array([.006, -.007, .008])


class AnalyticalCalculator(Calculator):
    """Raw anchor force is deliberately large; only the second atom is free."""
    implemented_properties = ["energy", "forces"]
    total_calls = 0

    def __init__(self, settings, event_log=None):
        super().__init__()
        self.settings = settings
        self.diagnostics = {}

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        type(self).total_calls += 1
        delta = atoms.positions[1] - REFERENCE[1]
        free_gradient = GRADIENT + HESSIAN @ delta
        self.results = {"energy": float(GRADIENT @ delta + .5 * delta @ HESSIAN @ delta),
                        "forces": np.array([[4., -5., 6.], -free_gradient])}
        self.diagnostics = {
            "settings": asdict(self.settings), "scf_initial_guess": self.settings.scf_initial_guess,
            "scf_converged": True, "gradient_completed": True,
            "call_id": f"analytical-l1-{type(self).total_calls}",
            "versions": {"fixture": "analytical Cartesian quadratic; not PySCF"},
        }


@pytest.fixture(params=["workflow", "h1"])
def produced(tmp_path, monkeypatch, request):
    atoms = Atoms("HeH", positions=REFERENCE)
    atoms.set_masses([5., 2.])
    atoms.set_constraint(FixAtoms(indices=[0]))
    settings = QuantumSettings(spin=1, dispersion=None, threads=1)
    write(tmp_path / "initial.xyz", atoms)
    write(tmp_path / "final.xyz", atoms)
    design = {"schema_version": 1, "length_unit": "angstrom", "initial": "initial.xyz",
              "final": "final.xyz", "fixed_indices": [0], "quantum": asdict(settings)}
    design_path = tmp_path / "design.json"
    design_path.write_text(json.dumps(design))
    # Source bytes are binary and contain a selected frame, unlike the snapshot.
    other = atoms.copy()
    other.positions[1, 0] += .1
    source = tmp_path / "frames.traj"
    write(source, [other, atoms])
    before_source = source.read_bytes()
    output = tmp_path / "output"
    if request.param == "workflow":
        monkeypatch.setattr(workflow, "PySCFCalculator", AnalyticalCalculator)
        result = workflow.run_characterization(design_path, source, output, image=1, fmax=.03)
    else:
        name = "l1_actual_h1_producer"
        spec = importlib.util.spec_from_file_location(name, ROOT / "research/characterization-resume/prototype.py")
        module = importlib.util.module_from_spec(spec)
        monkeypatch.setitem(sys.modules, name, module)
        spec.loader.exec_module(module)
        monkeypatch.setattr(module, "PySCFCalculator", AnalyticalCalculator)
        atoms.calc = AnalyticalCalculator(settings)
        result = module.run_characterization_checkpoint(atoms, output, settings=settings,
            input_context={"design": design, "input_hashes": {
                "structure_sha256": hashlib.sha256(before_source).hexdigest(),
                "design_sha256": hashlib.sha256(design_path.read_bytes()).hexdigest()},
                "fixture": "synthetic analytical forces"}, fmax=.03)
    assert result["status"] == "completed"
    return output / "result.json", source, before_source


def test_actual_producer_roundtrip_no_calls_or_changes(produced):
    result_path, source, original = produced
    files = {p: p.read_bytes() for p in result_path.parent.iterdir() if p.is_file()}
    calls = AnalyticalCalculator.total_calls
    report = diagnose_saved_result(result_path, **OPTIONS)
    assert report["status"] == "estimated", report
    assert report["saved_force_reconstruction"]["status"] == "passed"
    expected = -np.linalg.solve(HESSIAN, GRADIENT)
    estimate = report["estimate"]
    np.testing.assert_allclose(estimate["displacement_free_atoms_angstrom"], [expected], atol=1e-12)
    np.testing.assert_array_equal(estimate["displacement_all_atoms_angstrom"][0], [0., 0., 0.])
    assert estimate["signed_quadratic_energy_change_ev"] == pytest.approx(.5 * GRADIENT @ expected)
    np.testing.assert_array_equal(report["reference_positions_angstrom"], REFERENCE)
    assert report["geometry_binding"]["snapshot_relation"] == "reference_rounded_to_eight_decimal_angstrom"
    assert report["source"]["source_structure_sha256"] != report["source"]["input_snapshot_sha256"]
    assert not report["energy_error_bound_established"]
    assert not report["actual_relaxation_performed"]
    assert AnalyticalCalculator.total_calls == calls
    assert source.read_bytes() == original
    assert all(p.read_bytes() == data for p, data in files.items())


def test_snapshot_tamper_is_refused(produced):
    result_path, _, _ = produced
    snapshot = result_path.with_name("input.extxyz")
    snapshot.write_bytes(snapshot.read_bytes() + b"\n")
    report = diagnose_saved_result(result_path, **OPTIONS)
    assert report["status"] == "refused"
    assert report["estimate"] is None
    assert "snapshot hash" in report["findings"][0]["message"]


def test_force_tamper_is_refused(produced):
    result_path, _, _ = produced
    data = json.loads(result_path.read_bytes())
    data["stationary"]["finite_difference_evidence"]["baseline_forces_ev_per_angstrom"][1][0] += .001
    result_path.write_text(json.dumps(data))
    report = diagnose_saved_result(result_path, **OPTIONS)
    assert report["status"] == "refused"
    assert report["estimate"] is None
    assert report["saved_force_reconstruction"]["status"] == "failed"


def test_missing_historical_forces_remain_unavailable(produced):
    result_path, _, _ = produced
    data = json.loads(result_path.read_bytes())
    del data["stationary"]["finite_difference_evidence"]
    result_path.write_text(json.dumps(data))
    report = diagnose_saved_result(result_path, **OPTIONS)
    assert report["status"] == "unavailable"
    assert report["estimate"] is None


def test_cli_does_not_overwrite_source(produced, capsys):
    result_path, _, _ = produced
    before = result_path.read_bytes()
    args = [str(result_path), "--curvature-floor", ".001", "--max-condition", "100000",
            "--max-atom-shift", ".1", "--max-relative-asymmetry", ".001", "--expected-index", "0"]
    assert main(args) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "estimated"
    assert main(args + ["--output", str(result_path)]) == 2
    assert result_path.read_bytes() == before
