"""Real H1-to-N1 artifact integration using analytical forces, never quantum data.

Only the acquisition calculator is replaced. H1 writes its own result, snapshot,
and force evidence; these tests do not synthesize or repair its output schema.
"""
from copy import deepcopy
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
from ase.io import read
import numpy as np
import pytest

from nanodesign.quantum import QuantumSettings


HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from compare_modes import ComparisonError, compare_bound, load_saved_result


class AnalyticalAcquisition(Calculator):
    """Synthetic harmonic test backend; successful flags exercise H1's contract."""

    implemented_properties = ["forces"]

    def __init__(self, settings, reference):
        super().__init__()
        self.settings = settings
        self.reference = np.array(reference, dtype=float, copy=True)
        self.diagonal = np.arange(1., self.reference.size + 1.)
        self.evaluations = 0
        self.diagnostics = {"fixture": "analytical forces, no quantum calculation"}

    def calculate(self, atoms=None, properties=("forces",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations += 1
        delta = (atoms.positions - self.reference).ravel()
        self.results = {"forces": -(self.diagonal * delta).reshape(self.reference.shape)}
        self.diagnostics = {
            "settings": asdict(self.settings),
            "scf_initial_guess": self.settings.scf_initial_guess,
            "scf_converged": True,
            "gradient_completed": True,
            "call_id": f"synthetic-h1-n1-{id(self)}-{self.evaluations}",
            "versions": {"test_backend": "analytical harmonic fixture, not PySCF"},
        }


@pytest.fixture(scope="module")
def h1_module():
    name = "n1_actual_h1_producer_test"
    path = ROOT / "research/characterization-resume/prototype.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    yield module
    sys.modules.pop(name, None)


@pytest.fixture
def produced_pair(tmp_path, monkeypatch, h1_module):
    monkeypatch.setattr(h1_module, "PySCFCalculator", AnalyticalAcquisition)
    source = tmp_path / "original-precise-source.xyz"
    source.write_bytes(
        b"2\nSynthetic analytical fixture, no molecular prediction\n"
        b"H 0.123456789123 0.023456789123 -0.034567891234\n"
        b"H 0.864197532468 0.023456789123 -0.034567891234\n"
    )
    original_bytes = source.read_bytes()
    atoms = read(source, format="xyz")
    atoms.set_masses([1.5, 2.5])
    atoms.set_constraint(FixAtoms(indices=[0]))
    positions = atoms.positions.copy()
    settings = QuantumSettings(spin=0, dispersion=None, threads=1)
    design = {"name": "synthetic H1/N1 integration fixture",
              "quantum": asdict(settings), "fixed_indices": [0]}
    design_path = tmp_path / "source-design.json"
    design_path.write_text(json.dumps(design), encoding="utf-8")
    input_context = {
        "design": design,
        "input_hashes": {
            "structure_sha256": hashlib.sha256(original_bytes).hexdigest(),
            "design_sha256": hashlib.sha256(design_path.read_bytes()).hexdigest(),
        },
        "source_structure_path": str(source),
        "synthetic_fixture": True,
    }
    context_before = deepcopy(input_context)
    design_before = design_path.read_bytes()
    outputs = []
    for label, step in (("coarse", .005), ("fine", .0025)):
        atoms.calc = AnalyticalAcquisition(settings, positions)
        output = tmp_path / label
        result = h1_module.run_characterization_checkpoint(
            atoms, output, settings=settings, input_context=input_context, step=step,
        )
        assert result["status"] == "completed"
        assert atoms.calc.evaluations == 7  # baseline plus +/- on three free coordinates
        assert json.loads((output / "result.json").read_text()) == result
        outputs.append((output, result))
    assert source.read_bytes() == original_bytes
    assert design_path.read_bytes() == design_before
    assert input_context == context_before
    np.testing.assert_array_equal(atoms.positions, positions)
    return {"outputs": outputs, "positions": positions,
            "source": source, "source_bytes": original_bytes,
            "source_hash": input_context["input_hashes"]["structure_sha256"]}


def test_actual_h1_records_compare_with_distinct_source_snapshot_and_precise_geometry(produced_pair):
    loaded = []
    immutable = {}
    for output, result in produced_pair["outputs"]:
        snapshot = output / "input.extxyz"
        snapshot_hash = hashlib.sha256(snapshot.read_bytes()).hexdigest()
        assert result["input_hashes"]["structure_sha256"] == produced_pair["source_hash"]
        assert result["input_hashes"]["input_snapshot_sha256"] == snapshot_hash
        assert snapshot_hash != produced_pair["source_hash"]
        assert result["units"]["length"] == "angstrom"
        assert result["units"]["energy"] == "eV"
        assert result["units"]["frequency"] == "cm^-1"
        immutable.update({path: path.read_bytes() for path in output.rglob("*") if path.is_file()})
        bound = load_saved_result(output / "result.json")
        assert bound["source"]["snapshot_hash_binding"] == "input_snapshot_sha256"
        assert bound["source"]["source_structure_sha256"] == produced_pair["source_hash"]
        assert bound["geometry_binding"]["snapshot_relation"] == "reference_rounded_to_eight_decimal_angstrom"
        np.testing.assert_array_equal(bound["positions_angstrom"], produced_pair["positions"])
        assert not np.array_equal(read(snapshot).positions, produced_pair["positions"])
        loaded.append(bound)
    comparison = compare_bound(*loaded)
    assert comparison["status"] == "compared", comparison["findings"]
    assert comparison["steps_angstrom"] == [.005, .0025]
    assert comparison["distinct_steps"]
    assert comparison["resolved_negative_mode_counts"] == [0, 0]
    assert comparison["spectral_rank_frequency_changes_cm1"] == pytest.approx([0.] * 3, abs=1e-8)
    assert not comparison["step_convergence_certified"]
    assert not comparison["electronic_state_verified"]
    assert all(path.read_bytes() == content for path, content in immutable.items())
    assert produced_pair["source"].read_bytes() == produced_pair["source_bytes"]


def test_actual_h1_snapshot_tampering_is_rejected_without_rewriting_evidence(produced_pair):
    output, result = produced_pair["outputs"][0]
    snapshot = output / "input.extxyz"
    result_path = output / "result.json"
    saved_result = result_path.read_bytes()
    assert "input_snapshot_sha256" in result["input_hashes"]
    snapshot.write_bytes(snapshot.read_bytes() + b"\n")
    tampered = snapshot.read_bytes()
    with pytest.raises(ComparisonError, match="snapshot hash.*mismatched"):
        load_saved_result(result_path)
    assert snapshot.read_bytes() == tampered
    assert result_path.read_bytes() == saved_result
    assert produced_pair["source"].read_bytes() == produced_pair["source_bytes"]
