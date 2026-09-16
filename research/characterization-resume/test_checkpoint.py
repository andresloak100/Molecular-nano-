"""Checkpoint lifecycle/provenance tests; every force evaluation is synthetic.

H1 test ownership: internal mode_contract helper. No quantum backend is run.
The independent E1 consumer owns saved-force algebra verification.
"""

from copy import deepcopy
from dataclasses import asdict, replace
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
import numpy as np
import pytest

from nanodesign.quantum import QuantumSettings


class MockQuantum(Calculator):
    """Known harmonic forces and independently tracked successful evaluations."""

    implemented_properties = ["forces"]

    def __init__(self, settings, *, fail_at=None, failure_type=KeyboardInterrupt, free_force=0.0):
        super().__init__()
        self.settings = settings
        self.reference = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
        self.matrix = np.diag([1.0, 2, 3, 2, 3, 4])
        self.offset = np.array([[7, -8, 9], [free_force, 0, 0]], dtype=float)
        self.fail_at = fail_at
        self.failure_type = failure_type
        self.evaluations = []
        self.completed_positions = []
        self.diagnostics = {"before": True}

    def forces_at(self, positions):
        delta = (positions - self.reference).ravel()
        return self.offset - (self.matrix @ delta).reshape((-1, 3))

    def calculate(self, atoms=None, properties=("forces",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations.append(atoms.positions.copy())
        number = len(self.evaluations)
        self.diagnostics = {
            "settings": asdict(self.settings),
            "scf_converged": False,
            "gradient_completed": False,
            "call_id": f"mock-quantum-{id(self)}-{number}",
        }
        if number == self.fail_at:
            self.results = {}
            raise self.failure_type("synthetic force evaluation interrupted")
        self.diagnostics.update(scf_converged=True, gradient_completed=True)
        self.results = {"forces": self.forces_at(atoms.positions)}
        self.completed_positions.append(atoms.positions.copy())


@pytest.fixture(scope="module")
def prototype_module():
    name = "characterization_resume_prototype_tests"
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name("prototype.py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    yield module
    sys.modules.pop(name, None)


@pytest.fixture
def prototype(prototype_module, monkeypatch):
    monkeypatch.setattr(prototype_module, "PySCFCalculator", MockQuantum)
    return prototype_module


@pytest.fixture
def settings():
    return QuantumSettings(spin=0, dispersion=None, threads=1)


@pytest.fixture
def input_context():
    return {
        "input_hashes": {"design_sha256": "a" * 64, "structure_sha256": "b" * 64},
        "structure_image": 0,
    }


def make_atoms(settings, **calculator_options):
    atoms = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]], masses=[1.5, 2.5])
    atoms.set_constraint(FixAtoms(indices=[0]))
    atoms.calc = MockQuantum(settings, **calculator_options)
    return atoms


def read_checkpoint(directory):
    return json.loads((directory / "force_checkpoint.json").read_text())


def rewrite_checkpoint(directory, envelope, *, rehash=True):
    if rehash:
        encoded = json.dumps(envelope["payload"], sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        envelope["sha256"] = hashlib.sha256(encoded).hexdigest()
    (directory / "force_checkpoint.json").write_text(json.dumps(envelope, allow_nan=False))


def directory_bytes(directory):
    return {
        str(path.relative_to(directory)): path.read_bytes()
        for path in directory.rglob("*") if path.is_file()
    }


def snapshot_calculator(atoms):
    calculator = atoms.calc
    return {
        "positions": atoms.positions.copy(),
        "calculator": calculator,
        "calculator_positions": None if calculator.atoms is None else calculator.atoms.positions.copy(),
        "results": deepcopy(calculator.results),
        "diagnostics": deepcopy(calculator.diagnostics),
    }


def assert_calculator_restored(atoms, saved):
    calculator = saved["calculator"]
    assert atoms.calc is calculator
    np.testing.assert_array_equal(atoms.positions, saved["positions"])
    assert calculator.diagnostics == saved["diagnostics"]
    assert calculator.results.keys() == saved["results"].keys()
    for key, value in saved["results"].items():
        np.testing.assert_array_equal(calculator.results[key], value)
    if saved["calculator_positions"] is None:
        assert calculator.atoms is None
    else:
        np.testing.assert_array_equal(calculator.atoms.positions, saved["calculator_positions"])


def completed_source(prototype, tmp_path, settings, input_context):
    atoms = make_atoms(settings)
    output = tmp_path / "source"
    result = prototype.run_characterization_checkpoint(
        atoms, output, settings=settings, input_context=input_context,
    )
    assert result["status"] == "completed"
    return output


def test_interrupted_attempt_resumes_only_missing_records_without_changing_source(
    prototype, tmp_path, settings, input_context,
):
    original = make_atoms(settings, fail_at=3)
    original_state = snapshot_calculator(original)
    source = tmp_path / "interrupted"
    with pytest.raises(KeyboardInterrupt, match="synthetic"):
        prototype.run_characterization_checkpoint(
            original, source, settings=settings, input_context=input_context,
        )

    checkpoint = read_checkpoint(source)["payload"]
    assert checkpoint["status"] == "interrupted"
    assert len(original.calc.evaluations) == 3
    assert len(checkpoint["records"]) == len(original.calc.completed_positions) == 2
    assert all(record["status"] == "completed" for record in checkpoint["records"])
    assert_calculator_restored(original, original_state)
    source_before = directory_bytes(source)

    resumed = make_atoms(settings)
    resumed_state = snapshot_calculator(resumed)
    output = tmp_path / "resumed"
    result = prototype.run_characterization_checkpoint(
        resumed, output, settings=settings, input_context=input_context, resume_from=source,
    )

    assert result["status"] == "completed"
    assert result["checkpoint"]["reused_force_records"] == 2
    assert result["checkpoint"]["new_force_evaluations"] == len(resumed.calc.evaluations) == 5
    assert result["checkpoint"]["completed_force_records"] == 7
    assert result["quantum_settings"] == asdict(settings)
    assert all(result["input_hashes"][key] == value for key, value in input_context["input_hashes"].items())
    assert result["input_hashes"]["input_snapshot_sha256"] == hashlib.sha256((output / "input.extxyz").read_bytes()).hexdigest()
    assert result["stationary"]["transition_state_verified"] is False
    # Missing-work acquisition starts with minus, not baseline or saved plus.
    expected_first_missing = resumed.positions.copy()
    expected_first_missing[1, 0] -= 0.005
    np.testing.assert_array_equal(resumed.calc.evaluations[0], expected_first_missing)
    new_checkpoint = read_checkpoint(output)["payload"]
    assert new_checkpoint["attempt_id"] != checkpoint["attempt_id"]
    assert new_checkpoint["records"][:2] == checkpoint["records"]
    assert all(record["origin_attempt_id"] == checkpoint["attempt_id"] for record in new_checkpoint["records"][:2])
    assert all(record["origin_attempt_id"] == new_checkpoint["attempt_id"] for record in new_checkpoint["records"][2:])
    assert directory_bytes(source) == source_before
    assert_calculator_restored(resumed, resumed_state)


@pytest.mark.parametrize("change", ["geometry", "anchors", "masses", "settings", "guess", "step", "input_context"])
def test_incompatible_resume_rejected_before_calls_or_output(
    prototype, tmp_path, settings, input_context, change,
):
    source = completed_source(prototype, tmp_path, settings, input_context)
    source_before = directory_bytes(source)
    context = deepcopy(input_context)
    options = {}
    if change == "settings":
        settings = replace(settings, basis="sto-3g")
    if change == "guess":
        settings = replace(settings, scf_initial_guess="atom")
    atoms = make_atoms(settings)
    if change == "geometry":
        atoms.positions[1, 1] += 0.01
    elif change == "anchors":
        atoms.set_constraint([])
    elif change == "masses":
        atoms.set_masses([1.5, 3.0])
    elif change == "step":
        options["step"] = 0.01
    elif change == "input_context":
        context["input_hashes"]["structure_sha256"] = "c" * 64
    output = tmp_path / "incompatible"

    with pytest.raises(ValueError):
        prototype.run_characterization_checkpoint(
            atoms, output, settings=settings, input_context=context, resume_from=source, **options,
        )

    assert atoms.calc.evaluations == []
    assert not output.exists()
    assert directory_bytes(source) == source_before


@pytest.mark.parametrize("corruption", [
    "digest", "unknown_request", "force_shape", "nonfinite_force",
    "scf_failed", "gradient_failed", "record_failed", "duplicate_request",
])
def test_corrupt_or_semantically_invalid_resume_rejected_before_work(
    prototype, tmp_path, settings, input_context, corruption,
):
    source = completed_source(prototype, tmp_path, settings, input_context)
    envelope = read_checkpoint(source)
    records = envelope["payload"]["records"]
    if corruption == "digest":
        envelope["sha256"] = "0" * 64
    elif corruption == "unknown_request":
        records[1]["request"] = {"kind": "not-a-request"}
    elif corruption == "force_shape":
        records[1]["forces_ev_per_angstrom"] = [[0, 0, 0]]
    elif corruption == "nonfinite_force":
        records[1]["forces_ev_per_angstrom"][1][0] = "nan"
    elif corruption == "scf_failed":
        records[1]["quantum_diagnostics"]["scf_converged"] = False
    elif corruption == "gradient_failed":
        records[1]["quantum_diagnostics"]["gradient_completed"] = False
    elif corruption == "record_failed":
        records[1]["status"] = "failed"
    elif corruption == "duplicate_request":
        records[1] = deepcopy(records[0])
    rewrite_checkpoint(source, envelope, rehash=corruption != "digest")
    source_before = directory_bytes(source)
    atoms = make_atoms(settings)
    output = tmp_path / "rejected"

    with pytest.raises(ValueError):
        prototype.run_characterization_checkpoint(
            atoms, output, settings=settings, input_context=input_context, resume_from=source,
        )

    assert atoms.calc.evaluations == []
    assert not output.exists()
    assert directory_bytes(source) == source_before


def test_nonstationary_baseline_is_saved_with_raw_anchor_forces_before_guard_failure(
    prototype, tmp_path, settings, input_context,
):
    atoms = make_atoms(settings, free_force=0.2)
    output = tmp_path / "nonstationary"
    saved = snapshot_calculator(atoms)

    with pytest.raises(ValueError, match="stationary|force"):
        prototype.run_characterization_checkpoint(
            atoms, output, settings=settings, input_context=input_context,
        )

    payload = read_checkpoint(output)["payload"]
    assert payload["status"] == "failed"
    assert len(payload["records"]) == len(atoms.calc.evaluations) == 1
    record = payload["records"][0]
    assert record["status"] == "completed"
    np.testing.assert_array_equal(record["forces_ev_per_angstrom"], atoms.calc.forces_at(atoms.positions))
    assert np.linalg.norm(record["forces_ev_per_angstrom"][0]) > 1
    assert record["quantum_diagnostics"]["scf_converged"] is True
    assert record["quantum_diagnostics"]["gradient_completed"] is True
    assert_calculator_restored(atoms, saved)


@pytest.mark.parametrize("failing", [False, True])
@pytest.mark.parametrize("existing_cache", [False, True])
def test_calculator_geometry_cache_and_diagnostics_restore_on_success_and_failure(
    prototype, tmp_path, settings, input_context, failing, existing_cache,
):
    atoms = make_atoms(settings, failure_type=RuntimeError)
    calculator = atoms.calc
    if existing_cache:
        # Restore a different preexisting cached geometry, not just the input.
        cached_geometry = atoms.copy()
        cached_geometry.positions[1, 1] += 0.1
        cached_geometry.calc = calculator
        cached_geometry.get_forces(apply_constraint=False)
    if failing:
        calculator.fail_at = len(calculator.evaluations) + 3
    saved = snapshot_calculator(atoms)
    output = tmp_path / "calculation"

    if failing:
        with pytest.raises(RuntimeError, match="synthetic"):
            prototype.run_characterization_checkpoint(
                atoms, output, settings=settings, input_context=input_context,
            )
        assert read_checkpoint(output)["payload"]["status"] == "failed"
    else:
        result = prototype.run_characterization_checkpoint(
            atoms, output, settings=settings, input_context=input_context,
        )
        assert result["status"] == "completed"

    assert_calculator_restored(atoms, saved)


def test_snapshot_hash_is_separate_from_multiframe_source_and_full_precision_geometry(
    prototype, tmp_path, settings,
):
    from ase.io import read

    source = tmp_path / 'frames.xyz'
    source.write_text('2\nframe zero\nH 0 0 0\nH 1 0 0\n'
                      '2\nselected frame\nH 0 0 0\nH 1.123456789012345 0 0\n')
    atoms = read(source, index=1)
    atoms.set_constraint(FixAtoms(indices=[0]))
    atoms.calc = MockQuantum(settings)
    atoms.calc.reference = atoms.positions.copy()
    exact = atoms.positions.copy()
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    context = {'input_hashes': {'structure_sha256': source_hash}, 'structure_image': 1}
    original_context = deepcopy(context)
    output = tmp_path / 'first'
    result = prototype.run_characterization_checkpoint(
        atoms, output, settings=settings, input_context=context,
    )
    snapshot_hash = hashlib.sha256((output / 'input.extxyz').read_bytes()).hexdigest()
    assert snapshot_hash != source_hash
    assert result['input_hashes'] == {
        'structure_sha256': source_hash, 'input_snapshot_sha256': snapshot_hash,
    }
    assert result['units'] == {'length': 'angstrom', 'energy': 'eV', 'force': 'eV/angstrom', 'frequency': 'cm^-1'}
    assert context == original_context
    assert result['input_context'] == original_context
    assert result['checkpoint_manifest']['input_context'] == original_context
    np.testing.assert_array_equal(result['checkpoint_manifest']['geometry']['positions_angstrom'], exact)
    np.testing.assert_array_equal(result['stationary']['finite_difference_evidence']['reference_positions_angstrom'], exact)
    assert not np.array_equal(read(output / 'input.extxyz').positions, exact)
    assert json.loads((output / 'result.json').read_text())['input_hashes'] == result['input_hashes']

    original_bytes = directory_bytes(output)
    calls = len(atoms.calc.evaluations)
    resumed = prototype.run_characterization_checkpoint(
        atoms, tmp_path / 'resumed', settings=settings, input_context=context, resume_from=output,
    )
    assert len(atoms.calc.evaluations) == calls
    assert resumed['checkpoint']['new_force_evaluations'] == 0
    assert resumed['input_hashes']['structure_sha256'] == source_hash
    assert resumed['input_hashes']['input_snapshot_sha256'] == hashlib.sha256((tmp_path / 'resumed/input.extxyz').read_bytes()).hexdigest()
    assert directory_bytes(output) == original_bytes


def test_snapshot_hash_does_not_mutate_caller_snapshot_provenance(
    prototype, tmp_path, settings, input_context,
):
    context = deepcopy(input_context)
    context['input_hashes']['input_snapshot_sha256'] = 'c' * 64
    original = deepcopy(context)
    output = tmp_path / 'snapshot'
    result = prototype.run_characterization_checkpoint(
        make_atoms(settings), output, settings=settings, input_context=context,
    )
    assert context == original
    assert result['input_context'] == original
    assert result['checkpoint_manifest']['input_context'] == original
    assert result['input_hashes']['structure_sha256'] == context['input_hashes']['structure_sha256']
    assert result['input_hashes']['input_snapshot_sha256'] == hashlib.sha256((output / 'input.extxyz').read_bytes()).hexdigest()
    assert result['input_hashes']['input_snapshot_sha256'] != context['input_hashes']['input_snapshot_sha256']
