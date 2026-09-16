"""Independent capture/cache edge cases using synthetic solver boundaries only.

The first eight cases replace PySCF and serialization to isolate the lifecycle.
Two additional roundtrips use the real E3 serializer with a synthetic solver.
No SCF, gradient solver, or AO integrals are evaluated.
"""

from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.constraints import FixAtoms
from ase.data import atomic_numbers
from ase.units import Hartree

from nanodesign.quantum import PySCFCalculator, QuantumCalculationError, QuantumSettings
from nanodesign.stationary import characterize_stationary_point


@pytest.fixture
def synthetic_boundaries(monkeypatch):
    state = SimpleNamespace(
        solver_calls=[], captures=[], gradient_calls=0, fail_gradient_at=None,
        primary_error=None, fail_restore=False, threads=3, thread_requests=[],
    )

    class Molecule:
        def __init__(self, **values):
            self.symbols, positions = zip(*values["atom"])
            self.positions = np.asarray(positions, dtype=float)
            self.natm = len(self.symbols)
            self.charge = values["charge"]
            self.spin = values["spin"]
            self.nelectron = sum(atomic_numbers[symbol] for symbol in self.symbols) - values["charge"]
            self.nelec = ((self.nelectron + self.spin) // 2, (self.nelectron - self.spin) // 2)
            self.cart = False
            self.nbas = self.natm
            self.ecp = self.pseudo = False
            self._basis = {symbol: [[0, [1.0, 1.0]]] for symbol in self.symbols}

        def nao_nr(self):
            return 2

        def atom_symbol(self, index):
            return self.symbols[index]

        def atom_charges(self):
            return np.array([atomic_numbers[symbol] for symbol in self.symbols])

        def atom_coords(self, *, unit):
            assert unit == "Bohr"
            return self.positions / 0.52917721092

        def ao_labels(self, *, fmt):
            assert fmt is False
            return [(index, symbol, "1s", "") for index, symbol in enumerate(self.symbols)]

        def bas_atom(self, index):
            return index

        def bas_angular(self, index):
            return 0

        def bas_kappa(self, index):
            return 0

        def bas_exp(self, index):
            return np.array([1.0])

        def bas_ctr_coeff(self, index):
            return np.array([[1.0]])

        def _libcint_ctr_coeff(self, index):
            return np.array([[1.0]])

    class MeanField:
        def __init__(self, molecule):
            self.mol = molecule
            self.grids = SimpleNamespace(level=None)
            self.converged = True
            self.cycles = 1
            self.density_fitted = False
            self.mo_coeff = np.eye(2) if molecule.spin == 0 else np.array([np.eye(2), np.eye(2)])
            self.mo_occ = (np.array([2.0, 0.0]) if molecule.spin == 0
                           else np.array([[1.0, 0.0], [0.0, 0.0]]))

        def istype(self, name):
            return name == ("RKS" if self.mol.spin == 0 else "UKS")

        def get_ovlp(self):
            # Deliberately synthetic AO metric; no integral evaluation occurs.
            return np.eye(2)

        def density_fit(self):
            self.density_fitted = True
            return self

        def kernel(self):
            state.solver_calls.append(self)
            # E=-1 Ha + (1 eV/A^2)/2 * sum(x_i^2), a synthetic quadratic.
            self.e_tot = -1.0 + 0.5 * float(np.sum(self.mol.positions ** 2)) / Hartree
            return self.e_tot

        def spin_square(self):
            spin = self.mol.spin / 2
            return spin * (spin + 1), self.mol.spin + 1

        def nuc_grad_method(self):
            def gradient():
                state.gradient_calls += 1
                if state.gradient_calls == state.fail_gradient_at:
                    raise state.primary_error
                return self.mol.positions * 0.52917721092 / Hartree
            return SimpleNamespace(kernel=gradient)

    def num_threads(value=None):
        if value is not None:
            state.thread_requests.append(value)
            if value == 3 and state.fail_restore:
                raise OSError("Synthetic thread restoration failed")
            state.threads = value
        return state.threads

    modules = {name: ModuleType(name) for name in (
        "pyscf", "pyscf.dft", "pyscf.gto", "pyscf.gto.mole", "pyscf.lib", "pyscf.lib.param",
        "pyscf.scf", "pyscf.scf.dispersion",
    )}
    modules["pyscf"].__version__ = "synthetic-test-only"
    modules["pyscf"].__path__ = []
    modules["pyscf.scf"].__path__ = []
    modules["pyscf.dft"].RKS = MeanField
    modules["pyscf.dft"].UKS = MeanField
    modules["pyscf.dft"].libxc = SimpleNamespace(
        libxc_version=lambda: "synthetic-test-only",
        test_deriv_order=lambda *args, **kwargs: None, is_nlc=lambda xc: False,
    )
    modules["pyscf.gto"].M = Molecule
    modules["pyscf.gto"].mole = modules["pyscf.gto.mole"]
    modules["pyscf.gto.mole"].NORMALIZE_GTO = True
    modules["pyscf.lib"].num_threads = num_threads
    modules["pyscf.lib.param"].BOHR = 0.52917721092
    modules["pyscf.lib"].param = modules["pyscf.lib.param"]
    modules["pyscf.scf.dispersion"].parse_dft = lambda xc: (xc, None, None)
    for part in ("dft", "gto", "lib", "scf"):
        setattr(modules["pyscf"], part, modules["pyscf." + part])
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)

    # This boundary intentionally does not imitate or test E3's snapshot schema.
    capture_module = ModuleType("nanodesign.electronic_state")

    def capture_snapshot(mean_field, *, settings, call_id):
        assert mean_field is state.solver_calls[-1]
        snapshot = {
            "kind": "synthetic_capture_boundary_only", "call_id": call_id,
            "positions_angstrom": mean_field.mol.positions.tolist(),
            "settings": deepcopy(settings),
        }
        state.captures.append(deepcopy(snapshot))
        return snapshot

    def save_snapshot(path, snapshot):
        raw = json.dumps(snapshot, allow_nan=False).encode()
        with Path(path).open("xb") as stream:
            stream.write(raw)
        return {"path": str(Path(path).resolve()), "sha256": hashlib.sha256(raw).hexdigest(),
                "size_bytes": len(raw)}

    capture_module.capture_snapshot = capture_snapshot
    capture_module.save_snapshot = save_snapshot
    monkeypatch.setitem(sys.modules, "nanodesign.electronic_state", capture_module)
    return state


def make_calculator(tmp_path):
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.9]])
    atoms.set_constraint(FixAtoms(indices=[1]))
    atoms.calc = PySCFCalculator(
        QuantumSettings(spin=0, threads=1, basis="sto-3g", dispersion=None),
        event_log=tmp_path / "events.jsonl", electronic_state_directory=tmp_path / "snapshots",
    )
    return atoms, atoms.calc


@pytest.mark.parametrize("displacement_fails", [False, True], ids=["completed-hessian", "failed-displacement"])
def test_stationary_restores_baseline_capture_with_baseline_cache(
        tmp_path, synthetic_boundaries, displacement_fails):
    state = synthetic_boundaries
    atoms, calculator = make_calculator(tmp_path)
    baseline_forces = atoms.get_forces(apply_constraint=False)
    baseline_energy = atoms.get_potential_energy()
    baseline_diagnostics = deepcopy(calculator.diagnostics)
    baseline_reference = baseline_diagnostics["electronic_state_snapshot"]
    baseline_path = Path(baseline_reference["path"])
    baseline_bytes = baseline_path.read_bytes()
    baseline_positions = atoms.positions.copy()

    if displacement_fails:
        # Baseline call 1, positive x displacement call 2, negative x call 3.
        state.fail_gradient_at = 3
        state.primary_error = RuntimeError("Synthetic displaced gradient failed")
        with pytest.raises(QuantumCalculationError, match="Synthetic displaced gradient failed"):
            characterize_stationary_point(atoms, step_angstrom=0.005, max_free_coordinates=3)
        expected_calls = 3
    else:
        record = characterize_stationary_point(atoms, step_angstrom=0.005, max_free_coordinates=3)
        np.testing.assert_allclose(record["hessian_ev_per_angstrom2"], np.eye(3))
        assert record["finite_difference_evidence"]["baseline_calculation_call_id"] == baseline_diagnostics["call_id"]
        assert len({row["calculation_call_id"] for row in record["finite_difference_evidence"]["displacements"]}) == 6
        expected_calls = 7

    assert calculator.diagnostics == baseline_diagnostics
    assert calculator.diagnostics["electronic_state_snapshot"]["whole_evaluation_accepted"] is True
    np.testing.assert_array_equal(calculator.atoms.positions, baseline_positions)
    np.testing.assert_array_equal(atoms.positions, baseline_positions)
    np.testing.assert_array_equal(atoms.get_forces(apply_constraint=False), baseline_forces)
    assert atoms.get_potential_energy() == baseline_energy
    assert len(state.solver_calls) == len(state.captures) == state.gradient_calls == expected_calls
    assert state.threads == 3
    assert baseline_path.read_bytes() == baseline_bytes
    saved = [json.loads(path.read_bytes()) for path in (tmp_path / "snapshots").glob("*.json")]
    assert len(saved) == expected_calls
    assert len({item["call_id"] for item in saved}) == expected_calls
    for item in saved:
        if item["call_id"] != baseline_diagnostics["call_id"]:
            assert item["positions_angstrom"] != baseline_positions.tolist()
    if displacement_fails:
        events = [json.loads(line) for line in (tmp_path / "events.jsonl").read_text().splitlines()]
        failure = events[-1]
        assert failure["event"] == "calculation_failed"
        assert failure["call_id"] != baseline_diagnostics["call_id"]
        assert any(item["call_id"] == failure["call_id"] for item in saved)


@pytest.mark.parametrize("error_type", [RuntimeError, KeyboardInterrupt, SystemExit])
@pytest.mark.parametrize("failure_log_fails", [False, True], ids=["normal-failure-log", "failed-failure-log"])
def test_primary_failure_survives_thread_restoration_and_optional_log_failure(
        tmp_path, monkeypatch, synthetic_boundaries, error_type, failure_log_fails):
    state = synthetic_boundaries
    atoms, calculator = make_calculator(tmp_path)
    primary = error_type("Synthetic primary gradient failure")
    state.fail_gradient_at = 1
    state.primary_error = primary
    state.fail_restore = True
    original_event = calculator._event

    if failure_log_fails:
        def event(event_name, *args, **kwargs):
            if event_name in {"calculation_failed", "calculation_interrupted"}:
                raise OSError("Synthetic failure log unavailable")
            return original_event(event_name, *args, **kwargs)
        monkeypatch.setattr(calculator, "_event", event)

    expected_type = QuantumCalculationError if error_type is RuntimeError else error_type
    with pytest.raises(expected_type, match="Synthetic primary gradient failure") as observed:
        atoms.get_forces()
    if error_type is RuntimeError:
        assert observed.value.__cause__ is primary
    else:
        assert observed.value is primary
    assert calculator.results == {}
    diagnostic = calculator.diagnostics
    assert diagnostic["error"] == "Synthetic primary gradient failure"
    assert diagnostic["calculation_status"] == ("failed" if error_type is RuntimeError else "interrupted")
    assert diagnostic["scf_converged"] is True and diagnostic["gradient_completed"] is False
    assert diagnostic["thread_restoration_error"] == {
        "type": "OSError", "message": "Synthetic thread restoration failed",
    }
    assert ("event_log_error" in diagnostic) is failure_log_fails
    if failure_log_fails:
        assert diagnostic["event_log_error"] == "Synthetic failure log unavailable"
    reference = diagnostic["electronic_state_snapshot"]
    assert reference["phase"] == "converged_scf_only"
    assert reference["whole_evaluation_accepted"] is False
    raw = Path(reference["path"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == reference["sha256"]
    assert json.loads(raw)["call_id"] == reference["call_id"] == diagnostic["call_id"]
    assert len(state.solver_calls) == len(state.captures) == state.gradient_calls == 1
    assert state.thread_requests == [1, 3]
    json.dumps(diagnostic, allow_nan=False)


@pytest.mark.parametrize("spin,gradient_fails", [(0, False), (1, True)],
                         ids=["rks-completed", "uks-failed-gradient"])
def test_real_e3_serializer_roundtrip_through_production_hook(
        tmp_path, monkeypatch, synthetic_boundaries, spin, gradient_fails):
    # Load the actual E3 implementation; only the solver/integral side stays fake.
    source = Path(__file__).resolve().parents[2] / "nanodesign/electronic_state.py"
    spec = importlib.util.spec_from_file_location("nanodesign._capture_review_electronic_state", source)
    electronic_state = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(electronic_state)
    monkeypatch.setitem(sys.modules, "nanodesign.electronic_state", electronic_state)
    monkeypatch.setattr(electronic_state, "version", lambda package: "synthetic-solver-test-only")
    atoms = Atoms("H2", positions=[[0.1234567891, 0, 0], [0, 0, 0.9]])
    atoms.calc = calculator = PySCFCalculator(
        QuantumSettings(spin=spin, charge=spin, threads=1, basis="sto-3g", dispersion=None),
        event_log=tmp_path / "events.jsonl", electronic_state_directory=tmp_path / "snapshots",
    )
    state = synthetic_boundaries
    if gradient_fails:
        state.fail_gradient_at = 1
        state.primary_error = RuntimeError("Synthetic post-export gradient failure")
        with pytest.raises(QuantumCalculationError, match="Synthetic post-export gradient failure"):
            atoms.get_forces()
    else:
        atoms.get_forces()
    reference = calculator.diagnostics["electronic_state_snapshot"]
    raw = Path(reference["path"]).read_bytes()
    snapshot = electronic_state.read_snapshot(reference["path"])
    assert hashlib.sha256(raw).hexdigest() == reference["sha256"]
    assert len(raw) == reference["size_bytes"]
    assert snapshot["call_id"] == reference["call_id"] == calculator.diagnostics["call_id"]
    assert snapshot["phase"] == reference["phase"] == "converged_scf_only"
    assert reference["whole_evaluation_accepted"] is (not gradient_fails)
    assert (calculator.results == {}) is gradient_fails
    assert snapshot["settings"] == calculator.settings.to_dict()
    assert snapshot["reference"] == ("RKS" if spin == 0 else "UKS")
    np.testing.assert_allclose(
        np.array(snapshot["geometry"]["positions_bohr"]) * snapshot["geometry"]["bohr_angstrom"],
        atoms.positions, rtol=0, atol=1e-15,
    )
    assert snapshot["ground_state_verified"] is False
    assert snapshot["electronic_state_identity_verified"] is False
    assert "whole_evaluation_accepted" not in snapshot
    assert len(state.solver_calls) == state.gradient_calls == 1
    assert state.captures == []  # Confirms the boundary stub was bypassed.
    # Saved coefficients are detached from the solver's mutable buffers.
    state.solver_calls[0].mo_coeff[:] = 0
    assert Path(reference["path"]).read_bytes() == raw
    assert electronic_state.read_snapshot(reference["path"]) == snapshot
