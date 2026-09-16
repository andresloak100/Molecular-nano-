"""Real scan/capture persistence with installed imports and synthetic numerics.

The actual PySCF module layout and E3 capture/save/read remain in use. Molecule
and solver factories, gradients, and AO overlap values are synthetic; no SCF,
quantum gradients or AO integrals are computed.
"""

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
from uuid import uuid4

import numpy as np
import pytest
from ase import Atoms
from ase.io import write
from ase.units import Hartree
from pyscf import dft, gto, lib

from nanodesign import electronic_state, state_scan
from nanodesign.quantum import PySCFCalculator, QuantumSettings, SCF_INITIAL_GUESSES


@pytest.fixture
def synthetic_solver(monkeypatch):
    state = SimpleNamespace(calls=[], gradients=[], failed_guesses=set(), threads=3)

    class Molecule:
        def __init__(self, **values):
            self.symbols, positions = zip(*values["atom"])
            assert self.symbols == ("H", "H")
            self.positions = np.asarray(positions)
            self.natm = self.nbas = 2
            self.spin, self.charge = values["spin"], values["charge"]
            self.nelectron = 2 - self.charge
            self.nelec = ((self.nelectron + self.spin) // 2, (self.nelectron - self.spin) // 2)
            self.cart = False
            self.ecp = self.pseudo = False
            self._basis = {"H": [[0, [1.0, 1.0]]]}

        def nao_nr(self):
            return 2

        def atom_symbol(self, index):
            return self.symbols[index]

        def atom_charges(self):
            return np.ones(2, dtype=int)

        def atom_coords(self, *, unit):
            assert unit == "Bohr"
            return self.positions / lib.param.BOHR

        def ao_labels(self, *, fmt):
            assert fmt is False
            return [(i, "H", "1s", "") for i in range(2)]

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
            self.converged, self.cycles, self.fitted = True, 1, False
            self.mo_coeff = np.eye(2) if molecule.spin == 0 else np.array([np.eye(2), np.eye(2)])
            self.mo_occ = (np.array([2.0, 0.0]) if molecule.spin == 0
                           else np.array([[1.0, 0.0], [0.0, 0.0]]))

        def istype(self, name):
            return name == ("RKS" if self.mol.spin == 0 else "UKS")

        def get_ovlp(self):
            return np.eye(2)  # Synthetic AO metric: no integral evaluation.

        def density_fit(self):
            result = MeanField(self.mol)
            result.__dict__.update(self.__dict__)
            result.fitted = True
            return result

        def kernel(self):
            self.e_tot = {"minao": -1.0, "atom": -1.1, "1e": -1.1, "huckel": -1.0}[self.init_guess]
            state.calls.append(self)
            return self.e_tot

        def spin_square(self):
            s = self.mol.spin / 2
            return s * (s + 1), self.mol.spin + 1

        def nuc_grad_method(self):
            def gradient():
                state.gradients.append(self.init_guess)
                if self.init_guess in state.failed_guesses:
                    raise RuntimeError("Synthetic gradient failure")
                return np.zeros((2, 3))
            return SimpleNamespace(kernel=gradient)

    def num_threads(value=None):
        if value is not None:
            state.threads = value
        return state.threads

    monkeypatch.setattr(gto, "M", Molecule)
    monkeypatch.setattr(dft, "RKS", MeanField)
    monkeypatch.setattr(dft, "UKS", MeanField)
    monkeypatch.setattr(lib, "num_threads", num_threads)
    return state


def create_scan(tmp_path, *, guesses=SCF_INITIAL_GUESSES, capture=True, spin=1):
    source = tmp_path / "two frames.extxyz"
    frames = [Atoms("H2", positions=[[0.12345678, 0, 0], [0, 0, distance]]) for distance in (0.9, 1.2)]
    write(source, frames)
    root = tmp_path / "original scan"
    settings = QuantumSettings(charge=spin, spin=spin, basis="sto-3g", dispersion=None,
                               grid_level=1, threads=1, density_fit=True)
    state_scan.create_state_scan(source, root, settings, image=0, guesses=guesses,
                                 capture_electronic_state=capture)
    return root, source, settings


def record_path(root, guess="minao"):
    return root / "attempts" / guess / "result.json"


def read_record(root, guess="minao"):
    return json.loads(record_path(root, guess).read_bytes())


def snapshot_path(root, record, guess="minao"):
    call_id = record["quantum_diagnostics"]["electronic_state_snapshot"]["call_id"]
    return root / "attempts" / guess / "electronic-states" / f"{call_id}.json"


def replace_record(root, record, guess="minao"):
    path = record_path(root, guess)
    path.write_text(json.dumps(record, allow_nan=False))
    ledger_path = root / "scan.json"
    ledger = json.loads(ledger_path.read_bytes())
    ledger["attempts"][guess]["result_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    ledger_path.write_text(json.dumps(ledger, allow_nan=False))


def replace_snapshot_and_hashes(root, record, snapshot, guess="minao"):
    raw = json.dumps(snapshot, allow_nan=False).encode()
    snapshot_path(root, record, guess).write_bytes(raw)
    reference = record["quantum_diagnostics"]["electronic_state_snapshot"]
    reference.update(sha256=hashlib.sha256(raw).hexdigest(), size_bytes=len(raw))
    replace_record(root, record, guess)


def tree_bytes(root):
    return {str(path.relative_to(root)): path.read_bytes() for path in root.rglob("*") if path.is_file()}


@pytest.mark.parametrize("spin", [0, 1], ids=["RKS", "UKS"])
def test_all_guesses_capture_real_schema_settings_frame_and_call_binding(tmp_path, synthetic_solver, spin):
    root, _, settings = create_scan(tmp_path, spin=spin)
    plan = json.loads((root / "plan.json").read_bytes())
    assert plan["capture_electronic_state"] is True
    assert state_scan.state_scan_report(root)["capture_electronic_state"] is True
    assert synthetic_solver.calls == []
    report = state_scan.run_state_scan(root, max_jobs=4)
    assert report["converged_guesses"] == list(SCF_INITIAL_GUESSES)
    assert report["complete_four_guess_check"] is True
    assert report["ground_state_verified"] is report["electronic_state_identity_verified"] is False
    calls = []
    for row, solver in zip(report["rows"], synthetic_solver.calls):
        guess = row["guess"]
        record = read_record(root, guess)
        assert record["capture_electronic_state"] is True
        reference = record["quantum_diagnostics"]["electronic_state_snapshot"]
        path = root / row["electronic_state_snapshot"]["path"]
        raw = path.read_bytes()
        snapshot = electronic_state.read_snapshot(path)
        assert path == snapshot_path(root, record, guess)
        assert snapshot["call_id"] == reference["call_id"] == record["quantum_diagnostics"]["call_id"]
        assert snapshot["settings"] == {**asdict(settings), "scf_initial_guess": guess}
        assert snapshot["phase"] == reference["phase"] == "converged_scf_only"
        assert reference["whole_evaluation_accepted"] is True
        assert record["quantum_diagnostics"]["calculation_status"] == "completed"
        assert hashlib.sha256(raw).hexdigest() == reference["sha256"]
        assert len(raw) == reference["size_bytes"]
        assert snapshot["scf_energy_hartree"] == solver.e_tot
        assert record["energy_ev"] == pytest.approx(solver.e_tot * Hartree)
        np.testing.assert_allclose(np.array(snapshot["geometry"]["positions_bohr"]) * lib.param.BOHR,
                                   plan["geometry"]["positions_angstrom"], atol=1e-15, rtol=0)
        assert snapshot["ground_state_verified"] is snapshot["electronic_state_identity_verified"] is False
        assert solver.fitted is True and solver.init_guess == guess
        calls.append(snapshot["call_id"])
    assert len(set(calls)) == 4
    assert synthetic_solver.gradients == list(SCF_INITIAL_GUESSES)
    before = tree_bytes(root)
    assert state_scan.state_scan_report(root) == report
    assert tree_bytes(root) == before
    assert len(synthetic_solver.calls) == 4


def test_relocated_scan_uses_fixed_local_paths_and_can_continue(tmp_path, synthetic_solver):
    root, source, _ = create_scan(tmp_path)
    state_scan.run_state_scan(root)
    record = read_record(root)
    decoy = tmp_path / "untrusted saved absolute path.json"
    decoy.write_text("This is deliberately not electronic evidence")
    record["quantum_diagnostics"]["electronic_state_snapshot"]["path"] = str(decoy)
    replace_record(root, record)
    saved_result = record_path(root).read_bytes()
    moved = tmp_path / "moved scan"
    shutil.move(root, moved)
    source.unlink()
    before = tree_bytes(moved)
    report = state_scan.state_scan_report(moved)
    reference = report["rows"][0]["electronic_state_snapshot"]
    assert not Path(reference["path"]).is_absolute()
    assert (moved / reference["path"]).is_file()
    assert tree_bytes(moved) == before
    resumed = state_scan.run_state_scan(moved, max_jobs=3)
    assert resumed["complete_four_guess_check"] is True
    assert record_path(moved).read_bytes() == saved_result
    assert decoy.read_text() == "This is deliberately not electronic evidence"
    assert [solver.init_guess for solver in synthetic_solver.calls] == list(SCF_INITIAL_GUESSES)


@pytest.mark.parametrize("corruption", ["same-length-bytes", "digest", "size", "bool-size", "deleted"])
def test_exact_snapshot_bytes_and_size_are_required_before_continuation(tmp_path, synthetic_solver, corruption):
    root, _, _ = create_scan(tmp_path)
    state_scan.run_state_scan(root)
    record = read_record(root)
    path = snapshot_path(root, record)
    reference = record["quantum_diagnostics"]["electronic_state_snapshot"]
    if corruption == "same-length-bytes":
        raw = path.read_bytes()
        changed = raw.replace(b'"pbe0"', b'"m06l"', 1)
        assert changed != raw and len(changed) == len(raw)
        path.write_bytes(changed)
    elif corruption == "deleted":
        path.unlink()
    else:
        reference[{"digest": "sha256", "size": "size_bytes", "bool-size": "size_bytes"}[corruption]] = {
            "digest": "0" * 64, "size": reference["size_bytes"] + 1, "bool-size": True,
        }[corruption]
        replace_record(root, record)
    before = tree_bytes(root)
    for inspect in (state_scan.state_scan_report, state_scan.run_state_scan):
        with pytest.raises((ValueError, OSError)):
            inspect(root)
    assert len(synthetic_solver.calls) == 1
    assert tree_bytes(root) == before


@pytest.mark.parametrize("context", ["call", "guess", "method", "geometry", "energy", "bohr-conversion"])
def test_recomputed_hashes_do_not_accept_a_snapshot_from_different_context(tmp_path, synthetic_solver, context):
    root, _, _ = create_scan(tmp_path)
    state_scan.run_state_scan(root)
    record = read_record(root)
    snapshot = electronic_state.read_snapshot(snapshot_path(root, record))
    if context == "call":
        snapshot["call_id"] = str(uuid4())
    elif context == "guess":
        snapshot["settings"]["scf_initial_guess"] = "atom"
    elif context == "method":
        snapshot["settings"]["xc"] = "b3lyp"
    elif context == "geometry":
        snapshot["geometry"]["positions_bohr"][0][0] += 1e-6
    elif context == "energy":
        snapshot["scf_energy_hartree"] += 0.01
    else:
        snapshot["geometry"]["bohr_angstrom"] *= 1.01
        snapshot["geometry"]["positions_bohr"] = (np.array(snapshot["geometry"]["positions_bohr"]) / 1.01).tolist()
    electronic_state.validate_snapshot(snapshot)  # Still a valid standalone E3 record.
    replace_snapshot_and_hashes(root, record, snapshot)
    before = tree_bytes(root)
    with pytest.raises(ValueError, match="snapshot"):
        state_scan.state_scan_report(root)
    with pytest.raises(ValueError, match="snapshot"):
        state_scan.run_state_scan(root)
    assert len(synthetic_solver.calls) == 1 and tree_bytes(root) == before


def test_gradient_failure_retains_scf_only_evidence_and_remaining_guesses_continue(tmp_path, synthetic_solver):
    root, _, _ = create_scan(tmp_path, guesses=("minao", "atom"))
    synthetic_solver.failed_guesses = {"minao"}
    report = state_scan.run_state_scan(root)
    assert report["failed_guesses"] == ["minao"] and report["pending_guesses"] == ["atom"]
    record = read_record(root)
    diagnostic = record["quantum_diagnostics"]
    assert diagnostic["calculation_status"] == "failed"
    assert diagnostic["scf_converged"] is True and diagnostic["gradient_completed"] is False
    assert diagnostic["electronic_state_snapshot"]["whole_evaluation_accepted"] is False
    assert "energy_ev" not in record and "energy_ev" not in report["rows"][0]
    saved = snapshot_path(root, record).read_bytes()
    assert electronic_state.read_snapshot_bytes(saved)["phase"] == "converged_scf_only"
    assert state_scan.state_scan_report(root) == report
    synthetic_solver.failed_guesses.clear()
    resumed = state_scan.run_state_scan(root)
    assert resumed["converged_guesses"] == ["atom"] and resumed["failed_guesses"] == ["minao"]
    assert snapshot_path(root, record).read_bytes() == saved
    assert [solver.init_guess for solver in synthetic_solver.calls] == ["minao", "atom"]


def test_backend_capture_acceptance_is_distinct_from_downstream_survey_rejection(tmp_path, synthetic_solver, monkeypatch):
    class WrongReturnedEnergy(PySCFCalculator):
        def calculate(self, *args, **kwargs):
            super().calculate(*args, **kwargs)
            self.results["energy"] += 1.0  # S2 must reject inconsistent energy accounting.

    monkeypatch.setattr(state_scan, "PySCFCalculator", WrongReturnedEnergy)
    root, _, _ = create_scan(tmp_path, guesses=("minao",))
    report = state_scan.run_state_scan(root)
    assert report["failed_guesses"] == ["minao"]
    record = read_record(root)
    assert "energy units/accounting" in record["error"]["message"]
    assert record["quantum_diagnostics"]["calculation_status"] == "completed"
    assert record["quantum_diagnostics"]["electronic_state_snapshot"]["whole_evaluation_accepted"] is True
    assert report["rows"][0]["electronic_state_snapshot"]["whole_evaluation_accepted"] is True
    assert "energy_ev" not in report["rows"][0]
    assert state_scan.state_scan_report(root) == report


@pytest.mark.parametrize("context", ["reference", "basis_functions"])
def test_failed_gradient_snapshot_still_binds_available_scf_diagnostics(tmp_path, synthetic_solver, context):
    root, _, _ = create_scan(tmp_path, guesses=("minao",), spin=0)
    synthetic_solver.failed_guesses = {"minao"}
    state_scan.run_state_scan(root)
    record = read_record(root)
    if context == "reference":
        snapshot = electronic_state.read_snapshot(snapshot_path(root, record))
        snapshot["reference"] = "UKS"
        snapshot["source_occupations"] = {"alpha": [1.0, 0.0], "beta": [1.0, 0.0]}
        electronic_state.validate_snapshot(snapshot)
        replace_snapshot_and_hashes(root, record, snapshot)
    else:
        record["quantum_diagnostics"]["basis_functions"] = 3
        replace_record(root, record)
    with pytest.raises(ValueError, match="snapshot"):
        state_scan.state_scan_report(root)


@pytest.mark.parametrize("already_completed", [False, True], ids=["legacy-pending", "legacy-completed"])
def test_legacy_plan_without_capture_option_stays_solver_compatible(tmp_path, synthetic_solver, already_completed):
    root, _, _ = create_scan(tmp_path, guesses=("minao",), capture=False)
    if already_completed:
        state_scan.run_state_scan(root)
        record = read_record(root)
        record.pop("capture_electronic_state")
        replace_record(root, record)
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_bytes())
    plan.pop("capture_electronic_state")
    plan_path.write_text(json.dumps(plan))
    ledger_path = root / "scan.json"
    ledger = json.loads(ledger_path.read_bytes())
    ledger["plan_sha256"] = hashlib.sha256(plan_path.read_bytes()).hexdigest()
    ledger_path.write_text(json.dumps(ledger))
    assert state_scan.state_scan_report(root)["capture_electronic_state"] is False
    report = state_scan.run_state_scan(root)
    assert report["converged_guesses"] == ["minao"]
    assert len(synthetic_solver.calls) == 1
    assert "electronic_state_snapshot" not in report["rows"][0]
    assert "capture_electronic_state" not in read_record(root)
    assert not (root / "attempts/minao/electronic-states").exists()


def test_frozen_capture_choice_cannot_be_removed_from_an_existing_attempt(tmp_path, synthetic_solver):
    root, _, _ = create_scan(tmp_path, guesses=("minao",))
    state_scan.run_state_scan(root)
    record = read_record(root)
    record.pop("capture_electronic_state")
    replace_record(root, record)
    with pytest.raises(ValueError, match="capture_electronic_state"):
        state_scan.state_scan_report(root)
