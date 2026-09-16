"""Bounded fixed-geometry starting-guess scans, with synthetic calculators only."""

from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from ase.units import Hartree

from nanodesign import state_scan
from nanodesign.quantum import QuantumSettings, SCF_INITIAL_GUESSES


class SyntheticDFT(Calculator):
    """Different known solutions, without importing or invoking a quantum solver."""

    implemented_properties = ["energy", "forces"]
    instances = []
    calls = []
    failure_by_guess = {}
    energies = {"minao": -10.0, "atom": -10.4, "1e": -10.4, "huckel": -10.0}

    def __init__(self, settings=None, *, event_log=None):
        super().__init__()
        self.settings = settings
        self.event_log = Path(event_log) if event_log is not None else None
        self.diagnostics = self.base_diagnostics()
        type(self).instances.append(self)

    def base_diagnostics(self):
        return {
            "synthetic_test_only": True,
            "settings": asdict(self.settings),
            "scf_initial_guess": self.settings.scf_initial_guess,
            "scf_converged": False,
            "gradient_completed": False,
            "ground_state_verified": False,
            "electronic_state_identity_verified": False,
            "initial_guess_scan_performed": False,
            "stability_checked": False,
            "requested_pyscf_threads": self.settings.threads,
            "effective_pyscf_threads": 1,
            "threads_honored": self.settings.threads == 1,
        }

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        guess = self.settings.scf_initial_guess
        type(self).calls.append({
            "guess": guess,
            "settings": asdict(self.settings),
            "numbers": atoms.numbers.copy(),
            "positions": atoms.positions.copy(),
        })
        self.diagnostics = self.base_diagnostics()
        if self.event_log is not None:
            self.event_log.parent.mkdir(parents=True, exist_ok=True)
            with self.event_log.open("a") as stream:
                stream.write(json.dumps({"synthetic_test_only": True, "guess": guess}) + "\n")
        failure = type(self).failure_by_guess.get(guess)
        if failure is not None:
            self.results = {}
            self.diagnostics["error"] = "Synthetic starting-guess failure"
            raise failure("Synthetic starting-guess failure")
        energy = self.energies[guess]
        expected_s = self.settings.spin / 2
        expected_s2 = expected_s * (expected_s + 1)
        observed_s2 = expected_s2 + (0.05 if guess in {"minao", "huckel"} else 0.45)
        forces = np.zeros((len(atoms), 3))
        forces[:, 0] = 0.02
        self.results = {"energy": energy, "forces": forces}
        self.diagnostics.update(
            scf_converged=True,
            gradient_completed=True,
            s2=observed_s2,
            expected_s2=expected_s2,
            s2_deviation=observed_s2 - expected_s2,
            total_energy_hartree=energy / Hartree,
            hartree_eV=float(Hartree),
            elapsed_seconds=0.001,
        )


@pytest.fixture(autouse=True)
def synthetic_backend(monkeypatch):
    SyntheticDFT.instances = []
    SyntheticDFT.calls = []
    SyntheticDFT.failure_by_guess = {}
    monkeypatch.setattr(state_scan, "PySCFCalculator", SyntheticDFT)


def source_geometry(tmp_path, *, atoms=None):
    if atoms is None:
        atoms = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
    path = tmp_path / "source.extxyz"
    write(path, atoms)
    return path


def make_scan(tmp_path, *, guesses=SCF_INITIAL_GUESSES, settings=None):
    source = source_geometry(tmp_path)
    output = tmp_path / "scan"
    settings = settings or QuantumSettings(spin=0, threads=1)
    report = state_scan.create_state_scan(source, output, settings, guesses=guesses)
    return source, output, report


def result_file(root, row):
    return root / row["output"] / "result.json"


def assert_no_certification(report):
    assert report["ground_state_verified"] is False
    assert report["electronic_state_identity_verified"] is False


def test_create_is_solver_free_and_bounded_resumption_keeps_order(tmp_path):
    _, root, report = make_scan(tmp_path)
    guesses = list(SCF_INITIAL_GUESSES)
    assert report["status"] == "pending"
    assert report["requested_guesses"] == report["pending_guesses"] == guesses
    assert report["attempted_guesses"] == []
    assert not SyntheticDFT.instances
    assert_no_certification(report)

    first = state_scan.run_state_scan(root)
    assert first["status"] == "partial"
    assert first["attempted_guesses"] == first["converged_guesses"] == ["minao"]
    assert first["pending_guesses"] == guesses[1:]
    assert first["energy_spread"]["scope"] == "converged_subset"
    assert first["energy_spread"]["member_guesses"] == ["minao"]
    assert first["energy_spread"]["energy_spread_ev"] is None
    assert first["energy_spread"]["energy_spread_hartree"] is None
    assert first["all_requested_attempts_finished"] is False
    assert first["all_requested_guesses_converged"] is False
    assert first["complete_four_guess_check"] is False

    second = state_scan.run_state_scan(root, max_jobs=2)
    assert second["attempted_guesses"] == guesses[:3]
    assert second["pending_guesses"] == ["huckel"]
    final = state_scan.run_state_scan(root, max_jobs=1)
    saved = [result_file(root, row).read_bytes() for row in final["rows"]]
    assert state_scan.run_state_scan(root, max_jobs=4)["status"] == "finished"
    assert [call["guess"] for call in SyntheticDFT.calls] == guesses
    assert len(SyntheticDFT.instances) == 4
    assert [result_file(root, row).read_bytes() for row in final["rows"]] == saved


def test_all_four_preserve_exact_geometry_and_only_guess_setting_changes(tmp_path):
    settings = QuantumSettings(spin=0, threads=1, density_fit=True, basis="sto-3g", grid_level=1)
    source, root, _ = make_scan(tmp_path, settings=settings)
    expected = read(source)
    report = state_scan.run_state_scan(root, max_jobs=4)
    guesses = list(SCF_INITIAL_GUESSES)
    assert report["status"] == "finished"
    assert report["converged_guesses"] == guesses
    assert report["all_requested_attempts_finished"] is True
    assert report["all_requested_guesses_converged"] is True
    assert report["complete_four_guess_check"] is True
    assert report["pending_guesses"] == report["failed_guesses"] == []
    assert report["energy_spread"]["member_guesses"] == guesses
    assert report["energy_spread"]["energy_spread_ev"] == pytest.approx(0.4)
    assert report["energy_spread"]["energy_spread_hartree"] == pytest.approx(0.4 / Hartree)
    assert_no_certification(report)
    assert len(SyntheticDFT.instances) == len(SyntheticDFT.calls) == 4
    assert len({instance.event_log for instance in SyntheticDFT.instances}) == 4
    baseline = asdict(settings)
    baseline.pop("scf_initial_guess")
    for call, row in zip(SyntheticDFT.calls, report["rows"]):
        observed = deepcopy(call["settings"])
        assert observed.pop("scf_initial_guess") == row["guess"]
        assert observed == baseline
        np.testing.assert_array_equal(call["positions"], expected.positions)
        np.testing.assert_array_equal(call["numbers"], expected.numbers)
        assert row["energy_ev"] == SyntheticDFT.energies[row["guess"]]


def test_successful_subset_is_not_a_complete_four_guess_check(tmp_path):
    guesses = ["huckel", "atom"]
    _, root, _ = make_scan(tmp_path, guesses=guesses)
    report = state_scan.run_state_scan(root, max_jobs=2)
    assert report["requested_guesses"] == report["converged_guesses"] == guesses
    assert report["all_requested_guesses_converged"] is True
    assert report["complete_four_guess_check"] is False
    assert [call["guess"] for call in SyntheticDFT.calls] == guesses
    assert_no_certification(report)


def test_failure_is_preserved_while_remaining_guesses_continue(tmp_path):
    _, root, _ = make_scan(tmp_path)
    SyntheticDFT.failure_by_guess = {"atom": RuntimeError}
    report = state_scan.run_state_scan(root, max_jobs=4)
    assert report["status"] == "finished"
    assert report["attempted_guesses"] == list(SCF_INITIAL_GUESSES)
    assert report["failed_guesses"] == ["atom"]
    assert report["converged_guesses"] == ["minao", "1e", "huckel"]
    assert report["all_requested_attempts_finished"] is True
    assert report["all_requested_guesses_converged"] is False
    assert report["complete_four_guess_check"] is False
    assert report["energy_spread"]["member_guesses"] == ["minao", "1e", "huckel"]
    assert_no_certification(report)
    failed_row = next(row for row in report["rows"] if row["guess"] == "atom")
    failed_path = result_file(root, failed_row)
    preserved = failed_path.read_bytes()
    saved = json.loads(preserved)
    assert saved["status"] == "failed"
    assert "Synthetic starting-guess failure" in json.dumps(saved)
    assert "energy_ev" not in failed_row or failed_row["energy_ev"] is None
    SyntheticDFT.failure_by_guess = {}
    assert state_scan.run_state_scan(root, max_jobs=4)["failed_guesses"] == ["atom"]
    assert failed_path.read_bytes() == preserved
    assert len(SyntheticDFT.calls) == 4


def test_interrupt_preserves_attempt_and_leaves_remaining_guesses_pending(tmp_path):
    _, root, _ = make_scan(tmp_path)
    SyntheticDFT.failure_by_guess = {"atom": KeyboardInterrupt}
    with pytest.raises(KeyboardInterrupt):
        state_scan.run_state_scan(root, max_jobs=4)
    report = state_scan.state_scan_report(root)
    assert report["status"] == "partial"
    assert report["converged_guesses"] == ["minao"]
    assert report["interrupted_guesses"] == ["atom"]
    assert report["pending_guesses"] == ["1e", "huckel"]
    assert report["all_requested_attempts_finished"] is False
    assert report["complete_four_guess_check"] is False
    interrupted = next(row for row in report["rows"] if row["guess"] == "atom")
    interrupted_path = result_file(root, interrupted)
    saved = interrupted_path.read_bytes()
    assert json.loads(saved)["status"] == "interrupted"
    SyntheticDFT.failure_by_guess = {}
    resumed = state_scan.run_state_scan(root, max_jobs=4)
    assert resumed["all_requested_attempts_finished"] is True
    assert resumed["interrupted_guesses"] == ["atom"]
    assert resumed["all_requested_guesses_converged"] is False
    assert interrupted_path.read_bytes() == saved
    assert [call["guess"] for call in SyntheticDFT.calls] == list(SCF_INITIAL_GUESSES)


@pytest.mark.parametrize("artifact", ["plan", "captured_input"])
def test_tampered_plan_or_captured_input_rejected_before_computation(tmp_path, artifact):
    _, root, _ = make_scan(tmp_path)
    plan = json.loads((root / "plan.json").read_text())
    path = root / "plan.json" if artifact == "plan" else root / plan["input"]["path"]
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError):
        state_scan.state_scan_report(root)
    with pytest.raises(ValueError):
        state_scan.run_state_scan(root)
    assert not SyntheticDFT.instances


@pytest.mark.parametrize("change", ["modify", "remove"])
def test_terminal_result_integrity_checked_before_read_or_resume(tmp_path, change):
    _, root, _ = make_scan(tmp_path)
    report = state_scan.run_state_scan(root)
    path = result_file(root, report["rows"][0])
    if change == "modify":
        path.write_bytes(path.read_bytes() + b"\n")
    else:
        path.unlink()
    with pytest.raises(ValueError):
        state_scan.state_scan_report(root)
    with pytest.raises(ValueError):
        state_scan.run_state_scan(root, max_jobs=4)
    assert len(SyntheticDFT.calls) == 1


def test_source_changes_and_cached_old_method_do_not_change_frozen_scan(tmp_path):
    initial = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
    initial.calc = SinglePointCalculator(initial, energy=-999.0, forces=np.full((2, 3), 99.0))
    source = source_geometry(tmp_path, atoms=initial)
    root = tmp_path / "scan"
    state_scan.create_state_scan(source, root, QuantumSettings(spin=0), guesses=["atom", "minao"])
    initial.positions[1, 2] = 1.50
    initial.calc = None
    write(source, initial)
    report = state_scan.run_state_scan(root, max_jobs=2)
    assert report["converged_guesses"] == ["atom", "minao"]
    for call in SyntheticDFT.calls:
        assert call["positions"][1, 2] == 0.74
    assert [row["energy_ev"] for row in report["rows"]] == [-10.4, -10.0]


def test_explicit_frame_selection_and_single_atom_supported(tmp_path):
    source = tmp_path / "frames.extxyz"
    frames = [Atoms("H", positions=[[0.0, 0.0, 0.0]]), Atoms("H", positions=[[1.5, 0.0, 0.0]])]
    write(source, frames)
    root = tmp_path / "scan"
    state_scan.create_state_scan(source, root, QuantumSettings(spin=1), image=0, guesses=["atom"])
    report = state_scan.run_state_scan(root)
    assert report["converged_guesses"] == ["atom"]
    np.testing.assert_array_equal(SyntheticDFT.calls[0]["positions"], frames[0].positions)


@pytest.mark.parametrize("guesses", [[], ["atom", "atom"], ["invalid"], [True]])
def test_invalid_guess_plan_is_rejected_without_output_or_calculator(tmp_path, guesses):
    source = source_geometry(tmp_path)
    root = tmp_path / "scan"
    with pytest.raises((TypeError, ValueError)):
        state_scan.create_state_scan(source, root, QuantumSettings(spin=0), guesses=guesses)
    assert not root.exists()
    assert not SyntheticDFT.instances


@pytest.mark.parametrize("image", [True, 0.5, ":"])
def test_invalid_frame_selector_rejected_without_calculator(tmp_path, image):
    source = source_geometry(tmp_path)
    with pytest.raises((TypeError, ValueError)):
        state_scan.create_state_scan(source, tmp_path / "scan", QuantumSettings(spin=0), image=image)
    assert not SyntheticDFT.instances


@pytest.mark.parametrize("max_jobs", [0, -1, True, 1.5, 5])
def test_invalid_job_bound_rejected_before_any_attempt(tmp_path, max_jobs):
    _, root, _ = make_scan(tmp_path)
    with pytest.raises((TypeError, ValueError)):
        state_scan.run_state_scan(root, max_jobs=max_jobs)
    assert state_scan.state_scan_report(root)["attempted_guesses"] == []
    assert not SyntheticDFT.instances


@pytest.mark.parametrize("invalid", ["periodic", "nonfinite", "unsupported_element"])
def test_invalid_geometry_rejected_before_calculator(tmp_path, invalid):
    atoms = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
    if invalid == "periodic":
        atoms.pbc = True
        atoms.cell = [3.0, 3.0, 3.0]
    elif invalid == "nonfinite":
        atoms.positions[1, 0] = np.nan
    else:
        atoms = Atoms("Rb", positions=[[0.0, 0.0, 0.0]])
    source = source_geometry(tmp_path, atoms=atoms)
    with pytest.raises(ValueError):
        state_scan.create_state_scan(source, tmp_path / "scan", QuantumSettings(spin=0))
    assert not SyntheticDFT.instances


def test_existing_scan_cannot_be_overwritten(tmp_path):
    source, root, _ = make_scan(tmp_path)
    plan = (root / "plan.json").read_bytes()
    with pytest.raises(FileExistsError):
        state_scan.create_state_scan(source, root, QuantumSettings(spin=0))
    assert (root / "plan.json").read_bytes() == plan
    assert not SyntheticDFT.instances


@pytest.mark.parametrize("digest", ["absent", None, ""])
def test_completed_attempt_requires_valid_recorded_result_digest(tmp_path, digest):
    _, root, _ = make_scan(tmp_path, guesses=["minao"])
    state_scan.run_state_scan(root)
    ledger_path = root / "scan.json"
    ledger = json.loads(ledger_path.read_text())
    if digest == "absent":
        ledger["attempts"]["minao"].pop("result_sha256")
    else:
        ledger["attempts"]["minao"]["result_sha256"] = digest
    ledger_path.write_text(json.dumps(ledger))
    with pytest.raises(ValueError):
        state_scan.state_scan_report(root)
    with pytest.raises(ValueError):
        state_scan.run_state_scan(root)
    assert len(SyntheticDFT.calls) == 1


@pytest.mark.parametrize("field", ["image", "quantum_settings", "geometry"])
def test_saved_binding_rejects_booleans_aliasing_frozen_numeric_values(tmp_path, field):
    source = source_geometry(tmp_path)
    root = tmp_path / "scan"
    state_scan.create_state_scan(source, root, QuantumSettings(spin=0), image=0, guesses=["minao"])
    report = state_scan.run_state_scan(root)
    path = result_file(root, report["rows"][0])
    record = json.loads(path.read_text())
    if field == "image":
        assert record["image"] == 0
        record["image"] = False
    elif field == "quantum_settings":
        assert record["quantum_settings"]["charge"] == 0
        record["quantum_settings"]["charge"] = False
    else:
        assert record["geometry"]["positions_angstrom"][0][0] == 0.0
        record["geometry"]["positions_angstrom"][0][0] = False
    raw = json.dumps(record).encode()
    path.write_bytes(raw)
    ledger_path = root / "scan.json"
    ledger = json.loads(ledger_path.read_text())
    # Only update the file digest: this isolates a malformed writer's semantic
    # binding from ordinary file-tamper detection, using synthetic evidence.
    ledger["attempts"]["minao"]["result_sha256"] = hashlib.sha256(raw).hexdigest()
    ledger_path.write_text(json.dumps(ledger))
    with pytest.raises(ValueError):
        state_scan.state_scan_report(root)
    with pytest.raises(ValueError):
        state_scan.run_state_scan(root)
    assert len(SyntheticDFT.calls) == 1


def test_same_scan_lock_rejects_another_worker_before_calculator(tmp_path):
    _, root, _ = make_scan(tmp_path)
    with state_scan._lock(root):
        with pytest.raises(RuntimeError):
            state_scan.run_state_scan(root)
    assert state_scan.state_scan_report(root)["attempted_guesses"] == []
    assert not SyntheticDFT.instances
