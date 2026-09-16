"""Actual H1 producer/consumer checks with a synthetic energy/force backend."""
from copy import deepcopy
from dataclasses import asdict
import hashlib
import importlib.util
from pathlib import Path
from uuid import uuid4

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
from ase.units import Hartree

from nanodesign.quantum import QuantumSettings

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("f1_adapter_under_test", HERE / "consistency.py")
f1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f1)


class SyntheticTotalPotential(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, settings, reference, *, fail_on=None, omit_energy=False):
        super().__init__()
        self.settings, self.reference = settings, reference.copy()
        self.fail_on, self.omit_energy, self.calls = fail_on, omit_energy, 0

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls += 1
        self.diagnostics = {"settings": asdict(self.settings), "scf_initial_guess": self.settings.scf_initial_guess,
            "call_id": str(uuid4()), "scf_converged": False, "gradient_completed": False}
        if self.calls == self.fail_on:
            self.diagnostics["error"] = "deliberate synthetic failure"
            raise RuntimeError("deliberate synthetic failure")
        q = atoms.positions - self.reference
        dft = -2.0 + .5 * np.sum(q * q)
        dispersion = .125 + .2 * q[1, 0] + .25 * q[1, 0] ** 2 + .4 * q[0, 0]
        forces = -q.copy()
        forces[1, 0] -= .2 + .5 * q[1, 0]
        forces[0, 0] -= .4
        total = float(dft + dispersion)
        self.results = {"energy": total, "forces": forces}
        self.diagnostics.update(scf_converged=True, gradient_completed=True)
        if not self.omit_energy:
            self.diagnostics.update(total_energy_hartree=total / Hartree,
                dft_energy_hartree=float(dft) / Hartree, dispersion_energy_hartree=float(dispersion) / Hartree,
                hartree_eV=float(Hartree))


@pytest.fixture
def producer(monkeypatch, tmp_path):
    h1 = f1._h1_module()
    monkeypatch.setattr(h1, "PySCFCalculator", SyntheticTotalPotential)
    settings = QuantumSettings(spin=1, threads=1)
    atoms = Atoms("HeH", positions=[[0., 0., 0.], [1.23456789123, .314159265359, .271828182846]])
    atoms.set_constraint(FixAtoms(indices=[0]))
    def make(name, step=.005, fail_on=None, omit_energy=False):
        calc = SyntheticTotalPotential(settings, atoms.positions, fail_on=fail_on, omit_energy=omit_energy)
        atoms.calc = calc
        output = tmp_path / name
        kwargs = dict(settings=settings, input_context={"fixture": "synthetic total potential only"},
                      step=step, fmax=1., max_free_coordinates=3)
        if fail_on:
            with pytest.raises(RuntimeError, match="deliberate synthetic failure"):
                h1.run_characterization_checkpoint(atoms, output, **kwargs)
        else:
            h1.run_characterization_checkpoint(atoms, output, **kwargs)
        return output, calc
    return h1, make


def contents(path):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in path.iterdir() if p.is_file()}


def test_actual_h1_two_steps_total_energy_raw_forces_and_no_new_calculations(producer):
    h1, make = producer
    coarse, first = make("coarse")
    fine, second = make("fine", step=.0025)
    before = [contents(path) for path in (coarse, fine)]
    calls = [first.calls, second.calls]
    document = f1.from_h1_checkpoints([coarse, fine])
    report = f1.analyze_stencils(document, force_tolerance_ev_per_angstrom=1e-10)
    assert report["status"] == "comparison_produced", report
    assert len(report["comparisons"]) == 6
    assert report["all_comparisons_within_declared_tolerance"] is True
    assert all(r["distinct_stencils"] == 2 for r in report["step_sensitivity"])
    assert [first.calls, second.calls] == calls == [7, 7]
    assert [contents(path) for path in (coarse, fine)] == before
    reference = document["stencils"][0]["reference"]
    assert reference["positions_angstrom"][1][0] == 1.23456789123
    assert reference["forces_ev_per_angstrom"][0][0] == -.4
    x = report["comparisons"][0]
    assert x["analytic_force_ev_per_angstrom"] == -.2
    assert x["finite_difference_force_ev_per_angstrom"] == pytest.approx(-.2, abs=1e-10)
    assert x["evidence"]["reference"]["source_energy_conversion"]["from"] == "Hartree"
    assert not report["electronic_branch_continuity_verified"]


def test_valid_force_only_native_checkpoint_is_unavailable_not_invented(producer):
    _, make = producer
    output, calc = make("force-only", omit_energy=True)
    document = f1.from_h1_checkpoints([output])
    assert all(s["reference"]["energy_ev"] is None for s in document["stencils"])
    report = f1.analyze_stencils(document)
    assert report["status"] == "unavailable"
    assert len(report["unavailable_stencils"]) == 3
    assert calc.calls == 7


def test_missing_one_energy_conversion_is_unavailable_not_context_mismatch(producer):
    h1, make = producer
    output, _ = make("missing-conversion")
    payload, _ = h1.load_checkpoint(output)
    del payload["records"][1]["quantum_diagnostics"]["hartree_eV"]
    h1._save_checkpoint(output, payload)
    report = f1.analyze_stencils(f1.from_h1_checkpoints([output]))
    assert report["status"] == "partial", report
    assert len(report["comparisons"]) == 2
    assert len(report["unavailable_stencils"]) == 1


def test_scf_only_energy_does_not_replace_missing_total_energy(producer):
    h1, make = producer
    output, _ = make("scf-only")
    payload, _ = h1.load_checkpoint(output)
    for row in payload["records"]:
        del row["quantum_diagnostics"]["total_energy_hartree"]
    h1._save_checkpoint(output, payload)
    report = f1.analyze_stencils(f1.from_h1_checkpoints([output]))
    assert report["status"] == "unavailable"


def test_partial_failed_checkpoint_preserves_error_and_accepted_subset(producer):
    _, make = producer
    output, calc = make("partial", fail_on=4)
    before = contents(output)
    document = f1.from_h1_checkpoints([output])
    source = document["source_checkpoints"][0]
    assert source["status"] == "failed"
    assert source["error"]["message"] == "deliberate synthetic failure"
    assert source["failed_calculation"]["request"]["id"] == "displacement-0002"
    report = f1.analyze_stencils(document)
    assert report["status"] == "partial", report
    assert len(report["comparisons"]) == 1
    assert len(report["unavailable_stencils"]) == 2
    assert contents(output) == before
    assert calc.calls == 4


def test_no_baseline_native_failure_remains_unavailable(producer):
    _, make = producer
    output, _ = make("no-baseline", fail_on=1)
    report = f1.analyze_stencils(f1.from_h1_checkpoints([output]))
    assert report["status"] == "unavailable", report


def test_corrupted_native_checkpoint_rejected(producer):
    _, make = producer
    output, _ = make("corrupt")
    path = output / "force_checkpoint.json"
    path.write_text(path.read_text().replace('"status":"completed"', '"status":"failed"', 1))
    with pytest.raises(ValueError, match="digest mismatch"):
        f1.from_h1_checkpoints([output])


def test_inconsistent_known_conversions_rejected(producer):
    h1, make = producer
    output, _ = make("conversion-change")
    payload, _ = h1.load_checkpoint(output)
    payload["records"][1]["quantum_diagnostics"]["hartree_eV"] = 27.3
    h1._save_checkpoint(output, payload)
    with pytest.raises(ValueError, match="conversions differ"):
        f1.from_h1_checkpoints([output])


def test_duplicate_input_json_keys_and_nonfinite_values_rejected(tmp_path):
    path = tmp_path / "input.json"
    path.write_text('{"schema_version":1,"schema_version":2}')
    with pytest.raises(ValueError, match="Duplicate JSON"):
        f1.read_stencil_document(path)
    path.write_text('{"schema_version":1,"value":NaN}')
    with pytest.raises(ValueError):
        f1.read_stencil_document(path)


def test_offline_cli_refuses_output_overwrite_and_never_calls_backend(producer, tmp_path):
    _, make = producer
    output, calc = make("cli")
    report = tmp_path / "report.json"
    assert f1.main(["--checkpoints", str(output), "--output", str(report)]) == 0
    before = report.read_bytes()
    with pytest.raises(FileExistsError):
        f1.main(["--checkpoints", str(output), "--output", str(report)])
    assert report.read_bytes() == before
    assert calc.calls == 7
