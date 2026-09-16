"""Numerical checks of the actual quantum backend (no surrogate potentials)."""

from dataclasses import replace
from datetime import datetime
import json
from uuid import UUID

from ase import Atoms
import numpy as np
import pytest

from nanodesign.quantum import PySCFCalculator, QuantumCalculationError, QuantumSettings


def small_settings(**overrides):
    defaults = dict(spin=0, basis="sto-3g", dispersion=None, grid_level=3, threads=1, conv_tol=1e-11)
    defaults.update(overrides)
    return QuantumSettings(**defaults)


def evaluate(atoms, settings):
    atoms = atoms.copy()
    atoms.calc = PySCFCalculator(settings)
    forces = atoms.get_forces()
    return atoms.get_potential_energy(), forces, atoms.calc.diagnostics


def finite_difference_force(atoms, settings, atom, axis, step=1e-4):
    plus, minus = atoms.copy(), atoms.copy()
    plus.positions[atom, axis] += step
    minus.positions[atom, axis] -= step
    return -(evaluate(plus, settings)[0] - evaluate(minus, settings)[0]) / (2 * step)


@pytest.mark.parametrize("dispersion,density_fit", [(None, False), ("d3bj", False), (None, True)])
def test_closed_shell_gradient_matches_energy_derivative(dispersion, density_fit):
    """A displaced H2 bond tests the force sign and both energy/length units."""
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    settings = small_settings(dispersion=dispersion, density_fit=density_fit)
    energy, forces, diagnostics = evaluate(atoms, settings)
    finite_difference = finite_difference_force(atoms, settings, 1, 2)
    assert forces[1, 2] == pytest.approx(finite_difference, abs=3e-5)
    assert energy < -25
    assert np.linalg.norm(forces.sum(axis=0)) < 1e-8
    assert diagnostics["scf_converged"] is True
    assert diagnostics["gradient_completed"] is True
    assert diagnostics["grid_response"] is True
    assert diagnostics["reference"] == "RKS"
    assert diagnostics["s2"] == pytest.approx(0, abs=1e-10)
    assert diagnostics["expected_s2"] == 0
    assert diagnostics["stability_checked"] is False
    if dispersion:
        assert diagnostics["dispersion_energy_hartree"] < 0
        assert diagnostics["dispersion_three_body"] is False
    else:
        assert diagnostics["dispersion_energy_hartree"] == 0
    json.dumps(diagnostics, allow_nan=False)


def test_finite_grid_rotation_error_is_small():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    settings = small_settings()
    energy, forces, _ = evaluate(atoms, settings)
    rotated = atoms.copy()
    rotated.rotate(37, [1, 2, 3])
    rotated_energy, rotated_forces, _ = evaluate(rotated, settings)
    expected_forces = Atoms("H2", positions=forces)
    expected_forces.rotate(37, [1, 2, 3])
    # Atom-centred numerical quadrature has a small, nonzero orientation error.
    assert rotated_energy == pytest.approx(energy, abs=2e-6)
    np.testing.assert_allclose(rotated_forces, expected_forces.positions, atol=2e-5, rtol=0)


def test_ethynyl_radical_gradient_and_translation_invariance():
    atoms = Atoms("C2H", positions=[[0, 0, 0], [0, 0, 1.25], [0.08, 0, 2.34]])
    settings = small_settings(spin=1)
    energy, forces, diagnostics = evaluate(atoms, settings)
    assert forces[2, 2] == pytest.approx(finite_difference_force(atoms, settings, 2, 2), abs=1e-4)
    translated = atoms.copy()
    translated.translate([2.4, -1.8, 3.1])
    shifted_energy, shifted_forces, _ = evaluate(translated, settings)
    assert shifted_energy == pytest.approx(energy, abs=2e-7)
    np.testing.assert_allclose(shifted_forces, forces, atol=2e-5, rtol=0)
    assert np.linalg.norm(forces.sum(axis=0)) < 1e-7
    assert diagnostics["reference"] == "UKS"
    assert diagnostics["electron_count"] == 13
    assert diagnostics["alpha_electrons"] == 7
    assert diagnostics["beta_electrons"] == 6
    assert diagnostics["expected_s2"] == 0.75
    assert diagnostics["s2"] >= 0.75 - 1e-7
    assert diagnostics["s2_deviation"] == pytest.approx(diagnostics["s2"] - 0.75)


def test_nonconverged_scf_is_rejected():
    atoms = Atoms("C2H", positions=[[0, 0, 0], [0, 0, 1.25], [0.08, 0, 2.34]])
    atoms.calc = PySCFCalculator(small_settings(spin=1, max_cycle=1))
    with pytest.raises(QuantumCalculationError, match="SCF did not converge"):
        atoms.get_forces()
    assert atoms.calc.results == {}
    assert atoms.calc.diagnostics["scf_converged"] is False


def test_effective_thread_count_is_reported_not_assumed():
    """Builds without OpenMP ignore the thread request; the record must show it."""
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    _, _, diagnostics = evaluate(atoms, small_settings(threads=2))
    effective = diagnostics["effective_pyscf_threads"]
    assert diagnostics["requested_pyscf_threads"] == 2
    assert isinstance(effective, int) and effective >= 1
    assert diagnostics["threads_honored"] is (effective == 2)
    assert ("threading_note" in diagnostics) is (effective != 2)
    json.dumps(diagnostics, allow_nan=False)


def test_ase_cache_invalidation_after_geometry_change():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    atoms.calc = PySCFCalculator(small_settings())
    first = atoms.get_potential_energy()
    assert atoms.calc.results["forces"].shape == (2, 3)
    atoms.positions[1, 2] = 1.10
    second = atoms.get_potential_energy()
    assert abs(second - first) > 0.01


@pytest.mark.parametrize(
    "settings,error",
    [
        ({"spin": 0.5}, "integer"),
        ({"spin": -1}, "nonnegative"),
        ({"charge": True}, "integer"),
        ({"grid_level": 10}, "between"),
        ({"conv_tol": float("nan")}, "finite"),
        ({"max_cycle": 0}, "positive"),
        ({"dispersion": "fake"}, "dispersion"),
        ({"xc": ""}, "nonempty"),
        ({"basis": ""}, "nonempty"),
        ({"density_fit": "false"}, "boolean"),
    ],
)
def test_invalid_settings_are_rejected(settings, error):
    with pytest.raises(ValueError, match=error):
        QuantumSettings(**settings)


def test_invalid_geometry_and_state_are_rejected():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    atoms.calc = PySCFCalculator(small_settings(spin=1))
    with pytest.raises(ValueError, match="parity mismatch"):
        atoms.get_forces()
    atoms.calc = PySCFCalculator(small_settings())
    atoms.pbc = True
    with pytest.raises(ValueError, match="Periodic"):
        atoms.get_forces()
    atoms.pbc = False
    atoms.positions[1] = atoms.positions[0]
    with pytest.raises(ValueError, match="closer"):
        atoms.get_forces()
    assert atoms.calc.results == {}


def test_method_changes_cannot_silently_reuse_cached_results():
    calculator = PySCFCalculator(small_settings())
    with pytest.raises(AttributeError):
        calculator.settings = replace(calculator.settings, xc="pbe")
    with pytest.raises(ValueError, match="QuantumSettings"):
        calculator.set(xc="pbe")


def test_invalid_functional_does_not_fall_back():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    atoms.calc = PySCFCalculator(small_settings(xc="not_a_real_functional"))
    with pytest.raises((QuantumCalculationError, ValueError)):
        atoms.get_forces()
    assert atoms.calc.results == {}


def test_event_log_reports_real_calculation_stages_and_appends(tmp_path):
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    log_path = tmp_path / "run" / "electronic.jsonl"
    atoms.calc = PySCFCalculator(small_settings(), event_log=log_path)
    atoms.get_forces()
    text_before_cached_read = log_path.read_text()
    atoms.get_potential_energy()
    assert log_path.read_text() == text_before_cached_read

    records = [json.loads(line) for line in text_before_cached_read.splitlines()]
    assert [item["event"] for item in records[-3:]] == [
        "scf_completed", "gradient_started", "calculation_completed",
    ]
    cycles = records[:-3]
    assert cycles and all(item["event"] == "scf_cycle" for item in cycles)
    assert [item["cycle"] for item in cycles] == list(range(1, len(cycles) + 1))
    assert all(np.isfinite(item["e_tot"]) for item in cycles)
    assert records[-3]["converged"] is True
    assert records[-1]["energy_eV"] == atoms.get_potential_energy()
    call_id = atoms.calc.diagnostics["call_id"]
    UUID(call_id)
    assert {item["call_id"] for item in records} == {call_id}
    for item in records:
        assert datetime.fromisoformat(item["timestamp"]).utcoffset().total_seconds() == 0
        assert np.isfinite(item["elapsed_seconds"])
        json.dumps(item, allow_nan=False)

    # A failed later geometry appends an identified failure without overwriting
    # the successful calculation or treating an ASE cache read as a new call.
    atoms.positions[1] = atoms.positions[0]
    with pytest.raises(ValueError, match="closer"):
        atoms.get_forces()
    updated = log_path.read_text()
    assert updated.startswith(text_before_cached_read)
    last = json.loads(updated.splitlines()[-1])
    assert last["event"] == "calculation_failed"
    assert last["call_id"] != call_id
    assert last["error_type"] == "ValueError"
