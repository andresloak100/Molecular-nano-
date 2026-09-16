"""Analytical test potentials check Hessian numerics, not chemistry accuracy."""

from copy import deepcopy
import json

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms, FixBondLength
import numpy as np
import pytest
from scipy.constants import atomic_mass, electron_volt, speed_of_light

from nanodesign.stationary import characterize_stationary_point


class HarmonicTestCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, reference, matrix, fail_after=None):
        super().__init__()
        self.reference = np.array(reference, dtype=float)
        self.matrix = np.array(matrix, dtype=float)
        self.evaluations = []
        self.fail_after = fail_after
        self.diagnostics = {"before": True}

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations.append(atoms.positions.copy())
        self.diagnostics = {"evaluation": len(self.evaluations)}
        if self.fail_after is not None and len(self.evaluations) > self.fail_after:
            self.results = {"forces": np.full_like(atoms.positions, np.nan)}
            raise RuntimeError("deliberate displaced force failure")
        delta = (atoms.positions - self.reference).ravel()
        self.results = {
            "energy": float(delta @ self.matrix @ delta / 2),
            "forces": -(self.matrix @ delta).reshape((-1, 3)),
        }


class DiatomicBondTestCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, spring=2.0, length=0.8):
        super().__init__()
        self.spring, self.length = spring, length
        self.evaluations = []

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations.append(atoms.positions.copy())
        vector = atoms.positions[1] - atoms.positions[0]
        distance = np.linalg.norm(vector)
        bond_force = -self.spring * (distance - self.length) * vector / distance
        self.results = {
            "energy": self.spring * (distance - self.length) ** 2 / 2,
            "forces": np.array([-bond_force, bond_force]),
        }


class DoubleWellTestCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x, y, z = atoms.positions[1]
        self.results = {
            "energy": (x*x - 1)**2 + 2*y*y + 3*z*z,
            "forces": np.array([[0, 0, 0], [-4*x*(x*x - 1), -4*y, -6*z]]),
        }


def test_anchored_diatomic_frequency_has_physical_units_and_no_external_projection():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.8]], masses=[2.0, 1.5])
    atoms.set_constraint(FixAtoms(indices=[0]))
    calculator = DiatomicBondTestCalculator()
    atoms.calc = calculator
    original = atoms.positions.copy()
    result = characterize_stationary_point(atoms, step_angstrom=0.001)
    # k=2 eV/A^2 = 2*eV/1e-20 N/m; only the unfrozen isotope mass moves.
    expected = np.sqrt(2 * electron_volt / 1e-20 / (1.5 * atomic_mass)) / (2*np.pi * speed_of_light * 100)
    assert result["frequencies_cm1"][-1] == pytest.approx(expected, rel=1e-10)
    assert result["free_coordinate_count"] == 3
    assert result["free_masses_amu"] == [1.5]
    assert result["negative_mode_count"] == 0
    assert result["unresolved_near_zero_mode_count"] == 2
    assert result["external_modes_removed"] is False
    assert result["hessian_asymmetry_max_abs_ev_per_angstrom2"] < 1e-12
    assert result["force_requests"] == len(calculator.evaluations) == 7
    assert all(np.array_equal(positions[0], original[0]) for positions in calculator.evaluations)
    np.testing.assert_array_equal(atoms.positions, original)
    assert atoms.calc is calculator
    assert calculator.atoms is None
    assert calculator.results == {}
    modes = np.array(result["cartesian_modes_per_sqrt_amu"])
    np.testing.assert_allclose(np.sum(1.5 * modes**2, axis=(1, 2)), 1)
    json.dumps(result, allow_nan=False)


def test_doublewell_saddle_has_one_resolved_imaginary_mode_but_is_not_verified():
    atoms = Atoms("H2", positions=[[0, 0, -2], [0, 0, 0]], masses=[1, 1])
    atoms.set_constraint(FixAtoms(indices=[0]))
    atoms.calc = DoubleWellTestCalculator()
    step = 0.002
    result = characterize_stationary_point(atoms, step_angstrom=step)
    assert result["stationary_within_force_tolerance"] is True
    assert result["negative_mode_count"] == 1
    assert result["resolved_positive_mode_count"] == 2
    assert result["unresolved_near_zero_mode_count"] == 0
    assert result["classification"] == "first_order_saddle_candidate"
    assert result["frequencies_cm1"][0] < 0
    np.testing.assert_allclose(result["hessian_ev_per_angstrom2"], np.diag([-4 + 4*step**2, 4, 6]), atol=1e-12)
    assert result["irc_performed"] is False
    assert result["transition_state_verified"] is False


def test_free_diatomic_reports_all_3n_modes_with_unresolved_external_modes():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.8]], masses=[1, 1])
    atoms.calc = DiatomicBondTestCalculator()
    result = characterize_stationary_point(atoms, step_angstrom=0.001)
    assert len(result["frequencies_cm1"]) == 6
    assert result["resolved_positive_mode_count"] == 1
    assert result["unresolved_near_zero_mode_count"] == 5
    assert result["external_modes_removed"] is False
    assert "3N" in result["coordinate_space"]


def test_multiple_frozen_groups_omit_anchor_forces_from_stationarity_residual():
    reference = np.zeros((3, 3))
    atoms = Atoms("H3", positions=[[3, 0, 0], [0.02, 0, 0], [0, 0, 5]])
    atoms.set_constraint([FixAtoms(indices=[0]), FixAtoms(indices=[2])])
    atoms.calc = HarmonicTestCalculator(reference, np.eye(9))
    result = characterize_stationary_point(atoms)
    assert result["free_atom_indices"] == [1]
    assert result["frozen_atom_indices"] == [0, 2]
    assert result["free_force_max_ev_per_angstrom"] == pytest.approx(0.02)
    assert result["stationary_within_force_tolerance"] is True
    np.testing.assert_allclose(result["hessian_ev_per_angstrom2"], np.eye(3), atol=1e-12)
    for positions in atoms.calc.evaluations:
        np.testing.assert_array_equal(positions[[0, 2]], atoms.positions[[0, 2]])


def test_nonstationary_geometry_is_not_called_a_saddle_even_with_one_negative_mode():
    atoms = Atoms("H", positions=[[0.2, 0, 0]])
    atoms.calc = HarmonicTestCalculator([[0, 0, 0]], np.diag([-1, 2, 3]))
    result = characterize_stationary_point(atoms)
    assert result["negative_mode_count"] == 1
    assert result["free_force_max_ev_per_angstrom"] == pytest.approx(0.2)
    assert result["stationary_within_force_tolerance"] is False
    assert result["classification"] == "not_stationary_within_force_tolerance"


def test_asymmetric_force_derivatives_are_reported_before_symmetrizing():
    # Deliberately nonconservative forces test the diagnostic; these are not a
    # physically valid molecular energy/force pair.
    matrix = np.array([[3, 0.5, 0], [0.1, 2, 0.2], [0, 0.2, 1]])
    atoms = Atoms("H", positions=[[0, 0, 0]])
    atoms.calc = HarmonicTestCalculator(atoms.positions, matrix)
    result = characterize_stationary_point(atoms)
    np.testing.assert_allclose(result["hessian_ev_per_angstrom2"], (matrix + matrix.T)/2)
    assert result["hessian_asymmetry_max_abs_ev_per_angstrom2"] == pytest.approx(0.4)
    assert result["hessian_asymmetry_relative_frobenius"] == pytest.approx(np.linalg.norm(matrix-matrix.T)/np.linalg.norm(matrix))


@pytest.mark.parametrize("failing", [False, True])
def test_existing_calculator_cache_diagnostics_and_geometry_are_restored(failing):
    atoms = Atoms("H", positions=[[0.1, 0.2, 0.3]])
    calculator = HarmonicTestCalculator([[0, 0, 0]], np.eye(3), fail_after=3 if failing else None)
    atoms.calc = calculator
    before_forces = atoms.get_forces()
    before_positions = atoms.positions.copy()
    before_calc_positions = calculator.atoms.positions.copy()
    before_diagnostics = deepcopy(calculator.diagnostics)
    before_energy = atoms.get_potential_energy()
    if failing:
        with pytest.raises(RuntimeError, match="deliberate displaced"):
            characterize_stationary_point(atoms)
    else:
        characterize_stationary_point(atoms)
    np.testing.assert_array_equal(atoms.positions, before_positions)
    np.testing.assert_array_equal(calculator.atoms.positions, before_calc_positions)
    np.testing.assert_array_equal(calculator.results["forces"], before_forces)
    assert calculator.diagnostics == before_diagnostics
    assert calculator.results["energy"] == before_energy
    number_of_calculations = len(calculator.evaluations)
    np.testing.assert_array_equal(atoms.get_forces(), before_forces)
    assert len(calculator.evaluations) == number_of_calculations


def test_small_negative_eigenvalues_are_unresolved_not_resolved_imaginary_modes():
    atoms = Atoms("H", positions=[[0, 0, 0]], masses=[1])
    atoms.calc = HarmonicTestCalculator(atoms.positions, np.diag([-1e-4, 1, 2]))
    result = characterize_stationary_point(atoms)
    assert result["raw_negative_eigenvalue_count"] == 1
    assert result["negative_mode_count"] == 0
    assert result["unresolved_near_zero_mode_count"] == 1
    assert result["mode_classifications"][0] == "unresolved_near_zero"
    tighter = characterize_stationary_point(atoms, frequency_tolerance_cm1=1)
    assert tighter["negative_mode_count"] == 1


@pytest.mark.parametrize("options", [
    {"step_angstrom": 0}, {"step_angstrom": float("nan")},
    {"force_tolerance_ev_per_angstrom": -1}, {"frequency_tolerance_cm1": True},
    {"max_free_coordinates": 2.5}, {"max_free_coordinates": 2},
])
def test_invalid_settings_and_size_limits_fail_before_any_evaluation(options):
    atoms = Atoms("H", positions=[[0, 0, 0]])
    atoms.calc = HarmonicTestCalculator(atoms.positions, np.eye(3))
    with pytest.raises(ValueError):
        characterize_stationary_point(atoms, **options)
    assert atoms.calc.evaluations == []


def test_unsupported_constraints_and_fully_frozen_system_are_rejected():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.8]])
    atoms.calc = DiatomicBondTestCalculator()
    atoms.set_constraint(FixBondLength(0, 1))
    with pytest.raises(ValueError, match="Only FixAtoms"):
        characterize_stationary_point(atoms)
    atoms.set_constraint(FixAtoms(indices=[0, 1]))
    with pytest.raises(ValueError, match="must be free"):
        characterize_stationary_point(atoms)
    assert atoms.calc.evaluations == []
