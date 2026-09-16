"""Step resolution and reconstructible force evidence, using synthetic fields."""

from copy import deepcopy
import json

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
import numpy as np
import pytest

from nanodesign.stationary import characterize_stationary_point


class AnalyticForceField(Calculator):
    """A declared linear force field; nonsymmetric matrices are nonconservative."""

    implemented_properties = ["forces"]

    def __init__(self, reference, matrix, offset=None):
        super().__init__()
        self.reference = np.array(reference, dtype=float)
        self.matrix = np.array(matrix, dtype=float)
        self.offset = (
            np.zeros_like(self.reference) if offset is None
            else np.array(offset, dtype=float)
        )
        self.evaluations = []
        self.diagnostics = {"before": True}

    def forces_at(self, positions):
        delta = (positions - self.reference).ravel()
        return self.offset - (self.matrix @ delta).reshape((-1, 3))

    def calculate(self, atoms=None, properties=("forces",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations.append(atoms.positions.copy())
        self.diagnostics = {
            "evaluation": len(self.evaluations),
            "call_id": f"analytic-force-call-{len(self.evaluations)}",
        }
        self.results = {"forces": self.forces_at(atoms.positions)}


@pytest.mark.parametrize("coordinate,options", [
    pytest.param(0.0, {"step_angstrom": 1e-15}, id="cache-invisible-step"),
    pytest.param(1e16, {}, id="default-step-rounds-to-zero"),
    pytest.param(1e8, {"step_angstrom": 1e-8}, id="step-size-distorted-by-rounding"),
])
def test_unresolved_displacements_are_rejected_before_any_evaluation(coordinate, options):
    # Put the large coordinate last so a late coordinate also has to be
    # preflighted before evaluating earlier, otherwise valid, displacements.
    atoms = Atoms("H2", positions=[[0, 0, 0], [2, 0, coordinate]])
    calculator = AnalyticForceField(atoms.positions, np.eye(6))
    atoms.calc = calculator
    original_positions = atoms.positions.copy()

    with pytest.raises(ValueError, match="step|displacement|resolution"):
        characterize_stationary_point(atoms, **options)

    assert calculator.evaluations == []
    assert calculator.atoms is None
    assert calculator.results == {}
    assert calculator.diagnostics == {"before": True}
    np.testing.assert_array_equal(atoms.positions, original_positions)


def coupled_anchored_field(*, nonconservative=False):
    atoms = Atoms(
        "HCH", positions=[[0.25, 0, -0.5], [1.0, 0.2, 0.3], [2.0, -0.1, 0.75]],
        masses=[1.5, 13.0, 2.5],
    )
    atoms.set_constraint(FixAtoms(indices=[1]))
    matrix = np.diag([-2.0, 2, 3, 4, 5, 6, 7, 8, 9])
    for left, right, value in [(0, 6, 0.4), (1, 7, 0.6), (2, 8, -0.2), (3, 0, 0.5)]:
        matrix[left, right] = matrix[right, left] = value
    if nonconservative:
        matrix[0, 7] = 0.8
        matrix[7, 0] = -0.3
    # A large nonzero anchor load must survive the evidence archive even
    # though it is excluded from the stationarity residual and Hessian.
    offset = [[0.01, 0.005, 0], [5, -4, 7], [0, -0.006, 0.007]]
    atoms.calc = AnalyticForceField(atoms.positions, matrix, offset)
    return atoms, matrix


@pytest.mark.parametrize("cache_baseline", [False, True])
def test_default_step_preserves_coupled_unequal_mass_modes_and_cache(cache_baseline):
    atoms, matrix = coupled_anchored_field()
    calculator = atoms.calc
    positions = atoms.positions.copy()
    if cache_baseline:
        atoms.get_forces(apply_constraint=False)
    prior_calculator_positions = None if calculator.atoms is None else calculator.atoms.positions.copy()
    prior_results = deepcopy(calculator.results)
    prior_diagnostics = deepcopy(calculator.diagnostics)

    result = characterize_stationary_point(atoms)

    free_coordinates = [0, 1, 2, 6, 7, 8]
    expected_hessian = matrix[np.ix_(free_coordinates, free_coordinates)]
    np.testing.assert_allclose(result["hessian_ev_per_angstrom2"], expected_hessian, atol=1e-12)
    masses = np.repeat([1.5, 2.5], 3)
    expected_eigenvalues = np.linalg.eigvalsh(expected_hessian / np.sqrt(np.outer(masses, masses)))
    np.testing.assert_allclose(
        result["mass_weighted_eigenvalues_ev_per_angstrom2_amu"], expected_eigenvalues, atol=1e-12,
    )
    modes = np.asarray(result["cartesian_modes_per_sqrt_amu"]).reshape((6, 6)).T
    np.testing.assert_allclose(modes.T @ np.diag(masses) @ modes, np.eye(6), atol=1e-12)
    assert result["negative_mode_count"] == 1
    assert result["stationary_within_force_tolerance"] is True
    assert result["settings"]["step_angstrom"] == 0.005
    assert result["force_requests"] == len(calculator.evaluations) == 13
    assert result["finite_difference_evidence"]["baseline_calculation_call_id"] == "analytic-force-call-1"
    assert atoms.calc is calculator
    np.testing.assert_array_equal(atoms.positions, positions)
    assert calculator.diagnostics == prior_diagnostics
    if cache_baseline:
        np.testing.assert_array_equal(calculator.atoms.positions, prior_calculator_positions)
        np.testing.assert_array_equal(calculator.results["forces"], prior_results["forces"])
        atoms.get_forces(apply_constraint=False)
        assert len(calculator.evaluations) == 13
    else:
        assert calculator.atoms is None
        assert calculator.results == prior_results == {}
    for evaluated in calculator.evaluations:
        np.testing.assert_array_equal(evaluated[1], positions[1])


def test_saved_full_force_evidence_reconstructs_hessian_and_nonzero_asymmetry():
    atoms, matrix = coupled_anchored_field(nonconservative=True)
    calculator = atoms.calc
    result = characterize_stationary_point(atoms)
    # Reconstruction must also work after the data's actual JSON round trip.
    saved = json.loads(json.dumps(result, allow_nan=False))
    evidence = saved["finite_difference_evidence"]
    reference = np.asarray(evidence["reference_positions_angstrom"])
    baseline = np.asarray(evidence["baseline_forces_ev_per_angstrom"])
    free = saved["free_atom_indices"]
    step = saved["settings"]["step_angstrom"]
    records = evidence["displacements"]
    dimension = saved["free_coordinate_count"]

    np.testing.assert_array_equal(reference, atoms.positions)
    np.testing.assert_allclose(baseline, calculator.forces_at(reference), atol=0)
    assert baseline.shape == (len(atoms), 3)
    assert np.linalg.norm(baseline[1]) > 1
    np.testing.assert_allclose(saved["free_gradient_ev_per_angstrom"], -baseline[free], atol=0)
    assert len(records) == 2 * dimension
    assert evidence["baseline_calculation_call_id"] == "analytic-force-call-1"
    raw_hessian = np.empty((dimension, dimension))

    for column in range(dimension):
        atom_index = free[column // 3]
        axis = column % 3
        pair_forces = []
        for pair_index, sign in enumerate((1, -1)):
            record = records[2 * column + pair_index]
            assert record["calculation_call_id"] == f"analytic-force-call-{2 + 2 * column + pair_index}"
            assert record["atom_index"] == atom_index
            assert record["axis"] == axis
            assert record["requested_offset_angstrom"] == sign * step
            displaced = reference.copy()
            displaced[atom_index, axis] = record["displaced_coordinate_angstrom"]
            actual_offset = displaced[atom_index, axis] - reference[atom_index, axis]
            assert record["actual_offset_angstrom"] == actual_offset
            assert actual_offset == pytest.approx(sign * step, rel=1e-12, abs=0)
            np.testing.assert_array_equal(displaced, calculator.evaluations[1 + 2 * column + pair_index])
            forces = np.asarray(record["forces_ev_per_angstrom"])
            assert forces.shape == (len(atoms), 3)
            np.testing.assert_allclose(forces, calculator.forces_at(displaced), atol=1e-15)
            pair_forces.append(forces[free].ravel())
        raw_hessian[:, column] = -(pair_forces[0] - pair_forces[1]) / (2 * step)

    free_coordinates = [3 * index + axis for index in free for axis in range(3)]
    expected_raw = matrix[np.ix_(free_coordinates, free_coordinates)]
    np.testing.assert_allclose(raw_hessian, expected_raw, atol=1e-12)
    symmetrized = (raw_hessian + raw_hessian.T) / 2
    asymmetry = raw_hessian - raw_hessian.T
    np.testing.assert_allclose(saved["hessian_ev_per_angstrom2"], symmetrized, atol=1e-12)
    assert saved["hessian_asymmetry_max_abs_ev_per_angstrom2"] == pytest.approx(1.1, abs=1e-12)
    assert saved["hessian_asymmetry_max_abs_ev_per_angstrom2"] == pytest.approx(np.abs(asymmetry).max())
    assert saved["hessian_asymmetry_relative_frobenius"] == pytest.approx(
        np.linalg.norm(asymmetry) / np.linalg.norm(raw_hessian),
    )


@pytest.mark.parametrize("cache_baseline", [False, True])
def test_coarser_calculator_cache_rejects_before_displaced_work_and_restores(cache_baseline):
    class CoarseCacheField(AnalyticForceField):
        def check_state(self, atoms, tol=0.02):
            return super().check_state(atoms, tol=tol)

    atoms = Atoms("H", positions=[[0.1, 0.2, 0.3]])
    calculator = CoarseCacheField(atoms.positions, np.diag([-1, 2, 3]))
    atoms.calc = calculator
    positions = atoms.positions.copy()
    if cache_baseline:
        atoms.get_forces()
    prior_results = deepcopy(calculator.results)
    prior_diagnostics = deepcopy(calculator.diagnostics)

    with pytest.raises(ValueError, match="cache"):
        characterize_stationary_point(atoms)

    # Only the baseline was evaluated (before the call when already cached).
    assert len(calculator.evaluations) == 1
    np.testing.assert_array_equal(calculator.evaluations[0], positions)
    np.testing.assert_array_equal(atoms.positions, positions)
    assert calculator.diagnostics == prior_diagnostics
    if cache_baseline:
        np.testing.assert_array_equal(calculator.atoms.positions, positions)
        np.testing.assert_array_equal(calculator.results["forces"], prior_results["forces"])
        atoms.get_forces()
        assert len(calculator.evaluations) == 1
    else:
        assert calculator.atoms is None
        assert calculator.results == prior_results == {}


@pytest.mark.parametrize("call_id", [None, 123])
def test_absent_or_nonstring_calculation_identifiers_are_explicit_nulls(call_id):
    class NoStringCallIdField(AnalyticForceField):
        def calculate(self, atoms=None, properties=("forces",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            if call_id is None:
                del self.diagnostics["call_id"]
            else:
                self.diagnostics["call_id"] = call_id

    atoms = Atoms("H", positions=[[0, 0, 0]])
    atoms.calc = NoStringCallIdField(atoms.positions, np.eye(3))
    evidence = characterize_stationary_point(atoms)["finite_difference_evidence"]

    assert evidence["baseline_calculation_call_id"] is None
    assert len(evidence["displacements"]) == 6
    assert all(record["calculation_call_id"] is None for record in evidence["displacements"])
