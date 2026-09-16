"""Analytic cache/precision regression; not a molecular calculation."""
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from nanodesign.stationary import characterize_stationary_point


class Quadratic(Calculator):
    implemented_properties = ["forces"]

    def __init__(self, center=(0.0, 0.0, 0.0)):
        super().__init__()
        self.center = np.asarray(center)
        self.evaluations = 0

    def calculate(self, atoms=None, properties=("forces",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.evaluations += 1
        self.results = {"forces": -(atoms.positions - self.center) * [-1.0, 2.0, 3.0]}


def test_default_step_preserves_known_negative_mode():
    atoms = Atoms("H", positions=[[0, 0, 0]], masses=[1])
    atoms.calc = Quadratic()
    result = characterize_stationary_point(atoms)
    np.testing.assert_allclose(result["mass_weighted_eigenvalues_ev_per_angstrom2_amu"], [-1, 2, 3])
    assert result["negative_mode_count"] == 1
    assert result["transition_state_verified"] is False
    assert atoms.calc.evaluations == 7


@pytest.mark.parametrize("step", [1e-15, 1e-16])
def test_cache_invisible_step_is_rejected(step):
    atoms = Atoms("H", positions=[[0, 0, 0]], masses=[1])
    atoms.calc = Quadratic()
    with pytest.raises(ValueError, match="step"):
        characterize_stationary_point(atoms, step_angstrom=step)


@pytest.mark.parametrize("center,step", [(1e16, 0.005), (1e8, 1e-8)])
def test_unrepresentable_or_distorted_step_is_rejected(center, step):
    # First case rounds entirely to zero; second shifts by 1.49e-8 instead of 1e-8.
    atoms = Atoms("H", positions=[[center, 0, 0]], masses=[1])
    atoms.calc = Quadratic(center=(center, 0, 0))
    with pytest.raises(ValueError, match="step"):
        characterize_stationary_point(atoms, step_angstrom=step)
