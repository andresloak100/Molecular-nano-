"""Tiny real AO-integral checks with synthetic determinants; no SCF or gradients.

The E3 producer captures already-populated, metric-orthonormal fixture orbitals.
Every molecular object below has at most two atoms and ten AO functions. Neither
the fixture energies nor the occupied spaces are molecular prediction results.
"""
from math import exp
from pathlib import Path
import sys

import numpy as np
import pytest
from pyscf import gto

from nanodesign.electronic_state import capture_snapshot, save_snapshot
from nanodesign.quantum import QuantumSettings


HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))


class SyntheticMetricDeterminant:
    """Synthetic occupied spaces over real Gaussian AOs, never solved by SCF."""

    def __init__(self, molecule):
        self.mol = molecule
        self.converged = True  # Fixture flag needed to exercise the actual capture API.
        self.e_tot = 0.0  # Synthetic placeholder, never interpreted as a physical energy.
        self.reference = "RKS" if molecule.spin == 0 else "UKS"
        overlap = self.get_ovlp()
        eigenvalues, eigenvectors = np.linalg.eigh(overlap)
        coefficients = eigenvectors @ np.diag(1 / np.sqrt(eigenvalues)) @ eigenvectors.T
        alpha, beta = molecule.nelec
        assert max(alpha, beta) <= molecule.nao_nr()
        if self.reference == "RKS":
            self.mo_coeff = coefficients
            self.mo_occ = np.zeros(molecule.nao_nr())
            self.mo_occ[:alpha] = 2
        else:
            self.mo_coeff = np.stack([coefficients, coefficients])
            self.mo_occ = np.zeros((2, molecule.nao_nr()))
            self.mo_occ[0, :alpha] = 1
            self.mo_occ[1, :beta] = 1

    def istype(self, kind):
        return kind == self.reference

    def get_ovlp(self):
        return self.mol.intor_symmetric("int1e_ovlp")

    def kernel(self, *args, **kwargs):
        raise AssertionError("Integral fixture must never launch SCF")

    def nuc_grad_method(self, *args, **kwargs):
        raise AssertionError("Integral fixture must never launch gradients")


def molecule(positions, basis, *, unit="Bohr", cart=False):
    result = gto.Mole()
    result.atom = [("H", [float(x) for x in point]) for point in positions]
    result.unit = unit
    result.basis = basis
    result.charge = 0
    result.spin = len(positions) % 2
    result.cart = cart
    result.symmetry = False
    result.verbose = 0
    result.build(parse_arg=False, dump_input=False)
    assert result.natm <= 4 and result.nao_nr() <= 10
    return result


def captured_file(directory, name, mol):
    settings = QuantumSettings(spin=mol.spin, basis="synthetic-expanded-integral-fixture",
                               dispersion=None, threads=1).to_dict()
    snapshot = capture_snapshot(SyntheticMetricDeterminant(mol), settings=settings,
                                call_id=f"synthetic-determinant-integral-only-{name}")
    path = directory / f"{name}.json"
    save_snapshot(path, snapshot)
    return path


def bridge(left, right, atom_count):
    # Import the owner's actual adapter; no result or source metadata is patched.
    from cross_overlap import build_cross_overlap
    before = {path: path.read_bytes() for path in (left, right)}
    result = build_cross_overlap(left, right, common_frame_id="synthetic-common-Cartesian-frame",
                                 atom_mapping=list(range(atom_count)))
    assert result["status"] == "computed"
    assert all(path.read_bytes() == content for path, content in before.items())
    return np.asarray(result["cross_overlap"], dtype=float)


@pytest.mark.parametrize("alpha,displacement", [
    (1.37, [.4, 0., 0.]),
    (.4, [.2, -.3, .4]),
])
def test_normalized_single_gaussian_shift_matches_analytical_overlap(tmp_path, alpha, displacement):
    # For equal normalized s primitives, S_AB = exp(-alpha * |A-B|^2 / 2).
    basis = {"H": [[0, [alpha, 1.]]]}
    a = captured_file(tmp_path, "left", molecule([[0., 0., 0.]], basis))
    b = captured_file(tmp_path, "right", molecule([displacement], basis))
    overlap = bridge(a, b, 1)
    expected = exp(-alpha * float(np.dot(displacement, displacement)) / 2)
    assert overlap.shape == (1, 1)
    assert overlap[0, 0] == pytest.approx(expected, abs=2e-13, rel=0)


def test_identical_captured_pair_has_unit_single_ao_overlap(tmp_path):
    basis = {"H": [[0, [.75, 1.]]]}
    mol = molecule([[.123456789123, -.2, .4]], basis)
    a = captured_file(tmp_path, "same-left", mol)
    b = captured_file(tmp_path, "same-right", mol)
    assert bridge(a, b, 1) == pytest.approx(np.ones((1, 1)), abs=2e-13, rel=0)


def test_angstrom_input_is_captured_and_compared_in_bohr_without_double_conversion(tmp_path):
    alpha, shift_angstrom = 1.1, .25
    basis = {"H": [[0, [alpha, 1.]]]}
    left_mol = molecule([[0., 0., 0.]], basis, unit="Angstrom")
    right_mol = molecule([[shift_angstrom, 0., 0.]], basis, unit="Angstrom")
    a = captured_file(tmp_path, "angstrom-left", left_mol)
    b = captured_file(tmp_path, "angstrom-right", right_mol)
    delta_bohr = right_mol.atom_coords(unit="Bohr") - left_mol.atom_coords(unit="Bohr")
    expected = exp(-alpha * float(np.sum(delta_bohr**2)) / 2)
    assert bridge(a, b, 1)[0, 0] == pytest.approx(expected, abs=2e-13, rel=0)
    assert abs(expected - exp(-alpha * shift_angstrom**2 / 2)) > .05


def test_asymmetric_two_atom_pair_reverses_to_transpose(tmp_path):
    left_mol = molecule([[0., 0., 0.], [1.4, 0., 0.]], "sto-3g")
    right_mol = molecule([[.1, .2, 0.], [1.7, -.1, .05]], "sto-3g")
    a = captured_file(tmp_path, "h2-left", left_mol)
    b = captured_file(tmp_path, "h2-right", right_mol)
    forward = bridge(a, b, 2)
    reverse = bridge(b, a, 2)
    expected = gto.intor_cross("int1e_ovlp_sph", left_mol, right_mol)
    assert forward.shape == (2, 2)
    assert forward == pytest.approx(expected, abs=2e-13, rel=0)
    assert reverse == pytest.approx(forward.T, abs=2e-13, rel=0)
    assert np.max(np.abs(forward - forward.T)) > .01


@pytest.mark.parametrize("cart,expected_nao", [(False, 9), (True, 10)])
def test_angular_ao_reconstruction_preserves_spherical_or_cartesian_order(tmp_path, cart, expected_nao):
    # One s, one p and one d shell distinguish the five real-spherical d
    # functions from the six Cartesian d functions, within the ten-AO budget.
    basis = {"H": [[0, [1.1, 1.]], [1, [.8, 1.]], [2, [.6, 1.]]]}
    left_mol = molecule([[0., 0., 0.]], basis, cart=cart)
    right_mol = molecule([[.2, -.15, .4]], basis, cart=cart)
    assert left_mol.nao_nr() == right_mol.nao_nr() == expected_nao
    a = captured_file(tmp_path, "angular-left", left_mol)
    b = captured_file(tmp_path, "angular-right", right_mol)
    operator = "int1e_ovlp_cart" if cart else "int1e_ovlp_sph"
    expected = gto.intor_cross(operator, left_mol, right_mol)
    assert bridge(a, b, 1) == pytest.approx(expected, abs=2e-13, rel=0)
