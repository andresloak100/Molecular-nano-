"""Small actual quantum references plus failure-path checks for CCSD(T)."""

from dataclasses import replace
import json

from ase import Atoms
from ase.units import Hartree
import numpy as np
import pytest
from pyscf import cc, fci, gto, lib, scf

from nanodesign.highlevel import CCSettings, CoupledClusterCalculationError, coupled_cluster_energy
from nanodesign.quantum import _PYSCF_LOCK


def small_settings(**overrides):
    values = dict(basis="sto-3g", scf_conv_tol=1e-12, cc_conv_tol=1e-11, threads=1)
    values.update(overrides)
    return CCSettings(**values)


@pytest.mark.quantum
def test_two_electron_ccsd_matches_fci_in_same_finite_basis():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    before = atoms.positions.copy()
    result = coupled_cluster_energy(atoms, small_settings())
    with _PYSCF_LOCK:
        previous_threads = lib.num_threads()
        try:
            lib.num_threads(1)
            molecule = gto.M(atom="H 0 0 0; H 0 0 0.91", basis="sto-3g", verbose=0)
            reference = scf.RHF(molecule).run(conv_tol=1e-12)
            full_ci = fci.FCI(reference)
            full_ci.conv_tol = 1e-12
            exact_in_basis, _ = full_ci.kernel()
        finally:
            lib.num_threads(previous_threads)
    assert result["ccsd_total_energy_hartree"] == pytest.approx(exact_in_basis, abs=2e-9)
    assert result["triples_correction_hartree"] == pytest.approx(0, abs=1e-14)
    assert result["total_energy_ev"] == pytest.approx(exact_in_basis * Hartree, abs=6e-8)
    assert result["total_energy_hartree"] == pytest.approx(
        result["hf_energy_hartree"] + result["ccsd_correlation_energy_hartree"] + result["triples_correction_hartree"], abs=1e-12)
    assert result["reference"] == "RHF"
    assert result["variant"] == "RCCSD(T)"
    assert result["scf_converged"] and result["ccsd_converged"] and result["triples_completed"]
    assert result["hf_s2"] == pytest.approx(0)
    assert result["hf_expected_s2"] == 0
    assert result["electron_count"] == 2
    assert result["basis_functions"] == 2
    assert result["occupied_orbitals"] == {"alpha": 1, "beta": 1}
    assert result["virtual_orbitals"] == {"alpha": 1, "beta": 1}
    assert result["forces_computed"] is False
    assert result["geometry_optimized"] is False
    assert result["chemical_accuracy_validated"] is False
    assert result["cc_wavefunction_s2_evaluated"] is False
    assert result["dispersion"] is None
    np.testing.assert_array_equal(atoms.positions, before)
    assert atoms.calc is None
    json.dumps(result, allow_nan=False)


@pytest.mark.quantum
def test_open_shell_lih_cation_has_explicit_uhf_reference_and_finite_triples():
    atoms = Atoms("LiH", positions=[[0, 0, 0], [0, 0, 1.6]])
    result = coupled_cluster_energy(atoms, small_settings(charge=1, spin=1))
    assert result["reference"] == "UHF"
    assert result["variant"] == "UCCSD(T)"
    assert result["status"] == "completed"
    assert result["electron_count"] == 3
    assert result["alpha_electrons"] == 2
    assert result["beta_electrons"] == 1
    assert result["occupied_orbitals"] == {"alpha": 2, "beta": 1}
    assert result["hf_expected_s2"] == 0.75
    assert result["hf_s2"] >= 0.75 - 1e-10
    assert result["hf_s2_deviation"] == pytest.approx(result["hf_s2"] - 0.75)
    assert np.isfinite(result["triples_correction_hartree"])
    assert result["total_energy_hartree"] < result["hf_energy_hartree"]
    assert any("ROHF-based RCCSD(T)" in text for text in result["limitations"])
    json.dumps(result, allow_nan=False)


@pytest.mark.quantum
def test_nonconverged_scf_rejected_before_coupled_cluster(monkeypatch):
    def forbidden_cc(*args, **kwargs):
        pytest.fail("CCSD must not run after rejected SCF")
    monkeypatch.setattr(cc, "UCCSD", forbidden_cc)
    atoms = Atoms("LiH", positions=[[0, 0, 0], [0, 0, 1.6]])
    with pytest.raises(CoupledClusterCalculationError, match="SCF did not converge") as caught:
        coupled_cluster_energy(atoms, small_settings(charge=1, spin=1, max_cycle=1))
    assert caught.value.diagnostics["scf_converged"] is False
    assert "total_energy_hartree" not in caught.value.diagnostics
    json.dumps(caught.value.diagnostics, allow_nan=False)


@pytest.mark.quantum
def test_nonconverged_ccsd_rejected_without_triples(monkeypatch):
    real_factory = cc.RCCSD
    def short_cc(*args, **kwargs):
        solver = real_factory(*args, **kwargs)
        original_kernel = solver.kernel
        def stop_early(*args, **kwargs):
            solver.max_cycle = 1
            return original_kernel(*args, **kwargs)
        solver.kernel = stop_early
        def forbidden_triples(*args, **kwargs):
            pytest.fail("Triples must not run after rejected CCSD")
        solver.ccsd_t = forbidden_triples
        return solver
    monkeypatch.setattr(cc, "RCCSD", short_cc)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    with pytest.raises(CoupledClusterCalculationError, match="CCSD did not converge") as caught:
        coupled_cluster_energy(atoms, small_settings())
    diagnostics = caught.value.diagnostics
    assert diagnostics["scf_converged"] is True
    assert diagnostics["ccsd_converged"] is False
    assert diagnostics["triples_completed"] is False
    assert "total_energy_hartree" not in diagnostics


def test_basis_size_limit_rejects_before_scf_and_restores_thread_count(monkeypatch):
    def forbidden_scf(*args, **kwargs):
        pytest.fail("SCF must not run after basis-size rejection")
    monkeypatch.setattr(scf, "RHF", forbidden_scf)
    previous_threads = lib.num_threads()
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    with pytest.raises(CoupledClusterCalculationError, match="exceeding max_basis_functions=1") as caught:
        coupled_cluster_energy(atoms, small_settings(max_basis_functions=1))
    assert caught.value.diagnostics["basis_functions"] == 2
    assert caught.value.diagnostics["scf_converged"] is False
    assert lib.num_threads() == previous_threads


@pytest.mark.quantum
def test_nonfinite_triples_rejected_and_diagnostics_remain_json_safe(monkeypatch):
    real_factory = cc.RCCSD
    def invalid_triples(*args, **kwargs):
        solver = real_factory(*args, **kwargs)
        solver.ccsd_t = lambda: float("nan")
        return solver
    monkeypatch.setattr(cc, "RCCSD", invalid_triples)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    with pytest.raises(CoupledClusterCalculationError, match="Nonfinite perturbative triples") as caught:
        coupled_cluster_energy(atoms, small_settings())
    assert caught.value.diagnostics["status"] == "failed"
    assert "total_energy_hartree" not in caught.value.diagnostics
    json.dumps(caught.value.diagnostics, allow_nan=False)


@pytest.mark.parametrize("overrides,error", [
    ({"charge": True}, "integer"), ({"spin": 1.5}, "integer"),
    ({"spin": -1}, "nonnegative"), ({"threads": 0}, "positive"),
    ({"max_cycle": 0}, "positive"), ({"memory_mb": 0}, "positive"),
    ({"max_basis_functions": 0}, "positive"), ({"basis": ""}, "nonempty"),
    ({"scf_conv_tol": float("nan")}, "finite"), ({"cc_conv_tol": 0}, "finite"),
    ({"cc_conv_tol": True}, "finite"),
])
def test_invalid_settings_rejected(overrides, error):
    with pytest.raises(ValueError, match=error):
        CCSettings(**overrides)


def test_invalid_geometry_and_electronic_state_rejected():
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    settings = small_settings()
    with pytest.raises(ValueError, match="parity mismatch"):
        coupled_cluster_energy(atoms, replace(settings, spin=1))
    atoms.pbc = True
    with pytest.raises(ValueError, match="Periodic"):
        coupled_cluster_energy(atoms, settings)
    atoms.pbc = False
    atoms.positions[1] = atoms.positions[0]
    with pytest.raises(ValueError, match="closer"):
        coupled_cluster_energy(atoms, settings)
    atoms.positions[1, 0] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        coupled_cluster_energy(atoms, settings)
    with pytest.raises(ValueError, match="nonempty"):
        coupled_cluster_energy(Atoms(), settings)
    with pytest.raises(TypeError, match="CCSettings"):
        coupled_cluster_energy(atoms, {})
