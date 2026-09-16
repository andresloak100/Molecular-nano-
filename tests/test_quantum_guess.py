"""Explicit DFT starting guesses, with no implicit state selection or fallback."""

import json

from ase import Atoms
import numpy as np
import pytest
from pyscf.dft import uks
from pyscf.scf import hf

from nanodesign.quantum import PySCFCalculator, QuantumCalculationError, QuantumSettings


GUESSES = ("minao", "atom", "1e", "huckel")


def test_legacy_settings_keep_minao_default():
    # Documents written before the setting existed must retain PySCF's prior
    # starting guess; do not silently change their electronic-state search.
    legacy = {"spin": 0, "basis": "sto-3g", "dispersion": None}
    settings = QuantumSettings(**legacy)
    assert settings.scf_initial_guess == "minao"
    assert QuantumSettings(**json.loads(json.dumps(settings.to_dict()))) == settings


@pytest.mark.parametrize("guess", [None, True, 1, [], {}, "", "ATOM", "chk", "vsap", "invented"])
def test_unsupported_guess_is_rejected_before_calculation(guess):
    # PySCF itself can silently map an unknown key to minao. The public
    # settings boundary must reject it before any solver is constructed.
    with pytest.raises(ValueError, match="scf_initial_guess"):
        QuantumSettings(scf_initial_guess=guess)


@pytest.mark.parametrize("guess", GUESSES)
@pytest.mark.parametrize("spin", [0, 1], ids=["rks", "uks"])
@pytest.mark.parametrize("density_fit", [False, True], ids=["direct", "density-fit"])
def test_guess_reaches_final_solver_and_failure_is_preserved(monkeypatch, tmp_path, guess, spin, density_fit):
    calls = []

    def fail_at_kernel(mean_field, *args, **kwargs):
        calls.append((mean_field.init_guess, mean_field.mol.spin, bool(getattr(mean_field, "with_df", None))))
        raise RuntimeError("controlled SCF failure")

    monkeypatch.setattr(hf.SCF, "kernel", fail_at_kernel)
    settings = QuantumSettings(
        spin=spin, charge=spin, basis="sto-3g", dispersion=None,
        threads=1, density_fit=density_fit, scf_initial_guess=guess,
    )
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.91]])
    log = tmp_path / "electronic.jsonl"
    atoms.calc = calculator = PySCFCalculator(settings, event_log=log)
    with pytest.raises(QuantumCalculationError, match="controlled SCF failure"):
        atoms.get_forces()

    # One call proves there was no silent retry with a different guess. Test
    # the final density-fitted object as well as the unrestricted/restricted
    # factories, instead of merely testing that a setting was echoed.
    assert calls == [(guess, spin, density_fit)]
    assert calculator.results == {}
    diagnostics = calculator.diagnostics
    assert diagnostics["settings"]["scf_initial_guess"] == diagnostics["scf_initial_guess"] == guess
    for key in ("scf_converged", "gradient_completed", "initial_guess_scan_performed",
                "ground_state_verified", "electronic_state_identity_verified"):
        assert diagnostics[key] is False
    assert "controlled SCF failure" in diagnostics["error"]
    json.dumps(diagnostics, allow_nan=False)
    event = json.loads(log.read_text().splitlines()[-1])
    assert event["event"] == "calculation_failed"
    assert event["call_id"] == diagnostics["call_id"]


@pytest.mark.quantum
def test_real_open_shell_density_fitted_calculation_uses_atom_guess(monkeypatch):
    observed = []
    original = uks.UKS.get_init_guess

    def record_guess(mean_field, mol=None, key="minao", **kwargs):
        observed.append(key)
        return original(mean_field, mol, key, **kwargs)

    # Observe the molecular UKS dispatcher. Its atom guess performs separate
    # atomic HF solves internally; those auxiliary guesses are not molecular
    # retries or changes to the requested molecular starting guess.
    # LiH+ is a small three-electron UKS calculation, with analytic forces.
    monkeypatch.setattr(uks.UKS, "get_init_guess", record_guess)
    atoms = Atoms("LiH", positions=[[0, 0, 0], [0, 0, 1.6]])
    atoms.calc = calculator = PySCFCalculator(QuantumSettings(
        charge=1, spin=1, basis="sto-3g", dispersion=None, grid_level=1,
        threads=1, memory_mb=256, density_fit=True, scf_initial_guess="atom",
    ))
    forces = atoms.get_forces()
    assert observed == ["atom"]
    assert forces.shape == (2, 3) and np.isfinite(forces).all()
    assert np.isfinite(atoms.get_potential_energy())
    assert observed == ["atom"]  # ASE's cache does not launch a second guess.
    diagnostics = calculator.diagnostics
    assert diagnostics["reference"] == "UKS"
    assert diagnostics["electron_count"] == 3
    assert diagnostics["scf_converged"] is diagnostics["gradient_completed"] is True
    assert diagnostics["settings"]["scf_initial_guess"] == diagnostics["scf_initial_guess"] == "atom"
    assert diagnostics["initial_guess_scan_performed"] is False
    assert diagnostics["ground_state_verified"] is False
    assert diagnostics["electronic_state_identity_verified"] is False
    json.dumps(diagnostics, allow_nan=False)
