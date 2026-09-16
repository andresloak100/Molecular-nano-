"""Fixed-geometry, all-electron CCSD(T) references for small finite clusters.

Closed-shell states use RHF/RCCSD(T); open-shell states use UHF/UCCSD(T).
The latter is NOT the open-shell ROHF-based RCCSD(T) variant in Temelso 2006.
This module supplies energies only: it neither relaxes geometry nor establishes
a barrier, reaction connectivity, a complete-basis limit, or chemical accuracy.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import math
from numbers import Integral, Real
import time
from typing import Any

import numpy as np
from ase import Atoms
from ase.units import Hartree

from .quantum import QuantumCalculationError, _PYSCF_LOCK, _package_version


# Documented standalone HF guesses: https://pyscf.org/user/scf.html#initial-guess
# DFT-only 'vsap' and checkpoint guesses (no checkpoint input here) are excluded.
SCF_INITIAL_GUESSES = ("minao", "atom", "1e", "huckel")


@dataclass(frozen=True)
class CCSettings:
    """Explicit state and numerical controls; spin means N_alpha - N_beta.

    No orbitals are frozen. ``memory_mb`` is PySCF's advisory memory setting,
    not a hard process limit. The basis-function cap is checked before SCF and
    before construction of the expensive CC tensors; it is not a cost estimate.
    A restricted singlet cannot describe a broken-symmetry open-shell singlet.
    ``scf_initial_guess`` selects one explicit SCF starting guess. A converged
    solution from any one guess does not certify the electronic ground state.
    """

    charge: int = 0
    spin: int = 0
    basis: str = "cc-pvdz"
    scf_conv_tol: float = 1e-10
    cc_conv_tol: float = 1e-9
    max_cycle: int = 150
    threads: int = 1
    memory_mb: int = 2000
    max_basis_functions: int = 150
    scf_initial_guess: str = "minao"

    def __post_init__(self) -> None:
        for name in ("charge", "spin", "max_cycle", "threads", "memory_mb", "max_basis_functions"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise ValueError(f"{name} must be an integer")
            object.__setattr__(self, name, int(value))
        if self.spin < 0:
            raise ValueError("spin must be nonnegative and equal to N_alpha - N_beta")
        for name in ("max_cycle", "threads", "memory_mb", "max_basis_functions"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        for name in ("scf_conv_tol", "cc_conv_tol"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or not 0 < value < 1:
                raise ValueError(f"{name} must be finite, positive, and less than 1 Hartree")
            object.__setattr__(self, name, float(value))
        if not isinstance(self.basis, str) or not self.basis.strip():
            raise ValueError("basis must be a nonempty string")
        if not isinstance(self.scf_initial_guess, str) or self.scf_initial_guess not in SCF_INITIAL_GUESSES:
            raise ValueError(f"scf_initial_guess must be one of {SCF_INITIAL_GUESSES}; checkpoint and DFT-only guesses are unsupported")


class CoupledClusterCalculationError(QuantumCalculationError):
    """Rejected reference calculation with JSON-safe diagnostics attached."""

    def __init__(self, message: str, diagnostics: dict[str, Any]):
        super().__init__(message)
        self.diagnostics = diagnostics


def _validate_atoms(atoms: Atoms, settings: CCSettings) -> None:
    if not isinstance(atoms, Atoms) or len(atoms) == 0:
        raise ValueError("A nonempty ASE Atoms geometry is required")
    if np.any(atoms.pbc):
        raise ValueError("Periodic systems are unsupported; use a finite nonperiodic cluster")
    if not np.all(np.isfinite(atoms.positions)):
        raise ValueError("Atomic coordinates must be finite")
    if np.any(atoms.numbers < 1) or np.any(atoms.numbers > 36):
        raise ValueError("This all-electron nonrelativistic backend supports H through Kr; no ECPs")
    electrons = int(np.sum(atoms.numbers)) - settings.charge
    if electrons <= 0:
        raise ValueError("The requested charge leaves no electrons")
    if settings.spin > electrons or (electrons - settings.spin) % 2:
        raise ValueError(f"Electron/spin parity mismatch: {electrons} electrons and spin={settings.spin}")
    if len(atoms) > 1:
        distances = atoms.get_all_distances()
        distances[np.diag_indices_from(distances)] = np.inf
        if float(np.min(distances)) < 0.1:
            raise ValueError("Atoms are closer than 0.1 Angstrom; check the geometry and units")


def _finite(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite {label}")
    return number


def coupled_cluster_energy(atoms: Atoms, settings: CCSettings | None = None) -> dict[str, Any]:
    """Return JSON-safe CCSD(T) energy diagnostics on the unchanged geometry.

    Input errors raise ``ValueError``/``TypeError``. Rejected or failed quantum
    work raises ``CoupledClusterCalculationError`` with partial diagnostics;
    no accepted total energy is returned. Interrupts are not swallowed. No
    functional, basis, reference, or state fallback is attempted.
    """
    settings = CCSettings() if settings is None else settings
    if not isinstance(settings, CCSettings):
        raise TypeError("settings must be a CCSettings instance")
    _validate_atoms(atoms, settings)
    # Snapshot coordinates and atomic identities once; do not touch ASE's
    # attached calculator, constraints, results, or user geometry.
    atoms = atoms.copy()
    started = time.monotonic()
    unrestricted = settings.spin != 0
    diagnostics: dict[str, Any] = {
        "schema_version": 2,
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "settings": asdict(settings),
        "reference": "UHF" if unrestricted else "RHF",
        "variant": "UCCSD(T)" if unrestricted else "RCCSD(T)",
        "model": "nonperiodic all-electron nonrelativistic single-reference coupled cluster",
        "energy_scope": "Born-Oppenheimer electronic energy including nuclear repulsion; fixed geometry",
        "geometry": {"symbols": atoms.get_chemical_symbols(), "positions_angstrom": atoms.positions.tolist()},
        "units": {"energy": "Hartree", "secondary_energy": "eV", "length": "angstrom"},
        "hartree_to_ev": float(Hartree),
        "scf_converged": False,
        "ccsd_converged": False,
        "triples_completed": False,
        "geometry_optimized": False,
        "forces_computed": False,
        "scf_stability_checked": False,
        "scf_initial_guess": settings.scf_initial_guess,
        "initial_guess_scan_performed": False,
        "ground_state_verified": False,
        "electronic_state_identity_verified": False,
        "cc_wavefunction_s2_evaluated": False,
        "frozen_core": False,
        "density_fitting": False,
        "dispersion": None,
        "chemical_accuracy_validated": False,
        "limitations": [
            "CCSD(T) is a nonvariational single-reference approximation; convergence does not establish chemical accuracy or suitability for strong correlation.",
            "Basis convergence, multireference character, reference stability and electronic-state identity require separate checks.",
            "One explicitly selected SCF guess is used. Convergence, and even local orbital stability, do not certify the ground state or exclude other self-consistent solutions. No automatic state selection or initial-guess scan is performed.",
            "Reported spin-square values describe the Hartree-Fock determinant, not the correlated CC wavefunction.",
            "Open-shell UHF/UCCSD(T) is different from the ROHF-based RCCSD(T) used in the published Temelso Table 5 comparator.",
            "Fixed-geometry energies omit zero-point/thermal effects and do not verify a saddle, reaction path or molecular tool.",
            "The basis-function cap and advisory memory setting do not guarantee a bounded runtime or peak memory use.",
        ],
    }
    try:
        import pyscf
        from pyscf import cc, gto, lib, scf

        diagnostics["versions"] = {"pyscf": pyscf.__version__, "numpy": np.__version__, "ase": _package_version("ase")}
        with _PYSCF_LOCK:
            previous_threads = lib.num_threads()
            try:
                lib.num_threads(settings.threads)
                effective_threads = int(lib.num_threads())
                diagnostics.update(requested_pyscf_threads=settings.threads,
                                   effective_pyscf_threads=effective_threads,
                                   threads_honored=effective_threads == settings.threads)
                molecule = gto.M(
                    atom=list(zip(atoms.get_chemical_symbols(), atoms.positions.tolist())),
                    unit="Angstrom", basis=settings.basis, charge=settings.charge,
                    spin=settings.spin, symmetry=False, verbose=0, max_memory=settings.memory_mb,
                )
                nao = int(molecule.nao_nr())
                diagnostics.update(basis_functions=nao, electron_count=int(molecule.nelectron),
                                   alpha_electrons=int(molecule.nelec[0]), beta_electrons=int(molecule.nelec[1]))
                if nao > settings.max_basis_functions:
                    raise ValueError(f"Basis has {nao} functions, exceeding max_basis_functions={settings.max_basis_functions}; no SCF or CC calculation started")
                mean_field = scf.UHF(molecule) if unrestricted else scf.RHF(molecule)
                mean_field.init_guess = settings.scf_initial_guess
                mean_field.conv_tol = settings.scf_conv_tol
                mean_field.max_cycle = settings.max_cycle
                mean_field.max_memory = settings.memory_mb
                hf_energy = _finite(mean_field.kernel(), "Hartree-Fock energy")
                diagnostics.update(scf_converged=bool(mean_field.converged),
                                   scf_cycles=int(getattr(mean_field, "cycles", 0)),
                                   hf_energy_hartree=hf_energy)
                if not mean_field.converged:
                    raise ValueError(f"SCF did not converge within {settings.max_cycle} cycles")
                hf_s2, hf_multiplicity = mean_field.spin_square()
                hf_s2 = _finite(hf_s2, "Hartree-Fock spin square")
                expected_s2 = settings.spin / 2 * (settings.spin / 2 + 1)
                diagnostics.update(hf_s2=hf_s2, hf_expected_s2=expected_s2,
                                   hf_s2_deviation=hf_s2 - expected_s2,
                                   hf_effective_multiplicity=_finite(hf_multiplicity, "Hartree-Fock multiplicity"))
                solver = cc.UCCSD(mean_field, frozen=0) if unrestricted else cc.RCCSD(mean_field, frozen=0)
                solver.conv_tol = settings.cc_conv_tol
                solver.conv_tol_normt = math.sqrt(settings.cc_conv_tol)
                solver.max_cycle = settings.max_cycle
                solver.max_memory = settings.memory_mb
                nmo = tuple(int(v) for v in solver.nmo) if unrestricted else (int(solver.nmo),) * 2
                nocc = tuple(int(v) for v in solver.nocc) if unrestricted else (int(solver.nocc),) * 2
                diagnostics.update(
                    molecular_orbitals={"alpha": nmo[0], "beta": nmo[1]},
                    occupied_orbitals={"alpha": nocc[0], "beta": nocc[1]},
                    virtual_orbitals={"alpha": nmo[0] - nocc[0], "beta": nmo[1] - nocc[1]},
                    cc_amplitude_convergence_tolerance=float(solver.conv_tol_normt),
                    frozen_orbitals=0,
                )
                correlation, t1, t2 = solver.kernel()
                diagnostics["ccsd_converged"] = bool(solver.converged)
                if not solver.converged:
                    raise ValueError(f"CCSD did not converge within {settings.max_cycle} cycles")
                correlation = _finite(correlation, "CCSD correlation energy")
                amplitudes = (*t1, *t2) if unrestricted else (t1, t2)
                if any(not np.all(np.isfinite(value)) for value in amplitudes):
                    raise ValueError("Nonfinite CCSD amplitudes")
                ccsd_energy = _finite(solver.e_tot, "CCSD total energy")
                if not math.isclose(ccsd_energy, hf_energy + correlation, rel_tol=0, abs_tol=1e-8):
                    raise ValueError("Inconsistent Hartree-Fock and CCSD energy accounting")
                diagnostics.update(ccsd_correlation_energy_hartree=correlation,
                                   ccsd_total_energy_hartree=ccsd_energy)
                triples = _finite(solver.ccsd_t(), "perturbative triples energy")
                total = _finite(ccsd_energy + triples, "CCSD(T) total energy")
                total_ev = _finite(total * Hartree, "CCSD(T) total energy in eV")
                diagnostics.update(
                    status="completed", triples_completed=True,
                    triples_correction_hartree=triples,
                    total_correlation_energy_hartree=correlation + triples,
                    total_energy_hartree=total, total_energy_ev=total_ev,
                    nuclear_repulsion_energy_hartree=_finite(molecule.energy_nuc(), "nuclear repulsion"),
                )
            finally:
                lib.num_threads(previous_threads)
    except Exception as error:
        diagnostics.update(status="failed", error_type=type(error).__name__, error=str(error))
        # Even a late failure cannot leave a usable accepted total energy.
        diagnostics.pop("total_energy_hartree", None)
        diagnostics.pop("total_energy_ev", None)
        raise CoupledClusterCalculationError(f"Requested CCSD(T) calculation failed: {error}", diagnostics) from error
    finally:
        diagnostics["elapsed_seconds"] = time.monotonic() - started
    return diagnostics
