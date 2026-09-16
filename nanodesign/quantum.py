"""Finite-cluster Kohn--Sham DFT energies and analytical forces for ASE.

This is a nonrelativistic, all-electron electronic-structure calculator, not a
validated model of a complete molecular assembler. Energies are Born--Oppenheimer
electronic energies; nuclear zero-point energy, finite-temperature free energies,
environmental effects, and dynamical reaction probabilities are not included.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import json
import math
from numbers import Integral, Real
from pathlib import Path
import threading
import time
from typing import Any
from uuid import uuid4

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree


class QuantumCalculationError(RuntimeError):
    """The requested electronic-structure calculation did not produce a result."""


# Standalone guesses need no checkpoint input and match the supported HF set.
SCF_INITIAL_GUESSES = ("minao", "atom", "1e", "huckel")


@dataclass(frozen=True)
class QuantumSettings:
    """Explicit electronic state, approximation, and numerical controls.

    ``spin`` is N_alpha - N_beta (2S for the intended pure spin state), not the
    multiplicity. A restricted singlet cannot represent a broken-symmetry
    open-shell singlet. ``dispersion='d3bj'`` includes the pairwise D3(BJ)
    correction; the ATM three-body term is explicitly off in this version.
    ``scf_initial_guess`` selects one explicit starting guess. Convergence from
    that guess does not establish electronic-state identity or the ground state.
    """

    charge: int = 0
    spin: int = 1
    xc: str = "pbe0"
    basis: str = "def2-svp"
    dispersion: str | None = "d3bj"
    grid_level: int = 3
    conv_tol: float = 1e-9
    max_cycle: int = 150
    threads: int = 2
    memory_mb: int = 2000
    density_fit: bool = False
    scf_initial_guess: str = "minao"

    def __post_init__(self) -> None:
        for name in ("charge", "spin", "grid_level", "max_cycle", "threads", "memory_mb"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise ValueError(f"{name} must be an integer")
            object.__setattr__(self, name, int(value))
        if self.spin < 0:
            raise ValueError("spin must be nonnegative and equal to N_alpha - N_beta")
        if not 0 <= self.grid_level <= 9:
            raise ValueError("grid_level must be between 0 and 9")
        for name in ("max_cycle", "threads", "memory_mb"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        if (
            isinstance(self.conv_tol, bool)
            or not isinstance(self.conv_tol, Real)
            or not math.isfinite(self.conv_tol)
            or not 0 < self.conv_tol < 1
        ):
            raise ValueError("conv_tol must be finite, positive, and less than 1 Hartree")
        object.__setattr__(self, "conv_tol", float(self.conv_tol))
        for name in ("xc", "basis"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a nonempty string")
        if self.dispersion not in (None, "d3bj", "d3zero"):
            raise ValueError("dispersion must be null, 'd3bj', or 'd3zero'")
        if not isinstance(self.density_fit, bool):
            raise ValueError("density_fit must be a boolean")
        if not isinstance(self.scf_initial_guess, str) or self.scf_initial_guess not in SCF_INITIAL_GUESSES:
            raise ValueError(f"scf_initial_guess must be one of {SCF_INITIAL_GUESSES}; checkpoint guesses are unsupported")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# PySCF's thread count is process-global. Independent instances in the same
# process must not race while changing it. Use processes for parallel images.
_PYSCF_LOCK = threading.RLock()


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


class PySCFCalculator(Calculator):
    """ASE calculator returning energy in eV and forces in eV/Angstrom.

    Every evaluation computes both energy and force and caches them using ASE's
    geometry-change tracking. Failed or nonconverged calculations leave no
    usable energy/force result. ``diagnostics`` describes the latest attempt.
    No functional, spin, basis, or SCF fallback is performed automatically.

    ``event_log`` optionally appends JSON records for SCF iterations and major
    calculation stages. It is an output destination, not a physical setting.
    SCF cycle numbers are one-based; ``e_tot`` in those records is in Hartree
    before adding dispersion. A new call UUID distinguishes successive ASE
    geometries; reads served from ASE's cache do not create additional events.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        settings: QuantumSettings | None = None,
        *,
        event_log: str | Path | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        self._settings = settings if settings is not None else QuantumSettings()
        if not isinstance(self.settings, QuantumSettings):
            raise TypeError("settings must be a QuantumSettings instance")
        self.diagnostics: dict[str, Any] = {}
        self.event_log = Path(event_log) if event_log is not None else None

    @property
    def settings(self) -> QuantumSettings:
        """Read-only settings; construct a new calculator to change the method."""
        return self._settings

    def set(self, **kwargs: Any):
        if kwargs:
            raise ValueError("Set quantum parameters with QuantumSettings when constructing a new calculator")
        return {}

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {}
        started = time.monotonic()
        call_id = str(uuid4())
        self.diagnostics = {
            "call_id": call_id,
            "settings": self.settings.to_dict(),
            "scf_converged": False,
            "gradient_completed": False,
            "scf_initial_guess": self.settings.scf_initial_guess,
            "initial_guess_scan_performed": False,
            "ground_state_verified": False,
            "electronic_state_identity_verified": False,
            "initial_guess_scope": (
                "One explicitly selected SCF starting guess; no automatic scan or "
                "state selection. Convergence does not verify the ground state "
                "or exclude other self-consistent solutions."
            ),
            "model": "nonperiodic all-electron nonrelativistic Kohn-Sham DFT",
            "energy_scope": "Born-Oppenheimer electronic energy plus nuclear repulsion and specified D3",
        }
        try:
            self._validate_atoms()
            # Keep import failures explicit, including missing compiled libraries.
            import pyscf
            from pyscf import dft, gto, lib
            from pyscf.scf.dispersion import parse_dft

            self.diagnostics["versions"] = {
                "pyscf": pyscf.__version__,
                "libxc": dft.libxc.libxc_version(),
                "numpy": np.__version__,
                "ase": _package_version("ase"),
                "dftd3": _package_version("dftd3"),
            }
            parsed_xc, _, native_dispersion = parse_dft(self.settings.xc)
            if native_dispersion is not None:
                raise ValueError("Specify dispersion in the dispersion setting, not inside xc")
            dft.libxc.test_deriv_order(parsed_xc, 1, raise_error=True)
            if self.settings.dispersion and dft.libxc.is_nlc(parsed_xc):
                raise ValueError("Combining nonlocal-correlation xc with D3 requires a separately validated method")

            with _PYSCF_LOCK:
                previous_threads = lib.num_threads()
                lib.num_threads(self.settings.threads)
                # A PySCF build without OpenMP silently ignores the request, so
                # record what the process actually used rather than what was asked.
                effective_threads = int(lib.num_threads())
                self.diagnostics["requested_pyscf_threads"] = self.settings.threads
                self.diagnostics["effective_pyscf_threads"] = effective_threads
                self.diagnostics["threads_honored"] = effective_threads == self.settings.threads
                if effective_threads != self.settings.threads:
                    self.diagnostics["threading_note"] = (
                        f"This PySCF build ran on {effective_threads} thread(s) despite "
                        f"threads={self.settings.threads}; timings reflect the effective count. "
                        "Builds without OpenMP cannot use additional cores within one calculation."
                    )
                try:
                    molecule = gto.M(
                        atom=list(zip(self.atoms.get_chemical_symbols(), self.atoms.positions.tolist())),
                        unit="Angstrom",
                        basis=self.settings.basis,
                        charge=self.settings.charge,
                        spin=self.settings.spin,
                        symmetry=False,
                        verbose=0,
                        max_memory=self.settings.memory_mb,
                    )
                    mean_field = dft.RKS(molecule) if self.settings.spin == 0 else dft.UKS(molecule)
                    mean_field.xc = self.settings.xc
                    # D3 is added once, explicitly, through simple-dftd3 below.
                    mean_field.disp = False
                    mean_field.conv_tol = self.settings.conv_tol
                    mean_field.max_cycle = self.settings.max_cycle
                    mean_field.grids.level = self.settings.grid_level
                    if self.settings.density_fit:
                        mean_field = mean_field.density_fit()
                    mean_field.init_guess = self.settings.scf_initial_guess
                    if self.event_log is not None:
                        def log_scf_cycle(environment):
                            self._event(
                                "scf_cycle", call_id, started,
                                cycle=int(environment["cycle"]) + 1,
                                e_tot=self._finite_or_none(environment["e_tot"]),
                                energy_unit="Hartree",
                            )

                        mean_field.callback = log_scf_cycle
                    dft_energy = float(mean_field.kernel())
                    self.diagnostics.update(
                        scf_converged=bool(mean_field.converged),
                        scf_cycles=int(getattr(mean_field, "cycles", 0)),
                        electron_count=int(molecule.nelectron),
                        alpha_electrons=int(molecule.nelec[0]),
                        beta_electrons=int(molecule.nelec[1]),
                        basis_functions=int(molecule.nao_nr()),
                        reference="RKS" if self.settings.spin == 0 else "UKS",
                    )
                    self._event(
                        "scf_completed", call_id, started,
                        converged=bool(mean_field.converged),
                        cycles=self.diagnostics["scf_cycles"],
                        e_tot=self._finite_or_none(dft_energy),
                        energy_unit="Hartree",
                    )
                    if not mean_field.converged:
                        raise QuantumCalculationError(
                            f"SCF did not converge within {self.settings.max_cycle} cycles; no energy or forces accepted"
                        )
                    actual_s2, multiplicity = mean_field.spin_square()
                    if not math.isfinite(actual_s2) or not math.isfinite(multiplicity):
                        raise QuantumCalculationError("Electronic calculation returned nonfinite spin diagnostics")
                    expected_s = self.settings.spin / 2
                    expected_s2 = expected_s * (expected_s + 1)
                    self.diagnostics.update(
                        s2=float(actual_s2),
                        expected_s2=expected_s2,
                        s2_deviation=float(actual_s2 - expected_s2),
                        spin_contamination=float(actual_s2 - expected_s2),
                        effective_multiplicity=float(multiplicity),
                        stability_checked=False,
                    )
                    gradients = mean_field.nuc_grad_method()
                    gradients.grid_response = True
                    if self.settings.density_fit:
                        gradients.auxbasis_response = True
                    self._event("gradient_started", call_id, started)
                    gradient = np.asarray(gradients.kernel(), dtype=float)
                    dispersion_energy, dispersion_gradient = self._dispersion(molecule)
                    total_energy = dft_energy + dispersion_energy
                    full_gradient = gradient + dispersion_gradient
                    # Use PySCF's own length conversion for an exact derivative
                    # of the geometry supplied to its molecular constructor.
                    bohr_angstrom = float(lib.param.BOHR)
                    forces = -full_gradient * Hartree / bohr_angstrom
                    if not math.isfinite(total_energy) or not np.all(np.isfinite(forces)):
                        raise QuantumCalculationError("Electronic calculation returned a nonfinite energy or force")
                    if forces.shape != (len(self.atoms), 3):
                        raise QuantumCalculationError("Electronic calculation returned an invalid force array")
                    self.diagnostics.update(
                        gradient_completed=True,
                        grid_response=True,
                        auxiliary_basis_response=self.settings.density_fit,
                        dft_energy_hartree=dft_energy,
                        dispersion_energy_hartree=dispersion_energy,
                        total_energy_hartree=total_energy,
                        dispersion_three_body=False,
                        energy_unit="eV",
                        force_unit="eV/Angstrom",
                        hartree_eV=float(Hartree),
                        bohr_angstrom=bohr_angstrom,
                        net_force_eV_A=np.sum(forces, axis=0).tolist(),
                    )
                    self.results = {"energy": total_energy * Hartree, "forces": forces}
                    self._event(
                        "calculation_completed", call_id, started,
                        energy_eV=float(self.results["energy"]),
                        scf_converged=True,
                        gradient_completed=True,
                    )
                finally:
                    lib.num_threads(previous_threads)
        except Exception as error:
            self.results = {}
            self.diagnostics["error"] = str(error)
            try:
                self._event(
                    "calculation_failed", call_id, started,
                    error_type=type(error).__name__, error=str(error),
                )
            except OSError as log_error:
                # Preserve the original scientific or I/O failure if the log
                # destination itself cannot be written.
                self.diagnostics["event_log_error"] = str(log_error)
            if isinstance(error, (QuantumCalculationError, ValueError)):
                raise
            raise QuantumCalculationError(f"Requested quantum calculation failed: {error}") from error
        finally:
            self.diagnostics["elapsed_seconds"] = time.monotonic() - started

    @staticmethod
    def _finite_or_none(value: Real) -> float | None:
        number = float(value)
        return number if math.isfinite(number) else None

    def _event(self, event: str, call_id: str, started: float, **fields: Any) -> None:
        """Append one independently parseable record without retaining a buffer."""
        if self.event_log is None:
            return
        record = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "call_id": call_id,
            "event": event,
            "elapsed_seconds": time.monotonic() - started,
            **fields,
        }
        line = json.dumps(record, allow_nan=False, separators=(",", ":")) + "\n"
        self.event_log.parent.mkdir(parents=True, exist_ok=True)
        with self.event_log.open("a", encoding="utf-8") as stream:
            stream.write(line)

    def _validate_atoms(self) -> None:
        if self.atoms is None or len(self.atoms) == 0:
            raise ValueError("At least one atom is required")
        if np.any(self.atoms.pbc):
            raise ValueError("Periodic boundary conditions are unsupported; use a finite nonperiodic cluster")
        if not np.all(np.isfinite(self.atoms.positions)):
            raise ValueError("Atomic coordinates must be finite")
        numbers = self.atoms.numbers
        if np.any(numbers < 1) or np.any(numbers > 36):
            raise ValueError("This all-electron nonrelativistic backend supports H through Kr; ECPs are not implemented")
        electrons = int(np.sum(numbers)) - self.settings.charge
        if electrons <= 0:
            raise ValueError("The requested charge leaves no electrons")
        if self.settings.spin > electrons or (electrons - self.settings.spin) % 2:
            raise ValueError(
                f"Electron/spin parity mismatch: {electrons} electrons and spin={self.settings.spin}; "
                "spin is N_alpha - N_beta, not multiplicity"
            )
        if len(self.atoms) > 1:
            distances = self.atoms.get_all_distances()
            distances[np.diag_indices_from(distances)] = np.inf
            if float(np.min(distances)) < 0.1:
                raise ValueError("Atoms are closer than 0.1 Angstrom; check the input geometry and units")

    def _dispersion(self, molecule) -> tuple[float, np.ndarray]:
        if self.settings.dispersion is None:
            return 0.0, np.zeros((molecule.natm, 3))
        from dftd3.interface import DispersionModel, RationalDampingParam, ZeroDampingParam

        damping_class = RationalDampingParam if self.settings.dispersion == "d3bj" else ZeroDampingParam
        parameter = damping_class(method=self.settings.xc, atm=False)
        model = DispersionModel(np.asarray(molecule.atom_charges()), molecule.atom_coords())
        correction = model.get_dispersion(parameter, grad=True)
        return float(correction["energy"]), np.asarray(correction["gradient"], dtype=float)
