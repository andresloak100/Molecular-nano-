"""Paired DFT and CCSD(T) energies on exactly the same archived geometries.

These finite-basis method differences are not verified barriers, experimental
errors or calibrated uncertainty estimates. Reference structures and their
separate CC BY-NC license are archived by the existing benchmark workflow.
"""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from io import StringIO
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
from ase.io import read
from ase.units import Hartree

from . import __version__
from .benchmark import EV_PER_KCAL_PER_MOL, SPECIES, _write_json, run_benchmark
from .highlevel import CCSettings, coupled_cluster_energy
from .quantum import QuantumSettings, SCF_INITIAL_GUESSES


class MethodComparisonError(RuntimeError):
    """The paired calculation could not be accepted; inspect saved evidence."""


def _finite_number(value: Any, label: str) -> float:
    if not isinstance(value, (float, int)) or isinstance(value, bool) or not math.isfinite(value):
        raise MethodComparisonError(f"Nonfinite or invalid {label}")
    return float(value)


def _check_spin_diagnostics(record: dict[str, Any], spin: int, prefix: str, label: str) -> None:
    actual = _finite_number(record.get(f"{prefix}s2"), f"{label} determinant spin square")
    expected_key = "hf_expected_s2" if prefix else "expected_s2"
    deviation_key = "hf_s2_deviation" if prefix else "spin_contamination"
    expected = _finite_number(record.get(expected_key), f"{label} intended spin square")
    deviation = _finite_number(record.get(deviation_key), f"{label} determinant spin deviation")
    target = spin / 2 * (spin / 2 + 1)
    if not math.isclose(expected, target, rel_tol=0, abs_tol=1e-10) or not math.isclose(deviation, actual - target, rel_tol=0, abs_tol=1e-8):
        raise MethodComparisonError(f"Inconsistent {label} determinant spin diagnostics")


def _energy_differences(energies: dict[str, float]) -> dict[str, float]:
    reactants = energies["methane"] + energies["ethynyl_radical"]
    nominal = energies["methane_ethynyl_ts"] - reactants
    reaction = energies["acetylene"] + energies["methyl_radical"] - reactants
    if not all(math.isfinite(value) for value in (reactants, nominal, reaction)):
        raise MethodComparisonError("Nonfinite aggregate energy")
    values = {
        "nominal_ts_relative_energy_ev": nominal,
        "nominal_ts_relative_energy_kcal_per_mol": nominal / EV_PER_KCAL_PER_MOL,
        "reaction_energy_ev": reaction,
        "reaction_energy_kcal_per_mol": reaction / EV_PER_KCAL_PER_MOL,
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise MethodComparisonError("Nonfinite aggregate energy after unit conversion")
    return values


def _check_dft_guess_provenance(
    dft: dict[str, Any], record: dict[str, Any], name: str, expected_guess: str,
) -> None:
    """Check recorded guesses against the request without rewriting evidence.

    This backend used minao before exposing a guess setting. Only a missing
    field has that historical meaning; explicit invalid values are not defaults.
    Other settings are not filled from today's QuantumSettings defaults.
    """
    diagnostics = record.get("quantum_diagnostics", {})
    if not isinstance(diagnostics, dict):
        raise MethodComparisonError(f"Invalid DFT quantum_diagnostics for {name}")
    sources = (
        ("benchmark quantum_settings", dft.get("quantum_settings", {})),
        ("species quantum_settings", record.get("quantum_settings", {})),
        ("quantum_diagnostics.settings", diagnostics.get("settings", {})),
        ("quantum_diagnostics", diagnostics),
    )
    for label, settings in sources:
        if not isinstance(settings, dict):
            raise MethodComparisonError(f"Invalid DFT {label} for {name}")
        guess = settings.get("scf_initial_guess", "minao")
        if not isinstance(guess, str) or guess not in SCF_INITIAL_GUESSES:
            raise MethodComparisonError(f"Invalid DFT scf_initial_guess in {label} for {name}")
        if guess != expected_guess:
            raise MethodComparisonError(f"DFT scf_initial_guess mismatch in {label} for {name}: expected {expected_guess!r}, recorded {guess!r}")


def _paired_input(
    dft: dict[str, Any], output: Path, name: str, basis: str,
    *, expected_dft_guess: str = "minao",
):
    spec = SPECIES[name]
    path = output / "dft" / "inputs" / spec["file"]
    # Parse the exact bytes that were hashed, rather than reading a mutable
    # path a second time after verifying its digest.
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    record = dft["species"][name]
    if digest != dft["input_hashes"][spec["file"]] or digest != record["file_sha256"]:
        raise MethodComparisonError(f"Archived DFT geometry hash mismatch for {name}")
    if record.get("status") != "completed":
        raise MethodComparisonError(f"Incomplete DFT calculation for {name}")
    if record.get("charge") != spec["charge"] or record.get("spin_2s") != spec["spin"]:
        raise MethodComparisonError(f"DFT nominal charge/spin mismatch for {name}")
    _check_dft_guess_provenance(dft, record, name, expected_dft_guess)
    settings = record["quantum_settings"]
    if settings.get("basis") != basis or settings.get("xc") != "pbe0" or settings.get("dispersion") != "d3bj":
        raise MethodComparisonError(f"Unexpected DFT method for {name}")
    if settings.get("charge") != spec["charge"] or settings.get("spin") != spec["spin"]:
        raise MethodComparisonError(f"DFT nominal charge/spin settings mismatch for {name}")
    if record.get("scf_converged") is not True or record.get("quantum_diagnostics", {}).get("gradient_completed") is not True:
        raise MethodComparisonError(f"Incomplete DFT electronic calculation for {name}")
    energy = _finite_number(record.get("energy_ev"), f"DFT energy for {name}")
    hartree = _finite_number(record.get("energy_hartree"), f"DFT Hartree energy for {name}")
    if not math.isclose(energy, hartree * Hartree, abs_tol=1e-8, rel_tol=0):
        raise MethodComparisonError(f"DFT energy unit mismatch for {name}")
    _check_spin_diagnostics(record, spec["spin"], "", f"DFT {name}")
    frames = read(StringIO(raw.decode("utf-8")), format="xyz", index=":")
    if len(frames) != 1:
        raise MethodComparisonError(f"Expected one archived geometry for {name}")
    atoms = frames[0]
    if atoms.get_chemical_formula() != spec["formula"] or not np.all(np.isfinite(atoms.positions)):
        raise MethodComparisonError(f"Unexpected archived geometry for {name}")
    return atoms, digest


def run_method_comparison(output, directory=None, basis: str = "cc-pvdz", *, cc_initial_guess: str = "minao") -> dict[str, Any]:
    """Persist a same-geometry PBE0-D3(BJ) versus CCSD(T) comparison.

    ``output`` must be new; there is no overwrite or implicit resume. The
    initial implementation accepts cc-pVDZ only, rejecting larger campaigns
    before any DFT work. The generic coupled-cluster backend is configurable.
    Completed CC species are saved independently; interrupted/failed runs
    retain their partial evidence and never return a completed comparison.
    ``cc_initial_guess`` applies one explicitly selected guess to every HF
    reference. It does not scan guesses, select a state, or alter DFT settings.
    """
    if not isinstance(basis, str) or basis.strip().lower() != "cc-pvdz":
        raise ValueError("The paired comparison currently supports only cc-pvdz; no calculation started")
    basis = "cc-pvdz"
    dft_settings = QuantumSettings(basis=basis)
    cc_base = CCSettings(basis=basis, scf_initial_guess=cc_initial_guess)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    cc_directory = output / "ccsd_t"
    cc_directory.mkdir()
    started = time.monotonic()
    result: dict[str, Any] = {
        "schema_version": 2, "comparison": "same_geometry_same_basis_pbe0_d3bj_vs_ccsd_t",
        "status": "running", "stage": "dft", "started_utc": datetime.now(timezone.utc).isoformat(),
        "software_version": __version__, "basis": basis,
        "dft_settings": asdict(dft_settings), "cc_settings": asdict(cc_base),
        "dft_evidence": "dft/benchmark.json", "input_manifest": "dft/input_manifest.json",
        "reference_data_license": "CC BY-NC 4.0; separate from repository code license",
        "species": {}, "units": {"energy": "eV", "secondary_energy": "kcal/mol", "length": "angstrom"},
        "geometry_basis_pairing_verified": False,
        "electronic_state_identity_verified": False,
        "ground_state_verified": False, "initial_guess_scan_performed": False,
        "scientific_status": "electronic_state_identity_unverified",
        "method_accuracy_validated": False, "saddle_verified": False, "tool_validated": False,
        "published_barrier_used_as_acceptance_target": False,
        "interpretation": "DFT minus coupled-cluster energy differences at identical fixed geometries and basis with matching nominal charge/spin inputs. Matching those inputs does not verify the same electronic state or a ground-state solution. These differences are not a true error or uncertainty bound.",
        "limitations": [
            "The nominal source transition structure is not a verified first-order saddle; its SI/main-paper geometry and absolute-energy discrepancies remain unresolved.",
            "Matching geometry and basis removes those differences from this paired calculation, but does not establish basis convergence or accuracy of either method.",
            "Open-shell coupled cluster uses UHF/UCCSD(T), not the ROHF-based RCCSD(T) in the published Table 5 comparator.",
            "Single-reference applicability and electronic-state stability are unchecked. HF spin diagnostics do not give correlated CC spin purity; UKS and UHF determinants need not describe the same electronic state.",
            "One explicit CC reference initial guess is used for all species. No automatic guess scan or state selection is performed; a lower result from another guess does not alone certify the ground state.",
            "PBE0-D3(BJ) includes its empirical dispersion correction; this paired difference therefore includes dispersion-model differences.",
            "No geometry relaxation, reaction connectivity, zero-point energy, thermal correction or operational reliability is established by this comparison.",
            "The nested DFT benchmark contains a historical published-value comparison; it is not used to accept this paired calculation.",
        ],
    }
    summary_path = output / "method_comparison.json"
    _write_json(summary_path, result)
    try:
        dft = run_benchmark(dft_settings, directory, output / "dft")
        if dft.get("status") != "completed" or set(dft.get("species", {})) != set(SPECIES):
            raise MethodComparisonError("The DFT campaign did not complete all five reference species")
        result["stage"] = "coupled_cluster"
        _write_json(summary_path, result)
        for name, spec in SPECIES.items():
            result["active_species"] = name
            result["elapsed_seconds"] = time.monotonic() - started
            _write_json(summary_path, result)
            record: dict[str, Any] = {"species": name, "status": "running", "basis": basis,
                                      "charge": spec["charge"], "spin_2s": spec["spin"],
                                      "scf_initial_guess": cc_initial_guess,
                                      "electronic_state_identity_verified": False}
            record_path = cc_directory / f"{name}.json"
            _write_json(record_path, record, exclusive=True)
            species_started = time.monotonic()
            try:
                atoms, digest = _paired_input(dft, output, name, basis, expected_dft_guess=dft_settings.scf_initial_guess)
                record.update(geometry_file=f"dft/inputs/{spec['file']}", dft_xyz_sha256=digest, cc_xyz_sha256=digest)
                _write_json(record_path, record)
                settings = CCSettings(charge=spec["charge"], spin=spec["spin"], basis=basis,
                                      scf_initial_guess=cc_initial_guess)
                expected_geometry = {"symbols": atoms.get_chemical_symbols(), "positions_angstrom": atoms.positions.tolist()}
                quantum = coupled_cluster_energy(atoms, settings)
                if quantum.get("status") != "completed" or not all(quantum.get(flag) is True for flag in ("scf_converged", "ccsd_converged", "triples_completed")):
                    raise MethodComparisonError(f"Incomplete coupled-cluster calculation for {name}")
                if quantum.get("settings") != asdict(settings):
                    raise MethodComparisonError(f"Coupled-cluster settings mismatch for {name}")
                if quantum.get("geometry") != expected_geometry:
                    raise MethodComparisonError(f"Coupled-cluster geometry mismatch for {name}")
                reference = "UHF" if spec["spin"] else "RHF"
                variant = "UCCSD(T)" if spec["spin"] else "RCCSD(T)"
                if quantum.get("reference") != reference or quantum.get("variant") != variant:
                    raise MethodComparisonError(f"Unexpected coupled-cluster reference/variant for {name}")
                _check_spin_diagnostics(quantum, spec["spin"], "hf_", f"CC reference {name}")
                energy = _finite_number(quantum.get("total_energy_ev"), f"coupled-cluster energy for {name}")
                hartree = _finite_number(quantum.get("total_energy_hartree"), f"coupled-cluster Hartree energy for {name}")
                if not math.isclose(energy, hartree * Hartree, abs_tol=1e-8, rel_tol=0):
                    raise MethodComparisonError(f"Coupled-cluster energy unit mismatch for {name}")
                try:
                    json.dumps(quantum, allow_nan=False)
                except (TypeError, ValueError) as error:
                    raise MethodComparisonError(f"Nonserializable coupled-cluster evidence for {name}") from error
                dft_record = dft["species"][name]
                difference = _finite_number(dft_record["energy_ev"] - energy, f"paired energy difference for {name}")
                record.update(status="completed", energy_ev=energy, energy_hartree=hartree, quantum_diagnostics=quantum)
                result["species"][name] = {
                    "geometry_file": record["geometry_file"], "dft_xyz_sha256": digest, "cc_xyz_sha256": digest,
                    "basis": basis, "charge": spec["charge"], "spin_2s": spec["spin"],
                    "cc_scf_initial_guess": cc_initial_guess,
                    "nominal_charge_spin_inputs_matched": True, "electronic_state_identity_verified": False,
                    "dft_energy_ev": dft_record["energy_ev"], "cc_energy_ev": energy,
                    "dft_minus_cc_energy_ev": difference,
                    "dft_s2": dft_record.get("s2"), "dft_expected_s2": dft_record.get("expected_s2"),
                    "dft_s2_deviation": dft_record.get("spin_contamination"),
                    "cc_hf_s2": quantum.get("hf_s2"), "cc_hf_expected_s2": quantum.get("hf_expected_s2"),
                    "cc_hf_s2_deviation": quantum.get("hf_s2_deviation"),
                    "dft_spin_diagnostic_scope": "Kohn-Sham determinant; does not establish electronic-state identity",
                    "cc_spin_diagnostic_scope": "Hartree-Fock reference determinant, not correlated CC wavefunction",
                    "cc_reference": quantum["reference"], "cc_variant": quantum["variant"],
                    "dft_record": f"dft/species/{name}.json", "cc_record": f"ccsd_t/{name}.json",
                }
            except BaseException as error:
                record.update(status="interrupted" if isinstance(error, (KeyboardInterrupt, SystemExit)) else "failed",
                              error_type=type(error).__name__, error=str(error))
                if hasattr(error, "diagnostics"):
                    record["quantum_diagnostics"] = error.diagnostics
                raise
            finally:
                record["elapsed_seconds"] = time.monotonic() - species_started
                _write_json(record_path, record)
            _write_json(summary_path, result)
        # A later accidental change to archived inputs invalidates pairing.
        for name in SPECIES:
            _paired_input(dft, output, name, basis, expected_dft_guess=dft_settings.scf_initial_guess)
        dft_relative = _energy_differences({name: entry["dft_energy_ev"] for name, entry in result["species"].items()})
        cc_relative = _energy_differences({name: entry["cc_energy_ev"] for name, entry in result["species"].items()})
        differences = {key: _finite_number(dft_relative[key] - cc_relative[key], f"aggregate method difference {key}") for key in dft_relative}
        result["computed"] = {
            "dft": dft_relative, "ccsd_t": cc_relative,
            "dft_minus_ccsd_t": differences,
            "definition": "Nominal TS-relative = E(source transition structure)-E(CH4)-E(C2H); reaction = E(C2H2)+E(CH3)-E(CH4)-E(C2H).",
        }
        result.update(status="completed", stage="completed", geometry_basis_pairing_verified=True,
                      completed_utc=datetime.now(timezone.utc).isoformat())
        result.pop("active_species", None)
    except BaseException as error:
        result.update(status="interrupted" if isinstance(error, (KeyboardInterrupt, SystemExit)) else "failed",
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        result["elapsed_seconds"] = time.monotonic() - started
        _write_json(summary_path, result)
    return result
