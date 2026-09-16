"""Auditable fixed-geometry comparison for methane/ethynyl H abstraction.

Single-point energies on published UCCSD(T)/cc-pVDZ structures do not establish
this DFT method's saddle, barrier, chemical accuracy, or tool performance. The
comparison to the main paper combines method, basis, and geometry differences.
The SI's inconsistent absolute energies are never used as calibration targets.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, replace
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import platform
import re
import time
from typing import Any

import numpy as np
from ase import Atoms
from ase.io import read
from ase.units import Hartree

from . import __version__
from .design import sha256
from .quantum import PySCFCalculator, QuantumSettings

EV_PER_KCAL_PER_MOL = 0.0433641153087705
# Reference data have a separate CC BY-NC license and are not package data.
# A wheel installation must pass the reference directory from a source checkout.
REFERENCE_DIRECTORY = Path(__file__).resolve().parent.parent / "data" / "reference"
SOURCE_SHA256 = "ab1c05fb9add21e85d5e97e25f105c2abc59ffe3106b9fe7b181f892d951ff6e"
ANGSTROM_PER_BOHR = 0.529177210544

# spin = N_alpha - N_beta; charge, formula, source label and atom order are checked.
SPECIES: dict[str, dict[str, Any]] = {
    "methane": {"file": "methane.xyz", "charge": 0, "spin": 0, "formula": "CH4", "source_label": "CH4"},
    "ethynyl_radical": {"file": "ethynyl_radical.xyz", "charge": 0, "spin": 1, "formula": "C2H", "source_label": "CCH"},
    "methane_ethynyl_ts": {"file": "methane_ethynyl_ts.xyz", "charge": 0, "spin": 1, "formula": "C3H5", "source_label": "CH3-H-CCH"},
    "acetylene": {"file": "acetylene.xyz", "charge": 0, "spin": 0, "formula": "C2H2", "source_label": "HCCH"},
    "methyl_radical": {"file": "methyl_radical.xyz", "charge": 0, "spin": 1, "formula": "CH3", "source_label": "CH3"},
}

PUBLISHED_BARRIER = {
    "quantity": "bare electronic barrier relative to separated reactants",
    "reaction": "C2H + CH4 -> C2H2 + CH3 (neutral doublet reaction)",
    "level_of_theory": "RCCSD(T)/cc-pVTZ",
    "value_kcal_per_mol": 2.2,
    "source_doi": "10.1021/jp061821e",
    "source_url": "https://vergil.chemistry.gatech.edu/static/pdfs/temelso_2006_11160.pdf",
    "source_location": "Table 5, methane reaction, Delta E dagger, CCSD(T)/TZ",
    "geometry_caveat": (
        "The supplied geometries are UCCSD(T)/cc-pVDZ; they do not reproduce "
        "the RCCSD(T)/cc-pVTZ stationary-structure calculation. This comparison "
        "combines method, basis, and geometry differences, not an isolated functional error."
    ),
    "reference_structure_caveat": (
        "Table 1 reports three imaginary modes (259i, 50i, 50i cm^-1) for the "
        "nominal collinear UCCSD(T)/cc-pVDZ methane transition structure, not a "
        "first-order saddle. The supplied SI geometry has not been independently reverified."
    ),
    "source_geometry_discrepancy": (
        "Main-paper Table 2 donor-H/acceptor-H distances are 1.148/1.678 angstrom "
        "at UCCSD(T)/cc-pVDZ; the supplied SI converts to 1.149220/1.672376 angstrom. "
        "The geometry inconsistency remains unresolved."
    ),
    "excluded_quantity": (
        "The 1.7 kcal/mol 0 K enthalpic barrier uses RCCSD(T)/cc-pVDZ vibrational "
        "corrections (Table 5 footnote d) and is not the electronic-energy target."
    ),
}
SOURCE_ENERGY_CAVEAT = (
    "The SI's labeled absolute methane, ethynyl and TS energies imply about "
    "-14.76 kcal/mol, inconsistent with the main-paper positive barrier. They are "
    "unverified and unused as energy targets; see SOURCE_NOTES.md in the source checkout."
)


class BenchmarkError(RuntimeError):
    """The specified comparison could not be completed; inspect saved evidence."""


def _write_json(path: Path, result: dict[str, Any], *, exclusive: bool = False) -> None:
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if exclusive:
        with path.open("x", encoding="utf-8") as stream:
            stream.write(encoded)
    else:
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(encoded, encoding="utf-8")
        temporary.replace(path)


def _source_geometry(source: str, label: str) -> tuple[list[str], np.ndarray]:
    """Extract only atom coordinates, never the unverified source energies."""
    pattern = rf"(?m)^\s*{re.escape(label)}\s+-\s*UCCSD\(T\)/cc-pVDZ[^\n]*\n((?:[CH]\s+[^\n]+\n)+)"
    matches = re.findall(pattern, source)
    if len(matches) != 1:
        raise BenchmarkError(f"Expected exactly one SI geometry for {label}")
    rows = [line.split() for line in matches[0].strip().splitlines()]
    return [row[0] for row in rows], np.asarray([[float(v) for v in row[1:]] for row in rows]) * ANGSTROM_PER_BOHR


def _snapshot_and_validate(directory: Path, output: Path) -> tuple[dict[str, Atoms], dict[str, Any]]:
    snapshot = output / "inputs"
    snapshot.mkdir()
    files = ["provenance.json", "temelso_2006_si.txt", "README.md"] + [spec["file"] for spec in SPECIES.values()]
    hashes = {}
    for filename in files:
        original = directory / filename
        if not original.is_file():
            raise BenchmarkError(
                f"Missing reference input: {original}. Reference data are distributed in the "
                "source checkout; installed packages must supply its data/reference directory."
            )
        # Hash and calculate from the same archived bytes, not a mutable original.
        with (snapshot / filename).open("xb") as stream:
            stream.write(original.read_bytes())
        hashes[filename] = sha256(snapshot / filename)
    provenance = json.loads((snapshot / "provenance.json").read_text())
    if hashes["temelso_2006_si.txt"] != SOURCE_SHA256 or provenance.get("source_sha256") != SOURCE_SHA256:
        raise BenchmarkError("Published SI source hash mismatch")
    if provenance.get("article_doi") != PUBLISHED_BARRIER["source_doi"]:
        raise BenchmarkError("Unexpected reference article DOI")
    if provenance.get("license") != "CC BY-NC 4.0" or not provenance.get("authors"):
        raise BenchmarkError("Reference attribution or separate data license is missing")
    conversion = provenance.get("conversion", {})
    if (conversion.get("source_length_unit") != "bohr" or conversion.get("output_length_unit") != "angstrom"
            or conversion.get("angstrom_per_bohr") != ANGSTROM_PER_BOHR):
        raise BenchmarkError("Reference length-unit conversion does not match the published input conversion")
    source = (snapshot / "temelso_2006_si.txt").read_text()
    geometries = {}
    for name, spec in SPECIES.items():
        entries = [entry for entry in provenance.get("structures", []) if entry.get("file") == spec["file"]]
        if len(entries) != 1:
            raise BenchmarkError(f"Expected one provenance record for {name}")
        metadata = entries[0]
        if any(type(metadata.get(field)) is not int for field in ("charge", "multiplicity", "atoms")):
            raise BenchmarkError(f"Reference charge, multiplicity and atom count must be integers for {name}")
        if (metadata.get("charge") != spec["charge"] or metadata.get("multiplicity") != spec["spin"] + 1
                or metadata.get("source_label") != spec["source_label"]):
            raise BenchmarkError(f"Electronic state or source identity mismatch for {name}")
        frames = read(snapshot / spec["file"], index=":", format="xyz")
        if len(frames) != 1 or len(frames[0]) == 0:
            raise BenchmarkError(f"Expected one nonempty geometry for {name}")
        atoms = frames[0]
        if atoms.get_chemical_formula() != spec["formula"] or len(atoms) != metadata.get("atoms"):
            raise BenchmarkError(f"Species formula or atom count mismatch for {name}")
        symbols, positions = _source_geometry(source, spec["source_label"])
        if atoms.get_chemical_symbols() != symbols:
            raise BenchmarkError(f"Source atom ordering mismatch for {name}")
        if not np.all(np.isfinite(atoms.positions)) or not np.allclose(atoms.positions, positions, atol=5e-11, rtol=0):
            raise BenchmarkError(f"Coordinates for {name} do not match the published SI conversion")
        if np.any(atoms.pbc):
            raise BenchmarkError(f"Periodic reference geometry unsupported for {name}")
        geometries[name] = atoms
    reactants = ("methane", "ethynyl_radical")
    products = ("acetylene", "methyl_radical")
    transition = ("methane_ethynyl_ts",)
    for names in (reactants, products, transition):
        counts = sum((Counter(geometries[name].get_chemical_symbols()) for name in names), Counter())
        charge = sum(SPECIES[name]["charge"] for name in names)
        spin = sum(SPECIES[name]["spin"] for name in names)
        if counts != Counter({"C": 3, "H": 5}) or charge != 0 or spin != 1:
            raise BenchmarkError("Reaction does not conserve expected elements, charge and doublet state")
    manifest = {
        "source_directory": str(directory.resolve()), "archived_inputs": "inputs",
        "file_sha256": hashes, "provenance": provenance,
        "source_energy_caveat": SOURCE_ENERGY_CAVEAT,
        "geometry_validation": "Exact atom sequence and converted SI coordinates checked within 5e-11 angstrom",
        "reaction_conservation": {"elements": {"C": 3, "H": 5}, "charge": 0, "electrons": 23, "spin_2s": 1},
        "data_license": "CC BY-NC 4.0; independent of repository code license",
    }
    _write_json(output / "input_manifest.json", manifest, exclusive=True)
    return geometries, manifest


def _species_energy(name: str, atoms: Atoms, settings: QuantumSettings, output: Path) -> dict[str, Any]:
    spec = SPECIES[name]
    species_settings = replace(settings, charge=spec["charge"], spin=spec["spin"])
    event_path = output / "electronic" / f"{name}.jsonl"
    calculator = PySCFCalculator(species_settings, event_log=event_path)
    atoms = atoms.copy()
    atoms.calc = calculator
    record: dict[str, Any] = {
        "species": name, "file": f"inputs/{spec['file']}",
        "file_sha256": sha256(output / "inputs" / spec["file"]),
        "formula": atoms.get_chemical_formula(), "atoms": len(atoms),
        "charge": spec["charge"], "spin_2s": spec["spin"],
        "quantum_settings": asdict(species_settings),
        "event_log": str(event_path.relative_to(output)),
        "status": "running",
    }
    started = time.monotonic()
    try:
        energy_ev = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces())
        diagnostics = calculator.diagnostics
        if diagnostics.get("scf_converged") is not True or diagnostics.get("gradient_completed") is not True:
            raise BenchmarkError(f"Incomplete or unconverged electronic calculation for {name}")
        if not math.isfinite(energy_ev) or forces.shape != (len(atoms), 3) or not np.all(np.isfinite(forces)):
            raise BenchmarkError(f"Nonfinite energy or invalid forces for {name}")
        hartree = diagnostics.get("total_energy_hartree")
        if hartree is None or not math.isfinite(hartree) or not math.isclose(hartree * Hartree, energy_ev, abs_tol=1e-8, rel_tol=0):
            raise BenchmarkError(f"Inconsistent energy units for {name}")
        record.update(
            status="completed", energy_ev=energy_ev, energy_hartree=hartree,
            scf_converged=True, max_force_ev_per_angstrom=float(np.max(np.linalg.norm(forces, axis=1))),
            forces_ev_per_angstrom=forces.tolist(), s2=diagnostics.get("s2"),
            expected_s2=diagnostics.get("expected_s2"), spin_contamination=diagnostics.get("s2_deviation"),
        )
    except BaseException as error:
        record.update(status="interrupted" if isinstance(error, (KeyboardInterrupt, SystemExit)) else "failed",
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        record.update(elapsed_seconds=time.monotonic() - started, quantum_diagnostics=calculator.diagnostics)
        _write_json(output / "species" / f"{name}.json", record, exclusive=True)
    return record


def run_benchmark(
    settings: QuantumSettings | None = None,
    directory: Path | None = None,
    output: Path | None = None,
) -> dict[str, Any]:
    """Evaluate five published structures and save incremental, auditable evidence.

    ``output`` is required and must not exist. Inputs are copied and checked
    before expensive work; each completed/failed species is saved independently.
    ``benchmark.json`` records running, completed, failed or interrupted status.
    A hard process kill can leave ``running`` with completed species records;
    it cannot leave a falsely completed comparison. No automatic resume occurs.
    """
    if output is None:
        raise ValueError("output is required to preserve benchmark inputs and incremental evidence")
    settings = settings if settings is not None else QuantumSettings()
    if not isinstance(settings, QuantumSettings):
        raise TypeError("settings must be a QuantumSettings instance")
    directory = Path(directory) if directory is not None else REFERENCE_DIRECTORY
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "species").mkdir()
    (output / "electronic").mkdir()
    started = time.monotonic()
    result: dict[str, Any] = {
        "schema_version": 2,
        "benchmark": "temelso_2006_methane_ethynyl_fixed_geometry_comparison",
        "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
        "software_version": __version__, "python_version": platform.python_version(),
        "quantum_settings": asdict(settings),
        "state_policy": "Base charge/spin are replaced by explicit published states for each species.",
        "units": {"energy": "eV", "secondary_energy": "kcal/mol", "length": "angstrom"},
        "species": {}, "published_reference": PUBLISHED_BARRIER,
        "source_energy_caveat": SOURCE_ENERGY_CAVEAT,
        "method_validated": False, "saddle_verified": False, "tool_validated": False,
        "outstanding": [
            "Optimise stationary structures and verify saddle modes and reaction connectivity",
            "Converge basis and grid; compare functionals and consistent high-level calculations",
            "Verify electronic-state stability; inspect spin contamination and state switching",
            "Add zero-point/thermal effects and tunneling before rate comparisons",
            "Validate transfer to a mounted diamondoid tool and substrate, including competing chemistry",
        ],
    }
    _write_json(output / "benchmark.json", result)
    try:
        geometries, manifest = _snapshot_and_validate(directory, output)
        result["input_manifest"] = "input_manifest.json"
        result["input_hashes"] = manifest["file_sha256"]
        for name, atoms in geometries.items():
            result["active_species"] = name
            result["elapsed_seconds"] = time.monotonic() - started
            _write_json(output / "benchmark.json", result)
            result["species"][name] = _species_energy(name, atoms, settings, output)
            _write_json(output / "benchmark.json", result)
        species = result["species"]
        reactants = species["methane"]["energy_ev"] + species["ethynyl_radical"]["energy_ev"]
        ts_relative = species["methane_ethynyl_ts"]["energy_ev"] - reactants
        reaction = species["acetylene"]["energy_ev"] + species["methyl_radical"]["energy_ev"] - reactants
        discrepancy = ts_relative / EV_PER_KCAL_PER_MOL - PUBLISHED_BARRIER["value_kcal_per_mol"]
        if not all(math.isfinite(value) for value in (reactants, ts_relative, reaction, discrepancy)):
            raise BenchmarkError("Nonfinite aggregate energy; no comparison can be reported")
        result["computed"] = {
            "reference_geometry_energy_ev": ts_relative,
            "reference_geometry_energy_kcal_per_mol": ts_relative / EV_PER_KCAL_PER_MOL,
            "reaction_energy_ev": reaction,
            "reaction_energy_kcal_per_mol": reaction / EV_PER_KCAL_PER_MOL,
            "definition": "E(nominal reference TS geometry) - E(CH4) - E(C2H); reaction energy = E(C2H2) + E(CH3) - E(CH4) - E(C2H)",
            "scope": "Electronic energies at fixed reference geometries; no ZPE/thermal corrections, saddle search, or verified DFT barrier",
        }
        result["comparison"] = {
            "signed_discrepancy_kcal_per_mol": discrepancy,
            "absolute_discrepancy_kcal_per_mol": abs(discrepancy),
            "interpretation": (
                "Combined method, basis and geometry discrepancy against a rounded published bare barrier; "
                "not calibrated accuracy or an error bound. The source nominal collinear structure has "
                "three reported imaginary modes, and SI/main-paper geometry differences remain unresolved."
            ),
        }
        result.update(status="completed", completed_utc=datetime.now(timezone.utc).isoformat())
        result.pop("active_species", None)
    except BaseException as error:
        result.update(status="interrupted" if isinstance(error, (KeyboardInterrupt, SystemExit)) else "failed",
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        result["elapsed_seconds"] = time.monotonic() - started
        _write_json(output / "benchmark.json", result)
    return result
