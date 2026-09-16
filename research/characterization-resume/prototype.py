"""Research prototype: durable force checkpoints and read-only Hessian replay.

Not a production CLI or a validated molecular model. The acquisition backend
must be the exact supported PySCFCalculator, whose evaluations independently
construct a new electronic calculation. Arbitrary stateful calculators cannot
be resumed safely from a class name alone. Tests replace this type explicitly.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import math
import os
from pathlib import Path
import platform
import tempfile
import time
from uuid import uuid4

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
from ase.io import write

from nanodesign import quantum, stationary
from nanodesign.quantum import PySCFCalculator, QuantumSettings


class CheckpointError(ValueError):
    """Checkpoint cannot be reused without contradicting its recorded inputs."""


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _digest(raw):
    return hashlib.sha256(raw).hexdigest()


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _atomic_json(path, value):
    encoded = _canonical(value) + b"\n"
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".checkpoint-", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        # Preserve the rename on filesystems supporting directory fsync.
        descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _save_checkpoint(root, payload):
    _atomic_json(root / "force_checkpoint.json", {"payload": payload, "sha256": _digest(_canonical(payload))})


@contextmanager
def _lock(root, *, create):
    import fcntl
    path = root / ".checkpoint.lock"
    if path.is_symlink():
        raise CheckpointError("Checkpoint lock must not be a symlink")
    with path.open("xb" if create else "rb") as stream:
        operation = fcntl.LOCK_EX if create else fcntl.LOCK_SH
        try:
            fcntl.flock(stream, operation | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise CheckpointError("Source checkpoint has an active worker") from error
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def _reject_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise CheckpointError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _parse(raw):
    try:
        return json.loads(raw, object_pairs_hook=_reject_duplicates,
                          parse_constant=lambda value: (_ for _ in ()).throw(CheckpointError(f"Nonfinite JSON: {value}")))
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise CheckpointError("Unreadable checkpoint; original bytes preserved") from error


def _versions():
    packages = {}
    for name in ("numpy", "ase", "pyscf", "dftd3", "scipy"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    return {"python": platform.python_version(), "packages": packages}


def _geometry(atoms):
    if not isinstance(atoms, Atoms) or not len(atoms):
        raise CheckpointError("A nonempty ASE geometry is required")
    if atoms.pbc.any() or not np.isfinite(atoms.positions).all():
        raise CheckpointError("Geometry must be finite and nonperiodic")
    masses = atoms.get_masses()
    if not np.isfinite(masses).all() or np.any(masses <= 0):
        raise CheckpointError("Every atomic mass must be positive and finite")
    fixed = set()
    for constraint in atoms.constraints:
        if not isinstance(constraint, FixAtoms):
            raise CheckpointError("Only whole-atom FixAtoms constraints are supported")
        indices = constraint.get_indices().tolist()
        if any(type(index) is not int or not 0 <= index < len(atoms) for index in indices):
            raise CheckpointError("Invalid fixed atom index")
        fixed.update(indices)
    geometry = {
        "atomic_numbers": atoms.numbers.tolist(), "positions_angstrom": atoms.positions.tolist(),
        "masses_amu": masses.tolist(), "fixed_indices": sorted(fixed),
        "cell_angstrom": atoms.cell.tolist(), "pbc": atoms.pbc.tolist(),
        "initial_charges": atoms.get_initial_charges().tolist(),
        "initial_magnetic_moments": atoms.get_initial_magnetic_moments().tolist(),
    }
    _canonical(geometry)
    return geometry


def _manifest(atoms, settings, input_context, step, fmax, frequency_tolerance, max_free_coordinates):
    if type(atoms.calc) is not PySCFCalculator:
        raise CheckpointError("Automatic resume supports the exact PySCFCalculator only")
    if not isinstance(settings, QuantumSettings) or getattr(atoms.calc, "settings", None) != settings:
        raise CheckpointError("Effective quantum settings differ from the attached calculator")
    if not isinstance(input_context, dict) or not input_context:
        raise CheckpointError("Provide explicit nonempty input_context provenance")
    context = deepcopy(input_context)
    if "input_hashes" in context and not isinstance(context["input_hashes"], dict):
        raise CheckpointError("input_context.input_hashes must be an object")
    _canonical(context)
    geometry = _geometry(atoms)
    step = stationary._positive_real("step", step)
    fmax = stationary._positive_real("fmax", fmax)
    frequency_tolerance = stationary._positive_real("frequency_tolerance", frequency_tolerance)
    if type(max_free_coordinates) is not int or max_free_coordinates < 1:
        raise CheckpointError("max_free_coordinates must be a positive integer")
    free = [index for index in range(len(atoms)) if index not in geometry["fixed_indices"]]
    if not free or 3 * len(free) > max_free_coordinates:
        raise CheckpointError("Free-coordinate count exceeds the bound or is zero")
    coordinates = (3 * np.asarray(free)[:, None] + np.arange(3)).ravel()
    plan = stationary._displacement_plan(atoms, coordinates, step)
    requests = [{"id": "baseline", "kind": "baseline"}]
    for index, point in enumerate(plan):
        requests.append({"id": f"displacement-{index:04d}", "kind": "displacement", **point})
    return {
        "schema_version": 1, "algorithm": "central-force-difference-checkpoint-v1",
        "geometry": geometry, "quantum_settings": asdict(settings), "input_context": context,
        "characterization_settings": {"step_angstrom": step, "force_tolerance_ev_per_angstrom": fmax,
            "frequency_tolerance_cm1": frequency_tolerance, "max_free_coordinates": max_free_coordinates},
        "backend": {"module": type(atoms.calc).__module__, "class": type(atoms.calc).__qualname__},
        "software": _versions(),
        "implementation_sha256": {name: _digest(Path(module.__file__).read_bytes())
            for name, module in (("stationary", stationary), ("quantum", quantum))},
        "prototype_sha256": _digest(Path(__file__).read_bytes()),
        "units": {"length": "angstrom", "force": "eV/angstrom", "mass": "amu"},
        "requests": requests,
    }


def _forces(value, count):
    if (not isinstance(value, list) or len(value) != count
            or any(not isinstance(row, list) or len(row) != 3 for row in value)
            or any(type(component) not in (int, float) or not math.isfinite(component)
                   for row in value for component in row)):
        raise CheckpointError("Force evidence must contain finite numeric N by 3 values")
    return np.asarray(value, dtype=float)


def _validate_record(record, request, manifest):
    if not isinstance(record, dict) or record.get("status") != "completed":
        raise CheckpointError("Failed or incomplete calculation is not reusable force evidence")
    if _canonical(record.get("request")) != _canonical(request):
        raise CheckpointError("Checkpoint displacement/axis/request mismatch")
    _forces(record.get("forces_ev_per_angstrom"), len(manifest["geometry"]["atomic_numbers"]))
    diagnostics = record.get("quantum_diagnostics")
    if (not isinstance(diagnostics, dict) or diagnostics.get("scf_converged") is not True
            or diagnostics.get("gradient_completed") is not True
            or _canonical(diagnostics.get("settings")) != _canonical(manifest["quantum_settings"])):
        raise CheckpointError("Force record lacks consistent successful quantum diagnostics")
    identifier = diagnostics.get("call_id")
    expected_id = identifier if isinstance(identifier, str) else None
    if record.get("calculation_call_id") != expected_id:
        raise CheckpointError("Force record call ID disagrees with its diagnostics")
    if not isinstance(record.get("origin_attempt_id"), str) or not record["origin_attempt_id"]:
        raise CheckpointError("Force record has no originating attempt")


def _validate_manifest_structure(manifest):
    """Check a manifest's internal geometry/plan binding without a live solver."""
    try:
        if manifest.get("algorithm") != "central-force-difference-checkpoint-v1":
            raise CheckpointError("Unknown checkpoint algorithm")
        geometry = manifest["geometry"]
        atoms = Atoms(numbers=geometry["atomic_numbers"], positions=geometry["positions_angstrom"],
                      masses=geometry["masses_amu"], cell=geometry["cell_angstrom"], pbc=geometry["pbc"])
        fixed = geometry["fixed_indices"]
        if not isinstance(fixed, list) or any(type(index) is not int for index in fixed):
            raise CheckpointError("Invalid fixed-index types")
        atoms.set_constraint(FixAtoms(indices=fixed))
        atoms.set_initial_charges(geometry["initial_charges"])
        atoms.set_initial_magnetic_moments(geometry["initial_magnetic_moments"])
        if _canonical(_geometry(atoms)) != _canonical(geometry):
            raise CheckpointError("Geometry is not in the declared canonical form")
        settings = QuantumSettings(**manifest["quantum_settings"])
        if _canonical(asdict(settings)) != _canonical(manifest["quantum_settings"]):
            raise CheckpointError("Incomplete or noncanonical effective settings")
        controls = manifest["characterization_settings"]
        step = stationary._positive_real("step", controls["step_angstrom"])
        stationary._positive_real("force tolerance", controls["force_tolerance_ev_per_angstrom"])
        stationary._positive_real("frequency tolerance", controls["frequency_tolerance_cm1"])
        bound = controls["max_free_coordinates"]
        free = [index for index in range(len(atoms)) if index not in fixed]
        if type(bound) is not int or not 0 < 3 * len(free) <= bound:
            raise CheckpointError("Invalid free-coordinate budget")
        coordinates = (3 * np.asarray(free)[:, None] + np.arange(3)).ravel()
        plan = [{"id": "baseline", "kind": "baseline"}]
        plan += [{"id": f"displacement-{index:04d}", "kind": "displacement", **point}
                 for index, point in enumerate(stationary._displacement_plan(atoms, coordinates, step))]
        if _canonical(plan) != _canonical(manifest["requests"]):
            raise CheckpointError("Requests do not match the recorded geometry/step/free space")
        if manifest["units"] != {"length": "angstrom", "force": "eV/angstrom", "mass": "amu"}:
            raise CheckpointError("Unsupported force-evidence units")
    except (KeyError, TypeError, IndexError, ValueError) as error:
        if isinstance(error, CheckpointError):
            raise
        raise CheckpointError(f"Malformed checkpoint manifest: {error}") from error


def _validate_payload(payload, expected_manifest=None):
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise CheckpointError("Unsupported checkpoint schema")
    manifest = payload.get("manifest")
    if not isinstance(manifest, dict) or manifest.get("schema_version") != 1:
        raise CheckpointError("Missing checkpoint manifest")
    if expected_manifest is not None and _canonical(manifest) != _canonical(expected_manifest):
        raise CheckpointError("Checkpoint manifest mismatch: geometry, boundary, settings, controls or provenance changed")
    if payload.get("manifest_sha256") != _digest(_canonical(manifest)):
        raise CheckpointError("Manifest fingerprint mismatch")
    _validate_manifest_structure(manifest)
    if payload.get("status") not in {"running", "completed", "failed", "interrupted"}:
        raise CheckpointError("Unknown checkpoint status")
    if not isinstance(payload.get("attempt_id"), str) or not payload["attempt_id"]:
        raise CheckpointError("Missing checkpoint attempt identity")
    requests = manifest.get("requests")
    records = payload.get("records")
    if not isinstance(requests, list) or not requests or not isinstance(records, list) or len(records) > len(requests):
        raise CheckpointError("Invalid checkpoint record count")
    if any(not isinstance(point, dict) or not isinstance(point.get("id"), str) for point in requests):
        raise CheckpointError("Malformed checkpoint request identity")
    if len({point["id"] for point in requests}) != len(requests):
        raise CheckpointError("Duplicate checkpoint request identity")
    for record, request in zip(records, requests):
        _validate_record(record, request, manifest)
    if payload["status"] == "completed" and len(records) != len(requests):
        raise CheckpointError("Completed checkpoint lacks required force records")
    return payload


def load_checkpoint(root, *, expected_manifest=None):
    """Read a stable snapshot under a shared lock; never modify source bytes."""
    root = Path(root).resolve()
    with _lock(root, create=False):
        path = root / "force_checkpoint.json"
        if path.is_symlink():
            raise CheckpointError("Checkpoint file must not be a symlink")
        raw = path.read_bytes()
    document = _parse(raw)
    if not isinstance(document, dict) or set(document) != {"payload", "sha256"}:
        raise CheckpointError("Invalid checkpoint envelope")
    if document["sha256"] != _digest(_canonical(document["payload"])):
        raise CheckpointError("Checkpoint digest mismatch; source preserved")
    payload = _validate_payload(document["payload"], expected_manifest)
    return deepcopy(payload), _digest(raw)


def _point_atoms(atoms, request):
    point = atoms.copy()
    if request["kind"] == "displacement":
        point.positions[request["atom_index"], request["axis"]] = request["displaced_coordinate_angstrom"]
    return point


@contextmanager
def _preserve_backend(calculator):
    saved_atoms = None if calculator.atoms is None else calculator.atoms.copy()
    saved_results = deepcopy(calculator.results)
    had_diagnostics = hasattr(calculator, "diagnostics")
    saved_diagnostics = deepcopy(getattr(calculator, "diagnostics", None))
    try:
        yield
    finally:
        calculator.atoms = saved_atoms
        calculator.results = saved_results
        if had_diagnostics:
            calculator.diagnostics = saved_diagnostics
        elif hasattr(calculator, "diagnostics"):
            del calculator.diagnostics


class _ReplayForces(Calculator):
    """Analysis adapter over validated force evidence; never calls a solver."""
    implemented_properties = ["forces"]

    def __init__(self, reference, payload):
        super().__init__()
        self.reference = reference
        self.payload = payload
        self.diagnostics = {}

    def calculate(self, atoms=None, properties=("forces",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        for record in self.payload["records"]:
            expected = _point_atoms(self.reference, record["request"])
            if _canonical(_geometry(atoms)) == _canonical(_geometry(expected)):
                self.results = {"forces": np.asarray(record["forces_ev_per_angstrom"], dtype=float)}
                self.diagnostics = deepcopy(record["quantum_diagnostics"])
                return
        raise CheckpointError("Replay requested a geometry absent from validated force evidence")


def run_characterization_checkpoint(atoms, output, *, settings, input_context,
                                    resume_from=None, step=0.005, fmax=0.03,
                                    frequency_tolerance=20.0, max_free_coordinates=120):
    """Acquire missing forces, save every accepted record, then replay the math.

    ``input_context`` is explicit caller-supplied source provenance. Production
    integration must provide hashes from the same byte snapshots it parsed.
    Reuse is opt-in and always writes a new attempt. Matching metadata and
    hashes are provenance checks, not proof of electronic-state identity.
    """
    if not isinstance(atoms, Atoms):
        raise CheckpointError("An ASE geometry is required")
    # Pin all in-memory inputs before calculations, just as the workflow pins
    # coordinate-file bytes before parsing. The caller's Atoms is not modified.
    reference = atoms.copy()
    reference.calc = atoms.calc
    atoms = reference
    manifest = _manifest(atoms, settings, input_context, step, fmax, frequency_tolerance, max_free_coordinates)
    started = time.monotonic()
    source = None
    records = []
    if resume_from is not None:
        prior, raw_hash = load_checkpoint(resume_from, expected_manifest=manifest)
        records = deepcopy(prior["records"])
        source = {"path": str(Path(resume_from).resolve()), "checkpoint_sha256": raw_hash,
                  "attempt_id": prior["attempt_id"], "status": prior["status"],
                  "scope": "Original record/call IDs retained; source attempt is unchanged"}
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    attempt_id = str(uuid4())
    payload = {"schema_version": 1, "attempt_id": attempt_id, "started_utc": _utc(),
               "manifest": manifest, "manifest_sha256": _digest(_canonical(manifest)),
               "records": records, "status": "running", "resume_source": source}
    original_records = len(records)
    new_evaluations = 0
    current_request = None
    result = {"schema_version": 1, "stage": "characterize", "status": "running",
              "started_utc": payload["started_utc"],
              "units": {"length": "angstrom", "energy": "eV", "force": "eV/angstrom", "frequency": "cm^-1"},
              "quantum_settings": asdict(settings), "input_hashes": deepcopy(input_context.get("input_hashes", {})),
              "input_context": deepcopy(input_context), "checkpoint_manifest": manifest,
              "characterization_settings": manifest["characterization_settings"],
              "validation": {"design_validated": False, "transition_state_validated": False,
                  "electronic_state_identity_verified": False},
              "prototype_scope": "Research checkpoint/replay prototype; no production CLI integration"}
    if "design" in input_context:
        result["design"] = deepcopy(input_context["design"])
    with _lock(root, create=True):
        try:
            _save_checkpoint(root, payload)
            _atomic_json(root / "result.json", result)
            snapshot = atoms.copy()
            snapshot.calc = None
            write(root / "input.extxyz", snapshot)
            result["input_hashes"]["input_snapshot_sha256"] = _digest((root / "input.extxyz").read_bytes())
            result["input_snapshot_hash_scope"] = "input_hashes.input_snapshot_sha256 binds this attempt's input.extxyz bytes; structure_sha256 retains the caller's original source-file hash"
            backend = atoms.calc
            with _preserve_backend(backend):
                for index, request in enumerate(manifest["requests"]):
                    current_request = request
                    if index >= len(payload["records"]):
                        # Supported PySCF builds a fresh mean field per evaluation.
                        # Reset only the backend's local ASE cache; never inject
                        # recorded forces into a live electronic calculator.
                        backend.reset()
                        point = _point_atoms(atoms, request)
                        point.calc = backend
                        new_evaluations += 1
                        try:
                            forces = point.get_forces(apply_constraint=False)
                        except BaseException as error:
                            payload["failed_calculation"] = {
                                "request": deepcopy(request),
                                "quantum_diagnostics": deepcopy(backend.diagnostics),
                                "error": {"type": type(error).__name__, "message": str(error)},
                                "reusable": False,
                            }
                            raise
                        diagnostics = deepcopy(backend.diagnostics)
                        identifier = diagnostics.get("call_id")
                        record = {"status": "completed", "request": deepcopy(request),
                                  "forces_ev_per_angstrom": np.asarray(forces).tolist(),
                                  "quantum_diagnostics": diagnostics,
                                  "calculation_call_id": identifier if isinstance(identifier, str) else None,
                                  "origin_attempt_id": attempt_id}
                        _validate_record(record, request, manifest)
                        payload["records"].append(record)
                        _save_checkpoint(root, payload)
                    if index == 0:
                        baseline = payload["records"][0]
                        free = [i for i in range(len(atoms)) if i not in manifest["geometry"]["fixed_indices"]]
                        residual = float(np.linalg.norm(np.asarray(baseline["forces_ev_per_angstrom"])[free], axis=1).max())
                        result["initial_free_force_max_ev_per_angstrom"] = residual
                        result["initial_force_guard_passed"] = residual <= fmax
                        result["quantum_diagnostics"] = deepcopy(baseline["quantum_diagnostics"])
                        if residual > fmax:
                            raise CheckpointError("Baseline is not stationary; saved full forces, no displaced work")
            _validate_payload(payload, manifest)
            replay = atoms.copy()
            replay.calc = _ReplayForces(atoms.copy(), payload)
            result["stationary"] = stationary.characterize_stationary_point(
                replay, step_angstrom=step, force_tolerance_ev_per_angstrom=fmax,
                frequency_tolerance_cm1=frequency_tolerance, max_free_coordinates=max_free_coordinates)
            result["status"] = payload["status"] = "completed"
        except BaseException as error:
            status = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
            result["status"] = payload["status"] = status
            failure = {"type": type(error).__name__, "message": str(error), "request": current_request}
            result["error"] = payload["error"] = failure
            raise
        finally:
            result["elapsed_seconds"] = time.monotonic() - started
            result["checkpoint"] = {"attempt_id": attempt_id, "manifest_sha256": payload["manifest_sha256"],
                "source": source, "reused_force_records": original_records,
                "new_force_evaluations": new_evaluations, "completed_force_records": len(payload["records"]),
                "planned_force_records": len(manifest["requests"]),
                "analysis_scope": "Hessian analysis replays saved forces; replay creates no electronic calculations or log events"}
            payload["updated_utc"] = _utc()
            _save_checkpoint(root, payload)
            _atomic_json(root / "result.json", result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description="Inspect a resumable-force research checkpoint without running a solver")
    parser.add_argument("checkpoint", type=Path, help="Attempt directory containing force_checkpoint.json")
    args = parser.parse_args(argv)
    payload, digest = load_checkpoint(args.checkpoint)
    print(json.dumps({"status": payload["status"], "checkpoint_sha256": digest,
                      "completed_force_records": len(payload["records"]),
                      "planned_force_records": len(payload["manifest"]["requests"]),
                      "identity_binding": "Not compared against new calculation inputs by this inspection command",
                      "scientific_validation": False}, indent=2))


if __name__ == "__main__":
    main()
