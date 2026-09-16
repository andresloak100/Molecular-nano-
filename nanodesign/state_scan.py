"""Bounded DFT starting-guess surveys at one immutable geometry.

Only the SCF starting guess changes. A completed survey, energy agreement or a
lower energy does not establish electronic-state identity or a ground state.
Every requested guess gets at most one attempt; reruns require a new survey.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
import math
from numbers import Real
import os
from pathlib import Path
import time

import numpy as np
from ase.units import Hartree

from .design import read_coordinate_snapshot
from .quantum import PySCFCalculator, QuantumSettings, SCF_INITIAL_GUESSES
from .workflow import json_write


TERMINAL = {"completed", "failed", "interrupted", "abandoned"}
EV_PER_KCAL_PER_MOL = 4184 / (6.02214076e23 * 1.602176634e-19)


def _now():
    return datetime.now(timezone.utc).isoformat()


def _load_average():
    try:
        return list(os.getloadavg())
    except (AttributeError, OSError):
        return None


def _digest(raw):
    return hashlib.sha256(raw).hexdigest()


def _object_digest(value):
    return _digest(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode())


def _guesses(values):
    if not isinstance(values, (list, tuple)) or not 1 <= len(values) <= 4:
        raise ValueError("Supply one to four explicit starting guesses in an ordered list or tuple.")
    if any(type(value) is not str or value not in SCF_INITIAL_GUESSES for value in values):
        raise ValueError(f"Starting guesses must be drawn from {SCF_INITIAL_GUESSES}.")
    if len(set(values)) != len(values):
        raise ValueError("Starting guesses must be distinct.")
    return list(values)


def _local(root, relative):
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("Expected a relative survey evidence path.")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Survey evidence path escapes its directory.")
    return path


def _json(raw):
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError("Survey records must be JSON objects.")
    json.dumps(value, allow_nan=False)
    return value


def _geometry(atoms):
    return {"symbols": atoms.get_chemical_symbols(),
            "positions_angstrom": atoms.positions.tolist(),
            "pbc": atoms.pbc.tolist()}


def _preflight(atoms, settings):
    if len(atoms) == 0 or np.any(atoms.pbc) or not np.all(np.isfinite(atoms.positions)):
        raise ValueError("Supply a nonempty, finite, nonperiodic geometry.")
    if np.any(atoms.numbers < 1) or np.any(atoms.numbers > 36):
        raise ValueError("The all-electron backend supports H through Kr only.")
    electrons = int(atoms.numbers.sum()) - settings.charge
    if electrons <= 0 or settings.spin > electrons or (electrons - settings.spin) % 2:
        raise ValueError("Electron count and requested charge/spin are inconsistent.")
    if len(atoms) > 1:
        distances = atoms.get_all_distances()
        np.fill_diagonal(distances, np.inf)
        if np.min(distances) < 0.1:
            raise ValueError("Atoms closer than 0.1 Angstrom; inspect the input geometry.")


@contextmanager
def _lock(root):
    try:
        import fcntl
    except ImportError as error:
        raise RuntimeError("Survey execution currently requires macOS or Linux advisory locks.") from error
    # Retain the lock file: unlinking it can permit two simultaneous workers.
    with _local(root, ".scan.lock").open("a+") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("Another worker is running this survey.") from error
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def create_state_scan(structure_path, output, settings: QuantumSettings, *,
                      image=-1, guesses=SCF_INITIAL_GUESSES):
    """Freeze exact source bytes, selected frame and settings without calculating.

    ``settings`` explicitly supplies the charge and alpha-minus-beta electron
    count. The file's calculator/constraints do not override them. All nuclei
    stay fixed and reported forces are unconstrained analytical forces.
    """
    if not isinstance(settings, QuantumSettings):
        raise TypeError("settings must be an explicit QuantumSettings instance.")
    if type(image) is not int:
        raise ValueError("image must be an integer ASE frame index.")
    guesses = _guesses(guesses)
    source = Path(structure_path).resolve()
    raw = source.read_bytes()
    atoms = read_coordinate_snapshot(source, raw, index=image)
    _preflight(atoms, settings)
    geometry = _geometry(atoms)
    input_record = {"path": f"inputs/{source.name}", "sha256": _digest(raw),
                    "source_path": str(source), "image": image}
    plan = {
        "schema_version": 1, "kind": "fixed_geometry_dft_starting_guess_survey",
        "created_utc": _now(), "input": input_record, "geometry": geometry,
        "geometry_sha256": _object_digest(geometry),
        "frame_sha256": _object_digest({"image": image, "geometry": geometry}),
        "quantum_settings": asdict(settings), "settings_sha256": _object_digest(asdict(settings)),
        "guesses": guesses, "units": {"length": "angstrom", "energy": "eV", "force": "eV/angstrom"},
        "geometry_optimized": False, "ground_state_verified": False,
        "electronic_state_identity_verified": False,
        "constraint_policy": "All coordinates frozen; source constraints do not mask reported forces.",
    }
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root / "inputs").mkdir()
    _local(root, input_record["path"]).write_bytes(raw)
    json_write(root / "plan.json", plan)
    json_write(root / "scan.json", {"schema_version": 1,
               "plan_sha256": _digest((root / "plan.json").read_bytes()),
               "attempts": {guess: None for guess in guesses}})
    return state_scan_report(root)


def _expected(plan, guess):
    settings = replace(QuantumSettings(**plan["quantum_settings"]), scf_initial_guess=guess)
    return asdict(settings)


def _binding(plan, guess):
    settings = _expected(plan, guess)
    return {"guess": guess, "input_sha256": plan["input"]["sha256"],
            "image": plan["input"]["image"], "frame_sha256": plan["frame_sha256"],
            "geometry_sha256": plan["geometry_sha256"], "geometry": plan["geometry"],
            "quantum_settings": settings, "settings_sha256": _object_digest(settings)}


def _finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def _validate_completed(record, plan, guess):
    energy = record.get("energy_ev")
    raw_forces = record.get("forces_ev_per_angstrom")
    forces = np.asarray(raw_forces, dtype=float)
    diagnostics = record.get("quantum_diagnostics", {})
    if not isinstance(diagnostics, dict):
        raise ValueError("Missing electronic diagnostics.")
    if (not _finite(energy) or forces.shape != (len(plan["geometry"]["symbols"]), 3)
            or not np.all(np.isfinite(forces))
            or diagnostics.get("scf_converged") is not True
            or diagnostics.get("gradient_completed") is not True):
        raise ValueError("Completed evaluation requires finite energy/raw forces and SCF/gradient success.")
    if any(not _finite(value) for row in raw_forces for value in row):
        raise ValueError("Force components must be finite numbers, not booleans or strings.")
    if (diagnostics.get("settings") != _expected(plan, guess)
            or diagnostics.get("scf_initial_guess") != guess):
        raise ValueError("Electronic diagnostics do not match the frozen per-guess settings.")
    # Validate settings types as well as equality (Python treats True == 1).
    QuantumSettings(**diagnostics["settings"])
    hartree = diagnostics.get("total_energy_hartree")
    if (not _finite(hartree) or not _finite(diagnostics.get("hartree_eV"))
            or not math.isclose(diagnostics["hartree_eV"], Hartree, rel_tol=0, abs_tol=1e-12)
            or not math.isclose(energy, hartree * Hartree, rel_tol=0, abs_tol=1e-8)):
        raise ValueError("Electronic energy units/accounting are inconsistent.")
    spin = plan["quantum_settings"]["spin"]
    expected = spin / 2 * (spin / 2 + 1)
    s2, deviation = diagnostics.get("s2"), diagnostics.get("s2_deviation")
    if (not all(_finite(value) for value in (s2, deviation, diagnostics.get("expected_s2")))
            or not math.isclose(diagnostics["expected_s2"], expected, abs_tol=1e-10, rel_tol=0)
            or not math.isclose(deviation, s2 - expected, abs_tol=1e-8, rel_tol=0)):
        raise ValueError("Electronic spin diagnostics are inconsistent.")


def _read_attempt(root, plan, guess, attempt):
    if not isinstance(attempt, dict) or attempt.get("status") not in TERMINAL | {"running"}:
        raise ValueError("Invalid survey attempt ledger.")
    if attempt.get("output") != f"attempts/{guess}":
        raise ValueError("Unexpected attempt output path.")
    digest = attempt.get("result_sha256")
    valid_digest = isinstance(digest, str) and len(digest) == 64 and all(char in "0123456789abcdef" for char in digest)
    if (("result_sha256" in attempt and not valid_digest)
            or (attempt["status"] == "completed" and not valid_digest)):
        raise ValueError("Completed evaluations require a valid recorded result digest.")
    path = _local(root, attempt["output"] + "/result.json")
    if not path.exists():
        if attempt.get("result_sha256") or attempt["status"] == "completed":
            raise ValueError("Recorded result is missing.")
        return None
    raw = path.read_bytes()
    if attempt.get("result_sha256") and _digest(raw) != attempt["result_sha256"]:
        raise ValueError("Recorded result has changed.")
    try:
        record = _json(raw)
    except (ValueError, UnicodeError):
        if attempt["status"] in {"running", "abandoned"}:
            return {"status": "unreadable", "evidence_error": "Incomplete result bytes retained; not a converged evaluation."}
        raise
    if record.get("status") not in {"running", "completed", "failed", "interrupted"}:
        raise ValueError("Invalid saved evaluation status.")
    if attempt["status"] in {"completed", "failed", "interrupted"} and record["status"] != attempt["status"]:
        raise ValueError("Evaluation status disagrees with the attempt ledger.")
    for key, expected in _binding(plan, guess).items():
        # Canonical JSON comparison also distinguishes booleans from numbers
        # and integer frame indices from floating-point aliases.
        if _object_digest(record.get(key)) != _object_digest(expected):
            raise ValueError(f"Saved evaluation differs from frozen {key}.")
    if record.get("units") != plan["units"]:
        raise ValueError("Saved evaluation units differ from the plan.")
    if record.get("ground_state_verified") is not False or record.get("electronic_state_identity_verified") is not False:
        raise ValueError("A starting-guess survey cannot certify state identity or the ground state.")
    if record["status"] == "completed":
        _validate_completed(record, plan, guess)
    return record


def _load(root):
    ledger = _json(_local(root, "scan.json").read_bytes())
    raw_plan = _local(root, "plan.json").read_bytes()
    if ledger.get("schema_version") != 1 or ledger.get("plan_sha256") != _digest(raw_plan):
        raise ValueError("Survey plan integrity check failed.")
    plan = _json(raw_plan)
    if plan.get("schema_version") != 1:
        raise ValueError("Unsupported survey schema.")
    guesses = _guesses(plan["guesses"])
    if not isinstance(ledger.get("attempts"), dict) or set(ledger["attempts"]) != set(guesses):
        raise ValueError("Attempt ledger does not match requested guesses.")
    settings = QuantumSettings(**plan["quantum_settings"])
    if plan["settings_sha256"] != _object_digest(asdict(settings)):
        raise ValueError("Frozen settings hash mismatch.")
    source = _local(root, plan["input"]["path"])
    raw = source.read_bytes()
    if _digest(raw) != plan["input"]["sha256"] or type(plan["input"]["image"]) is not int:
        raise ValueError("Frozen input integrity check failed.")
    atoms = read_coordinate_snapshot(source, raw, index=plan["input"]["image"])
    _preflight(atoms, settings)
    geometry = _geometry(atoms)
    if (geometry != plan["geometry"] or _object_digest(geometry) != plan["geometry_sha256"]
            or _object_digest({"image": plan["input"]["image"], "geometry": geometry}) != plan["frame_sha256"]):
        raise ValueError("Frozen frame/geometry integrity check failed.")
    records = {}
    for guess in guesses:
        attempt = ledger["attempts"][guess]
        if attempt is None:
            if _local(root, f"attempts/{guess}").exists():
                raise ValueError("Unrecorded attempt directory exists; refusing to overwrite evidence.")
            records[guess] = None
        else:
            records[guess] = _read_attempt(root, plan, guess, attempt)
    return plan, ledger, atoms, records


def _failure_diagnostics(calculator):
    diagnostics = getattr(calculator, "diagnostics", {})
    try:
        return json.loads(json.dumps(diagnostics, allow_nan=False))
    except (TypeError, ValueError):
        return {"unserializable_diagnostics_repr": repr(diagnostics),
                "evidence_error": "Rejected diagnostic values preserved as text; not accepted numerical evidence."}


def _evaluate(root, plan, guess, atoms, attempt):
    output = _local(root, attempt["output"])
    result = {"schema_version": 1, "status": "running", "started_utc": _now(),
              **_binding(plan, guess), "units": plan["units"],
              "ground_state_verified": False, "electronic_state_identity_verified": False,
              "geometry_optimized": False, "event_log": "electronic.jsonl",
              "load_average_start": _load_average(),
              "timing_caveat": "Elapsed time reflects concurrent host load; this is not a controlled speed benchmark."}
    calculator = None
    created = False
    started = time.monotonic()
    try:
        output.mkdir(parents=True, exist_ok=False)
        created = True
        json_write(output / "result.json", result)
        work = atoms.copy()
        work.calc = None
        work.set_constraint()
        calculator = PySCFCalculator(QuantumSettings(**result["quantum_settings"]), event_log=output / "electronic.jsonl")
        work.calc = calculator
        energy = work.get_potential_energy()
        forces = np.asarray(work.get_forces(apply_constraint=False))
        if isinstance(energy, (bool, np.bool_)) or not isinstance(energy, Real) or forces.dtype.kind not in "if":
            raise ValueError("Solver energy and raw forces must be real numerical values.")
        energy = float(energy)
        if _geometry(work) != plan["geometry"]:
            raise ValueError("Calculator modified the frozen input geometry.")
        if getattr(calculator, "atoms", None) is not None and _geometry(calculator.atoms) != plan["geometry"]:
            raise ValueError("Calculator used a geometry different from the frozen frame.")
        candidate = {**result, "status": "completed", "energy_ev": energy,
                     "forces_ev_per_angstrom": forces.tolist(),
                     "quantum_diagnostics": calculator.diagnostics}
        _validate_completed(candidate, plan, guess)
        json.dumps(candidate, allow_nan=False)
        result = candidate
        attempt["status"] = "completed"
    except BaseException as error:
        interrupted = isinstance(error, (KeyboardInterrupt, SystemExit))
        result.update(status="interrupted" if interrupted else "failed",
                      error={"type": type(error).__name__, "message": str(error)},
                      quantum_diagnostics=_failure_diagnostics(calculator))
        attempt.update(status=result["status"], error=result["error"])
        if interrupted or not isinstance(error, Exception):
            raise
    finally:
        result["elapsed_seconds"] = time.monotonic() - started
        result["load_average_end"] = _load_average()
        attempt["finished_utc"] = _now()
        if created:
            try:
                json_write(output / "result.json", result)
                attempt["result_sha256"] = _digest((output / "result.json").read_bytes())
            except Exception as error:
                attempt.update(status="failed", persistence_error={"type": type(error).__name__, "message": str(error)})
                raise


def run_state_scan(directory, *, max_jobs=1):
    """Launch at most 1–4 pending guesses, serially; never retry an old attempt.

    The job bound is not a runtime/memory limit. Each job computes both an SCF
    energy and an analytical gradient with the fixed plan's numerical settings.
    Ordinary failures consume a job and remain in the survey. Interrupts are
    saved and re-raised. A new invocation can continue untouched pending guesses.
    """
    if type(max_jobs) is not int or not 1 <= max_jobs <= 4:
        raise ValueError("max_jobs must be an integer from 1 through 4.")
    root = Path(directory).resolve()
    with _lock(root):
        plan, ledger, atoms, records = _load(root)
        for guess, attempt in ledger["attempts"].items():
            if attempt is not None and attempt["status"] == "running":
                record = records[guess]
                status = record.get("status") if record else None
                attempt["status"] = status if status in {"completed", "failed", "interrupted"} else "abandoned"
                attempt["recovered_utc"] = _now()
                path = _local(root, attempt["output"] + "/result.json")
                if path.exists():
                    attempt["result_sha256"] = _digest(path.read_bytes())
                if attempt["status"] == "abandoned":
                    attempt["error"] = {"type": "IncompleteAttempt", "message": "Prior worker left no terminal evaluation; original bytes retained, no automatic retry."}
        json_write(root / "scan.json", ledger)
        launched = 0
        for guess in plan["guesses"]:
            if ledger["attempts"][guess] is not None:
                continue
            if launched >= max_jobs:
                break
            attempt = {"status": "running", "started_utc": _now(), "output": f"attempts/{guess}"}
            ledger["attempts"][guess] = attempt
            json_write(root / "scan.json", ledger)
            launched += 1
            try:
                _evaluate(root, plan, guess, atoms, attempt)
            finally:
                json_write(root / "scan.json", ledger)
    return state_scan_report(root)


def state_scan_report(directory):
    """Read and verify all saved attempts without running or repairing anything."""
    root = Path(directory).resolve()
    plan, ledger, _, records = _load(root)
    rows = []
    for guess in plan["guesses"]:
        attempt, record = ledger["attempts"][guess], records[guess]
        row = {"guess": guess, "status": attempt["status"] if attempt else "pending"}
        if attempt:
            row["output"] = attempt["output"]
            if attempt.get("error"):
                row["error"] = attempt["error"]
        if record:
            row["calculation_status"] = record["status"]
            if record.get("evidence_error"):
                row["evidence_error"] = record["evidence_error"]
            if row["status"] == record["status"] == "completed":
                row.update(energy_ev=record["energy_ev"],
                           s2=record["quantum_diagnostics"]["s2"],
                           expected_s2=record["quantum_diagnostics"]["expected_s2"],
                           s2_deviation=record["quantum_diagnostics"]["s2_deviation"])
        rows.append(row)
    completed = [row for row in rows if row["status"] == "completed"]
    span = max(row["energy_ev"] for row in completed) - min(row["energy_ev"] for row in completed) if len(completed) >= 2 else None
    if span is not None and not math.isfinite(span):
        raise ValueError("Nonfinite aggregate energy spread.")
    spread_values = {"energy_spread_ev": span,
                     "energy_spread_hartree": span / Hartree if span is not None else None,
                     "energy_spread_kcal_per_mol": span / EV_PER_KCAL_PER_MOL if span is not None else None}
    if any(value is not None and not math.isfinite(value) for value in spread_values.values()):
        raise ValueError("Nonfinite aggregate energy spread after unit conversion.")
    finished = all(row["status"] in TERMINAL for row in rows)
    converged = len(completed) == len(rows)
    return {"schema_version": 1, "scan": str(root), "status": "finished" if finished else "partial" if any(ledger["attempts"].values()) else "pending",
            "requested_guesses": plan["guesses"],
            "attempted_guesses": [row["guess"] for row in rows if row["status"] != "pending"],
            "converged_guesses": [row["guess"] for row in completed],
            **{status + "_guesses": [row["guess"] for row in rows if row["status"] == status]
               for status in ("pending", "failed", "interrupted", "abandoned")},
            "all_requested_attempts_finished": finished,
            "all_requested_guesses_converged": converged,
            "complete_four_guess_check": converged and set(plan["guesses"]) == set(SCF_INITIAL_GUESSES),
            "energy_spread": {"scope": "converged_subset", "member_guesses": [row["guess"] for row in completed],
                              **spread_values},
            "unit_conversion": {"hartree_to_ev": float(Hartree),
                                "ev_per_kcal_per_mol": EV_PER_KCAL_PER_MOL,
                                "kcal_definition": "thermochemical calorie; exact SI electron-volt and Avogadro constants"},
            "rows": rows, "input": plan["input"], "geometry_sha256": plan["geometry_sha256"],
            "frame_sha256": plan["frame_sha256"], "quantum_settings": plan["quantum_settings"],
            "ground_state_verified": False, "electronic_state_identity_verified": False,
            "method_accuracy_validated": False, "geometry_optimized": False,
            "interpretation": "Finite energy spread among accepted SCF-and-gradient evaluations at one fixed geometry only. Unconverged, failed and unattempted guesses are excluded. Equal energies/spins do not identify the same state; no ground state or preferred guess is selected.",
            "resource_note": "max_jobs bounds serial starting guesses, not runtime or memory. No retries, fallback, orbital-stability analysis, relaxation or state certification."}
