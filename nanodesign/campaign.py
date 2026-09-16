"""Bounded, resumable calculation campaigns with preserved inputs and attempts.

A campaign enumerates poses; it is not an inverse designer or calibrated ranking
model. Each attempt runs the existing actual quantum workflow. No failed attempt
is overwritten, and restarting unfinished work requires an explicit retry flag.
"""
from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from dataclasses import asdict
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np

from .candidates import create_design
from .design import load_design, sha256
from .workflow import run


def _save(path, value):
    """Readers see either the previous complete record or the new one."""
    text = json.dumps(value, indent=2, allow_nan=False) + "\n"
    fd, temporary = tempfile.mkstemp(prefix=".record-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


@contextmanager
def _lock(root):
    # Advisory locks are released by the OS even if a worker is killed. The
    # file is intentionally retained: deleting lock files introduces races.
    try:
        import fcntl
    except ImportError as error:
        raise RuntimeError("Campaign execution currently requires macOS or Linux file locking.") from error
    with (root / ".campaign.lock").open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("Another worker is already running this campaign.") from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def _local(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Campaign artifact path escapes its directory.")
    return path


def create_campaign(designs, output, *, stage="singlepoint", state="initial",
                    fmax=0.03, steps=200, images=7):
    """Snapshot comparable designs and controls; do not start any calculations."""
    paths = [Path(path).resolve() for path in designs]
    if not 1 <= len(paths) <= 100:
        raise ValueError("A campaign must contain 1 to 100 explicit designs.")
    if stage not in {"singlepoint", "relax", "path"} or state not in {"initial", "final"}:
        raise ValueError("Invalid calculation stage or state.")
    if isinstance(fmax, bool) or not isinstance(fmax, (int, float)) or not math.isfinite(fmax) or fmax <= 0:
        raise ValueError("fmax must be a positive finite number.")
    if type(steps) is not int or steps < 1 or type(images) is not int or images < 5:
        raise ValueError("steps must be a positive integer and images at least five.")
    loaded = [load_design(path) for path in paths]
    first_data, first_atoms, _, first_settings, _ = loaded[0]
    for data, initial, _, settings, _ in loaded:
        if (not np.array_equal(initial.numbers, first_atoms.numbers)
                or settings != first_settings
                or data.get("fixed_indices", []) != first_data.get("fixed_indices", [])
                or data.get("hydrogen_transfer") != first_data.get("hydrogen_transfer")):
            raise ValueError("Campaign designs must share element order, quantum settings, anchor indices and reaction identities.")
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    plan = {"schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
            "kind": "explicit_pose_campaign", "quantum_settings": asdict(first_settings),
            "controls": {"stage": stage, "state": state, "fmax": fmax, "steps": steps, "images": images},
            "designs": [], "design_validated": False, "automatic_ranking": False,
            "interpretation": "Comparable numerical evaluations at explicit poses. Absolute energies include pose-dependent deformation and interaction; no success probability or best tool is inferred."}
    for index, (path, loaded_design) in enumerate(zip(paths, loaded), 1):
        data, _, _, _, source_hashes = loaded_design
        identity = f"pose-{index:04d}"
        destination = root / "designs" / identity
        destination.mkdir(parents=True)
        # Copy exact source bytes and keep the original design document. The
        # runnable snapshot only normalizes coordinate paths into this folder.
        for endpoint in ("initial", "final"):
            source = (path.parent / data[endpoint]).resolve()
            suffix = source.suffix
            target = destination / f"{endpoint}{suffix}"
            target.write_bytes(source.read_bytes())
            if sha256(target) != source_hashes[f"{endpoint}_sha256"]:
                raise ValueError("Source coordinates changed while creating campaign.")
        original = path.read_bytes()
        (destination / "source-design.json").write_bytes(original)
        if sha256(destination / "source-design.json") != source_hashes["design_sha256"]:
            raise ValueError("Source design changed while creating campaign.")
        snapshot = dict(data)
        for endpoint in ("initial", "final"):
            snapshot[endpoint] = f"{endpoint}{Path(data[endpoint]).suffix}"
        _save(destination / "design.json", snapshot)
        _, _, _, _, hashes = load_design(destination / "design.json")
        metadata = data.get("metadata", {})
        plan["designs"].append({"id": identity, "design": f"designs/{identity}/design.json",
            "input_hashes": hashes, "source_input_hashes": source_hashes,
            "pose": {key: metadata.get(key) for key in ("separation_angstrom", "lateral_offset_angstrom")}})
    _save(root / "plan.json", plan)
    _save(root / "campaign.json", {"schema_version": 1, "plan_sha256": sha256(root / "plan.json"),
          "attempts": {entry["id"]: [] for entry in plan["designs"]}})
    return campaign_report(root)


def create_pose_campaign(output, separations, offsets, *, settings=None, **controls):
    """Enumerate a finite grid of the supported H-abstraction candidate."""
    separations, offsets = list(separations), list(offsets)
    if not 1 <= len(separations) * len(offsets) <= 100:
        raise ValueError("The pose grid must contain 1 to 100 candidates.")
    if len(set(separations)) != len(separations) or len(set(offsets)) != len(offsets):
        raise ValueError("Pose grid values must be unique.")
    with tempfile.TemporaryDirectory(prefix="nanodesign-poses-") as temporary:
        paths = [create_design(Path(temporary) / f"pose-{i}-{j}", separation, offset, settings)
                 for i, separation in enumerate(separations) for j, offset in enumerate(offsets)]
        return create_campaign(paths, output, **controls)


def _load(root):
    state = json.loads((root / "campaign.json").read_text())
    if state.get("schema_version") != 1 or state.get("plan_sha256") != sha256(root / "plan.json"):
        raise ValueError("Campaign plan integrity check failed.")
    plan = json.loads((root / "plan.json").read_text())
    if plan.get("schema_version") != 1:
        raise ValueError("Unsupported campaign schema.")
    for entry in plan["designs"]:
        _, _, _, _, hashes = load_design(_local(root, entry["design"]))
        if hashes != entry["input_hashes"]:
            raise ValueError(f"Input integrity check failed for {entry['id']}.")
        if entry["id"] not in state["attempts"]:
            raise ValueError("Campaign attempt ledger is incomplete.")
    return plan, state


def _result(root, entry, attempt, controls):
    path = _local(root, attempt["output"]) / "result.json"
    if not path.exists():
        if attempt["status"] == "completed" or attempt.get("result_sha256"):
            raise ValueError("A recorded result is missing; completed evidence cannot be discarded.")
        return None
    if attempt.get("result_sha256") and sha256(path) != attempt["result_sha256"]:
        raise ValueError("A completed attempt's result has changed.")
    try:
        result = json.loads(path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError):
        if attempt["status"] in {"running", "abandoned"} and not attempt.get("result_sha256"):
            return {"status": "unreadable", "evidence_error": "Incomplete result document preserved; explicit retry creates a new attempt."}
        raise ValueError("The recorded result document is unreadable.") from None
    # Python's JSON reader accepts NaN/Infinity extensions; evidence may not.
    json.dumps(result, allow_nan=False)
    if result.get("status") not in {"running", "completed", "failed", "interrupted", "not_converged"}:
        raise ValueError("Unknown calculation status in saved evidence.")
    if attempt["status"] == "completed" and result["status"] != "completed":
        raise ValueError("Completed ledger entry disagrees with its saved result.")
    if result.get("input_hashes") != entry["input_hashes"]:
        raise ValueError(f"Result input hashes differ for {entry['id']}.")
    if result.get("stage") != controls["stage"] or result.get("state") != controls["state"]:
        raise ValueError("Result stage/state differs from the campaign plan.")
    if result.get("optimization") != {"fmax_ev_per_angstrom": controls["fmax"],
            "max_steps_per_stage": controls["steps"], "images": controls["images"]}:
        raise ValueError("Result optimization controls differ from the campaign plan.")
    data, atoms, _, settings, _ = load_design(_local(root, entry["design"]))
    expected_settings = asdict(settings)
    if result.get("quantum_settings") != expected_settings:
        raise ValueError("Result quantum settings differ from the campaign plan.")
    if result["status"] == "completed":
        _completed_evidence(result, controls, len(atoms), data.get("fixed_indices", []))
    return result


def _completed_evidence(result, controls, atom_count, fixed):
    def finite(value):
        return type(value) in (int, float) and math.isfinite(value)

    def structure(summary):
        energy, residual = summary.get("energy_ev"), summary.get("free_force_max_ev_per_angstrom")
        forces = np.asarray(summary.get("forces_ev_per_angstrom"), dtype=float)
        diagnostics = summary.get("quantum_diagnostics", {})
        if (not finite(energy) or not finite(residual) or residual < 0
                or forces.shape != (atom_count, 3) or not np.all(np.isfinite(forces))
                or diagnostics.get("scf_converged") is not True
                or diagnostics.get("gradient_completed") is not True):
            raise ValueError("Completed result lacks finite energy/force or converged electronic evidence.")
        forces = forces.copy()
        forces[fixed] = 0
        if not math.isclose(residual, float(np.linalg.norm(forces, axis=1).max()), rel_tol=1e-7, abs_tol=1e-9):
            raise ValueError("Reported force residual disagrees with saved forces.")
        return residual

    if controls["stage"] in {"singlepoint", "relax"}:
        residual = structure(result.get("structure", {}))
        if controls["stage"] == "relax" and (result.get("geometry_converged") is not True or residual > controls["fmax"]):
            raise ValueError("Completed relaxation lacks geometry-convergence evidence.")
    else:
        if result.get("endpoints_converged") is not True or result.get("neb_converged") is not True:
            raise ValueError("Completed path lacks endpoint/NEB convergence evidence.")
        endpoints, images = result.get("endpoints", []), result.get("images", [])
        if len(endpoints) != 2 or len(images) != controls["images"]:
            raise ValueError("Completed path has an incomplete set of structures.")
        for endpoint in endpoints:
            residual = structure(endpoint)
            if (endpoint.get("geometry_converged") is not True or residual > controls["fmax"]
                    or endpoint.get("endpoint_identity_ok") is False
                    or endpoint.get("topology_screen", {}).get("preserved") is False):
                raise ValueError("Completed path has invalid endpoint evidence.")
        for image in images:
            structure(image)
        relative = np.asarray(result.get("relative_energies_ev"), dtype=float)
        expected = np.array([image["energy_ev"] - images[0]["energy_ev"] for image in images])
        if relative.shape != expected.shape or not np.allclose(relative, expected, atol=1e-8, rtol=0):
            raise ValueError("Path energy profile disagrees with image energies.")
        barrier = result.get("candidate_electronic_barrier_ev")
        if barrier is not None:
            peak = result.get("peak_image")
            if (not finite(barrier) or type(peak) is not int or not 0 < peak < len(images)-1
                    or peak != int(np.argmax(relative)) or not math.isclose(barrier, relative[peak], abs_tol=1e-8)
                    or not result.get("barrier_interpretation")):
                raise ValueError("Candidate barrier lacks matching path evidence and interpretation.")


def run_campaign(directory, *, max_jobs=1, retry_incomplete=False):
    """Run at most max_jobs serial attempts; completed jobs are never repeated.

    A retry starts a new calculation from the preserved input. It does not
    resume optimizer coordinates or SCF state from a partial attempt.
    """
    if type(max_jobs) is not int or not 1 <= max_jobs <= 100:
        raise ValueError("max_jobs must be an integer between 1 and 100.")
    if type(retry_incomplete) is not bool:
        raise ValueError("retry_incomplete must be boolean.")
    root = Path(directory).resolve()
    with _lock(root):
        plan, state = _load(root)
        started = 0
        for entry in plan["designs"]:
            attempts = state["attempts"][entry["id"]]
            if attempts:
                previous = attempts[-1]
                saved = _result(root, entry, previous, plan["controls"])
                if previous["status"] == "running":
                    terminal = saved and saved.get("status") in {"completed", "failed", "interrupted", "not_converged"}
                    previous["status"] = saved["status"] if terminal else "abandoned"
                    if saved and saved.get("evidence_error"):
                        previous["evidence_error"] = saved["evidence_error"]
                    if terminal:
                        previous["result_sha256"] = sha256(_local(root, previous["output"]) / "result.json")
                    _save(root / "campaign.json", state)
                if previous["status"] == "completed" or not retry_incomplete:
                    continue
            if started >= max_jobs:
                break
            attempt = {"status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
                       "output": f"runs/{entry['id']}/attempt-{len(attempts)+1:04d}"}
            attempts.append(attempt)
            _save(root / "campaign.json", state)
            started += 1
            output = _local(root, attempt["output"])
            try:
                result = run(_local(root, entry["design"]), output, **plan["controls"])
                attempt["status"] = result["status"]
                saved = _result(root, entry, attempt, plan["controls"])
                if saved is None or saved["status"] != result["status"]:
                    raise ValueError("Returned calculation status lacks matching saved evidence.")
            except KeyboardInterrupt:
                attempt.update(status="interrupted", error="KeyboardInterrupt")
                raise
            except Exception as error:
                attempt.update(status="failed", error=f"{type(error).__name__}: {error}")
            finally:
                try:
                    if (output / "result.json").exists():
                        _result(root, entry, attempt, plan["controls"])
                        attempt["result_sha256"] = sha256(output / "result.json")
                except Exception as error:
                    attempt.update(status="failed", evidence_integrity_error=f"{type(error).__name__}: {error}")
                    raise
                finally:
                    attempt["finished_utc"] = datetime.now(timezone.utc).isoformat()
                    _save(root / "campaign.json", state)
    return campaign_report(root)


def campaign_report(directory):
    """Read the evidence in plan order, without selecting a scientific winner."""
    root = Path(directory).resolve()
    plan, state = _load(root)
    rows = []
    for entry in plan["designs"]:
        attempts = state["attempts"][entry["id"]]
        latest = attempts[-1] if attempts else None
        row = {"id": entry["id"], "pose": entry["pose"], "attempts": len(attempts),
               "status": latest["status"] if latest else "pending", "design_validated": False,
               "electronic_state_identity_verified": False}
        if latest:
            row["result"] = latest["output"] + "/result.json"
            result = _result(root, entry, latest, plan["controls"])
            if result:
                row["calculation_status"] = result["status"]
                if result.get("evidence_error"):
                    row["evidence_error"] = result["evidence_error"]
                if result["status"] == "completed" and latest["status"] == "completed":
                    structure = result.get("structure", {})
                    row.update(energy_ev=structure.get("energy_ev"),
                        free_force_max_ev_per_angstrom=structure.get("free_force_max_ev_per_angstrom"),
                        hydrogen_basin_preserved=structure.get("endpoint_identity_ok"),
                        topology_screen=structure.get("topology_screen"),
                        anchor_force_groups=structure.get("anchor_force_groups"),
                        candidate_electronic_barrier_ev=result.get("candidate_electronic_barrier_ev"),
                        barrier_interpretation=result.get("barrier_interpretation"),
                        electronic_diagnostics={key: structure.get("quantum_diagnostics", {}).get(key)
                            for key in ("scf_converged", "s2", "expected_s2", "s2_deviation", "stability_checked")},
                        endpoints=[{key: endpoint.get(key) for key in ("state", "geometry_converged",
                            "endpoint_identity_ok", "topology_screen", "free_force_max_ev_per_angstrom", "quantum_diagnostics")}
                            for endpoint in result.get("endpoints", [])],
                        neb_converged=result.get("neb_converged"),
                        validation=result.get("validation"))
        rows.append(row)
    return {"schema_version": 1, "campaign": str(root), "controls": plan["controls"],
            "quantum_settings": plan["quantum_settings"], "electronic_state_identity_verified": False,
            "counts": dict(Counter(row["status"] for row in rows)), "rows": rows,
            "design_validated": False, "automatic_ranking": False,
            "interpretation": plan["interpretation"],
            "resource_note": "max_jobs bounds attempts, not runtime or the evaluations within each relaxation/path. Execution is serial. No automatic retry or method change."}
