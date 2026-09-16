"""Real calculations with immutable outputs and conservative result interpretation."""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import time

import numpy as np
from ase.io import write
from ase.mep import NEB
from ase.optimize import FIRE

from . import __version__
from .design import load_design, endpoint_identity_ok, transfer_distances, topology_screen
from .quantum import PySCFCalculator


def json_write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def force_max(atoms, constrained=True):
    return float(np.linalg.norm(atoms.get_forces(apply_constraint=constrained), axis=1).max())


def summarize_atoms(atoms, reaction, state=None, metadata=None):
    raw = atoms.get_forces(apply_constraint=False)
    fixed = np.concatenate([c.get_indices() for c in atoms.constraints]).astype(int) if atoms.constraints else np.array([], dtype=int)
    groups = {}
    for name in ("substrate", "tool"):
        indices = sorted(set((metadata or {}).get(f"{name}_indices", [])) & set(fixed.tolist()))
        if indices:
            internal = raw[indices].sum(axis=0)
            groups[name] = {"anchor_indices": indices, "internal_force_sum_ev_per_angstrom": internal.tolist(),
                            "external_holding_force_ev_per_angstrom": (-internal).tolist()}
    return {
        "energy_ev": float(atoms.get_potential_energy()),
        "free_force_max_ev_per_angstrom": force_max(atoms),
        "forces_ev_per_angstrom": raw.tolist(),
        "anchor_force_sum_ev_per_angstrom": raw[fixed].sum(axis=0).tolist(),
        "anchor_force_groups": groups,
        "anchor_force_interpretation": "Internal force; global sum can cancel opposing loads. Separate group holding forces are negatives of group internal forces.",
        "hydrogen_transfer": transfer_distances(atoms, reaction),
        "endpoint_identity_ok": endpoint_identity_ok(atoms, reaction, state) if state else None,
        "topology_screen": topology_screen(atoms, (metadata or {}).get("reactant_bonds" if state == "initial" else "product_bonds")) if state else {"assessed": False, "preserved": None},
        "quantum_diagnostics": atoms.calc.diagnostics,
    }


def relax(atoms, out, label, fmax, steps):
    opt = FIRE(atoms, trajectory=str(out / f"{label}.traj"), logfile=str(out / f"{label}.log"), maxstep=0.08)
    converged = bool(opt.run(fmax=fmax, steps=steps))
    write(out / f"{label}.extxyz", atoms)
    return converged


def audit_result(result):
    """Implemented evidence gates; never infer manufacture/reliability from an energy."""
    stage = result.get("stage")
    numerical = result.get("status") == "completed"
    if stage == "path":
        numerical = numerical and result.get("neb_converged", False) and result.get("endpoints_converged", False)
    elif stage == "relax":
        numerical = numerical and result.get("geometry_converged", False)
    else:
        numerical = False
    missing = [
        "Independent reaction-specific quantum benchmark and error estimate",
        "Larger basis, integration-grid, and cluster/boundary convergence",
        "Electronic-state stability and alternative spin-state investigation",
        "Transition-state modes and connectivity (or a separately verified barrierless path)",
        "Competing reactions, approach/retraction paths, and tip regeneration",
        "Finite-temperature free energies, tunnelling, and positioning uncertainty",
        "Experimental confirmation and whole-machine integration",
    ]
    return {"design_validated": False,
            "numerical_stage_converged": bool(numerical),
            "accuracy_claim": "Uncalibrated finite-cluster quantum model; no guaranteed design accuracy.",
            "missing_evidence": missing}


def run(design_path, output, stage="singlepoint", state="initial", fmax=0.03, steps=200, images=7):
    if stage not in {"singlepoint", "relax", "path"} or state not in {"initial", "final"}:
        raise ValueError("Invalid calculation stage or state.")
    if not np.isfinite(fmax) or fmax <= 0 or type(steps) is not int or steps < 1:
        raise ValueError("fmax must be positive and finite, and steps a positive integer.")
    if type(images) is not int or images < 5:
        raise ValueError("Use at least five images (including endpoints).")
    data, initial, final, settings, hashes = load_design(design_path)
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    result = {"schema_version": 1, "stage": stage, "state": state, "status": "running",
              "started_utc": datetime.now(timezone.utc).isoformat(), "software_version": __version__,
              "python_version": platform.python_version(), "input_hashes": hashes,
              "quantum_settings": asdict(settings), "design": data,
              "optimization": {"fmax_ev_per_angstrom": fmax, "max_steps_per_stage": steps, "images": images},
              "units": {"length": "angstrom", "energy": "eV", "force": "eV/angstrom"}}
    json_write(out / "result.json", result)
    write(out / "input-initial.extxyz", initial)
    write(out / "input-final.extxyz", final)
    reaction = data.get("hydrogen_transfer")
    metadata = data.get("metadata", {})
    for atoms in (initial, final):
        atoms.calc = PySCFCalculator(settings, event_log=out / "electronic.jsonl")
    try:
        if stage in {"singlepoint", "relax"}:
            atoms = initial if state == "initial" else final
            if stage == "relax":
                result["geometry_converged"] = relax(atoms, out, state, fmax, steps)
            result["structure"] = summarize_atoms(atoms, reaction, state, metadata)
            write(out / "structure.extxyz", atoms)
            result["status"] = "completed" if stage == "singlepoint" or result["geometry_converged"] else "not_converged"
        else:
            results = []
            for label, atoms in (("initial", initial), ("final", final)):
                ok = relax(atoms, out, label, fmax, steps)
                summary = summarize_atoms(atoms, reaction, label, metadata)
                results.append({"state": label, "geometry_converged": ok, **summary})
                result["endpoints"] = results
                json_write(out / "result.json", result)
                if not ok:
                    raise ValueError(f"{label} endpoint did not converge; no path barrier reported.")
                if summary["endpoint_identity_ok"] is False:
                    raise ValueError(f"{label} endpoint changed reaction basin; choose another tool pose. No barrier reported.")
                if summary["topology_screen"]["preserved"] is False:
                    raise ValueError(f"{label} endpoint changed nominal connectivity under geometric screening; inspect its structure. No barrier reported.")
            result["endpoints_converged"] = True
            if np.max(np.abs(initial.positions - final.positions)) < 1e-3:
                raise ValueError("Both endpoints relaxed to the same structure; no reaction path.")
            band = [initial] + [initial.copy() for _ in range(images - 2)] + [final]
            for atoms in band[1:-1]:
                atoms.calc = PySCFCalculator(settings, event_log=out / "electronic.jsonl")
            neb = NEB(band, k=0.1, climb=False, method="improvedtangent", remove_rotation_and_translation=False)
            neb.interpolate(method="idpp", apply_constraint=True)
            preliminary = FIRE(neb, trajectory=str(out / "band-pre.traj"), logfile=str(out / "band-pre.log"), maxstep=0.05)
            pre_ok = bool(preliminary.run(fmax=max(0.1, 3 * fmax), steps=steps))
            result["pre_neb_converged"] = pre_ok
            if not pre_ok:
                write(out / "path.extxyz", band)
                raise ValueError("Initial NEB did not converge; no barrier reported.")
            neb.climb = True
            opt = FIRE(neb, trajectory=str(out / "band.traj"), logfile=str(out / "band.log"), maxstep=0.03)
            ok = bool(opt.run(fmax=fmax, steps=steps))
            result["neb_converged"] = ok
            result["neb_force_max_ev_per_angstrom"] = float(np.linalg.norm(neb.get_forces(), axis=1).max())
            energies = [float(a.get_potential_energy()) for a in band]
            result["images"] = [summarize_atoms(a, reaction, metadata=metadata) for a in band]
            result["relative_energies_ev"] = [e - energies[0] for e in energies]
            result["reaction_energy_ev"] = energies[-1] - energies[0]
            write(out / "path.extxyz", band)
            peak = int(np.argmax(energies))
            if ok and 0 < peak < len(band) - 1:
                result["candidate_electronic_barrier_ev"] = energies[peak] - energies[0]
                result["peak_image"] = peak
                result["barrier_interpretation"] = "Converged CI-NEB estimate on a fixed-anchor DFT potential surface; not a validated transition state, free-energy barrier, or error probability."
            else:
                result["candidate_electronic_barrier_ev"] = None
                result["barrier_interpretation"] = "No converged interior saddle candidate; this does not establish barrierless assembly."
            result["status"] = "completed" if ok else "not_converged"
    except KeyboardInterrupt:
        result["status"] = "interrupted"
        result["error"] = {"type": "KeyboardInterrupt", "message": "Interrupted before completion; partial data are not converged results."}
        raise
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = {"type": type(exc).__name__, "message": str(exc)}
        raise
    finally:
        result["elapsed_seconds"] = time.monotonic() - started
        result["validation"] = audit_result(result)
        json_write(out / "result.json", result)
    return result
