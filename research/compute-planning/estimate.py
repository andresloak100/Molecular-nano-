"""Read-only serial cost scenarios from archived timings; standard library only.

No solver, optimizer, subprocess, network, or output-file writes are used.
Counts describe declared geometry-evaluation scenarios, not runtime guarantees.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ARCHIVES = {
    "direct": "data/validation/h-abstraction-direct-initial",
    "density_fitting": "data/validation/h-abstraction-df-initial",
}


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number")
    return float(value)


def path_counts(images=7, endpoint_steps=200, neb_steps=200, climbing_steps=200):
    """Count stage-start and moved-geometry batches, with/without cache credit.

    Steps are assumed executed updates; endpoint_steps applies to each endpoint.
    All stages are assumed reached. Repeated same-geometry force/energy requests
    within logging/convergence checks are not separate quantum evaluations.
    """
    integer(images, "images", 5)
    for name, value in (("endpoint_steps", endpoint_steps), ("neb_steps", neb_steps),
                        ("climbing_steps", climbing_steps)):
        integer(value, name)
    active = images - 2
    iteration_only = 2 * endpoint_steps + active * (neb_steps + climbing_steps)
    stages = [
        {"stage": "initial_endpoint", "assumed_steps": endpoint_steps, "active_geometries": 1},
        {"stage": "final_endpoint", "assumed_steps": endpoint_steps, "active_geometries": 1},
        {"stage": "ordinary_neb", "assumed_steps": neb_steps, "active_geometries": active},
        {"stage": "climbing_neb", "assumed_steps": climbing_steps, "active_geometries": active},
    ]
    for stage in stages:
        stage["stage_start_geometry_requests"] = stage["active_geometries"]
        stage["geometry_evaluations_before_cross_stage_cache_credit"] = (
            stage["assumed_steps"] + 1) * stage["active_geometries"]
    before_credit = sum(stage["geometry_evaluations_before_cross_stage_cache_credit"] for stage in stages)
    return {
        "images_including_endpoints": images,
        "active_band_images": active,
        "stages": stages,
        "iteration_only_evaluations": iteration_only,
        "stage_start_geometry_requests": 2 + 2 * active,
        "geometry_evaluations_without_cross_stage_cache_credit": before_credit,
        "ordinary_to_climbing_cache_credit": active,
        "geometry_evaluations_with_cross_stage_cache": before_credit - active,
        "scope": "Two endpoint relaxations then ordinary and climbing NEB; serial images. Counts are scenarios, not strict call caps or raw API-call counts.",
    }


def hessian_counts(atom_count, fixed_indices, max_free_coordinates=120):
    integer(atom_count, "atom_count", 1)
    integer(max_free_coordinates, "max_free_coordinates", 1)
    if not isinstance(fixed_indices, (list, tuple)):
        raise ValueError("fixed_indices must be an explicit list of whole-atom indices")
    for index in fixed_indices:
        integer(index, "fixed index")
        if index >= atom_count:
            raise ValueError("fixed index is outside the structure")
    if len(set(fixed_indices)) != len(fixed_indices):
        raise ValueError("fixed indices must be unique")
    free = 3 * (atom_count - len(fixed_indices))
    if not free:
        raise ValueError("at least one atom must be free")
    return {
        "atoms": atom_count,
        "fixed_atoms": len(fixed_indices),
        "fixed_indices": list(fixed_indices),
        "free_cartesian_coordinates": free,
        "full_stencil_force_requests": 1 + 2 * free,
        "workflow_force_requests_including_guard": 2 + 2 * free,
        "cold_start_quantum_evaluations_if_all_guards_pass": 1 + 2 * free,
        "displaced_geometry_evaluations": 2 * free,
        "max_free_coordinates": max_free_coordinates,
        "coordinate_guard_allows_run": free <= max_free_coordinates,
        "quantum_evaluations_if_coordinate_guard_rejects": 0,
        "quantum_evaluations_if_initial_force_guard_rejects": 1,
        "scope": "One baseline plus two displacements per free Cartesian coordinate. Characterization reuses its initial force-guard cache. A complete stencil does not verify a saddle or electronic state.",
    }


def duration(evaluations, seconds_per_evaluation):
    integer(evaluations, "evaluations")
    seconds_per_evaluation = positive(seconds_per_evaluation, "seconds_per_evaluation")
    try:
        seconds = evaluations * seconds_per_evaluation
    except OverflowError as error:
        raise ValueError("scenario duration exceeds finite numeric range") from error
    if not math.isfinite(seconds):
        raise ValueError("scenario duration exceeds finite numeric range")
    return {"seconds": seconds, "hours": seconds / 3600, "days": seconds / 86400}


def load_measurements(project_root=PROJECT_ROOT):
    """Read exact result/geometry bytes; fail rather than invent missing timing."""
    root = Path(project_root)
    measurements = {}
    for name, relative in ARCHIVES.items():
        directory = root / relative
        raw = (directory / "result.json").read_bytes()
        record = json.loads(raw)
        interpretation = json.loads((directory / "interpretation.json").read_bytes())
        digest = hashlib.sha256(raw).hexdigest()
        if interpretation["file_sha256"]["result.json"] != digest:
            raise ValueError(f"{name}: archived result hash differs from its interpretation record")
        diagnostics = record["structure"]["quantum_diagnostics"]
        if (record.get("status") != "completed" or record.get("stage") != "singlepoint"
                or record.get("state") != "initial"
                or diagnostics.get("scf_converged") is not True
                or diagnostics.get("gradient_completed") is not True):
            raise ValueError(f"{name}: a completed initial energy-and-force single point is required")
        coordinates = (directory / "input-initial.extxyz").read_bytes()
        geometry_hash = hashlib.sha256(coordinates).hexdigest()
        if geometry_hash != record["input_hashes"]["initial_sha256"]:
            raise ValueError(f"{name}: archived input geometry hash mismatch")
        lines = coordinates.decode("utf-8").splitlines()
        atoms = integer(int(lines[0]), "XYZ atom count", 1)
        if len(lines) < atoms + 2 or any(line.strip() for line in lines[atoms + 2:]):
            raise ValueError(f"{name}: expected one complete XYZ structure")
        elements = Counter(line.split()[0] for line in lines[2:atoms + 2])
        settings = record["quantum_settings"]
        if settings.get("density_fit") is not (name == "density_fitting"):
            raise ValueError(f"{name}: density-fitting setting disagrees with archive label")
        effective_threads = diagnostics.get("effective_pyscf_threads")
        if effective_threads is not None:
            integer(effective_threads, "effective_pyscf_threads", 1)
        fixed = record["design"]["fixed_indices"]
        hessian_counts(atoms, fixed)  # validate the counted whole-atom constraints
        measurements[name] = {
            "source": f"{relative}/result.json",
            "source_sha256": digest,
            "geometry_sha256": geometry_hash,
            "atoms": atoms,
            "elements": dict(sorted(elements.items())),
            "fixed_indices": fixed,
            "basis_functions": integer(diagnostics["basis_functions"], "basis_functions", 1),
            "quantum_settings_as_recorded": settings,
            "singlepoint_run_elapsed_seconds": positive(record["elapsed_seconds"], "singlepoint elapsed"),
            "calculator_elapsed_seconds": positive(diagnostics["elapsed_seconds"], "calculator elapsed"),
            "requested_pyscf_threads": settings["threads"],
            "effective_pyscf_threads_as_recorded": effective_threads,
            "effective_thread_evidence": "Recorded in this run" if effective_threads is not None else "Not recorded; direct run predates thread instrumentation. Later same-build checks observed one effective PySCF thread.",
            "timing_caveat": interpretation["timing_caveat"],
        }
    first, second = measurements.values()
    for key in ("geometry_sha256", "atoms", "elements", "fixed_indices", "basis_functions"):
        if first[key] != second[key]:
            raise ValueError(f"Direct/DF measurement mismatch: {key}")
    if dict(first["quantum_settings_as_recorded"], density_fit=True) != second["quantum_settings_as_recorded"]:
        raise ValueError("Direct/DF physical/numerical settings differ beyond density fitting")
    return measurements


def estimate(project_root=PROJECT_ROOT, *, images=7, steps=200, max_free_coordinates=120):
    path = path_counts(images, steps, steps, steps)
    measurements = load_measurements(project_root)
    candidate = measurements["direct"]
    hessian = hessian_counts(candidate["atoms"], candidate["fixed_indices"], max_free_coordinates)
    projections = {}
    for name, measured in measurements.items():
        elapsed = measured["singlepoint_run_elapsed_seconds"]
        projections[name] = {
            "timing_basis": "Full archived singlepoint-run elapsed, including SCF plus analytical gradient and small workflow overhead; used as a constant per-geometry proxy.",
            "seconds_per_geometry_proxy": elapsed,
            "path_with_cross_stage_cache": duration(path["geometry_evaluations_with_cross_stage_cache"], elapsed),
            "path_without_cross_stage_cache_credit": duration(path["geometry_evaluations_without_cross_stage_cache_credit"], elapsed),
            "full_hessian_if_guards_pass": duration(hessian["cold_start_quantum_evaluations_if_all_guards_pass"], elapsed),
        }
    return {
        "schema_version": 1,
        "kind": "serial_computation_scenario",
        "quantum_jobs_launched": 0,
        "runtime_guaranteed": False,
        "scientific_model_validated": False,
        "measurements": measurements,
        "path": path,
        "hessian": hessian,
        "projections": projections,
        "assumptions_and_limits": [
            "All four path stages are reached; every assumed update changes the evaluated geometry. Early convergence or rejection can reduce work.",
            "Ordinary and climbing NEB share image calculators: the climbing stage's initial unchanged geometries normally reuse cached energies/forces. Endpoint and summary reads also normally reuse caches.",
            "IDPP initialization uses its auxiliary distance objective, not the quantum backend. Its runtime and optimizer/I/O overhead are not separately modeled.",
            "Each new-geometry PySCF call evaluates energy AND forces. The current adapter restarts SCF without carrying a saved density between changed geometries.",
            "Full initial singlepoint elapsed is a proxy, not gradient-only time. Geometry-dependent SCF iterations, failures, retries, cache invalidations and concurrent workloads can change cost substantially.",
            "The direct/DF observations came from overlapping uncontrolled workloads. Their projected times are separate scenarios, not uncertainty bounds or a guaranteed speedup.",
            "Production's default coordinate limit is 120. See coordinate_guard_allows_run for this structure and planning limit. Raising a planning limit here launches nothing; production separately enforces its coordinate limit and initial-force guard.",
            "A second independent displacement-step check costs another full stencil under the same assumptions. Connectivity checks and alternative electronic-state calculations are additional work.",
            "Production uses CPU PySCF and serial NEB images. GPU execution and parallel image scheduling are not integrated or benchmarked; no accelerator speedup is assumed.",
            "These timings apply to the recorded model/method only. A smaller handle or reactive region is a different physical model requiring its own measurements and fidelity checks.",
        ],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--images", type=int, default=7)
    parser.add_argument("--steps", type=int, default=200, help="Assumed executed updates in EACH of the four stages; scenario only")
    parser.add_argument("--max-free-coordinates", type=int, default=120, help="Report the characterization guard; never launches a calculation")
    args = parser.parse_args(argv)
    try:
        report = estimate(args.project_root, images=args.images, steps=args.steps,
                          max_free_coordinates=args.max_free_coordinates)
    except (ValueError, OSError, KeyError, TypeError, IndexError) as error:
        parser.exit(2, f"Cannot estimate from supplied evidence: {error}\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
