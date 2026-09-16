"""Cost of the 53-atom reaction path, from measurements already in the repository.

This script computes no chemistry. It reads the two archived 53-atom records
and the structure of ``nanodesign.workflow.run(stage="path")``, and turns them
into a wall-clock projection with its assumptions written down.

What the archived records actually are: one energy-and-gradient evaluation
each, at the unrelaxed pose, PBE0-D3(BJ)/def2-SVP, 463 basis functions, one
effective thread (``threads_honored`` is false; this PySCF build has no
OpenMP). Both ran while the host was carrying other jobs, so both wall-clock
figures include contention of an amount nobody recorded at the time. They are
upper bounds on the uncontended cost, and the direct-over-DF ratio is not a
controlled benchmark.

What the path stage does, read from the code rather than assumed:

    two endpoint relaxations   FIRE, up to `steps` evaluations each
    preliminary NEB            FIRE, up to `steps` steps x (images - 2) images
    climbing-image NEB         FIRE, up to `steps` steps x (images - 2) images

Images are evaluated serially inside one process and no SCF density is carried
between changed geometries, so every evaluation pays the full SCF cost. At the
default 7 images and 200 steps the ceiling is 2*200 + 2*200*5 = 2400
evaluations.

The ceiling is not a prediction. A FIRE relaxation from a hand-built starting
guess rarely needs its step limit, and the preliminary NEB converges to a
loose tolerance. The projection therefore reports three scenarios with
explicitly stated step counts, and the report quotes the range, not a point.

Run from the repository root:

    python research/candidate-feasibility/cost_model.py
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

LANE = Path(__file__).resolve().parent
EVIDENCE = LANE / "evidence"
VALIDATION = ROOT / "data" / "validation"

HOURS = 3600.0
DAYS = 86400.0


def _load_record(directory: Path) -> dict:
    return json.loads((directory / "result.json").read_text())


def measured_evaluation_costs() -> dict:
    """Seconds per energy+gradient evaluation, as archived. No new compute."""
    direct = _load_record(VALIDATION / "h-abstraction-direct-initial")
    density_fit = _load_record(VALIDATION / "h-abstraction-df-initial")
    comparison = json.loads((VALIDATION / "h-abstraction-df-initial" / "comparison.json").read_text())

    costs = {}
    for label, record in (("direct", direct), ("density_fitting", density_fit)):
        diagnostics = record["structure"]["quantum_diagnostics"]
        # The direct run started before the thread diagnostics were added, so
        # its effective thread count was never recorded. Report that gap
        # instead of assuming the value the DF run happens to show.
        threads_recorded = "effective_pyscf_threads" in diagnostics
        costs[label] = {
            "seconds_per_energy_and_gradient": diagnostics["elapsed_seconds"],
            "basis_functions": diagnostics["basis_functions"],
            "scf_cycles": diagnostics["scf_cycles"],
            "effective_threads": diagnostics.get("effective_pyscf_threads"),
            "threads_honored": diagnostics.get("threads_honored"),
            "effective_threads_recorded": threads_recorded,
            "settings": diagnostics["settings"],
            "source": str(record["stage"]) + " record",
        }
        if not threads_recorded:
            costs[label]["threads_note"] = (
                "This run predates the thread diagnostics (added in 121c30f), so "
                "its effective thread count is not recorded. The requested value "
                f"was {diagnostics['settings'].get('threads')}. This PySCF build has "
                "no OpenMP and was measured at one effective thread in every "
                "instrumented run, so one thread is the expected value here, but "
                "it is an inference from the build, not a measurement of this run."
            )

    # The DF run has stage instrumentation; the direct run predates it.
    events = [
        json.loads(line)
        for line in (VALIDATION / "h-abstraction-df-initial" / "electronic.jsonl").read_text().splitlines()
        if line.strip()
    ]
    scf_completed = next(e for e in events if e["event"] == "scf_completed")
    finished = next(e for e in events if e["event"] == "calculation_completed")
    costs["density_fitting"]["scf_seconds"] = scf_completed["elapsed_seconds"]
    costs["density_fitting"]["gradient_seconds"] = finished["elapsed_seconds"] - scf_completed["elapsed_seconds"]
    costs["direct"]["scf_seconds"] = None
    costs["direct"]["gradient_seconds"] = None
    costs["direct"]["stage_breakdown_note"] = (
        "The direct run predates stage instrumentation; no SCF/gradient split exists for it."
    )

    costs["contention"] = {
        "both_runs_concurrent": True,
        "controlled_benchmark": False,
        "note": (
            "Both timings were taken with other jobs on the host and no load "
            "average was recorded. Each is an upper bound on its uncontended "
            "cost by an unknown factor, and the "
            f"{comparison['observed_walltime_ratio_direct_over_density_fitting']:.2f}x "
            "direct-over-DF ratio is an observation, not a speedup measurement."
        ),
    }
    costs["density_fitting_accuracy_cost"] = {
        "max_force_component_deviation_ev_per_angstrom": comparison["max_force_component_difference_ev_per_angstrom"],
        "rms_force_component_deviation_ev_per_angstrom": comparison["rms_force_component_difference_ev_per_angstrom"],
        "s2_difference": abs(comparison["s2_direct"] - comparison["s2_density_fitting"]),
        "optimization_tolerance_ev_per_angstrom": 0.03,
        "note": (
            "Measured at one geometry (the unrelaxed pose) only. A deviation "
            "far below the force tolerance there does not guarantee the same "
            "near a saddle, where the surface is flatter."
        ),
    }
    return costs


def path_evaluation_counts(images: int = 7, steps: int = 200) -> dict:
    """Evaluation counts for workflow.run(stage='path'), read from its code."""
    interior = images - 2
    ceiling = 2 * steps + 2 * steps * interior
    return {
        "images": images,
        "step_limit_per_stage": steps,
        "interior_images_evaluated_per_neb_step": interior,
        "stages": {
            "endpoint_relaxations": {"count": 2, "evaluations_each_at_limit": steps},
            "preliminary_neb": {"evaluations_per_step": interior, "step_limit": steps},
            "climbing_image_neb": {"evaluations_per_step": interior, "step_limit": steps},
        },
        "ceiling_evaluations": ceiling,
        "serial": True,
        "density_reuse_between_geometries": False,
        "scenarios": {
            # Step counts are declared assumptions, not measurements. Nothing in
            # this repository has ever run a NEB on this system to calibrate them.
            "optimistic": {
                "endpoint_steps": 40,
                "pre_neb_steps": 30,
                "ci_neb_steps": 60,
                "basis": "Well-behaved FIRE relaxations from a good guess; loose pre-NEB tolerance (3x fmax) converges quickly; CI-NEB converges well inside its limit.",
            },
            "central": {
                "endpoint_steps": 80,
                "pre_neb_steps": 60,
                "ci_neb_steps": 120,
                "basis": "Ordinary behaviour for a 53-atom radical system with fixed anchors and a hand-built starting band.",
            },
            "ceiling": {
                "endpoint_steps": steps,
                "pre_neb_steps": steps,
                "ci_neb_steps": steps,
                "basis": "Every stage exhausts its step limit. This is what the code permits, not what it is expected to need.",
            },
        },
    }


def project(costs: dict, counts: dict) -> dict:
    """Wall-clock projection per method and scenario."""
    interior = counts["interior_images_evaluated_per_neb_step"]
    projection = {}
    for method in ("direct", "density_fitting"):
        per_evaluation = costs[method]["seconds_per_energy_and_gradient"]
        rows = {}
        for name, scenario in counts["scenarios"].items():
            evaluations = (
                2 * scenario["endpoint_steps"]
                + scenario["pre_neb_steps"] * interior
                + scenario["ci_neb_steps"] * interior
            )
            seconds = evaluations * per_evaluation
            rows[name] = {
                "evaluations": evaluations,
                "seconds": seconds,
                "hours": seconds / HOURS,
                "days": seconds / DAYS,
                "assumed_steps": {
                    "endpoint": scenario["endpoint_steps"],
                    "pre_neb": scenario["pre_neb_steps"],
                    "ci_neb": scenario["ci_neb_steps"],
                },
            }
        projection[method] = {
            "seconds_per_evaluation_as_archived": per_evaluation,
            "scenarios": rows,
        }
    return projection


def process_parallel_note(counts: dict, costs: dict) -> dict:
    """What separate OS processes can and cannot recover.

    NEB images are independent within one optimizer step, so the five interior
    images could in principle be evaluated in five processes. The endpoint
    relaxations are inherently serial and are unaffected.
    """
    interior = counts["interior_images_evaluated_per_neb_step"]
    cores = 8
    usable = min(interior, cores)
    central = counts["scenarios"]["central"]
    per_evaluation = costs["density_fitting"]["seconds_per_energy_and_gradient"]
    serial_band = (central["pre_neb_steps"] + central["ci_neb_steps"]) * interior * per_evaluation
    parallel_band = (central["pre_neb_steps"] + central["ci_neb_steps"]) * per_evaluation  # one image per process
    endpoints = 2 * central["endpoint_steps"] * per_evaluation
    return {
        "mechanism": "ASE's NEB evaluates images serially in one process; a parallel band needs separate OS processes and is not implemented here.",
        "interior_images": interior,
        "logical_cores": cores,
        "max_useful_processes_for_the_band": usable,
        "central_scenario_density_fitting": {
            "band_serial_hours": serial_band / HOURS,
            "band_ideally_parallel_hours": parallel_band / HOURS,
            "endpoints_serial_hours": endpoints / HOURS,
            "total_serial_hours": (serial_band + endpoints) / HOURS,
            "total_with_parallel_band_hours": (parallel_band + endpoints) / HOURS,
        },
        "caveats": (
            "Ideal scaling assumes one image per core with no memory pressure "
            "and no contention. This host has eight logical cores already "
            "shared by other lanes, and the endpoint relaxations stay serial, "
            "so the parallel figure is a floor that the current code cannot "
            "reach without being modified."
        ),
    }


def main() -> int:
    costs = measured_evaluation_costs()
    counts = path_evaluation_counts()
    projection = project(costs, counts)
    report = {
        "schema_version": 1,
        "question": "Is the 53-atom candidate reaction path computable on this hardware?",
        "computes_chemistry": False,
        "inputs": "Archived 53-atom records under data/validation/; no new calculations.",
        "system": {
            "formula": "C22H31",
            "atoms": 53,
            "basis_functions": costs["direct"]["basis_functions"],
            "method": "PBE0-D3(BJ)/def2-SVP, neutral doublet, UKS",
        },
        "measured_evaluation_costs": costs,
        "path_evaluation_counts": counts,
        "projection": projection,
        "process_parallelism": process_parallel_note(counts, costs),
        "interpretation": (
            "A cost projection says nothing about whether the reaction works. "
            "An affordable path is not evidence of a viable tool, and an "
            "unaffordable one is not evidence against it."
        ),
    }
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    out = EVIDENCE / "cost-model.json"
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    for method in ("direct", "density_fitting"):
        per = projection[method]["seconds_per_evaluation_as_archived"]
        print(f"{method}: {per:.1f} s per energy+gradient (as archived, contended)")
        for name, row in projection[method]["scenarios"].items():
            print(f"    {name:10s} {row['evaluations']:5d} evaluations  {row['hours']:8.1f} h  {row['days']:6.1f} d")
    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
