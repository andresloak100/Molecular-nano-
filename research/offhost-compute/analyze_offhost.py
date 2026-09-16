"""Compare off-host 53-atom single points with the archived Mac records.

Read-only against the archive. Emits offhost-comparison.json with per-run
energies, forces, timings, cross-platform deltas, and a rerun of the
compute-planning projection arithmetic using the off-host per-evaluation
cost. The evaluation count is taken from research/compute-planning
(2,407 energy-and-force evaluations for the default 7-image, 200-step,
four-stage scenario); this script does not re-derive it.
"""

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

EVALUATIONS_DEFAULT_SCENARIO = 2407  # from research/compute-planning/README.md

PAIRS = {
    "density_fitting": {
        "archived": ROOT / "data/validation/h-abstraction-df-initial/result.json",
        "offhost": HERE / "runs/53atom-df-4t/result.json",
    },
    "direct": {
        "archived": ROOT / "data/validation/h-abstraction-direct-initial/result.json",
        "offhost": HERE / "runs/53atom-direct-4t/result.json",
    },
}


def summarize(path):
    data = json.loads(Path(path).read_text())
    if data.get("status") != "completed":
        raise SystemExit(f"{path}: status is {data.get('status')!r}, not completed")
    structure = data["structure"]
    return {
        "path": str(Path(path).relative_to(ROOT)),
        "energy_ev": structure["energy_ev"],
        "free_force_max_ev_per_angstrom": structure["free_force_max_ev_per_angstrom"],
        "elapsed_seconds": data["elapsed_seconds"],
        "requested_threads": data["quantum_settings"]["threads"],
        "diagnostics": {
            k: v
            for k, v in structure.items()
            if "thread" in k or k in ("scf_converged", "s_squared", "basis_functions")
        },
    }


def main():
    out = {
        "generated_by": "research/offhost-compute/analyze_offhost.py",
        "offhost_platform": {
            "os": os.uname().sysname,
            "machine": os.uname().machine,
            "logical_cpus": os.cpu_count(),
            "loadavg_at_analysis": os.getloadavg(),
        },
        "interpretation_limits": [
            "Different hardware, OS, arch and BLAS: timing ratios mix machine and contention effects and are not a controlled benchmark of either machine.",
            "The archived Mac runs executed under measured heavy contention (see CLAUDE.md); the off-host runs executed serially on an otherwise idle container, load averages logged in timing-log.txt.",
            "Energy/force agreement at one unrelaxed geometry does not bound errors along a reaction path or validate the electronic method.",
            "Projection scenarios reuse the compute-planning evaluation count (2,407) and inherit all of its stated assumptions; they are conditional arithmetic, not schedules.",
        ],
        "runs": {},
        "cross_platform_deltas": {},
        "projection_scenarios_days": {},
    }

    for label, pair in PAIRS.items():
        archived = summarize(pair["archived"])
        offhost = summarize(pair["offhost"])
        out["runs"][label] = {"archived_mac": archived, "offhost_container": offhost}
        out["cross_platform_deltas"][label] = {
            "energy_delta_ev_offhost_minus_mac": offhost["energy_ev"] - archived["energy_ev"],
            "free_force_max_delta_ev_per_angstrom": (
                offhost["free_force_max_ev_per_angstrom"]
                - archived["free_force_max_ev_per_angstrom"]
            ),
            "elapsed_ratio_mac_over_offhost": archived["elapsed_seconds"] / offhost["elapsed_seconds"],
        }
        for source, run in (("mac_archived", archived), ("offhost_container", offhost)):
            days = EVALUATIONS_DEFAULT_SCENARIO * run["elapsed_seconds"] / 86400.0
            out["projection_scenarios_days"][f"{label}_{source}"] = round(days, 2)

    dest = HERE / "offhost-comparison.json"
    dest.write_text(json.dumps(out, indent=1))
    print(json.dumps(out["cross_platform_deltas"], indent=1))
    print(json.dumps(out["projection_scenarios_days"], indent=1))
    print(f"written: {dest.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
