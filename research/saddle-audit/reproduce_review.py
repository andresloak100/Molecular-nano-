"""Cheap S1 review reproductions; synthetic mode probes are NOT computed modes."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read
from mode_evidence import transfer_projection

ROOT = Path(__file__).resolve().parents[2]


def source_hash(path):
    return {"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def load_s1():
    sys.path.insert(0, str(ROOT))
    old = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec = importlib.util.spec_from_file_location("s1_read_only_review", ROOT / "research/reference-saddle/saddle_search.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.dont_write_bytecode = old


def mode_probe(module):
    path = ROOT / "data/reference/methane_ethynyl_ts.xyz"
    atoms = read(path)
    masses = atoms.get_masses()
    # A proper rotation with no dependence on the molecule's axis orientation.
    angle = .37
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    outputs = []
    for label, vector in (("transverse", [0., 1., 0.]),
                          ("transverse_with_axial_noise", [1e-14, 1., 0.]),
                          ("axial", [1., 0., 0.])):
        for variant, sign, matrix in (("original", 1, np.eye(3)),
                                      ("sign_reversed", -1, np.eye(3)),
                                      ("rotated", 1, rotation)):
            mode = np.zeros((len(atoms), 3))
            mode[3] = np.array(vector) * sign
            mode /= np.sqrt(np.sum(masses[:, None] * mode**2))
            geometry = atoms.positions @ matrix.T
            mode = mode @ matrix.T
            fixture = {"free_atom_indices": list(range(len(atoms))),
                       "frequencies_cm1": [-100.],
                       "cartesian_modes_per_sqrt_amu": [mode.tolist()],
                       "free_masses_amu": masses.tolist(),
                       "settings": {"frequency_tolerance_cm1": 20.},
                       "geometry_angstrom": geometry.tolist()}
            result = module.analyze_transfer_mode(fixture, atoms.get_chemical_symbols())
            independent = transfer_projection(fixture, {"acceptor_carbon": 2, "transferring_hydrogen": 3, "donor_carbon": 4}, atoms.get_chemical_symbols())
            outputs.append({"synthetic_direction": label, "transform": variant,
                            "reported_transfer_like": result["single_transfer_like_mode"],
                            "reported_mode_details": result["negative_modes"],
                            "independent_projection": independent})
    return {"warning": "Synthetic direction probes at real seed coordinates. Frequencies and displacements here are test inputs, never quantum results.",
            "geometry_source": source_hash(path), "cases": outputs}


def evaluation_counts():
    results = []
    for path in sorted((ROOT / "research/reference-saddle").rglob("optimization.json")):
        data = json.loads(path.read_text())
        log_path = path.parent / "electronic.jsonl"
        if not log_path.exists():
            continue
        completed = set()
        incomplete_lines = 0
        for line in log_path.read_text().splitlines():
            try:
                item = json.loads(line)
            except ValueError:
                incomplete_lines += 1
                continue
            if item.get("event") == "calculation_completed" and item.get("call_id"):
                completed.add(item["call_id"])
        results.append({"optimization_source": source_hash(path), "events_source": source_hash(log_path),
                        "optimization_status": data.get("status"),
                        "reported_electronic_evaluations": data.get("electronic_evaluations"),
                        "unique_completed_calculator_calls": len(completed),
                        "unreadable_event_lines": incomplete_lines})
    return results


def main():
    module = load_s1()
    result = {"quantum_jobs_started": 0,
              "reviewed_source": source_hash(ROOT / "research/reference-saddle/saddle_search.py"),
              "independent_projection_source": source_hash(Path(__file__).with_name("mode_evidence.py")),
              "mode_heuristic_probe": mode_probe(module),
              "saved_evaluation_counts": evaluation_counts()}
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
