"""Inspect GPU prerequisites or prepare exact-input parity plans. No solver runs."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    "https://pyscf.org/user/gpu.html",
    "https://github.com/pyscf/gpu4pyscf",
    "https://github.com/pyscf/gpu4pyscf/blob/master/requirements.txt",
]
TOLERANCES = {"energy_hartree": 1e-6, "max_force_component_ev_per_angstrom": 1e-4, "s2": 1e-5}


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def json_text(value):
    return json.dumps(value, indent=2, allow_nan=False) + "\n"


def write_new_json(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        stream.write(json_text(value))


def package_versions():
    names = ("pyscf", "gpu4pyscf", "cupy", "cutensor", "numpy", "scipy", "ase", "dftd3")
    found = {}
    for distribution in metadata.distributions():
        name = distribution.metadata.get("Name", "").lower().replace("_", "-")
        if any(name == prefix or name.startswith(prefix + "-") for prefix in names):
            found[name] = distribution.version
    return found


def run_command(args, timeout=15):
    try:
        result = subprocess.run(args, capture_output=True, text=True, timeout=timeout, check=False)
        return {"returncode": result.returncode, "stdout": result.stdout.strip(), "stderr": result.stderr.strip()}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"returncode": None, "error": f"{type(error).__name__}: {error}"}


def memory_bytes():
    if platform.system() == "Darwin":
        result = run_command(["sysctl", "-n", "hw.memsize"])
        if result.get("returncode") == 0:
            return int(result["stdout"])
    try:
        return int(os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE"))
    except (ValueError, OSError, AttributeError):
        return None


RUNTIME_PROBE = r'''
import json
report = {}
try:
    from pyscf import lib
    previous = lib.num_threads()
    try:
        lib.num_threads(2)
        effective = int(lib.num_threads())
        report['pyscf_thread_probe'] = {'requested': 2, 'effective': effective, 'honored': effective == 2}
    finally:
        lib.num_threads(previous)
except Exception as error:
    report['pyscf_thread_probe'] = {'error': type(error).__name__ + ': ' + str(error)}
try:
    import cupy as cp
    import gpu4pyscf
    devices = []
    for index in range(cp.cuda.runtime.getDeviceCount()):
        prop = cp.cuda.runtime.getDeviceProperties(index)
        name = prop['name']
        with cp.cuda.Device(index):
            free_memory, total_memory = cp.cuda.runtime.memGetInfo()
        devices.append({'index': index, 'name': name.decode() if isinstance(name, bytes) else str(name),
                        'compute_capability': [int(prop['major']), int(prop['minor'])],
                        'free_memory_bytes': int(free_memory), 'total_memory_bytes': int(total_memory)})
    report['cuda'] = {'status': 'available' if devices else 'no_visible_devices', 'devices': devices,
                      'runtime_version': int(cp.cuda.runtime.runtimeGetVersion()),
                      'driver_version': int(cp.cuda.runtime.driverGetVersion())}
except Exception as error:
    report['cuda'] = {'status': 'unavailable', 'error': type(error).__name__ + ': ' + str(error)}
print('G1_JSON=' + json.dumps(report, allow_nan=False))
'''


def environment_report(*, probe_runtime=False):
    versions = package_versions()
    smi = shutil.which("nvidia-smi")
    if smi:
        nvidia = run_command([smi, "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"])
    else:
        nvidia = {"status": "command_not_found", "interpretation": "No nvidia-smi in PATH; this alone does not establish GPU absence."}
    runtime = {"status": "not_requested", "cuda": {"status": "not_probed"}}
    if probe_runtime:
        child = run_command([sys.executable, "-c", RUNTIME_PROBE], timeout=20)
        lines = child.get("stdout", "").splitlines()
        payload = next((line[len("G1_JSON="):] for line in reversed(lines) if line.startswith("G1_JSON=")), None)
        if payload is not None and child.get("returncode") == 0:
            runtime = {"status": "completed", **json.loads(payload)}
        else:
            runtime = {"status": "failed", "probe": child, "cuda": {"status": "not_established"}}
    devices = runtime.get("cuda", {}).get("devices", [])
    supported_device = any(tuple(device["compute_capability"]) >= (7, 0) for device in devices)
    reasons = []
    for prefix in ("pyscf", "gpu4pyscf", "cupy", "numpy", "ase", "dftd3"):
        if not any(name == prefix or name.startswith(prefix + "-cuda") for name in versions):
            reasons.append(f"{prefix} distribution not recorded in this environment")
    if runtime.get("cuda", {}).get("status") != "available":
        reasons.append("working CUDA runtime with GPU4PySCF imports not established")
    elif not supported_device:
        reasons.append("no visible device meets the documented compiled-package compute capability >=7.0")
    for prefix in ("cupy", "gpu4pyscf"):
        variants = [name for name in versions if name == prefix or name.startswith(prefix + "-cuda")]
        if len(variants) > 1:
            reasons.append(f"multiple {prefix} distributions installed; isolate a compatible environment")
    return {
        "schema_version": 1, "kind": "gpu_prerequisite_inventory",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "host": {"system": platform.system(), "machine": platform.machine(),
                 "processor": platform.processor(), "logical_cpus": os.cpu_count(),
                 "memory_bytes": memory_bytes(), "python": platform.python_version()},
        "packages": versions, "nvidia_smi": nvidia, "runtime_probe": runtime,
        "prerequisites_observed": not reasons, "remaining_prerequisites": reasons,
        "production_gpu_adapter_implemented": False, "gpu_numerical_parity_measured": False,
        "gpu_speedup_measured": False, "design_validated": False,
        "quantum_evaluations": 0,
        "interpretation": "Inventory/import/device visibility only; exact method/gradient compatibility and memory fit still require the parity protocol. No production GPU adapter is supplied.",
        "sources": SOURCES,
    }


def read_inside(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Protocol input path escapes the selected repository")
    return path.read_bytes()


def build_plan(repo=ROOT, *, include_candidate=False):
    """Read exact recorded coordinates and state explicitly every planned change."""
    from io import StringIO
    from ase.io import read
    sys.path.insert(0, str(ROOT))
    from nanodesign.quantum import QuantumSettings

    root = Path(repo).resolve()
    cases = [
        ("h2-rks-direct", "data/validation/h2-integration/inputs/h2.xyz", "data/validation/h2-integration/inputs/design.json", "quantum", {"density_fit": False, "dispersion": None}),
        ("methane-rks-df", "data/validation/pbe0-svp/inputs/methane.xyz", "data/validation/pbe0-svp/species/methane.json", "quantum_settings", {"density_fit": True}),
        ("ethynyl-uks-df-minao", "data/validation/pbe0-svp/inputs/ethynyl_radical.xyz", "data/validation/pbe0-svp/species/ethynyl_radical.json", "quantum_settings", {"density_fit": True, "scf_initial_guess": "minao"}),
        ("ethynyl-uks-df-atom", "data/validation/pbe0-svp/inputs/ethynyl_radical.xyz", "data/validation/pbe0-svp/species/ethynyl_radical.json", "quantum_settings", {"density_fit": True, "scf_initial_guess": "atom"}),
        ("nominal-ts-uks-df", "data/validation/pbe0-svp/inputs/methane_ethynyl_ts.xyz", "data/validation/pbe0-svp/species/methane_ethynyl_ts.json", "quantum_settings", {"density_fit": True}),
    ]
    if include_candidate:
        cases.append(("candidate-53-uks-df", "data/validation/h-abstraction-df-initial/input-initial.extxyz", "data/validation/h-abstraction-df-initial/result.json", "quantum_settings", {"density_fit": True}))
    result = {
        "schema_version": 1, "kind": "cpu_gpu_equivalence_protocol", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "planned_not_executed", "quantum_evaluations_completed": 0,
        "planned_energy_force_evaluations": 2 * len(cases),
        "execution": {"order": "serial, one CPU and one GPU evaluation per case", "automatic_retries": 0,
                      "endpoint_relaxations": 0, "path_searches": 0,
                      "timing_scope": "Cold paired checks include initialization; record synchronized total energy-plus-force wall time and load. No speed guarantee."},
        "tolerances": dict(TOLERANCES),
        "tolerance_scope": "Provisional engineering parity thresholds, not chemical-accuracy targets; investigate mismatches, never loosen to obtain a pass without review.",
        "production_gpu_adapter_implemented": False, "design_validated": False,
        "electronic_state_identity_verified": False, "sources": SOURCES, "cases": [],
        "reference_data_license": "Published reference coordinates retain CC BY-NC 4.0; bundled source notes identify origin and conversion. Code license does not relicense those coordinates.",
    }
    payloads = {}
    for identity, coordinate_path, record_path, key, overrides in cases:
        raw = read_inside(root, coordinate_path)
        record_raw = read_inside(root, record_path)
        source_record = json.loads(record_raw)
        source_settings = source_record[key]
        effective = QuantumSettings(**{**source_settings, **overrides}).to_dict()
        frames = read(StringIO(raw.decode()), format="extxyz", index=":")
        if len(frames) != 1:
            raise ValueError("Each protocol case must contain one geometry")
        atoms = frames[0]
        if not all(math.isfinite(float(value)) for row in atoms.positions for value in row):
            raise ValueError("Nonfinite protocol coordinates")
        target = f"inputs/{identity}{Path(coordinate_path).suffix}"
        metadata_target = f"source-records/{identity}.json"
        payloads[target] = raw
        payloads[metadata_target] = record_raw
        changes = {name: {"recorded": source_settings.get(name), "recorded_key_present": name in source_settings,
                          "planned": value} for name, value in effective.items() if name not in source_settings or source_settings[name] != value}
        item = {
            "id": identity, "geometry": target, "geometry_sha256": digest(raw),
            "coordinate_source": coordinate_path, "source_record": metadata_target,
            "source_record_original": record_path, "source_record_sha256": digest(record_raw),
            "source_settings_key": key, "settings": effective, "planned_changes_from_record": changes,
            "atom_count": len(atoms), "symbols": atoms.get_chemical_symbols(), "formula": atoms.get_chemical_formula(),
            "length_unit": "angstrom", "force_scope": "raw forces on every atom, including any anchors; no constraint zeroing",
            "expected_reference": "RKS" if effective["spin"] == 0 else "UKS",
            "scientific_status": "fixed geometry only; no verified saddle, path or assembly operation",
        }
        if identity.startswith("candidate"):
            item["fixed_indices_context"] = source_record["design"]["fixed_indices"]
            item["hydrogen_transfer_context"] = source_record["design"]["hydrogen_transfer"]
        result["cases"].append(item)
    for filename in ("README.md", "provenance.json"):
        source = "data/validation/pbe0-svp/inputs/" + filename
        target = "source-notes/reference-" + filename
        payloads[target] = read_inside(root, source)
    result["bundled_files_sha256"] = {name: digest(raw) for name, raw in payloads.items()}
    return result, payloads


def save_plan(output, plan, payloads):
    root = Path(output)
    root.mkdir(parents=True, exist_ok=False)
    for relative, raw in payloads.items():
        target = root / relative
        target.parent.mkdir(exist_ok=True, parents=True)
        with target.open("xb") as stream:
            stream.write(raw)
    write_new_json(root / "plan.json", plan)
    return root / "plan.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    inspect = commands.add_parser("inspect", help="Inventory libraries/devices; never calculate energies")
    inspect.add_argument("--probe-runtime", action="store_true", help="Import CPU/GPU libraries and query CUDA in a timeout-bounded child process")
    inspect.add_argument("--output", type=Path, help="New JSON file; omission prints to stdout")
    plan = commands.add_parser("plan", help="Copy real archived inputs into a new unexecuted protocol bundle")
    plan.add_argument("--repo", type=Path, default=ROOT)
    plan.add_argument("--output", type=Path, required=True)
    plan.add_argument("--include-candidate", action="store_true", help="Add the expensive 53-atom case; default only plans 2–8 atom checks")
    args = parser.parse_args()
    if args.command == "inspect":
        result = environment_report(probe_runtime=args.probe_runtime)
        if args.output:
            write_new_json(args.output, result)
        else:
            print(json_text(result), end="")
    else:
        result, payloads = build_plan(args.repo, include_candidate=args.include_candidate)
        path = save_plan(args.output, result, payloads)
        print(json_text({"plan": str(path), "status": result["status"], "planned_energy_force_evaluations": result["planned_energy_force_evaluations"]}), end="")


if __name__ == "__main__":
    main()
