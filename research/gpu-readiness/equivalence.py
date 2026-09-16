"""Compare recorded CPU/GPU numerical evidence; never run quantum calculations.

Passing this checker means the supplied records satisfy a numerical protocol.
It does not independently authenticate execution or validate electronic states,
chemical accuracy, reaction connectivity, or a molecular-machine design.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import sys


SETTING_KEYS = {
    "charge", "spin", "xc", "basis", "dispersion", "grid_level", "conv_tol",
    "max_cycle", "threads", "memory_mb", "density_fit", "scf_initial_guess",
}
TOLERANCE_KEYS = {"energy_hartree", "max_force_component_ev_per_angstrom", "s2"}
GUESSES = {"minao", "atom", "1e", "huckel"}


def _number(value):
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(value)
    except (OverflowError, ValueError):
        return False


def _json_value(value):
    """Reject non-JSON objects/nonfinite values even for direct Python callers."""
    if value is None or type(value) in (str, bool):
        return True
    if type(value) in (int, float):
        return _number(value)
    if type(value) is list:
        return all(_json_value(item) for item in value)
    if type(value) is dict:
        return all(type(key) is str and _json_value(item) for key, item in value.items())
    return False


def _same_json(left, right):
    # Python equality otherwise considers True == 1 and 1 == 1.0.
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(right, sort_keys=True, allow_nan=False)


def _settings_errors(settings, label):
    errors = []
    if type(settings) is not dict or set(settings) != SETTING_KEYS:
        return [f"{label}: settings must contain exactly the complete QuantumSettings fields"]
    for key in ("charge", "spin", "grid_level", "max_cycle", "threads", "memory_mb"):
        value = settings[key]
        if type(value) is not int:
            errors.append(f"{label}: settings.{key} must be an integer")
        elif key == "spin" and value < 0:
            errors.append(f"{label}: settings.spin must be nonnegative")
        elif key == "grid_level" and not 0 <= value <= 9:
            errors.append(f"{label}: settings.grid_level must be between 0 and 9")
        elif key in ("max_cycle", "threads", "memory_mb") and value <= 0:
            errors.append(f"{label}: settings.{key} must be positive")
    for key in ("xc", "basis"):
        if type(settings[key]) is not str or not settings[key].strip():
            errors.append(f"{label}: settings.{key} must be a nonempty string")
    if not _number(settings["conv_tol"]) or not 0 < settings["conv_tol"] < 1:
        errors.append(f"{label}: settings.conv_tol must be finite and between 0 and 1")
    if settings["dispersion"] not in (None, "d3bj", "d3zero"):
        errors.append(f"{label}: settings.dispersion is unsupported")
    if type(settings["density_fit"]) is not bool:
        errors.append(f"{label}: settings.density_fit must be a boolean")
    if type(settings["scf_initial_guess"]) is not str or settings["scf_initial_guess"] not in GUESSES:
        errors.append(f"{label}: settings.scf_initial_guess is unsupported")
    return errors


def _plan_errors(plan):
    if type(plan) is not dict:
        return ["plan: expected an object"]
    errors = []
    if type(plan.get("schema_version")) is not int or plan["schema_version"] != 1:
        errors.append("plan: schema_version must be integer 1")
    if plan.get("kind") != "cpu_gpu_equivalence_protocol":
        errors.append("plan: kind must be cpu_gpu_equivalence_protocol")
    tolerances = plan.get("tolerances")
    if type(tolerances) is not dict or set(tolerances) != TOLERANCE_KEYS:
        errors.append("plan: all three named tolerances must be explicitly provided")
    else:
        for key, value in tolerances.items():
            if not _number(value) or value < 0:
                errors.append(f"plan: tolerance {key} must be finite and nonnegative")
    cases = plan.get("cases")
    if type(cases) is not list or not cases:
        errors.append("plan: cases must be a nonempty list")
        return errors
    ids = set()
    for index, case in enumerate(cases):
        label = f"plan case {index}"
        if type(case) is not dict:
            errors.append(f"{label}: expected an object")
            continue
        case_id = case.get("id")
        if type(case_id) is not str or not case_id.strip():
            errors.append(f"{label}: id must be a nonempty string")
        elif case_id in ids:
            errors.append(f"{label}: duplicate case id {case_id}")
        else:
            ids.add(case_id)
        digest = case.get("geometry_sha256")
        if type(digest) is not str or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            errors.append(f"{label}: geometry_sha256 must be a lowercase SHA-256 digest")
        if type(case.get("atom_count")) is not int or case["atom_count"] <= 0:
            errors.append(f"{label}: atom_count must be a positive integer")
        errors.extend(_settings_errors(case.get("settings"), label))
    return errors


def _record_errors(record, case, gpu):
    label = "GPU record" if gpu else "CPU record"
    if type(record) is not dict:
        return [f"{label}: expected an object"]
    errors = []
    if type(record.get("schema_version")) is not int or record["schema_version"] != 1:
        errors.append(f"{label}: schema_version must be integer 1")
    for key, expected in (
        ("case_id", case["id"]),
        ("geometry_sha256", case["geometry_sha256"]),
        ("backend", "gpu4pyscf" if gpu else "pyscf_cpu"),
        ("execution_device", "nvidia_cuda" if gpu else "cpu"),
        ("scf_initial_guess", case["settings"]["scf_initial_guess"]),
    ):
        if record.get(key) != expected:
            errors.append(f"{label}: {key} must match {expected!r}")
    if not _same_json(record.get("settings"), case["settings"]):
        errors.append(f"{label}: settings differ from the complete planned settings")
    for key in ("scf_converged", "gradient_completed", "grid_response"):
        if record.get(key) is not True:
            errors.append(f"{label}: {key} must be true")
    if record.get("auxiliary_basis_response") is not case["settings"]["density_fit"]:
        errors.append(f"{label}: auxiliary_basis_response must match density_fit")
    expected_d3 = 0 if case["settings"]["dispersion"] is None else 1
    if type(record.get("dispersion_evaluations")) is not int or record["dispersion_evaluations"] != expected_d3:
        errors.append(f"{label}: dispersion_evaluations must be {expected_d3}")
    for key in ("energy_hartree", "s2"):
        if not _number(record.get(key)):
            errors.append(f"{label}: {key} must be a finite number")
    forces = record.get("forces_ev_per_angstrom")
    if (type(forces) is not list or len(forces) != case["atom_count"]
            or any(type(row) is not list or len(row) != 3 or not all(_number(value) for value in row) for row in forces)):
        errors.append(f"{label}: forces_ev_per_angstrom must be a finite atom_count x 3 array")
    versions = record.get("versions")
    required_versions = {"pyscf", "gpu4pyscf", "cupy"} if gpu else {"pyscf"}
    if (type(versions) is not dict or not required_versions.issubset(versions)
            or any(type(value) is not str or not value.strip() for value in versions.values())):
        errors.append(f"{label}: versions must include nonempty strings for {sorted(required_versions)}")
    if gpu and record.get("gpu_synchronized") is not True:
        errors.append(f"{label}: gpu_synchronized must be true")
    return errors


def compare_records(plan, cpu, gpu):
    """Return conservative results and explicit mismatches for one planned case."""
    result = {
        "schema_version": 1,
        "kind": "cpu_gpu_equivalence_result",
        "case_id": None,
        "record_validation_passed": False,
        "numerical_parity_passed": False,
        "input_files_verified": False,
        "resolved_numerics_verified": False,
        "parity_evidence_accepted": False,
        "scientific_validated": False,
        "electronic_state_identity_verified": False,
        "ground_state_verified": False,
        "differences": {key: None for key in sorted(TOLERANCE_KEYS)},
        "mismatches": [],
        "scope": "Provisional numerical acceptance of supplied evidence; no independent execution authentication, chemical-accuracy, electronic-state, or speedup claim.",
    }
    errors = result["mismatches"]
    for label, value in (("plan", plan), ("CPU record", cpu), ("GPU record", gpu)):
        if not _json_value(value):
            errors.append(f"{label}: malformed JSON value or nonfinite number")
    if errors:
        return result
    errors.extend(_plan_errors(plan))
    if errors:
        return result
    result["tolerances"] = dict(plan["tolerances"])
    case_id = cpu.get("case_id") if type(cpu) is dict else None
    case = next((entry for entry in plan["cases"] if entry["id"] == case_id), None)
    if case is None:
        errors.append("CPU record: case_id does not identify a planned case")
        return result
    result["case_id"] = case_id
    errors.extend(_record_errors(cpu, case, False))
    errors.extend(_record_errors(gpu, case, True))
    if errors:
        return result
    differences = {
        "energy_hartree": abs(cpu["energy_hartree"] - gpu["energy_hartree"]),
        "max_force_component_ev_per_angstrom": max(
            abs(left - right)
            for cpu_row, gpu_row in zip(cpu["forces_ev_per_angstrom"], gpu["forces_ev_per_angstrom"])
            for left, right in zip(cpu_row, gpu_row)
        ),
        "s2": abs(cpu["s2"] - gpu["s2"]),
    }
    if not all(_number(value) for value in differences.values()):
        errors.append("comparison: finite inputs overflowed while computing differences")
        return result
    result["record_validation_passed"] = True
    result["differences"] = differences
    for key, difference in differences.items():
        if difference > plan["tolerances"][key]:
            errors.append(f"{key}: absolute difference {difference:.12g} exceeds tolerance {plan['tolerances'][key]:.12g}")
    result["numerical_parity_passed"] = not errors
    return result


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path):
    def reject_constant(value):
        raise ValueError(f"Nonfinite JSON constant: {value}")
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=_unique_object, parse_constant=reject_constant)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    compare = commands.add_parser("compare", help="Compare existing records only; no calculations are launched")
    compare.add_argument("--plan", required=True)
    compare.add_argument("--cpu", required=True)
    compare.add_argument("--gpu", required=True)
    compare.add_argument("--records-only", action="store_true",
                         help="Historical v1 numerical comparison only; cannot accept parity evidence without actual files/resolved numerics")
    args = parser.parse_args(argv)
    try:
        if args.records_only:
            result = compare_records(read_json(args.plan), read_json(args.cpu), read_json(args.gpu))
        else:
            from protocol import compare_protocol_files
            result = compare_protocol_files(args.plan, args.cpu, args.gpu)
    except (OSError, UnicodeError, ValueError) as error:
        result = {"schema_version": 1, "kind": "cpu_gpu_equivalence_result",
                  "numerical_parity_passed": False, "scientific_validated": False,
                  "mismatches": [f"Unable to read evidence: {error}"]}
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["numerical_parity_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
