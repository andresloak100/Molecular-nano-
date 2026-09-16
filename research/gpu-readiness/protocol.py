"""Verify local protocol bytes and resolved numerical evidence without a solver."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import stat

import equivalence

POLICY = "resolved-grid-basis-float64-v1"
MAX_FILES = 128
MAX_FILE_BYTES = 16 * 1024 * 1024
MAX_TOTAL_BYTES = 64 * 1024 * 1024


def _digest(raw):
    return hashlib.sha256(raw).hexdigest()


def _parse(raw):
    def reject(value):
        raise ValueError(f"Nonfinite JSON constant: {value}")
    return json.loads(raw.decode("utf-8"), object_pairs_hook=equivalence._unique_object, parse_constant=reject)


def _relative(name):
    if (not isinstance(name, str) or not name or "\\" in name or ":" in name
            or name.startswith("/") or "\x00" in name
            or any(part in ("", ".", "..") for part in name.split("/"))
            or PurePosixPath(name).as_posix() != name):
        raise ValueError("Protocol member must be a canonical relative path")
    return name.split("/")


def _read_member(root_fd, name):
    """Open each path component relative to a pinned directory; follow no links."""
    if not hasattr(os, "O_NOFOLLOW") or os.open not in os.supports_dir_fd:
        raise ValueError("No-follow protocol verification requires POSIX descriptor-relative opens")
    components = _relative(name)
    directory = os.dup(root_fd)
    try:
        for component in components[:-1]:
            child = os.open(component, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=directory)
            os.close(directory)
            directory = child
        descriptor = os.open(components[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=directory)
        try:
            info = os.fstat(descriptor)
            if not stat.S_ISREG(info.st_mode) or info.st_size > MAX_FILE_BYTES:
                raise ValueError("Protocol member must be a bounded regular file")
            chunks, total = [], 0
            while True:
                block = os.read(descriptor, min(65536, MAX_FILE_BYTES + 1 - total))
                if not block:
                    break
                total += len(block)
                if total > MAX_FILE_BYTES:
                    raise ValueError("Protocol member exceeded byte limit while reading")
                chunks.append(block)
            return b"".join(chunks)
        finally:
            os.close(descriptor)
    finally:
        os.close(directory)


def read_explicit(path):
    path = Path(path).absolute()
    descriptor = os.open(path.parent.resolve(), os.O_RDONLY | os.O_DIRECTORY)
    try:
        return _read_member(descriptor, path.name)
    finally:
        os.close(descriptor)


def _geometry(raw, case):
    lines = raw.decode("utf-8").splitlines()
    if len(lines) < 2 or not re.fullmatch(r"[0-9]+", lines[0].strip()):
        raise ValueError("Geometry must contain one declared XYZ frame")
    count = int(lines[0])
    if count != case["atom_count"] or len(lines) < count + 2 or any(line.strip() for line in lines[count + 2:]):
        raise ValueError("Geometry atom count/frame count disagrees with plan")
    properties = [token.split("=", 1)[1] for token in shlex.split(lines[1]) if token.startswith("Properties=")]
    if properties and (len(properties) != 1 or properties[0].split(":")[:6] != ["species", "S", "1", "pos", "R", "3"]):
        raise ValueError("Unsupported extended-XYZ property order")
    symbols = []
    for line in lines[2:count + 2]:
        fields = line.split()
        if len(fields) < 4 or re.fullmatch(r"[A-Z][a-z]?", fields[0]) is None:
            raise ValueError("Malformed atom row")
        if not all(math.isfinite(float(value)) for value in fields[1:4]):
            raise ValueError("Nonfinite atomic coordinates")
        symbols.append(fields[0])
    if symbols != case.get("symbols") or case.get("length_unit") != "angstrom":
        raise ValueError("Geometry element order or declared units disagree with plan")


def read_verified_protocol(path):
    """Return a parsed plan bound to the exact bytes read, or raise explicitly.

    Original-source paths are provenance strings only and are never followed.
    All declared snapshot members are read and hashed, including source notes.
    """
    path = Path(path).absolute()
    descriptor = os.open(path.parent.resolve(), os.O_RDONLY | os.O_DIRECTORY)
    try:
        raw = _read_member(descriptor, path.name)
        plan = _parse(raw)
        if not equivalence._json_value(plan):
            raise ValueError("Invalid/nonfinite protocol JSON")
        if type(plan) is not dict or type(plan.get("schema_version")) is not int or plan["schema_version"] != 2:
            raise ValueError("Version 2 protocol required; version 1 lacks resolved-numerics requirements")
        compatible = deepcopy(plan)
        compatible["schema_version"] = 1
        errors = equivalence._plan_errors(compatible)
        if errors:
            raise ValueError("; ".join(errors))
        if plan.get("numerical_metadata_policy") != POLICY:
            raise ValueError("Explicit resolved numerical metadata policy required")
        manifest = plan.get("bundled_files_sha256")
        if type(manifest) is not dict or not 1 <= len(manifest) <= MAX_FILES:
            raise ValueError("Protocol must list a bounded nonempty file manifest")
        members, total = {}, len(raw)
        for name, expected in manifest.items():
            if type(expected) is not str or re.fullmatch(r"[0-9a-f]{64}", expected) is None:
                raise ValueError("Malformed member SHA-256")
            if name == path.name:
                raise ValueError("Protocol cannot include its own digest recursively")
            content = _read_member(descriptor, name)
            total += len(content)
            if total > MAX_TOTAL_BYTES:
                raise ValueError("Protocol exceeds total byte limit")
            if _digest(content) != expected:
                raise ValueError(f"Snapshot hash mismatch: {name}")
            members[name] = content
        for case in plan["cases"]:
            geometry, record_name = case.get("geometry"), case.get("source_record")
            if (not isinstance(geometry, str) or geometry not in members
                    or not isinstance(record_name, str) or record_name not in members):
                raise ValueError("Case references an unlisted snapshot member")
            if case["geometry_sha256"] != manifest[geometry] or case.get("source_record_sha256") != manifest[record_name]:
                raise ValueError("Case digest disagrees with manifest")
            _geometry(members[geometry], case)
            record = _parse(members[record_name])
            key = case.get("source_settings_key")
            source = record.get(key) if type(record) is dict and isinstance(key, str) else None
            if type(source) is not dict or set(source) - equivalence.SETTING_KEYS:
                raise ValueError("Source settings must be a supported recorded settings object")
            changes = {name: {"recorded": source.get(name), "recorded_key_present": name in source,
                              "planned": value} for name, value in case["settings"].items()
                       if name not in source or not equivalence._same_json(source[name], value)}
            if not equivalence._same_json(changes, case.get("planned_changes_from_record")):
                raise ValueError("Requested settings changes disagree with exact source record")
            if case.get("resolved_numerics_status") != "unavailable_before_execution" or case.get("resolved_numerics", "absent") is not None:
                raise ValueError("Planning artifact must retain explicit unavailable resolved numerics")
        return plan, _digest(raw), {"verified_files": len(members), "verified_bytes": total}
    finally:
        os.close(descriptor)


def compare_protocol_files(plan_path, cpu_path, gpu_path):
    """Strict CLI path: actual snapshots and resolved numerics precede parity."""
    from resolved_numerics import validate_resolved_numerics
    result = {"schema_version": 2, "kind": "cpu_gpu_verified_equivalence_result",
              "input_files_verified": False, "resolved_numerics_verified": False,
              "record_validation_passed": False, "numerical_parity_passed": False,
              "parity_evidence_accepted": False, "scientific_validated": False,
              "electronic_state_identity_verified": False, "ground_state_verified": False,
              "mismatches": [], "scope": "Exact snapshot and recorded numerical-configuration checks; no execution attestation or chemical validation."}
    try:
        plan, plan_hash, inventory = read_verified_protocol(plan_path)
        result.update(input_files_verified=True, plan_sha256=plan_hash, **inventory)
        cpu_raw, gpu_raw = read_explicit(cpu_path), read_explicit(gpu_path)
        cpu, gpu = _parse(cpu_raw), _parse(gpu_raw)
        result["record_sha256"] = {"cpu": _digest(cpu_raw), "gpu": _digest(gpu_raw)}
        for label, record in (("CPU", cpu), ("GPU", gpu)):
            if not equivalence._json_value(record) or type(record) is not dict:
                raise ValueError(f"{label}: invalid record JSON")
            if type(record.get("schema_version")) is not int or record["schema_version"] != 2:
                raise ValueError(f"{label}: version 2 execution record required")
            if record.get("plan_sha256") != plan_hash:
                raise ValueError(f"{label}: record is not bound to these exact plan bytes")
        case = next((item for item in plan["cases"] if item["id"] == cpu.get("case_id")), None)
        if case is None:
            raise ValueError("CPU case ID is not in protocol")
        result["case_id"] = case["id"]
        numerical_errors = validate_resolved_numerics(cpu, gpu, case)
        if numerical_errors:
            result["mismatches"].extend(numerical_errors)
            return result
        result["resolved_numerics_verified"] = True
        compatible = [deepcopy(item) for item in (plan, cpu, gpu)]
        for item in compatible:
            item["schema_version"] = 1
        numerical = equivalence.compare_records(*compatible)
        for key in ("record_validation_passed", "numerical_parity_passed", "differences", "tolerances", "mismatches"):
            if key in numerical:
                result[key] = numerical[key]
        result["parity_evidence_accepted"] = result["numerical_parity_passed"]
    except (OSError, UnicodeError, ValueError, TypeError) as error:
        result["mismatches"].append(f"Unable to verify evidence: {type(error).__name__}: {error}")
    return result
