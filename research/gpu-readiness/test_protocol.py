"""Actual snapshot-boundary tests; no quantum results are computed."""
import hashlib
import json

import pytest

import gpu_readiness
import protocol


@pytest.fixture
def bundle(tmp_path):
    plan, payloads = gpu_readiness.build_plan(protocol_version=2)
    path = gpu_readiness.save_plan(tmp_path / "protocol", plan, payloads)
    return path, plan


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*") if p.is_file()}


def rehash_case_member(path, plan, key):
    case = plan["cases"][0]
    relative = case[key]
    checksum = hashlib.sha256((path.parent / relative).read_bytes()).hexdigest()
    plan["bundled_files_sha256"][relative] = checksum
    case[key + "_sha256"] = checksum
    save(path, plan)


def test_real_protocol_is_verified_without_resolving_unknown_numerics(bundle):
    path, original = bundle
    before = hashes(path.parent)
    plan, checksum, inventory = protocol.read_verified_protocol(path)
    assert plan == original
    assert checksum == before["plan.json"]
    assert inventory["verified_files"] == len(plan["bundled_files_sha256"])
    assert all(case["resolved_numerics"] is None for case in plan["cases"])
    assert before == hashes(path.parent)


@pytest.mark.parametrize("kind", ["geometry", "source_record", "source_note"])
def test_changed_declared_bytes_are_rejected(bundle, kind):
    path, plan = bundle
    relative = ("source-notes/reference-README.md" if kind == "source_note" else plan["cases"][0][kind])
    member = path.parent / relative
    member.write_bytes(member.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="Snapshot hash mismatch"):
        protocol.read_verified_protocol(path)


@pytest.mark.parametrize("kind", ["geometry", "source_record"])
def test_missing_member_is_not_a_declared_hash_pass(bundle, kind):
    path, plan = bundle
    (path.parent / plan["cases"][0][kind]).unlink()
    with pytest.raises(OSError):
        protocol.read_verified_protocol(path)


@pytest.mark.parametrize("where", ["plan", "member", "directory"])
def test_no_follow_boundary_rejects_symlinks_even_with_matching_contents(bundle, tmp_path, where):
    path, plan = bundle
    if where == "plan":
        copy = tmp_path / "external-plan.json"
        path.rename(copy)
        path.symlink_to(copy)
    elif where == "member":
        target = path.parent / plan["cases"][0]["geometry"]
        copy = tmp_path / "external.xyz"
        target.rename(copy)
        target.symlink_to(copy)
    else:
        target = path.parent / "inputs"
        copy = tmp_path / "external-inputs"
        target.rename(copy)
        target.symlink_to(copy, target_is_directory=True)
    with pytest.raises(OSError):
        protocol.read_verified_protocol(path)


@pytest.mark.parametrize("name", ["../escape", "/tmp/absolute", "a//b", "a/./b", "a/../b", "C:drive", "a\\b"])
def test_manifest_paths_are_canonical_and_relative(bundle, name):
    path, plan = bundle
    plan["bundled_files_sha256"][name] = "a" * 64
    save(path, plan)
    with pytest.raises(ValueError, match="canonical relative"):
        protocol.read_verified_protocol(path)


@pytest.mark.parametrize("change", ["unlisted", "wrong_case_digest", "unknown_source_key", "source_change_omitted",
                                   "resolved_invented", "resolved_missing", "policy_missing"])
def test_plan_metadata_must_match_the_snapshot_and_stay_honest(bundle, change):
    path, plan = bundle
    case = plan["cases"][0]
    if change == "unlisted":
        del plan["bundled_files_sha256"][case["source_record"]]
    elif change == "wrong_case_digest":
        case["geometry_sha256"] = "b" * 64
    elif change == "unknown_source_key":
        case["source_settings_key"] = "invented"
    elif change == "source_change_omitted":
        case["settings"]["basis"] = "def2-tzvp"
    elif change == "resolved_invented":
        case["resolved_numerics"] = {"precision": "float64"}
    elif change == "resolved_missing":
        del case["resolved_numerics"]
    else:
        del plan["numerical_metadata_policy"]
    save(path, plan)
    with pytest.raises(ValueError):
        protocol.read_verified_protocol(path)


@pytest.mark.parametrize("change", ["wrong_symbols", "nan_position", "two_frames", "wrong_property_order"])
def test_even_self_consistent_hashes_do_not_hide_invalid_geometry(bundle, change):
    path, plan = bundle
    case = plan["cases"][0]
    member = path.parent / case["geometry"]
    rows = member.read_text().splitlines()
    if change == "wrong_symbols":
        case["symbols"][0] = "C"
    elif change == "nan_position":
        fields = rows[2].split()
        fields[1] = "NaN"
        rows[2] = " ".join(fields)
    elif change == "two_frames":
        rows.extend(rows[:])
    else:
        rows[1] = "Properties=pos:R:3:species:S:1"
    member.write_text("\n".join(rows) + "\n")
    rehash_case_member(path, plan, "geometry")
    with pytest.raises(ValueError):
        protocol.read_verified_protocol(path)


def test_source_original_location_is_provenance_only_not_followed(bundle):
    path, plan = bundle
    for case in plan["cases"]:
        case["source_record_original"] = "/unavailable/original/computer/source.json"
        case["coordinate_source"] = "/unavailable/original/computer/source.xyz"
    save(path, plan)
    assert protocol.read_verified_protocol(path)[0] == plan


def test_v1_remains_preserved_but_cannot_pass_v2_verification(tmp_path):
    plan, payloads = gpu_readiness.build_plan(protocol_version=1)
    path = gpu_readiness.save_plan(tmp_path / "legacy", plan, payloads)
    before = hashes(path.parent)
    with pytest.raises(ValueError, match="Version 2"):
        protocol.read_verified_protocol(path)
    assert before == hashes(path.parent)


def test_file_count_limit_precedes_member_reads(bundle, monkeypatch):
    path, plan = bundle
    monkeypatch.setattr(protocol, "MAX_FILES", 1)
    with pytest.raises(ValueError, match="bounded nonempty"):
        protocol.read_verified_protocol(path)


def test_byte_limit_rejects_oversized_input_without_loading_it(bundle, monkeypatch):
    path, _ = bundle
    monkeypatch.setattr(protocol, "MAX_FILE_BYTES", 10)
    with pytest.raises(ValueError, match="bounded regular"):
        protocol.read_verified_protocol(path)


def test_total_bytes_are_bounded(bundle, monkeypatch):
    path, _ = bundle
    monkeypatch.setattr(protocol, "MAX_TOTAL_BYTES", path.stat().st_size + 1)
    with pytest.raises(ValueError, match="total byte"):
        protocol.read_verified_protocol(path)


def test_fifo_member_fails_instead_of_blocking(bundle):
    import os
    path, plan = bundle
    member = path.parent / plan["cases"][0]["geometry"]
    member.unlink()
    os.mkfifo(member)
    with pytest.raises(ValueError, match="regular file"):
        protocol.read_verified_protocol(path)


def synthetic_records(path, plan):
    """Fixtures deliberately do not represent an actual grid, energy or GPU run."""
    from copy import deepcopy
    from resolved_numerics import canonical_grid_signature, GRID_CANONICALIZATION
    case = plan["cases"][0]
    metadata = {
        "precision": {"energy": "float64", "forces": "float64", "density": "float64"},
        "orbital_basis": {"sha256": "c" * 64},
        "auxiliary_basis": {"status": "not_applicable"},
        "quadrature": {"level": case["settings"]["grid_level"],
            "atom_grid": {element: {"radial": 1, "angular": 1} for element in case["symbols"]},
            "pruning": "synthetic-pruning", "radial_method": "synthetic-radial",
            "becke_scheme": "synthetic-becke", "point_count": 1,
            "canonicalization": GRID_CANONICALIZATION,
            "signature_sha256": canonical_grid_signature([[0, 0, 0, 1]])},
    }
    cpu = {"schema_version": 2, "synthetic_fixture": True,
        "plan_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "case_id": case["id"], "geometry_sha256": case["geometry_sha256"],
        "settings": case["settings"], "backend": "pyscf_cpu", "execution_device": "cpu",
        "scf_initial_guess": case["settings"]["scf_initial_guess"],
        "scf_converged": True, "gradient_completed": True, "energy_hartree": -1.0,
        "forces_ev_per_angstrom": [[0.0, 0.0, 0.0] for _ in case["symbols"]], "s2": 0.0,
        "grid_response": True, "auxiliary_basis_response": False, "dispersion_evaluations": 0,
        "versions": {"pyscf": "synthetic"}, "resolved_numerics": metadata}
    gpu = deepcopy(cpu)
    gpu.update(backend="gpu4pyscf", execution_device="nvidia_cuda", gpu_synchronized=True,
               versions={"pyscf": "synthetic", "gpu4pyscf": "synthetic", "cupy": "synthetic"})
    cpu_path, gpu_path = path.parent / "cpu.json", path.parent / "gpu.json"
    save(cpu_path, cpu)
    save(gpu_path, gpu)
    return cpu_path, gpu_path, cpu, gpu


def test_strict_comparison_links_all_three_gates(bundle):
    path, plan = bundle
    cpu_path, gpu_path, _, _ = synthetic_records(path, plan)
    before = hashes(path.parent)
    result = protocol.compare_protocol_files(path, cpu_path, gpu_path)
    assert result["input_files_verified"] is True
    assert result["resolved_numerics_verified"] is True
    assert result["record_validation_passed"] is True
    assert result["numerical_parity_passed"] is True
    assert result["parity_evidence_accepted"] is True
    assert result["scientific_validated"] is False
    assert before == hashes(path.parent)


@pytest.mark.parametrize("change", ["missing_resolved", "missing_precision", "changed_grid", "wrong_plan_hash",
                                   "plan_bytes_changed", "input_changed", "source_changed", "changed_settings",
                                   "cpu_fallback", "legacy_record", "energy_mismatch"])
def test_each_evidence_layer_can_prevent_false_parity(bundle, change):
    path, plan = bundle
    cpu_path, gpu_path, cpu, gpu = synthetic_records(path, plan)
    if change == "missing_resolved":
        del gpu["resolved_numerics"]
    elif change == "missing_precision":
        gpu["resolved_numerics"]["precision"] = None
    elif change == "changed_grid":
        gpu["resolved_numerics"]["quadrature"]["signature_sha256"] = "e" * 64
    elif change == "wrong_plan_hash":
        gpu["plan_sha256"] = "b" * 64
    elif change == "plan_bytes_changed":
        path.write_text(path.read_text() + "\n")
    elif change in ("input_changed", "source_changed"):
        key = "geometry" if change == "input_changed" else "source_record"
        member = path.parent / plan["cases"][0][key]
        member.write_bytes(member.read_bytes() + b"\n")
    elif change == "changed_settings":
        gpu["settings"]["scf_initial_guess"] = "atom"
    elif change == "cpu_fallback":
        gpu["execution_device"] = "cpu"
    elif change == "legacy_record":
        gpu["schema_version"] = 1
    else:
        gpu["energy_hartree"] += 0.01
    save(gpu_path, gpu)
    result = protocol.compare_protocol_files(path, cpu_path, gpu_path)
    assert result["numerical_parity_passed"] is False
    assert result["parity_evidence_accepted"] is False
    assert result["mismatches"]
    assert result["scientific_validated"] is False


def test_default_cli_requires_verified_v2_inputs(bundle):
    import subprocess
    import sys
    from pathlib import Path
    path, plan = bundle
    cpu_path, gpu_path, _, _ = synthetic_records(path, plan)
    command = [sys.executable, str(Path(__file__).with_name("equivalence.py")), "compare",
               "--plan", str(path), "--cpu", str(cpu_path), "--gpu", str(gpu_path)]
    passed = subprocess.run(command, text=True, capture_output=True)
    assert passed.returncode == 0
    assert json.loads(passed.stdout)["input_files_verified"] is True
    (path.parent / plan["cases"][0]["geometry"]).unlink()
    failed = subprocess.run(command, text=True, capture_output=True)
    assert failed.returncode == 1
    assert json.loads(failed.stdout)["numerical_parity_passed"] is False
