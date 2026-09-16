"""Portable evidence is a byte-preserving snapshot, never a science verdict.

All fixtures are synthetic files. These tests do not import an electronic solver
or touch archived project evidence.
"""
from dataclasses import FrozenInstanceError, asdict
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

import nanodesign.bundle as bundle
from nanodesign.bundle import BundleError, BundleLimits, export_bundle, verify_bundle


def make_source(root, files=None, directories=()):
    root.mkdir()
    files = {"result.json": b'{"status":"failed"}\n'} if files is None else files
    for relative in directories:
        (root / relative).mkdir(parents=True, exist_ok=True)
    for relative, raw in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    return root


def file_bytes(root):
    return {path.relative_to(root).as_posix(): path.read_bytes()
            for path in root.rglob("*") if path.is_file()}


def read_manifest(root):
    return json.loads((root / "manifest.json").read_bytes())


def write_manifest(root, manifest):
    (root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def assert_unassessed(record):
    assert record["scientific_validation"] == "not_assessed"
    assert record["reference_closure_verified"] is False
    assert "design_validated" not in record
    assert "scientific_status" not in record


def assert_rejected_report(report):
    assert report["integrity_verified"] is False
    assert report["issues"]
    assert_unassessed(report)
    for issue in report["issues"]:
        assert isinstance(issue["code"], str) and issue["code"]
        assert isinstance(issue["path"], str)
        assert isinstance(issue["message"], str) and issue["message"]


@pytest.fixture
def simple_bundle(tmp_path):
    source = make_source(tmp_path / "source", {"result.json": b"original\n"},
                         directories=("empty",))
    destination = tmp_path / "bundle"
    export_bundle(source, destination)
    return source, destination


def test_roundtrip_preserves_raw_failures_hidden_files_and_all_attempts(tmp_path):
    originals = {
        "plan.json": b'{ "design_validated": true, "note": "source claim" }\r\n',
        "campaign.json": b'{"attempts": ["failed", "completed"]}\n',
        "runs/pose-0001/attempt-0001/result.json": b'{"status":"failed",',
        "runs/pose-0001/attempt-0002/result.json": b'{"status":"completed"}\n',
        "runs/pose-0001/attempt-0001/electronic.jsonl": b'{"event":"start"}\n{',
        ".campaign.lock": b"",
        ".hidden/evidence.bin": bytes(range(256)),
        "unicode α/input.traj.gz": b"\x1f\x8b\x00opaque compressed evidence\xff",
        "manifest.json": b"This is source evidence, not the export manifest.\n",
    }
    source = make_source(tmp_path / "source", originals, directories=("empty/nested",))
    before = file_bytes(source)
    destination = tmp_path / "bundle"
    manifest = export_bundle(source, destination)

    assert file_bytes(source) == before == originals
    assert file_bytes(destination / "payload") == originals
    assert manifest == read_manifest(destination)
    assert manifest["format"] == "nanodesign-evidence-bundle"
    assert manifest["format_version"] == 1
    assert manifest["file_count"] == len(originals)
    assert manifest["total_bytes"] == sum(map(len, originals.values()))
    assert manifest["limits"] == asdict(BundleLimits())
    assert manifest["created_utc"]
    assert manifest["snapshot"]["atomic"] is False
    assert manifest["snapshot"]["source_rechecked"] is True
    assert_unassessed(manifest)
    assert set(manifest["directories"]) == {
        path.relative_to(source).as_posix() for path in source.rglob("*") if path.is_dir()
    }
    entries = {entry["path"]: entry for entry in manifest["files"]}
    assert len(entries) == len(manifest["files"]) == len(originals)
    for relative, raw in originals.items():
        assert entries[relative]["size_bytes"] == len(raw)
        assert entries[relative]["sha256"] == hashlib.sha256(raw).hexdigest()

    captured = file_bytes(destination)
    report = verify_bundle(destination)
    assert report["integrity_verified"] is True
    assert report["issues"] == []
    assert_unassessed(report)
    assert file_bytes(destination) == captured
    assert file_bytes(source) == originals


def test_empty_source_is_a_valid_snapshot_with_empty_directories(tmp_path):
    source = make_source(tmp_path / "source", {}, directories=("empty/deeper",))
    destination = tmp_path / "bundle"
    manifest = export_bundle(source, destination)
    assert manifest["files"] == []
    assert manifest["file_count"] == manifest["total_bytes"] == 0
    assert set(manifest["directories"]) == {"empty", "empty/deeper"}
    assert (destination / "payload" / "empty" / "deeper").is_dir()
    assert verify_bundle(destination)["integrity_verified"] is True


def test_bundle_is_portable_after_rename_and_source_removal(simple_bundle, tmp_path):
    source, destination = simple_bundle
    (source / "result.json").unlink()
    (source / "empty").rmdir()
    source.rmdir()
    moved = tmp_path / "relocated evidence"
    destination.rename(moved)
    assert verify_bundle(moved)["integrity_verified"] is True


def test_bundle_limits_are_frozen_and_have_documented_defaults():
    limits = BundleLimits()
    assert asdict(limits) == {
        "max_files": 10000, "max_file_bytes": 32 * 1024**2,
        "max_total_bytes": 256 * 1024**2, "max_entries": 20000, "max_depth": 32,
    }
    with pytest.raises(FrozenInstanceError):
        limits.max_files = 1


@pytest.mark.parametrize("field", list(asdict(BundleLimits())))
@pytest.mark.parametrize("value", [0, -1, True, False, 1.5, "2", None])
def test_limits_require_positive_nonboolean_integers(field, value):
    with pytest.raises((TypeError, ValueError)):
        BundleLimits(**{field: value})


@pytest.mark.parametrize("kind", ["max_files", "max_file_bytes", "max_total_bytes", "max_entries", "max_depth"])
def test_export_limit_failure_never_publishes_manifest(tmp_path, kind):
    source = make_source(tmp_path / "source", {
        "nested/deeper/result.json": b"1234", "other.txt": b"5678",
    })
    before = file_bytes(source)
    destination = tmp_path / "bundle"
    limits = BundleLimits(**{kind: 1})
    with pytest.raises(BundleError):
        export_bundle(source, destination, limits=limits)
    assert not (destination / "manifest.json").exists()
    assert file_bytes(source) == before


@pytest.mark.parametrize("kind", ["max_files", "max_file_bytes", "max_total_bytes", "max_entries", "max_depth"])
def test_verify_applies_callers_limits(tmp_path, kind):
    source = make_source(tmp_path / "source", {
        "nested/deeper/result.json": b"1234", "other.txt": b"5678",
    })
    destination = tmp_path / "bundle"
    export_bundle(source, destination)
    assert_rejected_report(verify_bundle(destination, limits=BundleLimits(**{kind: 1})))


@pytest.mark.parametrize("existing_kind", ["empty_directory", "nonempty_directory", "file"])
def test_export_refuses_existing_destination_without_altering_it(tmp_path, existing_kind):
    source = make_source(tmp_path / "source")
    destination = tmp_path / "bundle"
    if existing_kind == "file":
        destination.write_bytes(b"existing file")
    else:
        destination.mkdir()
        if existing_kind == "nonempty_directory":
            (destination / "keep.txt").write_bytes(b"existing evidence")
    before = destination.read_bytes() if destination.is_file() else file_bytes(destination)
    with pytest.raises(BundleError):
        export_bundle(source, destination)
    after = destination.read_bytes() if destination.is_file() else file_bytes(destination)
    assert before == after


@pytest.mark.parametrize("location", ["same", "child", "grandchild"])
def test_export_refuses_destination_inside_source(tmp_path, location):
    source = make_source(tmp_path / "source")
    destination = {"same": source, "child": source / "bundle",
                   "grandchild": source / "nested" / "bundle"}[location]
    before = file_bytes(source)
    with pytest.raises(BundleError):
        export_bundle(source, destination)
    assert file_bytes(source) == before
    assert not (destination / "manifest.json").exists()


def test_destination_parent_replacement_cannot_create_output_inside_source(tmp_path, monkeypatch):
    source = make_source(tmp_path / "source")
    before = file_bytes(source)
    parent = tmp_path / "output-parent"
    parent.mkdir()
    moved_parent = tmp_path / "original-output-parent"
    destination = parent / "bundle"
    original_mkdir = os.mkdir
    redirected = False

    def replace_parent_before_mkdir(path, mode=0o777, *, dir_fd=None):
        nonlocal redirected
        if not redirected and Path(path).name == destination.name:
            redirected = True
            parent.rename(moved_parent)
            parent.symlink_to(source, target_is_directory=True)
        return original_mkdir(path, mode, dir_fd=dir_fd)

    monkeypatch.setattr(os, "mkdir", replace_parent_before_mkdir)
    with pytest.raises(BundleError):
        export_bundle(source, destination)

    assert redirected, "The destination-creation race must be exercised"
    assert file_bytes(source) == before
    assert not (source / destination.name).exists()
    assert not (destination / "manifest.json").exists()
    assert not (moved_parent / destination.name / "manifest.json").exists()


def test_modified_written_payload_prevents_manifest_publication(tmp_path, monkeypatch):
    source = make_source(tmp_path / "source", {
        "first.bin": b"AAAA", "second.bin": b"BBBB",
    })
    before = file_bytes(source)
    destination = tmp_path / "bundle"
    original_write = bundle._write
    modified = False

    def write_then_modify_payload(root, relative, raw):
        nonlocal modified
        result = original_write(root, relative, raw)
        if not modified and relative.startswith("payload/"):
            modified = True
            (root.path / relative).write_bytes(b"X" * len(raw))
        return result

    monkeypatch.setattr(bundle, "_write", write_then_modify_payload)
    with pytest.raises(BundleError):
        export_bundle(source, destination)

    assert modified, "The payload mutation must occur after its initial write"
    assert not (destination / "manifest.json").exists()
    assert file_bytes(source) == before


@pytest.mark.parametrize("which", ["source", "destination"])
def test_export_rejects_parent_components_even_when_they_resolve_locally(tmp_path, which):
    source = make_source(tmp_path / "source")
    (tmp_path / "intermediate").mkdir()
    destination = tmp_path / "bundle"
    if which == "source":
        source = tmp_path / "intermediate" / ".." / "source"
    else:
        destination = tmp_path / "intermediate" / ".." / "bundle"
    with pytest.raises(BundleError):
        export_bundle(source, destination)
    assert not (tmp_path / "bundle" / "manifest.json").exists()


@pytest.mark.parametrize("source_kind", ["missing", "file", "symlink", "descendant_symlink"])
def test_export_rejects_invalid_or_linked_sources(tmp_path, source_kind):
    source = tmp_path / "source"
    if source_kind == "file":
        source.write_bytes(b"not a directory")
    elif source_kind == "symlink":
        actual = make_source(tmp_path / "actual")
        source.symlink_to(actual, target_is_directory=True)
    elif source_kind == "descendant_symlink":
        make_source(source)
        (source / "linked.json").symlink_to(source / "result.json")
    destination = tmp_path / "bundle"
    with pytest.raises(BundleError):
        export_bundle(source, destination)
    assert not (destination / "manifest.json").exists()


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="Platform lacks named pipes")
def test_export_rejects_special_file_without_opening_it(tmp_path):
    source = make_source(tmp_path / "source")
    os.mkfifo(source / "pipe")
    destination = tmp_path / "bundle"
    with pytest.raises(BundleError):
        export_bundle(source, destination)
    assert not (destination / "manifest.json").exists()


@pytest.mark.parametrize("mutation", ["changed", "missing", "extra_file", "extra_directory", "missing_directory", "root_file", "root_directory"])
def test_verify_detects_changed_missing_extra_payload_and_root_entries(simple_bundle, mutation):
    _, destination = simple_bundle
    payload = destination / "payload"
    if mutation == "changed":
        (payload / "result.json").write_bytes(b"tampered\n")
    elif mutation == "missing":
        (payload / "result.json").unlink()
    elif mutation == "extra_file":
        (payload / "extra.json").write_bytes(b"{}")
    elif mutation == "extra_directory":
        (payload / "extra").mkdir()
    elif mutation == "missing_directory":
        (payload / "empty").rmdir()
    elif mutation == "root_file":
        (destination / "extra.json").write_bytes(b"{}")
    else:
        (destination / "extra").mkdir()
    before = file_bytes(destination)
    assert_rejected_report(verify_bundle(destination))
    assert file_bytes(destination) == before


@pytest.mark.parametrize("raw", [b"{", b"[]", b"null", b"\xff", b"{}"])
def test_unreadable_or_malformed_manifest_returns_report_not_exception(simple_bundle, raw):
    _, destination = simple_bundle
    (destination / "manifest.json").write_bytes(raw)
    assert_rejected_report(verify_bundle(destination))
    assert (destination / "manifest.json").read_bytes() == raw


def test_missing_manifest_returns_report_and_does_not_create_one(simple_bundle):
    _, destination = simple_bundle
    (destination / "manifest.json").unlink()
    assert_rejected_report(verify_bundle(destination))
    assert not (destination / "manifest.json").exists()


@pytest.mark.parametrize("field,value", [
    ("format", "other-format"), ("format_version", 2),
    ("file_count", 999), ("total_bytes", -1),
    ("files", {}), ("directories", {}),
])
def test_invalid_manifest_shape_or_totals_return_report(simple_bundle, field, value):
    _, destination = simple_bundle
    manifest = read_manifest(destination)
    manifest[field] = value
    write_manifest(destination, manifest)
    assert_rejected_report(verify_bundle(destination))


@pytest.mark.parametrize("field,value", [
    ("size_bytes", -1), ("size_bytes", True), ("size_bytes", 1.5),
    ("size_bytes", "9"), ("sha256", "x" * 64),
    ("sha256", "0" * 63), ("sha256", None),
])
def test_malformed_file_sizes_and_hashes_return_report(simple_bundle, field, value):
    _, destination = simple_bundle
    manifest = read_manifest(destination)
    manifest["files"][0][field] = value
    write_manifest(destination, manifest)
    assert_rejected_report(verify_bundle(destination))


@pytest.mark.parametrize("section", ["files", "directories"])
@pytest.mark.parametrize("unsafe", ["../outside.txt", "/absolute.txt", "nested/../result.json", "", "."])
def test_unsafe_manifest_paths_return_report(simple_bundle, section, unsafe):
    _, destination = simple_bundle
    manifest = read_manifest(destination)
    if section == "files":
        manifest["files"][0]["path"] = unsafe
    else:
        manifest["directories"][0] = unsafe
    write_manifest(destination, manifest)
    assert_rejected_report(verify_bundle(destination))


@pytest.mark.parametrize("section", ["files", "directories"])
def test_duplicate_manifest_paths_return_report(simple_bundle, section):
    _, destination = simple_bundle
    manifest = read_manifest(destination)
    manifest[section].append(manifest[section][0])
    write_manifest(destination, manifest)
    assert_rejected_report(verify_bundle(destination))


@pytest.mark.parametrize("target", ["manifest", "payload_file", "payload_directory", "bundle_root"])
def test_verify_rejects_symlinks_and_preserves_external_target(simple_bundle, tmp_path, target):
    _, destination = simple_bundle
    external = tmp_path / "external"
    external.mkdir()
    (external / "evidence").write_bytes(b"external bytes")
    if target == "manifest":
        original = (destination / "manifest.json").read_bytes()
        (external / "manifest.json").write_bytes(original)
        (destination / "manifest.json").unlink()
        (destination / "manifest.json").symlink_to(external / "manifest.json")
    elif target == "payload_file":
        (destination / "payload" / "result.json").unlink()
        (destination / "payload" / "result.json").symlink_to(external / "evidence")
    elif target == "payload_directory":
        (destination / "payload" / "empty").rmdir()
        (destination / "payload" / "empty").symlink_to(external, target_is_directory=True)
    else:
        alias = tmp_path / "linked-bundle"
        alias.symlink_to(destination, target_is_directory=True)
        destination = alias
    before = file_bytes(external)
    assert_rejected_report(verify_bundle(destination))
    assert file_bytes(external) == before


def test_oversized_manifest_is_rejected_before_json_parsing(simple_bundle, monkeypatch):
    _, destination = simple_bundle
    monkeypatch.setattr(bundle, "MAX_MANIFEST_BYTES", 128)
    (destination / "manifest.json").write_bytes(b" " * 129)

    def no_parse(*args, **kwargs):
        pytest.fail("Manifest bytes exceeding the configured bound must not be parsed")

    monkeypatch.setattr(json, "loads", no_parse)
    assert_rejected_report(verify_bundle(destination))


def test_bundle_module_import_needs_no_quantum_dependencies():
    project = str(Path(__file__).resolve().parents[1])
    program = f"""
import sys
sys.path.insert(0, {project!r})
class RejectScientificImports:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {{'numpy', 'scipy', 'ase', 'pyscf', 'dftd3'}}:
            raise RuntimeError('Unexpected scientific import: ' + fullname)
sys.meta_path.insert(0, RejectScientificImports())
import nanodesign.bundle
"""
    completed = subprocess.run([sys.executable, "-c", program], capture_output=True,
                               text=True, timeout=15, check=False)
    assert completed.returncode == 0, completed.stderr
