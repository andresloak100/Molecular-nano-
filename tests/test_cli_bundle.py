"""User-facing bundle commands preserve evidence and report transport integrity.

Only synthetic temporary files are used. No quantum calculation is requested.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from nanodesign.cli import main


LIMIT_FLAGS = (
    "--max-files", "--max-file-bytes", "--max-total-bytes", "--max-entries", "--max-depth",
)


def make_tree(directory, contents):
    directory.mkdir()
    for relative, raw in contents.items():
        path = directory / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    return directory


def tree_bytes(directory):
    return {path.relative_to(directory).as_posix(): path.read_bytes()
            for path in directory.rglob("*") if path.is_file()}


def cli_json(capsys, *arguments):
    code = main([str(argument) for argument in arguments])
    captured = capsys.readouterr()
    # Every operational bundle result is machine-readable, including failures.
    report = json.loads(captured.out)
    assert "Traceback" not in captured.err
    return code, report


def assert_summary(report, directory, *, verified):
    assert {
        "bundle", "integrity_verified", "issues", "file_count", "total_bytes",
        "manifest_sha256", "scientific_validation", "reference_closure_verified",
    } <= report.keys()
    assert report["bundle"] == str(directory.resolve())
    assert report["integrity_verified"] is verified
    assert report["scientific_validation"] == "not_assessed"
    assert report["reference_closure_verified"] is False
    assert type(report["file_count"]) is int and report["file_count"] >= 0
    assert type(report["total_bytes"]) is int and report["total_bytes"] >= 0
    if verified:
        assert report["issues"] == []
        assert report["manifest_sha256"] == hashlib.sha256((directory / "manifest.json").read_bytes()).hexdigest()
    else:
        assert report["issues"]
    assert "design_validated" not in report
    assert "method_validated" not in report


def test_create_and_verify_keep_failed_and_malformed_source_records_unchanged(tmp_path, capsys):
    originals = {
        "attempts/failed/result.json": b'{ "status": "failed", "error": "SCF did not converge" }\r\n',
        "attempts/interrupted/result.json": b'{"status":"running","unfinished":',
        "attempts/interrupted/electronic.jsonl": b'{"event":"scf_cycle"}\n{',
        "source-claim.json": b'{"design_validated":true,"reference":"../missing.xyz"}\n',
        ".opaque.dat": b"\x00\xff\xfe\n",
    }
    source = make_tree(tmp_path / "input evidence", originals)
    destination = tmp_path / "portable evidence"
    code, created = cli_json(capsys, "bundle-create", source, "--out", destination)
    assert code == 0
    assert_summary(created, destination, verified=True)
    assert created["file_count"] == len(originals)
    assert created["total_bytes"] == sum(map(len, originals.values()))
    assert tree_bytes(source) == originals
    assert tree_bytes(destination / "payload") == originals

    before_verification = tree_bytes(destination)
    code, verified = cli_json(capsys, "bundle-verify", destination)
    assert code == 0
    assert_summary(verified, destination, verified=True)
    assert verified["file_count"] == created["file_count"]
    assert verified["total_bytes"] == created["total_bytes"]
    assert verified["manifest_sha256"] == created["manifest_sha256"]
    assert tree_bytes(destination) == before_verification
    assert tree_bytes(source) == originals


@pytest.mark.parametrize("corruption", ["same_length_payload", "missing_payload", "malformed_manifest"])
def test_verify_corruption_returns_nonzero_json_and_does_not_repair(tmp_path, capsys, corruption):
    source = make_tree(tmp_path / "source", {"result.json": b"original\n"})
    destination = tmp_path / "bundle"
    assert cli_json(capsys, "bundle-create", source, "--out", destination)[0] == 0
    if corruption == "same_length_payload":
        (destination / "payload" / "result.json").write_bytes(b"tampered\n")
    elif corruption == "missing_payload":
        (destination / "payload" / "result.json").unlink()
    else:
        (destination / "manifest.json").write_bytes(b'{"files": [')
    before = tree_bytes(destination)

    code, report = cli_json(capsys, "bundle-verify", destination)
    assert code == 2
    assert_summary(report, destination, verified=False)
    assert tree_bytes(destination) == before
    assert (source / "result.json").read_bytes() == b"original\n"
    if corruption == "same_length_payload":
        assert any(issue["code"] == "hash_mismatch" for issue in report["issues"])


@pytest.mark.parametrize("limit_flag", LIMIT_FLAGS)
def test_each_limit_reaches_create_and_verify(tmp_path, capsys, limit_flag):
    # Every limit is exceeded independently by the same small source tree:
    # two four-byte files, three entries, and a nested path of depth two.
    originals = {"root.dat": b"1234", "nested/leaf.dat": b"5678"}
    source = make_tree(tmp_path / "source", originals)
    rejected = tmp_path / "rejected"
    code, report = cli_json(capsys, "bundle-create", source, "--out", rejected, limit_flag, "1")
    assert code == 2
    assert_summary(report, rejected, verified=False)
    assert not (rejected / "manifest.json").exists()
    assert tree_bytes(source) == originals

    destination = tmp_path / "valid"
    assert cli_json(capsys, "bundle-create", source, "--out", destination)[0] == 0
    before = tree_bytes(destination)
    code, report = cli_json(capsys, "bundle-verify", destination, limit_flag, "1")
    assert code == 2
    assert_summary(report, destination, verified=False)
    assert tree_bytes(destination) == before


@pytest.mark.parametrize("limit_flag", LIMIT_FLAGS)
def test_nonpositive_limit_reports_json_without_creating_output(tmp_path, capsys, limit_flag):
    source = make_tree(tmp_path / "source", {"result.json": b"evidence\n"})
    destination = tmp_path / "invalid-limit"
    code, report = cli_json(capsys, "bundle-create", source, "--out", destination, limit_flag, "0")
    assert code == 2
    assert_summary(report, destination, verified=False)
    assert "positive integer" in json.dumps(report["issues"])
    assert not destination.exists()


@pytest.mark.parametrize("existing", ["empty_directory", "completed_bundle"])
def test_create_refuses_existing_output_without_changing_any_bytes(tmp_path, capsys, existing):
    source = make_tree(tmp_path / "source", {"result.json": b'{"status":"failed"}\n'})
    destination = tmp_path / "existing"
    if existing == "empty_directory":
        destination.mkdir()
    else:
        assert cli_json(capsys, "bundle-create", source, "--out", destination)[0] == 0
    before = tree_bytes(destination)
    source_before = tree_bytes(source)
    code, report = cli_json(capsys, "bundle-create", source, "--out", destination)
    assert code == 2
    assert_summary(report, destination, verified=False)
    assert destination.is_dir() and tree_bytes(destination) == before
    assert tree_bytes(source) == source_before


def test_create_rejects_destination_inside_source_before_publication(tmp_path, capsys):
    source = make_tree(tmp_path / "source", {"result.json": b"evidence\n"})
    destination = source / "nested bundle"
    before = tree_bytes(source)
    code, report = cli_json(capsys, "bundle-create", source, "--out", destination)
    assert code == 2
    assert_summary(report, destination, verified=False)
    assert tree_bytes(source) == before
    assert not (destination / "manifest.json").exists()


def module_environment():
    # Exercise the command from unrelated working directories in both editable
    # installations and source-tree test environments, without relying on cwd.
    environment = os.environ.copy()
    project = str(Path(__file__).resolve().parents[1])
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, (project, environment.get("PYTHONPATH"))))
    return environment


def test_module_entry_point_verifies_relocated_bundle_after_source_removal(tmp_path):
    originals = {"failed.json": b'{"status":"failed"}\r\n', "binary.dat": bytes(range(256))}
    source = make_tree(tmp_path / "source α with spaces", originals)
    destination = tmp_path / "export β with spaces"
    environment = module_environment()
    created = subprocess.run(
        [sys.executable, "-m", "nanodesign", "bundle-create", source.name, "--out", destination.name],
        cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=30,
    )
    assert created.returncode == 0, created.stderr
    assert_summary(json.loads(created.stdout), destination, verified=True)
    moved = tmp_path / "relocated γ"
    destination.rename(moved)
    shutil.rmtree(source)
    other_cwd = tmp_path / "unrelated working directory"
    other_cwd.mkdir()
    before = tree_bytes(moved)
    verified = subprocess.run(
        [sys.executable, "-m", "nanodesign", "bundle-verify", str(moved)],
        cwd=other_cwd, env=environment, capture_output=True, text=True, timeout=30,
    )
    assert verified.returncode == 0, verified.stderr
    assert_summary(json.loads(verified.stdout), moved, verified=True)
    assert tree_bytes(moved / "payload") == originals
    assert tree_bytes(moved) == before


def test_installed_console_entry_point_reports_corruption_nonzero(tmp_path, capsys):
    console = Path(sys.executable).parent / "nanodesign"
    if not console.is_file() or not os.access(console, os.X_OK):
        pytest.skip("Console entry point is absent; the module entry point is tested separately")
    source = make_tree(tmp_path / "source", {"record.dat": b"original"})
    destination = tmp_path / "bundle"
    assert cli_json(capsys, "bundle-create", source, "--out", destination)[0] == 0
    (destination / "payload" / "record.dat").write_bytes(b"modified")
    result = subprocess.run(
        [str(console), "bundle-verify", str(destination)], cwd=tmp_path,
        env=module_environment(), capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 2, result.stderr
    assert_summary(json.loads(result.stdout), destination, verified=False)
