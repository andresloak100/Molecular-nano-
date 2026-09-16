"""Result records remain parseable when an update cannot be completed."""

import json
from pathlib import Path

import pytest

from nanodesign import workflow


@pytest.mark.parametrize("bad_data,exception", [
    ({"nonfinite": float("nan")}, ValueError),
    ({"unserializable": object()}, TypeError),
])
def test_serialization_failure_preserves_previous_json_without_temporary_file(tmp_path, bad_data, exception):
    path = tmp_path / "result.json"
    previous = '{"status": "running", "step": 2}\n'
    path.write_text(previous)
    with pytest.raises(exception):
        workflow.json_write(path, bad_data)
    assert path.read_text() == previous
    assert json.loads(path.read_text())["step"] == 2
    assert list(tmp_path.iterdir()) == [path]


def test_replacement_exposes_complete_new_record_after_fsync(tmp_path, monkeypatch):
    path = tmp_path / "result.json"
    path.write_text('{"status": "running"}\n')
    data = {"status": "completed", "values": list(range(100)), "unit": "Å"}
    real_fsync = workflow.os.fsync
    real_replace = workflow.os.replace
    synced = []

    def fsync(fd):
        real_fsync(fd)
        synced.append(fd)

    def replace(source, destination):
        source = Path(source)
        assert synced
        assert source.parent == path.parent
        assert destination == path
        assert json.loads(source.read_text()) == data
        assert json.loads(path.read_text()) == {"status": "running"}
        real_replace(source, destination)

    monkeypatch.setattr(workflow.os, "fsync", fsync)
    monkeypatch.setattr(workflow.os, "replace", replace)
    workflow.json_write(path, data)
    assert json.loads(path.read_text()) == data
    assert path.read_text().endswith("\n")
    assert list(tmp_path.iterdir()) == [path]


def test_replace_failure_preserves_original_and_cleans_temporary_file(tmp_path, monkeypatch):
    path = tmp_path / "result.json"
    previous = '{"status": "running", "step": 2}\n'
    path.write_text(previous)
    temporary_paths = []

    def broken_replace(source, destination):
        source = Path(source)
        temporary_paths.append(source)
        assert source.parent == path.parent
        assert json.loads(source.read_text()) == {"status": "completed"}
        raise OSError("replacement failed")

    monkeypatch.setattr(workflow.os, "replace", broken_replace)
    with pytest.raises(OSError, match="replacement failed"):
        workflow.json_write(path, {"status": "completed"})
    assert path.read_text() == previous
    assert temporary_paths and not temporary_paths[0].exists()
    assert list(tmp_path.iterdir()) == [path]
