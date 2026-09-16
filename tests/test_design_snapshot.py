"""Input provenance must describe the bytes actually used, despite later edits."""
from dataclasses import asdict
import hashlib
import json

from ase import Atoms
from ase.io import write
import pytest

from nanodesign import design
from nanodesign.quantum import QuantumSettings


def make_inputs(directory, suffix=".xyz"):
    directory.mkdir()
    initial = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    final = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.76]])
    for state, atoms in (("initial", initial), ("final", final)):
        write(directory / (state + suffix), atoms, format="extxyz" if not suffix else None)
    document = {
        "schema_version": 1, "length_unit": "angstrom",
        "initial": "initial" + suffix, "final": "final" + suffix,
        "fixed_indices": [0], "quantum": asdict(QuantumSettings(spin=0)),
    }
    path = directory / "design.json"
    path.write_text(json.dumps(document))
    return path


def input_hashes(path):
    document = json.loads(path.read_text())
    return {key + "_sha256": design.sha256(source) for key, source in (
        ("design", path), ("initial", path.parent / document["initial"]),
        ("final", path.parent / document["final"]))}


@pytest.mark.parametrize("changed_source", ["design", "initial", "final"])
def test_source_edit_after_first_parse_cannot_change_returned_provenance(tmp_path, monkeypatch, changed_source):
    path = make_inputs(tmp_path / "input")
    expected_hashes = input_hashes(path)
    actual_read = design.read_coordinate_snapshot
    changed = False

    def read_then_edit(*args, **kwargs):
        nonlocal changed
        atoms = actual_read(*args, **kwargs)
        if not changed:
            changed = True
            if changed_source == "design":
                document = json.loads(path.read_text())
                document["quantum"]["basis"] = "sto-3g"
                path.write_text(json.dumps(document))
            else:
                write(path.parent / f"{changed_source}.xyz",
                      Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.9]]))
        return atoms

    monkeypatch.setattr(design, "read_coordinate_snapshot", read_then_edit)
    data, initial, final, settings, hashes = design.load_design(path)
    assert changed
    assert hashes == expected_hashes
    assert input_hashes(path)[changed_source + "_sha256"] != hashes[changed_source + "_sha256"]
    assert settings.basis == data["quantum"]["basis"] == "def2-svp"
    assert initial.get_distance(0, 1) == pytest.approx(0.74)
    assert final.get_distance(0, 1) == pytest.approx(0.76)


def test_captured_inputs_survive_source_removal_during_parse(tmp_path, monkeypatch):
    path = make_inputs(tmp_path / "input")
    expected_hashes = input_hashes(path)
    actual_read = design.read_coordinate_snapshot
    removed = False

    def read_then_remove(*args, **kwargs):
        nonlocal removed
        atoms = actual_read(*args, **kwargs)
        if not removed:
            removed = True
            for source in (path, path.parent / "initial.xyz", path.parent / "final.xyz"):
                source.unlink()
        return atoms

    monkeypatch.setattr(design, "read_coordinate_snapshot", read_then_remove)
    _, initial, final, _, hashes = design.load_design(path)
    assert removed
    assert hashes == expected_hashes
    assert initial.get_distance(0, 1) == pytest.approx(0.74)
    assert final.get_distance(0, 1) == pytest.approx(0.76)


@pytest.mark.parametrize("suffix", [".xyz", ".extxyz", ".xyz.gz", ".xyz.bz2", ".xyz.xz", ".traj", ".traj.gz", ".json", ".db", ".pdb", ""])
def test_snapshots_preserve_text_binary_compressed_and_filename_reader_formats(tmp_path, suffix):
    path = make_inputs(tmp_path / "input", suffix)
    expected_hashes = input_hashes(path)
    _, initial, final, _, hashes = design.load_design(path)
    assert hashes == expected_hashes
    assert initial.get_chemical_symbols() == final.get_chemical_symbols() == ["H", "H"]
    assert initial.get_distance(0, 1) == pytest.approx(0.74)
    assert final.get_distance(0, 1) == pytest.approx(0.76)
    assert initial.constraints[0].get_indices().tolist() == [0]
    assert final.constraints[0].get_indices().tolist() == [0]


def test_relative_symlink_uses_resolved_format_and_original_compressed_bytes(tmp_path):
    path = make_inputs(tmp_path / "input", ".xyz.gz")
    target = path.parent / "initial.xyz.gz"
    alias = path.parent / "initial.not-the-format"
    alias.symlink_to(target.name)
    data = json.loads(path.read_text())
    data["initial"] = alias.name
    path.write_text(json.dumps(data))
    expected = hashlib.sha256(target.read_bytes()).hexdigest()
    _, initial, _, _, hashes = design.load_design(path)
    assert initial.get_distance(0, 1) == pytest.approx(0.74)
    assert hashes["initial_sha256"] == expected


@pytest.mark.parametrize("suffix", [".extxyz", ".extxyz.gz", ".traj", ".db"])
def test_exact_frame_selection_rejects_out_of_range_indices(tmp_path, suffix):
    source = tmp_path / ("frames" + suffix)
    frames = [Atoms("H2", positions=[[0, 0, 0], [0, 0, distance]]) for distance in (0.74, 1.2)]
    write(source, frames)
    raw = source.read_bytes()
    for index, distance in ((0, 0.74), (1, 1.2), (-1, 1.2), (-2, 0.74)):
        atoms = design.read_coordinate_snapshot(source, raw, index=index)
        assert atoms.get_distance(0, 1) == pytest.approx(distance)
    for index in (2, 99, -3, -99):
        with pytest.raises(ValueError, match="frame .*unavailable"):
            design.read_coordinate_snapshot(source, raw, index=index)
