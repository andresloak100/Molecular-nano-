"""CLI electronic-state comparison against real snapshot I/O and synthetic MOs.

The fixture producer contains already populated matrices. Its SCF, gradient,
and stability methods deliberately fail; no electronic calculation is run.
"""

from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

from nanodesign import cli
import nanodesign.electronic_state as electronic
from nanodesign.quantum import PySCFCalculator


# Load the E3-owned synthetic producer by its explicit file, avoiding imports
# from an unrelated installed package named "tests" and avoiding a second copy
# of the snapshot schema or a fake comparison implementation.
_fixture_spec = importlib.util.spec_from_file_location(
    "state_compare_cli_e3_fixtures", Path(__file__).with_name("test_electronic_state.py")
)
_fixtures = importlib.util.module_from_spec(_fixture_spec)
_fixture_spec.loader.exec_module(_fixtures)
FakeMeanField = _fixtures.FakeMeanField
settings_for = _fixtures.settings_for


@pytest.fixture(autouse=True)
def no_quantum_jobs(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("state-compare must inspect saved evidence without quantum jobs")

    monkeypatch.setattr(PySCFCalculator, "calculate", forbidden)
    monkeypatch.setattr(cli, "run", forbidden)
    monkeypatch.setattr(cli, "run_characterization", forbidden)


def snapshot(mean_field=None, *, call_id="synthetic-cli"):
    mean_field = FakeMeanField() if mean_field is None else mean_field
    return electronic.capture_snapshot(mean_field, settings=settings_for(mean_field), call_id=call_id)


def write_pair(tmp_path, left=None, right=None):
    left = snapshot() if left is None else left
    right = deepcopy(left) if right is None else right
    paths = tmp_path / "left evidence.json", tmp_path / "right evidence.json"
    for path, value in zip(paths, (left, right)):
        electronic.save_snapshot(path, value)
    return paths


def read_only_compare(capsys, left, right):
    before = {path: path.read_bytes() for path in (Path(left), Path(right)) if path.is_file()}
    code = cli.main(["state-compare", str(left), str(right)])
    captured = capsys.readouterr()
    for path, raw in before.items():
        assert path.read_bytes() == raw
    assert "Traceback" not in captured.err
    return code, captured


def assert_uncertified(report):
    assert report["electronic_state_identity_verified"] is False
    assert report["ground_state_verified"] is False


def assert_raw_sources(report, left, right):
    for label, path in (("left", Path(left)), ("right", Path(right))):
        raw = path.read_bytes()
        assert report["source_files"][label] == {
            "path": str(path.resolve()), "sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw),
        }


@pytest.mark.parametrize("reference", ["RKS", "UKS"])
@pytest.mark.parametrize("representation", ["sign", "occupied_rotation"])
def test_saved_orbital_signs_and_occupied_rotations_do_not_change_subspace(
        tmp_path, capsys, reference, representation):
    left_mf, right_mf = FakeMeanField(reference), FakeMeanField(reference)
    if representation == "sign":
        transform = np.diag([-1., 1.])
    else:
        angle = .41
        transform = np.array([[math.cos(angle), -math.sin(angle)],
                              [math.sin(angle), math.cos(angle)]])
    if reference == "RKS":
        right_mf.mo_coeff[:, :2] = right_mf.mo_coeff[:, :2] @ transform
    else:
        right_mf.mo_coeff[0, :, :2] = right_mf.mo_coeff[0, :, :2] @ transform
        right_mf.mo_coeff[1, :, 0] *= -1
    left_record, right_record = snapshot(left_mf, call_id="left-call"), snapshot(right_mf, call_id="right-call")
    left, right = write_pair(tmp_path, left_record, right_record)
    names_before = sorted(path.name for path in tmp_path.iterdir())

    code, captured = read_only_compare(capsys, left, right)

    assert code == 0 and captured.err == ""
    report = json.loads(captured.out)
    assert report["comparison_status"] == "compared"
    assert_uncertified(report)
    for channel in ("alpha", "beta"):
        metric = report["metrics"][channel]
        assert metric["singular_values"] == pytest.approx([1.] * metric["n_occupied"], abs=1e-12)
        assert metric["density_distance_squared"] == pytest.approx(0., abs=1e-12)
    assert {key: value for key, value in report.items() if key != "source_files"} == electronic.compare_snapshots(left_record, right_record)
    assert_raw_sources(report, left, right)
    assert sorted(path.name for path in tmp_path.iterdir()) == names_before


def test_distinct_occupied_subspace_is_a_successful_comparison_not_state_identity(tmp_path, capsys):
    right_mf = FakeMeanField()
    angle = math.pi / 3
    rotation = np.eye(4)
    rotation[np.ix_([0, 2], [0, 2])] = [[math.cos(angle), -math.sin(angle)],
                                       [math.sin(angle), math.cos(angle)]]
    right_mf.mo_coeff[1] = right_mf.mo_coeff[1] @ rotation
    left, right = write_pair(tmp_path, snapshot(), snapshot(right_mf))

    code, captured = read_only_compare(capsys, left, right)

    assert code == 0 and captured.err == ""
    report = json.loads(captured.out)
    assert report["comparison_status"] == "compared"
    assert report["metrics"]["alpha"]["density_distance_squared"] == pytest.approx(0., abs=1e-12)
    assert report["metrics"]["beta"]["singular_values"] == pytest.approx([.5], abs=1e-12)
    assert report["metrics"]["beta"]["principal_angles_radians"] == pytest.approx([angle], abs=1e-12)
    assert report["metrics"]["beta"]["density_distance_squared"] == pytest.approx(1.5, abs=1e-12)
    assert_uncertified(report)
    assert_raw_sources(report, left, right)


@pytest.mark.parametrize("context", ["geometry", "frame", "basis", "method", "grid", "producer"])
def test_incompatible_context_exits_two_with_no_similarity_metrics(tmp_path, capsys, context):
    left_record = snapshot()
    right_record = deepcopy(left_record)
    if context == "geometry":
        right_record["geometry"]["positions_bohr"][1][2] += .1
    elif context == "frame":
        right_record["geometry"]["coordinate_frame_id"] = "another-recorded-frame"
    elif context == "basis":
        right_record["basis"]["expanded_basis"]["He"][0][1][0] += .1
    elif context == "method":
        right_record["settings"]["xc"] = "b3lyp"
    elif context == "grid":
        right_record["settings"]["grid_level"] += 1
    else:
        right_record["producer"]["version"] = "different-recorded-version"
    left, right = write_pair(tmp_path, left_record, right_record)

    code, captured = read_only_compare(capsys, left, right)

    assert code == 2 and captured.err == ""
    report = json.loads(captured.out)
    assert report["comparison_status"] == "incompatible_context"
    assert report["metrics"] is None
    assert report["reasons"]
    assert_uncertified(report)
    assert_raw_sources(report, left, right)


def test_equivalent_json_formatting_preserves_exact_raw_hashes_and_absolute_paths(tmp_path, capsys, monkeypatch):
    record = snapshot()
    left, right = write_pair(tmp_path, record, record)
    right.write_bytes((json.dumps(record, indent=4, sort_keys=False) + "\n\n").encode())
    monkeypatch.chdir(tmp_path)

    code, captured = read_only_compare(capsys, Path(left.name), Path(right.name))

    assert code == 0 and captured.err == ""
    report = json.loads(captured.out)
    assert_raw_sources(report, left, right)
    assert report["source_files"]["left"]["sha256"] != report["source_files"]["right"]["sha256"]
    assert report["input_digests"]["left"] == report["input_digests"]["right"]
    assert report["source_files"]["left"]["path"] == str(left)
    assert report["source_files"]["right"]["path"] == str(right)


@pytest.mark.parametrize("raw", [b"", b"{", b"\xff\xfe", b"[]", b'{"schema_version":1}', b'{"energy":NaN}'])
@pytest.mark.parametrize("side", ["left", "right"])
def test_malformed_or_incomplete_input_is_rejected_without_rewriting_either_file(tmp_path, capsys, raw, side):
    left, right = write_pair(tmp_path)
    damaged = left if side == "left" else right
    damaged.write_bytes(raw)

    code, captured = read_only_compare(capsys, left, right)

    assert code == 2
    assert captured.err or captured.out
    if captured.out:
        report = json.loads(captured.out)
        assert report["comparison_status"] != "compared"
        assert report.get("metrics") is None
    assert sorted(path.name for path in tmp_path.iterdir()) == sorted([left.name, right.name])


def test_duplicate_json_key_is_not_silently_accepted(tmp_path, capsys):
    left, right = write_pair(tmp_path)
    raw = right.read_bytes()
    right.write_bytes(raw[:-1] + b',"schema_version":1}')

    code, captured = read_only_compare(capsys, left, right)

    assert code == 2
    assert "duplicate" in (captured.err + captured.out).lower()


@pytest.mark.parametrize("side", ["left", "right"])
def test_each_input_has_a_hard_sixteen_mib_bound(tmp_path, capsys, side):
    left, right = write_pair(tmp_path)
    oversized = left if side == "left" else right
    # Trailing JSON whitespace preserves document validity, isolating the bound.
    raw = oversized.read_bytes()
    oversized.write_bytes(raw + b" " * (16 * 1024 * 1024 + 1 - len(raw)))

    code, captured = read_only_compare(capsys, left, right)

    assert code == 2
    assert captured.err or captured.out
    assert "Traceback" not in captured.err


def test_exact_bound_is_per_file_and_accepts_two_valid_sixteen_mib_inputs(tmp_path, capsys):
    left, right = write_pair(tmp_path)
    for path in (left, right):
        raw = path.read_bytes()
        path.write_bytes(raw + b" " * (16 * 1024 * 1024 - len(raw)))

    code, captured = read_only_compare(capsys, left, right)

    assert code == 0 and captured.err == ""
    report = json.loads(captured.out)
    assert report["comparison_status"] == "compared"
    assert_raw_sources(report, left, right)
    assert report["source_files"]["left"]["size_bytes"] == 16 * 1024 * 1024
    assert report["source_files"]["right"]["size_bytes"] == 16 * 1024 * 1024


@pytest.mark.parametrize("side", ["left", "right"])
def test_missing_source_is_a_clean_read_only_error(tmp_path, capsys, side):
    left, right = write_pair(tmp_path)
    (left if side == "left" else right).unlink()

    code, captured = read_only_compare(capsys, left, right)

    assert code == 2 and captured.err
    assert "Traceback" not in captured.err
    assert len(list(tmp_path.iterdir())) == 1


def test_output_option_is_not_available_for_this_read_only_command(tmp_path, capsys):
    left, right = write_pair(tmp_path)
    forbidden_output = tmp_path / "unexpected output"

    with pytest.raises(SystemExit) as error:
        cli.main(["state-compare", str(left), str(right), "--out", str(forbidden_output)])

    assert error.value.code == 2
    assert "--out" in capsys.readouterr().err
    assert not forbidden_output.exists()
