"""Adversarial provenance/reconstruction checks on tiny E3 captures."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from cross_overlap import BridgeError, build_cross_overlap
from test_integrals import captured_file, molecule


@pytest.fixture
def pair(tmp_path):
    basis = {"H": [[0, [1., 1.]]]}
    left = captured_file(tmp_path, "left", molecule([[0.,0.,0.]], basis))
    right = captured_file(tmp_path, "right", molecule([[0.,0.,.2]], basis))
    return left, right


def compare(pair, **kwargs):
    options = {"common_frame_id":"fixture-frame","atom_mapping":[0]}
    options.update(kwargs)
    return build_cross_overlap(*pair, **options)


def change(path, mutate):
    value = json.loads(path.read_bytes())
    mutate(value)
    path.write_text(json.dumps(value))


def test_exact_file_hashes_and_direction_not_canonical_hash_only(pair):
    right = pair[1]
    right.write_text(json.dumps(json.loads(right.read_bytes()), indent=3)+"\n")
    result = compare(pair)
    for side, path in zip(("left","right"), pair):
        assert result["sources"][side]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        assert result["sources"][side]["size_bytes"] == len(path.read_bytes())
    assert result["sources"]["right"]["sha256"] != result["sources"]["right"]["canonical_content_sha256"]
    assert result["pair_binding"]["direction"].startswith("rows=left/bra")
    assert not result["branch_continuity_verified"]
    assert not result["electronic_state_identity_verified"]
    assert not result["ground_state_verified"]


@pytest.mark.parametrize("frame", [None,"", "  ", 12])
def test_explicit_common_frame_required(pair, frame):
    with pytest.raises(BridgeError, match="common_frame_id"):
        compare(pair, common_frame_id=frame)


@pytest.mark.parametrize("mapping", [None,[],[True],[1],[0,1]])
def test_atom_mapping_must_be_explicit_identity(pair, mapping):
    with pytest.raises(BridgeError, match="atom_mapping"):
        compare(pair, atom_mapping=mapping)


def test_captured_frame_disagreement_rejected(pair):
    change(pair[1], lambda v:v["geometry"].update(coordinate_frame_id="different-frame"))
    with pytest.raises(BridgeError, match="frame conflicts"):
        compare(pair)


def test_explicit_captured_frame_retained(pair):
    for path in pair:
        change(path, lambda v:v["geometry"].update(coordinate_frame_id="fixture-frame"))
    assert compare(pair)["pair_binding"]["recorded_frame_ids"] == ["fixture-frame"]*2


def test_guess_change_is_reported_not_rejected(pair):
    change(pair[1], lambda v:v["settings"].update(scf_initial_guess="atom"))
    assert compare(pair)["changed_comparison_axes"]["scf_initial_guess"] == {"left":"minao","right":"atom"}


def test_physical_method_change_rejected(pair):
    change(pair[1], lambda v:v["settings"].update(xc="b3lyp"))
    with pytest.raises(BridgeError, match="method/settings differ"):
        compare(pair)


def test_different_bases_with_identical_self_overlap_rejected(pair):
    path = pair[1]
    # Both normalized one-s-function captures have S=1. The complete captured
    # basis must still differ and the pair must be rejected.
    other = captured_file(path.parent, "other-exponent", molecule([[0.,0.,.2]], {"H":[[0,[2.,1.]]]}))
    with pytest.raises(BridgeError, match="Pair actual basis"):
        compare((pair[0],other))


def test_stale_expanded_constructor_cannot_pass_same_self_overlap(pair):
    # Change only constructor exponent; actual captured shell fingerprint stays.
    for path in pair:
        change(path, lambda v:v["basis"]["expanded_basis"]["H"][0][1].__setitem__(0,2.))
    with pytest.raises(BridgeError, match="actual captured functions"):
        compare(pair)


def test_actual_libcint_coefficients_are_checked(pair):
    for path in pair:
        change(path, lambda v:v["basis"]["shells"][0]["libcint_coefficients"][0].__setitem__(0,2.))
    with pytest.raises(BridgeError, match="libcint_coefficients"):
        compare(pair)


@pytest.mark.parametrize("basis_input", ["sto-3g", [[0,"1.0 1.0"]], [[0,[1.,0.]]], [[0,1,[1.,1.]]]])
def test_no_basis_names_expressions_or_zero_contractions(pair,basis_input):
    for path in pair:
        change(path, lambda v:v["basis"]["expanded_basis"].update(H=basis_input))
    with pytest.raises(BridgeError):
        compare(pair)


def test_self_overlap_must_reproduce_capture(pair):
    # Keep occupied-metric normalization internally consistent so E3 validation
    # passes; the real one-Gaussian self-overlap remains 1, exposing bad evidence.
    def mutate(v):
        v["ao_overlap"]=[[2.]]
        v["channels"]["alpha"]["occupied_coefficients"]=[[2**-.5]]
    for path in pair:
        change(path,mutate)
    with pytest.raises(BridgeError, match="self-overlap differs"):
        compare(pair)


def test_runtime_normalization_mismatch_refused(pair):
    for path in pair:
        change(path,lambda v:v["basis"].update(normalize_gto=False))
    with pytest.raises(BridgeError,match="normalization convention"):
        compare(pair)


def test_runtime_conversion_mismatch_refused(pair):
    for path in pair:
        change(path,lambda v:v["geometry"].update(bohr_angstrom=1.))
    with pytest.raises(BridgeError,match="Bohr conversion"):
        compare(pair)


def test_conflicting_call_identifier_refused(pair):
    identifier=json.loads(pair[0].read_bytes())["call_id"]
    change(pair[1],lambda v:v.update(call_id=identifier))
    with pytest.raises(BridgeError,match="Same call_id"):
        compare(pair)


def test_symlink_snapshot_refused(pair,tmp_path):
    alias=tmp_path/"alias.json"
    alias.symlink_to(pair[0])
    with pytest.raises(BridgeError):
        compare((alias,pair[1]))


def test_duplicate_json_key_refused(pair):
    pair[0].write_text('{"schema_version":1,"schema_version":1}')
    with pytest.raises(BridgeError,match="Duplicate JSON"):
        compare(pair)


def test_cli_preserves_sources_and_refuses_output_overwrite(pair,tmp_path):
    output=tmp_path/"bridge.json"
    args=[sys.executable,str(HERE/"cross_overlap.py"),*[str(p) for p in pair],
          "--common-frame-id","fixture-frame","--atom-mapping","0","--output",str(output)]
    first=subprocess.run(args,capture_output=True,text=True)
    assert first.returncode == 0,first.stderr
    raw=output.read_bytes()
    assert json.loads(raw)["status"]=="computed"
    second=subprocess.run(args,capture_output=True,text=True)
    assert second.returncode != 0
    assert output.read_bytes()==raw


def test_small_absolute_metric_error_cannot_hide_large_relative_error():
    from cross_overlap import _check_reconstructed_overlap, OVERLAP_ATOL
    captured=np.diag([1.,1e-11,1.])
    actual=np.diag([1.,2e-11,1.])
    assert np.max(np.abs(actual-captured)) < OVERLAP_ATOL
    with pytest.raises(BridgeError,match="captured AO metric"):
        _check_reconstructed_overlap(actual,captured)


def test_metric_relative_gate_accepts_small_relative_error():
    from cross_overlap import _check_reconstructed_overlap, OVERLAP_METRIC_TOL
    captured=np.diag([1.,1e-11,1.])
    actual=np.diag([1.,1e-11*(1+1e-10),1.])
    absolute,relative=_check_reconstructed_overlap(actual,captured)
    assert absolute < 1e-20
    assert relative < OVERLAP_METRIC_TOL
