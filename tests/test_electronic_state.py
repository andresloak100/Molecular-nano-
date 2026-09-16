"""Synthetic determinant evidence checks; no SCF, gradient or stability jobs."""
from copy import deepcopy
import hashlib
import json
import math
import os

import numpy as np
import pytest

import nanodesign.electronic_state as electronic


class FakeMol:
    """Four synthetic s AOs with the accessors used from an existing Mole."""

    cart = False
    ecp = {}
    pseudo = {}
    natm = 2
    nbas = 4

    def __init__(self, symbols, spin, charge=0):
        self.symbols = list(symbols)
        self.spin, self.charge = spin, charge
        self.positions = np.array([[0., 0., 0.], [0., 0., 2.]])
        self._basis = {symbol: [[0, [1.5, 1.0]], [0, [0.4, 1.0]]]
                       for symbol in symbols}
        electrons = sum(self.atom_charges()) - charge
        self.nelec = ((electrons + spin) // 2, (electrons - spin) // 2)

    def nao_nr(self):
        return 4

    def atom_symbol(self, index):
        return self.symbols[index]

    def atom_charges(self):
        return np.array([{"H": 1, "He": 2}[symbol] for symbol in self.symbols])

    def atom_coords(self, unit="Bohr"):
        assert unit == "Bohr"
        return self.positions.copy()

    def ao_labels(self, fmt=False):
        assert fmt is False
        return [(i // 2, self.symbols[i // 2], f"{i % 2 + 1}s", "")
                for i in range(4)]

    def bas_atom(self, index):
        return index // 2

    def bas_angular(self, index):
        return 0

    def bas_kappa(self, index):
        return 0

    def bas_exp(self, index):
        return np.array([1.5 if index % 2 == 0 else 0.4])

    def bas_ctr_coeff(self, index):
        return np.ones((1, 1))

    def _libcint_ctr_coeff(self, index):
        exponent = self.bas_exp(index)[0]
        return np.array([[2 * (2 * exponent) ** 0.75 / math.pi ** 0.25]])


class FakeMeanField:
    """An already populated interface; computational methods deliberately fail."""

    converged = True
    e_tot = -3.0

    def __init__(self, reference="UKS", *, symbols=None, spin=None, charge=0):
        self.reference = reference
        symbols = symbols or (("He", "He") if reference == "RKS" else ("He", "H"))
        spin = (0 if reference == "RKS" else 1) if spin is None else spin
        self.mol = FakeMol(symbols, spin, charge)
        self.overlap = np.array([
            [1.2, 0.2, 0.0, 0.0], [0.2, 1.1, 0.15, 0.0],
            [0.0, 0.15, 1.3, 0.1], [0.0, 0.0, 0.1, 0.9],
        ])
        coefficients = np.linalg.inv(np.linalg.cholesky(self.overlap).T)
        self.mo_coeff = coefficients if reference == "RKS" else np.stack([coefficients, coefficients])
        if reference == "RKS":
            self.mo_occ = np.zeros(4)
            self.mo_occ[:self.mol.nelec[0]] = 2
        else:
            self.mo_occ = np.zeros((2, 4))
            for sector, count in enumerate(self.mol.nelec):
                self.mo_occ[sector, :count] = 1
        self.overlap_reads = 0

    def istype(self, name):
        return self.reference == name

    def get_ovlp(self):
        self.overlap_reads += 1
        return self.overlap.copy()

    def kernel(self, *args, **kwargs):
        pytest.fail("Electronic evidence must not start an SCF solve")

    def nuc_grad_method(self, *args, **kwargs):
        pytest.fail("Electronic evidence must not start a gradient")

    def stability(self, *args, **kwargs):
        pytest.fail("Electronic evidence must not start orbital stability analysis")


def settings_for(mean_field):
    return {"charge": mean_field.mol.charge, "spin": mean_field.mol.spin,
            "xc": "pbe0", "basis": "synthetic-s", "dispersion": "d3bj",
            "grid_level": 3, "conv_tol": 1e-9, "max_cycle": 150,
            "threads": 1, "memory_mb": 2000, "density_fit": False,
            "scf_initial_guess": "minao"}


def capture(mean_field=None, **kwargs):
    mean_field = FakeMeanField() if mean_field is None else mean_field
    return electronic.capture_snapshot(mean_field, settings=settings_for(mean_field),
                                       call_id="synthetic-call-1", **kwargs)


def assert_not_certified(record):
    assert record["electronic_state_identity_verified"] is False
    assert record["ground_state_verified"] is False


def assert_no_comparison(report, status):
    assert report["comparison_status"] == status
    assert report["metrics"] is None
    assert report["reasons"]
    assert_not_certified(report)


def test_capture_save_read_compare_roundtrip_is_immutable_and_bounded(tmp_path):
    mean_field = FakeMeanField()
    settings = settings_for(mean_field)
    snapshot = electronic.capture_snapshot(mean_field, settings=settings,
                                           call_id="synthetic-call", coordinate_frame_id="lab")
    captured = deepcopy(snapshot)
    mean_field.mo_coeff[:] = 0
    mean_field.mol.positions[:] = 99
    mean_field.mol._basis["He"][0][1][0] = 88
    settings["xc"] = "changed-after-capture"
    assert snapshot == captured
    assert mean_field.overlap_reads == 1
    assert snapshot["phase"] == "converged_scf_only"
    assert snapshot["stability_assessed_by_capture"] is False
    assert_not_certified(snapshot)

    path = tmp_path / "snapshot.json"
    saved = electronic.save_snapshot(path, snapshot)
    raw = path.read_bytes()
    assert saved["sha256"] == hashlib.sha256(raw).hexdigest()
    assert saved["size_bytes"] == len(raw)
    assert electronic.read_snapshot(path, max_bytes=len(raw)) == captured
    report = electronic.compare_snapshots(snapshot, electronic.read_snapshot(path))
    assert report["comparison_status"] == "compared"
    assert_not_certified(report)
    assert report["input_digests"]["left"] == report["input_digests"]["right"]
    assert snapshot == captured and path.read_bytes() == raw
    for name, count in (("alpha", 2), ("beta", 1)):
        assert report["metrics"][name]["n_occupied"] == count
        assert report["metrics"][name]["singular_values"] == pytest.approx([1.] * count)
        assert report["metrics"][name]["density_distance_squared"] == pytest.approx(0., abs=1e-12)


def test_rks_maps_spatial_occupations_to_identical_spin_channels():
    snapshot = capture(FakeMeanField("RKS"))
    assert snapshot["source_occupations"] == {"spatial": [2., 2., 0., 0.]}
    assert snapshot["channels"]["alpha"] == snapshot["channels"]["beta"]
    diagnostics = electronic.validate_snapshot(snapshot)
    assert diagnostics["alpha"]["density_electron_trace"] == pytest.approx(2.)
    assert diagnostics["beta"]["density_electron_trace"] == pytest.approx(2.)
    assert diagnostics["ao_overlap_condition_number"] > 1


def test_uks_detects_beta_difference_using_nonorthogonal_metric():
    left_mf = FakeMeanField()
    right_mf = FakeMeanField()
    angle = math.pi / 3
    rotation = np.eye(4)
    rotation[np.ix_([0, 2], [0, 2])] = [[math.cos(angle), -math.sin(angle)],
                                        [math.sin(angle), math.cos(angle)]]
    right_mf.mo_coeff[1] = right_mf.mo_coeff[1] @ rotation
    report = electronic.compare_snapshots(capture(left_mf), capture(right_mf))
    assert report["comparison_status"] == "compared"
    assert report["metrics"]["alpha"]["density_distance_squared"] == pytest.approx(0., abs=1e-12)
    assert report["metrics"]["beta"]["singular_values"] == pytest.approx([math.cos(angle)])
    assert report["metrics"]["beta"]["principal_angles_radians"] == pytest.approx([angle])
    assert report["metrics"]["beta"]["density_distance_squared"] == pytest.approx(2 * math.sin(angle)**2)


def test_relative_soft_metric_change_is_declined_in_both_comparison_directions():
    snapshots = []
    for soft_eigenvalue in (1e-9, 2e-9):
        mean_field = FakeMeanField()
        mean_field.overlap = np.diag([1., soft_eigenvalue, 1., 1.])
        coefficients = np.diag([1., 1 / math.sqrt(soft_eigenvalue), 1., 1.])
        mean_field.mo_coeff = np.stack([coefficients, coefficients])
        snapshots.append(capture(mean_field))
    left, right = snapshots
    assert np.allclose(left["ao_overlap"], right["ao_overlap"],
                       rtol=0, atol=electronic.NUMERICAL_TOLERANCE)
    originals = deepcopy(snapshots)
    for first, second in ((left, right), (right, left)):
        report = electronic.compare_snapshots(first, second)
        assert_no_comparison(report, "incompatible_context")
    assert snapshots == originals


def test_small_well_conditioned_metric_perturbation_has_no_similarity_metrics():
    left_mf, right_mf = FakeMeanField(), FakeMeanField()
    right_mf.overlap += np.eye(4) * 1e-10
    coefficients = np.linalg.inv(np.linalg.cholesky(right_mf.overlap).T)
    right_mf.mo_coeff = np.stack([coefficients, coefficients])
    left, right = capture(left_mf), capture(right_mf)
    assert not np.array_equal(left["ao_overlap"], right["ao_overlap"])
    assert np.allclose(left["ao_overlap"], right["ao_overlap"],
                       rtol=0, atol=electronic.NUMERICAL_TOLERANCE)
    assert_no_comparison(electronic.compare_snapshots(left, right), "incompatible_context")


@pytest.mark.parametrize("reference", ["RKS", "UKS"])
@pytest.mark.parametrize("change", ["sign", "permutation", "rotation"])
def test_occupied_subspace_is_invariant_to_orbital_representation(reference, change):
    left_mf, right_mf = FakeMeanField(reference), FakeMeanField(reference)
    if change == "sign":
        transform = np.diag([-1., 1.])
    elif change == "permutation":
        transform = np.array([[0., 1.], [1., 0.]])
    else:
        theta = 0.37
        transform = np.array([[math.cos(theta), -math.sin(theta)],
                              [math.sin(theta), math.cos(theta)]])
    if reference == "RKS":
        right_mf.mo_coeff[:, :2] = right_mf.mo_coeff[:, :2] @ transform
    else:
        right_mf.mo_coeff[0, :, :2] = right_mf.mo_coeff[0, :, :2] @ transform
        right_mf.mo_coeff[1, :, 0] *= -1
    left, right = capture(left_mf), capture(right_mf)
    originals = deepcopy((left, right))
    report = electronic.compare_snapshots(left, right)
    assert report["comparison_status"] == "compared"
    for metric in report["metrics"].values():
        assert metric["singular_values"] == pytest.approx([1.] * metric["n_occupied"])
        assert metric["density_distance_squared"] == pytest.approx(0., abs=1e-12)
    assert (left, right) == originals
    assert_not_certified(report)


def test_empty_beta_sector_remains_explicit():
    snapshot = capture(FakeMeanField(charge=2))
    assert snapshot["channels"]["beta"]["n_electrons"] == 0
    report = electronic.compare_snapshots(snapshot, snapshot)
    assert report["comparison_status"] == "compared"
    assert report["metrics"]["beta"]["singular_values"] == []
    assert report["metrics"]["beta"]["density_distance_squared"] == 0


@pytest.mark.parametrize("reference", ["RKS", "UKS"])
@pytest.mark.parametrize("fault", ["fractional", "complex_coefficients", "complex_occupations", "boolean", "nonfinite"])
def test_capture_rejects_unsupported_occupations_and_numeric_types(reference, fault):
    mean_field = FakeMeanField(reference)
    if fault == "fractional":
        mean_field.mo_occ.flat[0] = 0.5
    elif fault == "complex_coefficients":
        mean_field.mo_coeff = mean_field.mo_coeff.astype(complex) + 0.01j
    elif fault == "complex_occupations":
        mean_field.mo_occ = mean_field.mo_occ.astype(complex)
    elif fault == "boolean":
        mean_field.mo_coeff = mean_field.mo_coeff.astype(bool)
    else:
        mean_field.mo_coeff.flat[0] = np.nan
    with pytest.raises(electronic.ElectronicStateError):
        capture(mean_field)


@pytest.mark.parametrize("fault", ["unconverged", "unsupported_reference", "ecp", "periodic", "missing_coefficients", "inconsistent_settings"])
def test_capture_requires_supported_converged_molecular_evidence(fault):
    mean_field = FakeMeanField()
    settings = settings_for(mean_field)
    if fault == "unconverged":
        mean_field.converged = False
    elif fault == "unsupported_reference":
        mean_field.reference = "GHF"
    elif fault == "ecp":
        mean_field.mol.ecp = {"He": "unsupported"}
    elif fault == "periodic":
        mean_field.mol.lattice_vectors = lambda: np.eye(3)
    elif fault == "missing_coefficients":
        del mean_field.mo_coeff
    else:
        settings["spin"] = 0
    with pytest.raises(electronic.ElectronicStateError):
        electronic.capture_snapshot(mean_field, settings=settings, call_id="synthetic")


@pytest.mark.parametrize("overlap", [
    [[1., 0.3, 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]],
    np.diag([1., 1., 1., -1.]), np.diag([1., 1., 1., 0.]),
    np.diag([1., 1., 1., 1e-14]),
])
def test_overlap_must_be_symmetric_positive_definite_and_resolved(overlap):
    snapshot = capture()
    snapshot["ao_overlap"] = np.asarray(overlap).tolist()
    with pytest.raises(electronic.ElectronicStateError):
        electronic.validate_snapshot(snapshot)
    assert_no_comparison(electronic.compare_snapshots(snapshot, snapshot), "insufficient_evidence")


@pytest.mark.parametrize("fault", ["nonorthonormal", "electron_count", "occupied_index", "occupation_count", "basis_dimension", "negative_exponent", "element_order", "false_identity_claim"])
def test_saved_record_consistency_is_checked_before_comparison(fault):
    snapshot = capture()
    if fault == "nonorthonormal":
        snapshot["channels"]["alpha"]["occupied_coefficients"][0][0] += 0.2
    elif fault == "electron_count":
        snapshot["channels"]["alpha"]["n_electrons"] = 1
    elif fault == "occupied_index":
        snapshot["channels"]["alpha"]["occupied_indices"] = [0, 2]
    elif fault == "occupation_count":
        snapshot["source_occupations"]["alpha"] = [1., 0., 0., 0.]
    elif fault == "basis_dimension":
        snapshot["basis"]["n_ao"] = 5
    elif fault == "negative_exponent":
        snapshot["basis"]["shells"][0]["exponents"] = [-1.]
    elif fault == "element_order":
        snapshot["geometry"]["symbols"].reverse()
    else:
        snapshot["electronic_state_identity_verified"] = True
    with pytest.raises(electronic.ElectronicStateError):
        electronic.validate_snapshot(snapshot)
    assert_no_comparison(electronic.compare_snapshots(snapshot, snapshot), "insufficient_evidence")


@pytest.mark.parametrize("key", ["geometry", "basis", "settings", "producer", "ao_overlap", "source_occupations", "channels", "call_id", "scf_energy_hartree"])
def test_missing_required_evidence_is_not_synthesized(key):
    snapshot = capture()
    snapshot.pop(key)
    assert_no_comparison(electronic.compare_snapshots(snapshot, snapshot), "insufficient_evidence")


@pytest.mark.parametrize("key", ["grid_level", "max_cycle", "scf_initial_guess"])
def test_incomplete_settings_cannot_acquire_current_defaults(key):
    snapshot = capture()
    snapshot["settings"].pop(key)
    with pytest.raises(electronic.ElectronicStateError):
        electronic.validate_snapshot(snapshot)
    assert_no_comparison(electronic.compare_snapshots(snapshot, snapshot), "insufficient_evidence")


@pytest.mark.parametrize("changed", ["geometry", "frame", "basis", "ao_order", "method", "grid", "density_fit", "version", "class"])
def test_valid_but_incompatible_context_produces_no_similarity(changed):
    left = capture()
    right = deepcopy(left)
    if changed == "geometry":
        right["geometry"]["positions_bohr"][1][2] += 0.1
    elif changed == "frame":
        right["geometry"]["coordinate_frame_id"] = "different-frame"
    elif changed == "basis":
        right["basis"]["expanded_basis"]["He"][0][1][0] += 0.1
    elif changed == "ao_order":
        right["basis"]["ao_labels"][0], right["basis"]["ao_labels"][1] = right["basis"]["ao_labels"][1], right["basis"]["ao_labels"][0]
    elif changed == "method":
        right["settings"]["xc"] = "b3lyp"
    elif changed == "grid":
        right["settings"]["grid_level"] = 4
    elif changed == "density_fit":
        right["settings"]["density_fit"] = True
    else:
        right["producer"][changed] += "-different"
    electronic.validate_snapshot(right)
    assert_no_comparison(electronic.compare_snapshots(left, right), "incompatible_context")


def test_reference_type_and_valid_atom_order_changes_are_incompatible():
    restricted = capture(FakeMeanField("RKS"))
    unrestricted = capture(FakeMeanField("UKS", symbols=("He", "He"), spin=0))
    assert_no_comparison(electronic.compare_snapshots(restricted, unrestricted), "incompatible_context")
    left = capture(FakeMeanField("UKS", symbols=("He", "H")))
    right = capture(FakeMeanField("UKS", symbols=("H", "He")))
    assert_no_comparison(electronic.compare_snapshots(left, right), "incompatible_context")


def test_explicit_starting_guess_and_requested_resources_are_reported_comparison_axes():
    left = capture()
    right = deepcopy(left)
    right["settings"].update(scf_initial_guess="atom", threads=8, memory_mb=4000)
    right["call_id"] = "synthetic-second-call"
    right["scf_energy_hartree"] = -2.9
    report = electronic.compare_snapshots(left, right)
    assert report["comparison_status"] == "compared"
    assert report["comparison_axis"] == {
        "scf_initial_guess": {"left": "minao", "right": "atom"},
        "threads": {"left": 1, "right": 8},
        "memory_mb": {"left": 2000, "right": 4000},
    }
    assert_not_certified(report)


def test_dimension_bound_is_checked_before_large_capture(monkeypatch):
    monkeypatch.setattr(electronic, "MAX_AO", 3)
    mean_field = FakeMeanField()
    with pytest.raises(electronic.ElectronicStateError):
        capture(mean_field)
    assert mean_field.overlap_reads == 0


@pytest.mark.parametrize("max_bytes", [0, -1, True, 1.5, "1000", electronic.DEFAULT_MAX_BYTES + 1])
def test_invalid_io_bounds_are_rejected(tmp_path, max_bytes):
    path = tmp_path / "snapshot.json"
    with pytest.raises(electronic.ElectronicStateError):
        electronic.save_snapshot(path, capture(), max_bytes=max_bytes)
    assert not path.exists()
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot(path, max_bytes=max_bytes)


def test_save_does_not_overwrite_or_create_missing_parents(tmp_path):
    snapshot = capture()
    path = tmp_path / "snapshot.json"
    path.write_bytes(b"existing evidence")
    with pytest.raises(electronic.ElectronicStateError):
        electronic.save_snapshot(path, snapshot)
    assert path.read_bytes() == b"existing evidence"
    with pytest.raises(electronic.ElectronicStateError):
        electronic.save_snapshot(tmp_path / "missing" / "snapshot.json", snapshot)
    assert not (tmp_path / "missing").exists()


def test_oversize_snapshot_is_rejected_without_creating_or_parsing(tmp_path, monkeypatch):
    path = tmp_path / "snapshot.json"
    with pytest.raises(electronic.ElectronicStateError):
        electronic.save_snapshot(path, capture(), max_bytes=1)
    assert not path.exists()
    path.write_bytes(b" " * 129)

    def no_parse(*args, **kwargs):
        pytest.fail("Oversized snapshot must be rejected before JSON parsing")

    monkeypatch.setattr(json, "loads", no_parse)
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot(path, max_bytes=128)
    assert path.read_bytes() == b" " * 129


@pytest.mark.parametrize("raw", [b'{"schema_version":1,"schema_version":1}', b'{"value":NaN}', b"{}", b"[]", b"null", b"{", b"\xff"])
def test_read_rejects_duplicate_json_nonfinite_missing_or_malformed_data(tmp_path, raw):
    path = tmp_path / "snapshot.json"
    path.write_bytes(raw)
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot(path)
    assert path.read_bytes() == raw


@pytest.mark.parametrize("field", ["scf_energy_hartree", "ao_overlap"])
def test_nonfinite_in_otherwise_valid_saved_record_is_rejected(tmp_path, field):
    snapshot = capture()
    if field == "scf_energy_hartree":
        snapshot[field] = float("nan")
    else:
        snapshot[field][0][0] = float("inf")
    path = tmp_path / "snapshot.json"
    path.write_text(json.dumps(snapshot), encoding="utf-8")
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot(path)


def test_exact_bytes_reader_uses_captured_bytes_after_source_changes(tmp_path):
    snapshot = capture()
    snapshot["call_id"] = "synthetic-α"
    path = tmp_path / "snapshot.json"
    raw = (json.dumps(snapshot, ensure_ascii=False, indent=3) + "\r\n").encode("utf-8")
    path.write_bytes(raw)
    captured = path.read_bytes()
    path.write_bytes(b"changed after caller captured and hashed the file")
    assert electronic.read_snapshot_bytes(captured, max_bytes=len(captured)) == snapshot
    assert captured == raw
    assert path.read_bytes() == b"changed after caller captured and hashed the file"


@pytest.mark.parametrize("raw", ["{}", bytearray(b"{}"), memoryview(b"{}"), None, b""])
def test_exact_bytes_reader_requires_nonempty_immutable_bytes(raw):
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot_bytes(raw)


def test_exact_bytes_reader_checks_bound_before_parsing(monkeypatch):
    def no_parse(*args, **kwargs):
        pytest.fail("Byte bound must be checked before parsing")

    monkeypatch.setattr(json, "loads", no_parse)
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot_bytes(b" " * 129, max_bytes=128)


@pytest.mark.parametrize("where", ["file", "parent"])
def test_snapshot_io_rejects_symlinks_without_altering_target(tmp_path, where):
    actual = tmp_path / "actual"
    actual.mkdir()
    target = actual / "snapshot.json"
    electronic.save_snapshot(target, capture())
    before = target.read_bytes()
    if where == "file":
        path = tmp_path / "linked.json"
        path.symlink_to(target)
    else:
        link = tmp_path / "linked-parent"
        link.symlink_to(actual, target_is_directory=True)
        path = link / "snapshot.json"
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot(path)
    with pytest.raises(electronic.ElectronicStateError):
        electronic.save_snapshot(path, capture())
    assert target.read_bytes() == before


def test_missing_or_directory_snapshot_is_rejected(tmp_path):
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot(tmp_path / "missing.json")
    directory = tmp_path / "directory"
    directory.mkdir()
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot(directory)


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="Platform lacks named pipes")
def test_nonregular_snapshot_is_rejected_without_blocking(tmp_path):
    path = tmp_path / "pipe"
    os.mkfifo(path)
    with pytest.raises(electronic.ElectronicStateError):
        electronic.read_snapshot(path)
