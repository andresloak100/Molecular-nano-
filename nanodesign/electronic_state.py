"""Optional determinant evidence; similarity is never physical state verification.

No SCF, gradient, stability analysis, molecule deserialization or state selection
is performed here. Capture accepts an already converged real RKS/UKS object.
"""
from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
from dataclasses import fields
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import stat

import numpy as np


DEFAULT_MAX_BYTES = 16 * 1024 * 1024
MAX_AO = 1024
NUMERICAL_TOLERANCE = 1e-8
MIN_METRIC_RELATIVE_EIGENVALUE = 1e-12


class ElectronicStateError(ValueError):
    """Missing, unsupported, inconsistent or unwritable electronic evidence."""


def _require(condition, message):
    if not condition:
        raise ElectronicStateError(message)


def _integer(value, name, minimum=0):
    _require(type(value) is int and value >= minimum, f"Invalid {name}")
    return value


def _text(value, name):
    _require(isinstance(value, str) and 0 < len(value) <= 4096, f"Invalid {name}")
    return value


def _array(value, name, shape=None):
    try:
        if isinstance(value, (list, tuple)):
            pending = list(value)
            while pending:
                item = pending.pop()
                if isinstance(item, (list, tuple)):
                    pending.extend(item)
                else:
                    _require(type(item) in (int, float), f"{name} must contain only real numeric values")
        raw = np.asarray(value)
        _require(raw.dtype.kind in "fiu", f"{name} must contain real numbers, not strings, booleans or complex values")
        arr = np.asarray(raw, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ElectronicStateError(f"Invalid {name}: {error}") from error
    _require(np.all(np.isfinite(arr)), f"Nonfinite {name}")
    if shape is not None:
        _require(arr.shape == shape, f"Invalid {name} shape: expected {shape}, got {arr.shape}")
    return arr


def _json_bytes(value):
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    except (ValueError, TypeError, OverflowError, RecursionError) as error:
        raise ElectronicStateError(f"Evidence must be finite JSON data: {error}") from error


def _metric(overlap, nao):
    s = _array(overlap, "AO overlap", (nao, nao))
    _require(np.max(np.abs(s - s.T)) <= NUMERICAL_TOLERANCE, "AO overlap is not symmetric")
    try:
        eigenvalues = np.linalg.eigvalsh((s + s.T) / 2)
    except np.linalg.LinAlgError as error:
        raise ElectronicStateError("AO metric eigenvalue evaluation failed") from error
    _require(eigenvalues[-1] > 0 and eigenvalues[0] > MIN_METRIC_RELATIVE_EIGENVALUE * eigenvalues[-1],
             "AO overlap must be positive definite and numerically well conditioned for this implementation")
    return s, eigenvalues


def _validate(snapshot):
    """Validate raw evidence and derive matrices; never trust saved summaries."""
    _require(isinstance(snapshot, dict), "Expected an electronic snapshot object")
    _require(type(snapshot.get("schema_version")) is int and snapshot["schema_version"] == 1,
             "Unsupported electronic snapshot schema")
    _require(snapshot.get("kind") == "mean_field_electronic_snapshot", "Not an electronic snapshot")
    _require(snapshot.get("phase") == "converged_scf_only", "Snapshot phase must be converged_scf_only")
    _require(snapshot.get("reference") in ("RKS", "UKS"), "Only RKS and UKS are supported")
    _require(snapshot.get("scf_converged") is True, "Snapshot must describe converged SCF")
    _require(snapshot.get("electronic_state_identity_verified") is False and
             snapshot.get("ground_state_verified") is False, "Snapshot cannot certify physical state identity")
    _text(snapshot.get("call_id"), "call_id")
    geometry, basis, settings = snapshot.get("geometry"), snapshot.get("basis"), snapshot.get("settings")
    _require(all(isinstance(x, dict) for x in (geometry, basis, settings)), "Missing geometry, basis or settings")
    from .quantum import QuantumSettings
    _require(set(settings) == {field.name for field in fields(QuantumSettings)}, "Complete QuantumSettings fields are required")
    try:
        QuantumSettings(**settings)
    except (TypeError, ValueError) as error:
        raise ElectronicStateError(f"Invalid requested quantum settings: {error}") from error
    symbols = geometry.get("symbols")
    _require(isinstance(symbols, list) and 0 < len(symbols) <= MAX_AO, "Invalid atom order")
    for symbol in symbols:
        _text(symbol, "element symbol")
    from ase.data import atomic_numbers
    try:
        numbers = [atomic_numbers[s] for s in symbols]
    except KeyError as error:
        raise ElectronicStateError("Unknown element symbol") from error
    saved_numbers = geometry.get("atomic_numbers")
    _require(isinstance(saved_numbers, list) and all(type(z) is int for z in saved_numbers) and
             all(z > 0 for z in numbers) and saved_numbers == numbers, "Atomic numbers/symbols mismatch")
    _array(geometry.get("positions_bohr"), "positions_bohr", (len(symbols), 3))
    conversion = geometry.get("bohr_angstrom")
    _require(type(conversion) is float and math.isfinite(conversion) and conversion > 0, "Missing recorded Bohr conversion")
    _require(geometry.get("coordinate_frame_id") is None or isinstance(geometry["coordinate_frame_id"], str), "Invalid coordinate frame ID")
    pbc = geometry.get("pbc")
    _require(geometry.get("frame") == "solver_cartesian" and isinstance(pbc, list) and
             len(pbc) == 3 and all(x is False for x in pbc),
             "Only the original nonperiodic solver frame is supported")
    charge = _integer(geometry.get("charge"), "charge", minimum=-MAX_AO)
    spin = _integer(geometry.get("spin"), "spin")
    nelec = sum(numbers) - charge
    _require(nelec > 0 and spin <= nelec and (nelec - spin) % 2 == 0, "Inconsistent charge/spin/electron count")
    _require(type(settings.get("charge")) is int and settings["charge"] == charge and
             type(settings.get("spin")) is int and settings["spin"] == spin, "Settings and molecular charge/spin disagree")
    for key in ("xc", "basis", "scf_initial_guess"):
        _text(settings.get(key), f"settings.{key}")
    _require(snapshot["reference"] != "RKS" or spin == 0, "RKS requires spin zero")
    producer = snapshot.get("producer")
    _require(isinstance(producer, dict) and producer.get("engine") == "PySCF", "Missing producer identity")
    _text(producer.get("version"), "producer.version")
    _text(producer.get("class"), "producer.class")
    energy = snapshot.get("scf_energy_hartree")
    _require(type(energy) in (int, float) and math.isfinite(energy), "Invalid SCF energy")
    nao = _integer(basis.get("n_ao"), "basis.n_ao", 1)
    _require(nao <= MAX_AO, f"AO dimension exceeds {MAX_AO}")
    _require(basis.get("representation") in ("spherical", "cartesian"), "Unknown AO convention")
    _require(type(basis.get("normalize_gto")) is bool, "Missing Gaussian normalization convention")
    _require(isinstance(basis.get("expanded_basis"), dict) and bool(basis["expanded_basis"]), "Missing expanded basis construction input")
    _require(basis.get("ecp") is False and basis.get("pseudo") is False, "Only all-electron basis evidence is supported")
    labels = basis.get("ao_labels")
    _require(isinstance(labels, list) and len(labels) == nao, "Missing ordered AO labels")
    for label in labels:
        _require(isinstance(label, list) and len(label) == 4, "Invalid AO label")
        atom = _integer(label[0], "AO atom index")
        _require(atom < len(symbols) and label[1] == symbols[atom] and
                 all(isinstance(x, str) for x in label[1:]), "AO labels do not match atom order")
    shells = basis.get("shells")
    _require(isinstance(shells, list) and 0 < len(shells) <= nao, "Missing expanded basis shells")
    expanded_nao = 0
    for shell in shells:
        _require(isinstance(shell, dict), "Invalid basis shell")
        atom = _integer(shell.get("atom_index"), "shell atom index")
        angular = _integer(shell.get("angular_momentum"), "angular momentum")
        _require(atom < len(symbols) and angular <= 12, "Unsupported basis shell")
        exponents = _array(shell.get("exponents"), "basis exponents")
        coeff = _array(shell.get("contraction_coefficients"), "contraction coefficients")
        libcint = _array(shell.get("libcint_coefficients"), "actual libcint coefficients")
        _require(exponents.ndim == 1 and 0 < len(exponents) <= 256 and np.all(exponents > 0), "Invalid primitive exponents")
        _require(coeff.ndim == 2 and coeff.shape[0] == len(exponents) and 0 < coeff.shape[1] <= nao,
                 "Invalid contraction dimensions")
        _require(libcint.shape == coeff.shape, "Actual basis coefficient dimensions mismatch")
        _integer(shell.get("kappa"), "basis kappa", minimum=-13)
        degeneracy = 2 * angular + 1 if basis["representation"] == "spherical" else (angular + 1) * (angular + 2) // 2
        expanded_nao += degeneracy * coeff.shape[1]
    _require(expanded_nao == nao, "Expanded shell dimensions do not match AO dimension")
    s, eigenvalues = _metric(snapshot.get("ao_overlap"), nao)
    source = snapshot.get("source_occupations")
    _require(isinstance(source, dict), "Missing source occupations")
    if snapshot["reference"] == "RKS":
        occ = _array(source.get("spatial"), "RKS occupations")
        _require(occ.ndim == 1 and 0 < len(occ) <= nao and np.all((occ == 0) | (occ == 2)),
                 "RKS requires exact integer 0/2 occupations")
        source_occ = {"alpha": occ / 2, "beta": occ / 2}
    else:
        source_occ = {name: _array(source.get(name), f"{name} occupations") for name in ("alpha", "beta")}
        for occ in source_occ.values():
            _require(occ.ndim == 1 and 0 < len(occ) <= nao and np.all((occ == 0) | (occ == 1)),
                     "UKS requires exact integer 0/1 occupations")
    channels = snapshot.get("channels")
    _require(isinstance(channels, dict) and set(channels) == {"alpha", "beta"}, "Missing explicit spin channels")
    parsed, diagnostics = {}, {}
    for name, electrons in (("alpha", (nelec + spin) // 2), ("beta", (nelec - spin) // 2)):
        channel = channels[name]
        _require(isinstance(channel, dict), "Invalid spin channel")
        _require(type(channel.get("n_electrons")) is int and channel["n_electrons"] == electrons, "Spin electron count mismatch")
        occupied = np.flatnonzero(source_occ[name] == 1).tolist()
        saved_indices = channel.get("occupied_indices")
        _require(isinstance(saved_indices, list) and all(type(x) is int for x in saved_indices) and
                 len(occupied) == electrons and saved_indices == occupied, "Occupied index/electron mismatch")
        c = _array(channel.get("occupied_coefficients"), f"{name} occupied coefficients", (nao, electrons))
        residual = float(np.max(np.abs(c.T @ s @ c - np.eye(electrons)))) if electrons else 0.0
        _require(residual <= NUMERICAL_TOLERANCE, "Occupied coefficients are not AO-metric orthonormal")
        density = c @ c.T
        trace = float(np.trace(density @ s))
        _require(abs(trace - electrons) <= NUMERICAL_TOLERANCE * max(1, electrons), "Density electron trace mismatch")
        parsed[name] = (c, density)
        diagnostics[name] = {"metric_orthonormality_max_abs_residual": residual, "density_electron_trace": trace}
    if snapshot["reference"] == "RKS":
        _require(channels["alpha"] == channels["beta"], "RKS alpha/beta mapping must share the occupied spatial orbitals")
    diagnostics["ao_overlap_min_eigenvalue"] = float(eigenvalues[0])
    diagnostics["ao_overlap_condition_number"] = float(eigenvalues[-1] / eigenvalues[0])
    _json_bytes(snapshot)
    return s, parsed, diagnostics


def validate_snapshot(snapshot):
    """Return arithmetic diagnostics, raising on unsupported/missing evidence."""
    return _validate(snapshot)[2]


def capture_snapshot(mean_field, *, settings: dict, call_id: str, coordinate_frame_id=None):
    """Copy a converged real molecular RKS/UKS determinant, without solving it.

    The AO overlap may require one-electron integrals. Basis shells are captured
    using Mole's public per-shell accessors; no executable molecule dump is used.
    """
    try:
        from pyscf.gto import mole as mole_module
        from pyscf import lib
        _require(isinstance(mean_field.converged, (bool, np.bool_)) and bool(mean_field.converged), "Cannot capture unconverged SCF")
        reference = "UKS" if mean_field.istype("UKS") else "RKS" if mean_field.istype("RKS") else None
        _require(reference is not None, "Only molecular RKS/UKS capture is supported")
        mol = mean_field.mol
        _require(not getattr(mol, "ecp", None) and not getattr(mol, "pseudo", None), "Only all-electron molecular capture is supported")
        _require(not hasattr(mol, "lattice_vectors"), "Periodic capture is unsupported")
        nao = int(mol.nao_nr())
        _require(0 < nao <= MAX_AO, f"AO dimension exceeds {MAX_AO}")
        coeff, occ = _array(mean_field.mo_coeff, "MO coefficients"), _array(mean_field.mo_occ, "MO occupations")
        if reference == "RKS":
            _require(coeff.ndim == 2 and occ.ndim == 1 and coeff.shape == (nao, len(occ)), "RKS MO dimensions mismatch")
            source = {"spatial": occ.tolist()}
            pairs = {"alpha": (coeff, occ / 2), "beta": (coeff, occ / 2)}
        else:
            _require(coeff.ndim == 3 and coeff.shape[0] == 2 and occ.ndim == 2 and occ.shape[0] == 2 and
                     coeff.shape[1:] == (nao, occ.shape[1]), "UKS MO dimensions mismatch")
            source = {"alpha": occ[0].tolist(), "beta": occ[1].tolist()}
            pairs = {"alpha": (coeff[0], occ[0]), "beta": (coeff[1], occ[1])}
        channels = {}
        for i, (name, (c, o)) in enumerate(pairs.items()):
            indices = np.flatnonzero(o == 1)
            channels[name] = {"n_electrons": int(mol.nelec[i]), "occupied_indices": indices.tolist(),
                              "occupied_coefficients": c[:, indices].tolist()}
        symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
        raw_overlap = mean_field.get_ovlp()
        snapshot = {
            "schema_version": 1, "kind": "mean_field_electronic_snapshot", "reference": reference,
            "phase": "converged_scf_only",
            "producer": {"engine": "PySCF", "version": version("pyscf"),
                         "class": f"{type(mean_field).__module__}.{type(mean_field).__qualname__}"},
            "call_id": call_id, "settings": deepcopy(settings), "scf_converged": True,
            "scf_energy_hartree": float(mean_field.e_tot),
            "geometry": {"symbols": symbols, "atomic_numbers": [int(x) for x in mol.atom_charges()],
                         "positions_bohr": _array(mol.atom_coords(unit="Bohr"), "solver coordinates").tolist(),
                         "bohr_angstrom": float(lib.param.BOHR), "coordinate_frame_id": coordinate_frame_id,
                         "charge": int(mol.charge), "spin": int(mol.spin),
                         "frame": "solver_cartesian", "pbc": [False, False, False]},
            "basis": {"n_ao": nao, "representation": "cartesian" if mol.cart else "spherical",
                      "expanded_basis": deepcopy(mol._basis), "normalize_gto": bool(mole_module.NORMALIZE_GTO),
                      "ecp": False, "pseudo": False,
                      "ao_labels": [list(x) for x in mol.ao_labels(fmt=False)],
                      "shells": [{"atom_index": int(mol.bas_atom(i)), "angular_momentum": int(mol.bas_angular(i)),
                                  "kappa": int(mol.bas_kappa(i)),
                                  "exponents": _array(mol.bas_exp(i), "basis exponents").tolist(),
                                  "contraction_coefficients": _array(mol.bas_ctr_coeff(i), "basis contractions").tolist(),
                                  "libcint_coefficients": _array(mol._libcint_ctr_coeff(i), "actual basis contractions").tolist()}
                                 for i in range(mol.nbas)]},
            "ao_overlap": _array(raw_overlap, "AO overlap").tolist(),
            "source_array_dtypes": {"mo_coeff": str(np.asarray(mean_field.mo_coeff).dtype),
                                    "mo_occ": str(np.asarray(mean_field.mo_occ).dtype),
                                    "ao_overlap": str(np.asarray(raw_overlap).dtype)},
            "source_occupations": source, "channels": channels,
            "electronic_state_identity_verified": False, "ground_state_verified": False,
            "stability_assessed_by_capture": False,
        }
        validate_snapshot(snapshot)
        return snapshot
    except ElectronicStateError:
        raise
    except (AttributeError, TypeError, ValueError, OverflowError, IndexError, np.linalg.LinAlgError) as error:
        raise ElectronicStateError(f"Cannot capture electronic snapshot: {error}") from error


@contextmanager
def _parent(path):
    """Pin existing non-symlink parent directories; never create parents."""
    p = Path(os.path.abspath(os.fspath(path)))
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
    fd = os.open(p.anchor, flags)
    try:
        for part in p.parts[1:-1]:
            child = os.open(part, flags, dir_fd=fd)
            os.close(fd)
            fd = child
        yield fd, p.name, p
    finally:
        os.close(fd)


def _bound(max_bytes):
    _require(type(max_bytes) is int and 0 < max_bytes <= DEFAULT_MAX_BYTES, "max_bytes must be positive and at most 16 MiB")


def save_snapshot(path, snapshot, *, max_bytes=DEFAULT_MAX_BYTES):
    """Exclusively create one bounded snapshot file. Never overwrite or repair."""
    _bound(max_bytes)
    validate_snapshot(snapshot)
    raw = _json_bytes(snapshot)
    _require(len(raw) <= max_bytes, "Electronic snapshot exceeds byte bound")
    try:
        with _parent(path) as (parent, name, absolute):
            fd = os.open(name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600, dir_fd=parent)
            with os.fdopen(fd, "wb") as stream:
                stream.write(raw)
                stream.flush()
                os.fsync(stream.fileno())
    except OSError as error:
        raise ElectronicStateError(f"Cannot save snapshot: {error}") from error
    return {"path": str(absolute), "sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result, f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def read_snapshot(path, *, max_bytes=DEFAULT_MAX_BYTES):
    """Read bounded JSON only. Symlinks and nonregular files are rejected."""
    _bound(max_bytes)
    try:
        with _parent(path) as (parent, name, _):
            fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=parent)
            with os.fdopen(fd, "rb") as stream:
                info = os.fstat(stream.fileno())
                _require(stat.S_ISREG(info.st_mode) and info.st_size <= max_bytes, "Snapshot is nonregular or exceeds byte bound")
                raw = stream.read(max_bytes + 1)
                _require(len(raw) <= max_bytes, "Snapshot exceeds byte bound")
        return read_snapshot_bytes(raw, max_bytes=max_bytes)
    except ElectronicStateError:
        raise
    except (OSError, ValueError, UnicodeError, RecursionError) as error:
        raise ElectronicStateError(f"Cannot read snapshot: {error}") from error


def read_snapshot_bytes(raw, *, max_bytes=DEFAULT_MAX_BYTES):
    """Validate the exact already-captured bytes used for a caller's checksum."""
    _bound(max_bytes)
    _require(type(raw) is bytes, "Snapshot input must be immutable bytes")
    _require(0 < len(raw) <= max_bytes, "Snapshot is empty or exceeds byte bound")
    try:
        value = json.loads(raw, object_pairs_hook=_unique_object)
        validate_snapshot(value)
        return value
    except ElectronicStateError:
        raise
    except (ValueError, UnicodeError, RecursionError) as error:
        raise ElectronicStateError(f"Cannot parse snapshot bytes: {error}") from error


def compare_snapshots(left, right):
    """Compare same-context determinants; unsupported cases return no metrics.

    Only initial guess, requested threads and memory are ignored in the strict
    settings comparison. No physical identity threshold is chosen.
    """
    report = {"comparison_status": "insufficient_evidence", "reasons": [], "metrics": None,
              "electronic_state_identity_verified": False, "ground_state_verified": False,
              "comparison_scope": "same_geometry_real_integer_occupied_RKS_UKS"}
    try:
        sl, pl, dl = _validate(left)
        sr, pr, dr = _validate(right)
    except ElectronicStateError as error:
        report["reasons"] = [str(error)]
        return report
    reasons = []
    if left["geometry"] != right["geometry"]:
        reasons.append("Geometry, atom order, charge/spin or frame differs; cross-geometry AO overlap is not implemented")
    if left["basis"] != right["basis"]:
        reasons.append("Expanded basis or ordered AO convention differs")
    if left["reference"] != right["reference"]:
        reasons.append("Reference type differs")
    if left["producer"] != right["producer"]:
        reasons.append("Producer engine/version/class differs")
    ignored = {"scf_initial_guess", "threads", "memory_mb"}
    settings_l = {k: v for k, v in left["settings"].items() if k not in ignored}
    settings_r = {k: v for k, v in right["settings"].items() if k not in ignored}
    if settings_l != settings_r:
        reasons.append("Calculation settings differ beyond initial guess/requested resources")
    if not np.array_equal(sl, sr):
        reasons.append("Same-geometry AO overlap matrices differ; this version requires exact metric equality")
    if reasons:
        report.update(comparison_status="incompatible_context", reasons=reasons)
        return report
    metrics = {}
    try:
        for name in ("alpha", "beta"):
            cl, density_l = pl[name]
            cr, density_r = pr[name]
            overlap = cl.T @ sl @ cr
            singular = np.linalg.svd(overlap, compute_uv=False)
            _require(np.all(singular <= 1 + 4 * NUMERICAL_TOLERANCE), "Subspace overlap exceeds numerical allowance")
            bounded = np.clip(singular, 0, 1)
            diff = density_l - density_r
            density_distance_squared = float(np.trace(diff @ sl @ diff @ sl))
            _require(math.isfinite(density_distance_squared) and density_distance_squared >= -NUMERICAL_TOLERANCE,
                     "Invalid metric density distance")
            metrics[name] = {"n_occupied": cl.shape[1], "singular_values": singular.tolist(),
                             "principal_angles_radians": np.arccos(bounded).tolist(),
                             "singular_value_clipping_max_abs": float(np.max(np.abs(singular - bounded))) if len(singular) else 0.0,
                             "density_distance_squared": max(0.0, density_distance_squared)}
    except (ElectronicStateError, np.linalg.LinAlgError) as error:
        report["reasons"] = [str(error)]
        return report
    report.update(comparison_status="compared", metrics=metrics,
                  call_ids={"left": left["call_id"], "right": right["call_id"]},
                  input_digests={"left": hashlib.sha256(_json_bytes(left)).hexdigest(),
                                 "right": hashlib.sha256(_json_bytes(right)).hexdigest()},
                  validation={"left": dl, "right": dr},
                  numerical_tolerance=NUMERICAL_TOLERANCE,
                  ao_overlap_max_abs_difference=float(np.max(np.abs(sl - sr))),
                  settings_fields_excluded_from_compatibility=sorted(ignored),
                  comparison_axis={key: {"left": left["settings"].get(key), "right": right["settings"].get(key)}
                                   for key in sorted(ignored) if left["settings"].get(key) != right["settings"].get(key)})
    return report
