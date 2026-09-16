"""Reconstruct bound E3 basis functions and compute cross-AO overlaps, without SCF."""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import stat
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from nanodesign.electronic_state import ElectronicStateError, validate_snapshot

MAX_BYTES = 16 * 1024 * 1024
MAX_ATOMS = 256
MAX_AO = 512
OVERLAP_ATOL = 1e-10
OVERLAP_METRIC_TOL = 1e-8
BASIS_ATOL = 1e-12
BASIS_RTOL = 1e-12


class BridgeError(ValueError):
    """Unsupported or inconsistently bound snapshot pair; no matrix is accepted."""


def require(condition, message):
    if not condition:
        raise BridgeError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def _unique(pairs):
    obj = {}
    for key, value in pairs:
        require(key not in obj, f"Duplicate JSON key: {key}")
        obj[key] = value
    return obj


def _read(path):
    """Hash and validate the same bounded bytes, with pinned non-symlink parents."""
    path = Path(os.path.abspath(os.fspath(path)))
    fd = os.open(path.anchor, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        for part in path.parts[1:-1]:
            child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
            os.close(fd)
            fd = child
        payload_fd = os.open(path.name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
        with os.fdopen(payload_fd, "rb") as stream:
            info = os.fstat(stream.fileno())
            require(stat.S_ISREG(info.st_mode) and info.st_size <= MAX_BYTES, "Snapshot is nonregular or exceeds 16 MiB")
            raw = stream.read(MAX_BYTES + 1)
            require(len(raw) <= MAX_BYTES, "Snapshot exceeds 16 MiB")
    finally:
        os.close(fd)
    snapshot = json.loads(raw, object_pairs_hook=_unique)
    require(isinstance(snapshot, dict), "Snapshot must be an object")
    basis = snapshot.get("basis", {})
    geometry = snapshot.get("geometry", {})
    require(isinstance(basis, dict) and isinstance(geometry, dict), "Missing basis/geometry")
    require(type(basis.get("n_ao")) is int and 0 < basis["n_ao"] <= MAX_AO, "X1 supports 1..512 AOs")
    require(isinstance(geometry.get("symbols"), list) and 0 < len(geometry["symbols"]) <= MAX_ATOMS,
            "X1 supports 1..256 atoms")
    validate_snapshot(snapshot)
    return snapshot, {"path": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
                      "size_bytes": len(raw), "canonical_content_sha256": digest(snapshot),
                      "call_id": snapshot["call_id"], "geometry_sha256": digest(snapshot["geometry"]),
                      "basis_sha256": digest(snapshot["basis"])}


def _numeric(value, name):
    a = np.asarray(value)
    require(a.dtype.kind in "iuf" and np.isfinite(a).all(), f"{name} must contain finite real numbers")
    return a.astype(float)


def _check_expanded(basis, symbols):
    """Reject string aliases/expressions and bound the actual constructor input."""
    expanded = basis["expanded_basis"]
    require(isinstance(expanded, dict) and set(expanded) == set(symbols),
            "Expanded basis must explicitly cover exactly the untagged element symbols")
    predicted_nao = 0
    for symbol in symbols:
        shells = expanded[symbol]
        require(isinstance(shells, list) and 0 < len(shells) <= MAX_AO, "Invalid expanded basis shell list")
        for shell in shells:
            require(isinstance(shell, list) and 2 <= len(shell) <= 257, "Invalid expanded shell")
            angular = shell[0]
            require(type(angular) is int and 0 <= angular <= 6, "X1 supports nonrelativistic shells l=0..6")
            require(isinstance(shell[1], list), "Kappa/spinor shell construction is unsupported")
            rows = _numeric(shell[1:], "expanded primitive rows")
            require(rows.ndim == 2 and 2 <= rows.shape[1] <= MAX_AO + 1 and np.all(rows[:, 0] > 0),
                    "Invalid primitive exponents or contractions")
            require(np.all(np.any(rows[:, 1:] != 0, axis=0)), "Zero basis contraction")
            degeneracy = 2 * angular + 1 if basis["representation"] == "spherical" else (angular+1)*(angular+2)//2
            predicted_nao += degeneracy * (rows.shape[1]-1)
    require(predicted_nao == basis["n_ao"], "Expanded basis constructor dimensions differ from capture")


def _check_reconstructed_overlap(actual, recorded):
    """Absolute and relative quadratic-form checks; do not hide soft AO errors."""
    actual, recorded = np.asarray(actual), np.asarray(recorded)
    require(actual.shape == recorded.shape and np.isfinite(actual).all(), "Invalid reconstructed self-overlap")
    residual = float(np.max(np.abs(actual - recorded)))
    require(residual <= OVERLAP_ATOL, "Reconstructed self-overlap differs from capture")
    values, vectors = np.linalg.eigh((recorded + recorded.T) / 2)
    require(np.all(values > 0), "Captured AO metric is not positive definite")
    whitener = vectors / np.sqrt(values)
    relative = whitener.T @ (actual - recorded) @ whitener
    relative_residual = float(np.linalg.norm(relative, ord=2))
    require(np.isfinite(relative_residual) and relative_residual <= OVERLAP_METRIC_TOL,
            "Reconstructed self-overlap differs in the captured AO metric")
    return residual, relative_residual


def _reconstruct(snapshot):
    import pyscf
    from pyscf import gto, lib

    geo, basis = snapshot["geometry"], snapshot["basis"]
    require(snapshot["producer"]["version"] == pyscf.__version__, "Captured/runtime PySCF versions differ")
    require(geo["bohr_angstrom"] == float(lib.param.BOHR), "Captured/runtime Bohr conversion differs")
    require(basis["normalize_gto"] == bool(gto.mole.NORMALIZE_GTO), "Gaussian normalization convention differs")
    require(basis["ecp"] is False and basis["pseudo"] is False, "ECP/pseudo unsupported")
    require(all(shell["kappa"] == 0 for shell in basis["shells"]), "Kappa/spinor shell evidence unsupported")
    _check_expanded(basis, geo["symbols"])
    # All strings here are validated chemical symbols; all coordinates and shell
    # primitives are finite numbers. No arbitrary object dump or kwargs are used.
    mol = gto.Mole()
    mol.atom = list(zip(geo["symbols"], deepcopy(geo["positions_bohr"])))
    mol.unit = "Bohr"
    mol.basis = deepcopy(basis["expanded_basis"])
    mol.charge, mol.spin = geo["charge"], geo["spin"]
    mol.cart = basis["representation"] == "cartesian"
    mol.symmetry, mol.verbose = False, 0
    mol.build(parse_arg=False, dump_input=False)
    require(mol.nao_nr() == basis["n_ao"] and mol.nbas == len(basis["shells"]), "Rebuilt AO/shell count differs")
    require([list(label) for label in mol.ao_labels(fmt=False)] == basis["ao_labels"], "Rebuilt ordered AO labels differ")
    require(np.array_equal(mol.atom_coords(unit="Bohr"), np.asarray(geo["positions_bohr"])), "Rebuilt solver coordinates differ")
    require(mol.atom_charges().tolist() == geo["atomic_numbers"], "Rebuilt ordered nuclear charges differ")
    expected_nelec = tuple(snapshot["channels"][s]["n_electrons"] for s in ("alpha", "beta"))
    require(tuple(mol.nelec) == expected_nelec, "Rebuilt spin electron counts differ")
    max_abs, max_relative = 0.0, 0.0
    for i, shell in enumerate(basis["shells"]):
        require((mol.bas_atom(i), mol.bas_angular(i), mol.bas_kappa(i)) ==
                (shell["atom_index"], shell["angular_momentum"], shell["kappa"]), "Rebuilt shell order/type differs")
        for name, actual in (("exponents", mol.bas_exp(i)), ("contraction_coefficients", mol.bas_ctr_coeff(i)),
                             ("libcint_coefficients", mol._libcint_ctr_coeff(i))):
            recorded = np.asarray(shell[name])
            require(actual.shape == recorded.shape and np.allclose(actual, recorded, atol=BASIS_ATOL, rtol=BASIS_RTOL),
                    f"Rebuilt shell {i} {name} differs from actual captured functions")
            residual = np.abs(actual - recorded)
            max_abs = max(max_abs, float(residual.max()))
            max_relative = max(max_relative, float(np.max(residual / np.maximum(np.abs(recorded), BASIS_ATOL))))
    operator = "int1e_ovlp_cart" if mol.cart else "int1e_ovlp_sph"
    overlap = mol.intor_symmetric(operator)
    residual, metric_residual = _check_reconstructed_overlap(overlap, snapshot["ao_overlap"])
    return mol, overlap, {"self_overlap_max_abs_difference": residual,
                          "self_overlap_metric_relative_norm": metric_residual,
                          "shell_values_max_abs_difference": max_abs,
                          "shell_values_max_relative_difference": max_relative}


def build_cross_overlap(left_path, right_path, *, common_frame_id, atom_mapping):
    """Compute S_AB for two bound E3 artifacts; caller explicitly declares frame.

    atom_mapping lists each left atom's right index. Only identity order is
    supported. A frame label is a provenance assertion, not a measured alignment.
    """
    try:
        require(isinstance(common_frame_id, str) and bool(common_frame_id.strip()) and len(common_frame_id) <= 4096,
                "An explicit nonempty common_frame_id is required")
        left, left_source = _read(left_path)
        right, right_source = _read(right_path)
        require(left["call_id"] != right["call_id"] or left_source["canonical_content_sha256"] == right_source["canonical_content_sha256"],
                "Same call_id is attached to contradictory snapshot content")
        a, b = left["geometry"], right["geometry"]
        require(isinstance(atom_mapping, list) and all(type(i) is int for i in atom_mapping) and
                atom_mapping == list(range(len(a["symbols"]))), "Explicit identity ordered atom_mapping is required")
        for key in ("symbols", "atomic_numbers", "charge", "spin", "frame", "pbc", "bohr_angstrom"):
            require(a.get(key) == b.get(key), f"Pair geometry context differs: {key}")
        for geo in (a, b):
            require(geo.get("coordinate_frame_id") in (None, common_frame_id), "Captured coordinate frame conflicts with declaration")
        require(left["basis"] == right["basis"], "Pair actual basis or AO convention differs")
        require(left["reference"] == right["reference"] and left["producer"] == right["producer"], "Pair reference/producer differs")
        allowed_axes = {"scf_initial_guess", "threads", "memory_mb"}
        physical = lambda x: {k:v for k,v in x["settings"].items() if k not in allowed_axes}
        require(physical(left) == physical(right), "Pair method/settings differ beyond initial guess/requested resources")
        lm, ls, lv = _reconstruct(left)
        rm, rs, rv = _reconstruct(right)
        from pyscf import gto, __version__ as pyscf_version
        operator = "int1e_ovlp_cart" if lm.cart else "int1e_ovlp_sph"
        cross = gto.intor_cross(operator, lm, rm)
        require(cross.shape == (lm.nao_nr(), rm.nao_nr()) and np.isfinite(cross).all(), "Invalid computed cross-overlap")
        if np.array_equal(a["positions_bohr"], b["positions_bohr"]):
            require(np.allclose(cross, ls, rtol=0, atol=OVERLAP_ATOL) and np.allclose(cross, rs, rtol=0, atol=OVERLAP_ATOL),
                    "Same-geometry cross-overlap does not reduce to the self-overlap")
        displacement = np.asarray(b["positions_bohr"]) - np.asarray(a["positions_bohr"])
        changed_axes = {key: {"left":left["settings"].get(key),"right":right["settings"].get(key)}
                        for key in sorted(allowed_axes) if left["settings"].get(key) != right["settings"].get(key)}
        return {"schema_version":1, "kind":"cross_ao_overlap_bridge", "status":"computed",
                "sources":{"left":left_source,"right":right_source},
                "pair_binding":{"common_frame_id":common_frame_id,"atom_mapping":atom_mapping,
                                "frame_evidence":"caller_declared; captured IDs checked when present",
                                "direction":"rows=left/bra AO order; columns=right/ket AO order",
                                "geometry_alignment_applied":False,"maximum_atom_displacement_bohr":float(np.linalg.norm(displacement,axis=1).max()),
                                "recorded_frame_ids":[a.get("coordinate_frame_id"),b.get("coordinate_frame_id")]},
                "operator":operator,"matrix_units":"dimensionless","cross_overlap":cross.tolist(),
                "cross_overlap_sha256":digest(cross.tolist()),
                "reconstruction":{"left":lv,"right":rv,"overlap_atol":OVERLAP_ATOL,
                                  "overlap_metric_relative_tolerance":OVERLAP_METRIC_TOL,
                                  "basis_atol":BASIS_ATOL,"basis_rtol":BASIS_RTOL,"pyscf_version":pyscf_version},
                "changed_comparison_axes":changed_axes,"scf_jobs_started":0,"gradient_jobs_started":0,
                "electronic_state_identity_verified":False,"ground_state_verified":False,"branch_continuity_verified":False,
                "limitations":["Reconstruction checks consistency of captured basis data, not its external authenticity.",
                               "Shared frame and ordered atom correspondence are caller declarations, not physical measurements.",
                               "A cross-AO overlap does not identify a physical electronic state or certify branch continuity.",
                               "Captured converged SCF does not imply force or whole-calculation success."]}
    except BridgeError:
        raise
    except (ElectronicStateError, OSError, ValueError, TypeError, KeyError, IndexError, OverflowError, RecursionError) as error:
        raise BridgeError(str(error)) from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", type=Path); parser.add_argument("right", type=Path)
    parser.add_argument("--common-frame-id", required=True)
    parser.add_argument("--atom-mapping", type=int, nargs="+", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = build_cross_overlap(args.left,args.right,common_frame_id=args.common_frame_id,atom_mapping=args.atom_mapping)
    except BridgeError as error:
        report = {"schema_version":1,"kind":"cross_ao_overlap_bridge","status":"not_computed",
                  "reason":str(error),"cross_overlap":None,"electronic_state_identity_verified":False,
                  "ground_state_verified":False,"branch_continuity_verified":False}
    report["bridge_implementation_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    encoded = json.dumps(report,indent=2,allow_nan=False)+"\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(encoded)
    else:
        print(encoded,end="")
    return 0 if report["status"] == "computed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
