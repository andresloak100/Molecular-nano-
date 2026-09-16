"""Offline checks of explicit numerical choices and exact paired-grid identity.

All fields are producer-supplied evidence. Matching hashes do not authenticate
expanded basis content, absent raw grids, hardware execution, or chemical states.
No quantum backend is imported or executed here.
"""

from __future__ import annotations

import hashlib
import math
import re
import struct


GRID_CANONICALIZATION = "paired-xyz-weight-f64le-sorted-v1"
PLACEHOLDERS = {"", "default", "auto", "unknown", "unresolved", "not_recorded"}
ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split())
_MISSING = object()


def canonical_grid_signature(rows):
    """Hash an actual grid without decoupling positions from their weights.

    ``rows`` must be a nonempty list/tuple of (x, y, z, weight) rows containing
    finite JSON numbers. Coordinates must be in Bohr and weights in Bohr**3.
    Values are converted to IEEE-754 binary64, every signed zero becomes +0.0,
    and full paired rows are sorted lexicographically. The byte stream is one
    little-endian unsigned 64-bit row count followed by little-endian binary64
    rows (four values each), with no header or separators. Duplicate rows remain.
    There is no rounding, quantization, tolerance, or independent weight sorting.

    Order alone cannot change this signature. Any changed binary64 grid value
    can change it, including harmless roundoff; such grids do not pass this
    exact-identity protocol even if their energies are close.
    """
    if type(rows) not in (list, tuple) or not rows:
        raise ValueError("grid rows must be a nonempty list or tuple")
    canonical = []
    for index, row in enumerate(rows):
        if type(row) not in (list, tuple) or len(row) != 4:
            raise ValueError(f"grid row {index} must contain x, y, z, weight")
        converted = []
        for value in row:
            if type(value) not in (int, float):
                raise ValueError(f"grid row {index} must contain finite numbers, not booleans")
            try:
                number = float(value)
            except (OverflowError, ValueError) as error:
                raise ValueError(f"grid row {index} cannot be represented as float64") from error
            if not math.isfinite(number):
                raise ValueError(f"grid row {index} must contain finite numbers")
            converted.append(0.0 if number == 0 else number)
        canonical.append(tuple(converted))
    canonical.sort()
    digest = hashlib.sha256(struct.pack("<Q", len(canonical)))
    for row in canonical:
        digest.update(struct.pack("<dddd", *row))
    return digest.hexdigest()


def _issue(errors, kind, field, message):
    errors.append(f"{kind}: {field}: {message}")


def _object(value, field, errors):
    if value is _MISSING or value is None:
        _issue(errors, "unavailable", field, "resolved object is required")
        return None
    if type(value) is not dict:
        _issue(errors, "malformed", field, "must be an object")
        return None
    return value


def _identity(value, field, errors):
    if value is _MISSING or value is None:
        _issue(errors, "unavailable", field, "explicit resolved identity is required")
        return None
    if type(value) is not str:
        _issue(errors, "malformed", field, "must be a string")
        return None
    normalized = value.strip()
    if normalized.casefold().replace("-", "_").replace(" ", "_") in PLACEHOLDERS:
        _issue(errors, "unavailable", field, "placeholder is not a resolved identity")
        return None
    return normalized


def _hash(value, field, errors):
    normalized = _identity(value, field, errors)
    if normalized is None:
        return None
    if re.fullmatch(r"[0-9a-fA-F]{64}", normalized) is None:
        _issue(errors, "malformed", field, "must be a 64-digit SHA-256 hex string")
        return None
    return normalized.lower()


def _integer(value, field, errors, minimum=1):
    if value is _MISSING or value is None:
        _issue(errors, "unavailable", field, "explicit integer is required")
        return None
    if type(value) is not int or value < minimum:
        _issue(errors, "malformed", field, f"must be an integer >= {minimum}; booleans are invalid")
        return None
    return value


def _resolved(record, label, case, errors):
    record = _object(record, label, errors)
    if record is None:
        return None
    root = _object(record.get("resolved_numerics", _MISSING), f"{label}.resolved_numerics", errors)
    if root is None:
        return None
    result = {}
    precision = _object(root.get("precision", _MISSING), f"{label}.precision", errors)
    if precision is not None:
        result["precision"] = {}
        for key in ("energy", "forces", "density"):
            value = _identity(precision.get(key, _MISSING), f"{label}.precision.{key}", errors)
            result["precision"][key] = value
            if value is not None and value != "float64":
                _issue(errors, "mismatch", f"{label}.precision.{key}", "this protocol requires float64")
    basis = _object(root.get("orbital_basis", _MISSING), f"{label}.orbital_basis", errors)
    if basis is not None:
        result["orbital_basis"] = _hash(basis.get("sha256", _MISSING), f"{label}.orbital_basis.sha256", errors)
    auxiliary = _object(root.get("auxiliary_basis", _MISSING), f"{label}.auxiliary_basis", errors)
    if auxiliary is not None:
        status = _identity(auxiliary.get("status", _MISSING), f"{label}.auxiliary_basis.status", errors)
        expected = "resolved" if case["settings"]["density_fit"] else "not_applicable"
        if status is not None and status != expected:
            _issue(errors, "mismatch", f"{label}.auxiliary_basis.status", f"must be {expected!r} for requested density_fit")
        result["auxiliary_basis"] = {"status": status}
        if case["settings"]["density_fit"]:
            result["auxiliary_basis"]["sha256"] = _hash(auxiliary.get("sha256", _MISSING), f"{label}.auxiliary_basis.sha256", errors)
        elif "sha256" in auxiliary:
            _issue(errors, "malformed", f"{label}.auxiliary_basis.sha256", "must be absent when density fitting is not applicable")
    quadrature = _object(root.get("quadrature", _MISSING), f"{label}.quadrature", errors)
    if quadrature is not None:
        normalized = {}
        result["quadrature"] = normalized
        level = _integer(quadrature.get("level", _MISSING), f"{label}.quadrature.level", errors, minimum=0)
        normalized["level"] = level
        if level is not None and level != case["settings"]["grid_level"]:
            _issue(errors, "mismatch", f"{label}.quadrature.level", "does not match requested grid_level")
        normalized["point_count"] = _integer(quadrature.get("point_count", _MISSING), f"{label}.quadrature.point_count", errors)
        for key in ("pruning", "radial_method", "becke_scheme", "canonicalization"):
            normalized[key] = _identity(quadrature.get(key, _MISSING), f"{label}.quadrature.{key}", errors)
        if normalized["canonicalization"] is not None and normalized["canonicalization"] != GRID_CANONICALIZATION:
            _issue(errors, "mismatch", f"{label}.quadrature.canonicalization", f"requires {GRID_CANONICALIZATION!r}")
        normalized["signature_sha256"] = _hash(quadrature.get("signature_sha256", _MISSING), f"{label}.quadrature.signature_sha256", errors)
        for key in ("grid_points_sha256", "grid_weights_sha256"):
            if key in quadrature:
                _hash(quadrature[key], f"{label}.quadrature.{key}", errors)
        grids = _object(quadrature.get("atom_grid", _MISSING), f"{label}.quadrature.atom_grid", errors)
        if grids is not None:
            normalized["atom_grid"] = {}
            if not grids:
                _issue(errors, "unavailable", f"{label}.quadrature.atom_grid", "resolved per-element counts are required")
            for element, counts in grids.items():
                if type(element) is not str or element not in ELEMENTS:
                    _issue(errors, "malformed", f"{label}.quadrature.atom_grid", "keys must be chemical element symbols")
                    continue
                fields = _object(counts, f"{label}.quadrature.atom_grid.{element}", errors)
                if fields is None:
                    continue
                normalized["atom_grid"][element] = {
                    key: _integer(fields.get(key, _MISSING), f"{label}.quadrature.atom_grid.{element}.{key}", errors)
                    for key in ("radial", "angular")
                }
            missing = sorted(set(case["symbols"]) - set(grids))
            if missing:
                _issue(errors, "unavailable", f"{label}.quadrature.atom_grid", f"missing case elements {missing}")
    return result


def validate_resolved_numerics(cpu, gpu, case):
    """Return unavailable/malformed/mismatch strings; an empty list passes.

    Required record.resolved_numerics schema:
      precision: {energy: float64, forces: float64, density: float64}
      orbital_basis: {sha256: <actual expanded orbital basis digest>}
      auxiliary_basis: {status: not_applicable} without density fitting, or
                       {status: resolved, sha256: <actual expanded basis digest>}
      quadrature: {level, atom_grid: {element: {radial, angular}}, pruning,
                   radial_method, becke_scheme, point_count, canonicalization,
                   signature_sha256}
    Integer counts are positive; level matches the requested case setting.
    Every case element must occur in atom_grid. Canonicalization is exactly
    GRID_CANONICALIZATION and signature_sha256 comes from canonical_grid_signature.
    Optional grid_points_sha256/grid_weights_sha256 must be valid digests when
    supplied but are not compared: paired canonical rows determine grid identity.
    Other descriptive fields are permitted and are not identity evidence.

    The caller separately validates geometry/plan binding and all energy/force
    evidence. Case requires symbols, atom_count and settings.grid_level/density_fit.
    No missing field is filled from a default or the other backend's record.
    """
    errors = []
    case = _object(case, "case", errors)
    if case is None:
        return errors
    symbols = case.get("symbols")
    if type(symbols) is not list or not symbols:
        _issue(errors, "unavailable", "case.symbols", "ordered element symbols are required")
    elif any(type(symbol) is not str or symbol not in ELEMENTS for symbol in symbols):
        _issue(errors, "malformed", "case.symbols", "must contain chemical element symbols")
    atom_count = _integer(case.get("atom_count", _MISSING), "case.atom_count", errors)
    if type(symbols) is list and atom_count is not None and len(symbols) != atom_count:
        _issue(errors, "mismatch", "case.atom_count", "does not equal the number of symbols")
    settings = _object(case.get("settings", _MISSING), "case.settings", errors)
    if settings is not None:
        level = _integer(settings.get("grid_level", _MISSING), "case.settings.grid_level", errors, minimum=0)
        if level is not None and level > 9:
            _issue(errors, "malformed", "case.settings.grid_level", "must be at most 9")
        density_fit = settings.get("density_fit", _MISSING)
        if density_fit is _MISSING or density_fit is None:
            _issue(errors, "unavailable", "case.settings.density_fit", "explicit boolean is required")
        elif type(density_fit) is not bool:
            _issue(errors, "malformed", "case.settings.density_fit", "must be a boolean")
    if errors:
        return errors
    cpu_resolved = _resolved(cpu, "CPU", case, errors)
    gpu_resolved = _resolved(gpu, "GPU", case, errors)
    if errors:
        return errors
    for key in ("precision", "orbital_basis", "auxiliary_basis"):
        if cpu_resolved[key] != gpu_resolved[key]:
            _issue(errors, "mismatch", key, "CPU and GPU resolved choices differ")
    for key in ("level", "atom_grid", "pruning", "radial_method", "becke_scheme",
                "point_count", "canonicalization", "signature_sha256"):
        if cpu_resolved["quadrature"][key] != gpu_resolved["quadrature"][key]:
            note = "exact grid identity is unverified; canonical paired-grid digests differ" if key == "signature_sha256" else "CPU and GPU resolved choices differ"
            _issue(errors, "mismatch", f"quadrature.{key}", note)
    return errors
