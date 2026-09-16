"""Independently reconstruct a saved finite-difference characterization.

Reads explicit JSON only. No production/ASE/quantum imports or calculations.
Passing means internal numerical consistency, never physical validation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

WAVENUMBER_FACTOR = math.sqrt(1.602176634e-19 / (1e-20 * 1.66053906892e-27)) / (2 * math.pi * 299792458 * 100)
MAX_COORDINATES = 600
MAX_INPUT_BYTES = 64 * 1024 * 1024


class EvidenceError(ValueError):
    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


def demand(condition, code, message):
    if not condition:
        raise EvidenceError(code, message)


def number(value, label, *, positive=False):
    demand(type(value) in (float, int) and math.isfinite(value), "invalid_number", f"{label} must be a finite real number")
    demand(not positive or value > 0, "invalid_number", f"{label} must be positive")
    return float(value)


def array(value, shape, label):
    def numeric_leaves(item):
        if isinstance(item, list):
            return all(numeric_leaves(child) for child in item)
        return type(item) in (int, float) and math.isfinite(item)
    demand(isinstance(value, list) and numeric_leaves(value), "invalid_array", f"{label} must contain only finite real numbers (no booleans or strings)")
    result = np.asarray(value, dtype=float)
    demand(result.shape == shape, "invalid_shape", f"{label} expected shape {shape}, got {result.shape}")
    return result


def near(actual, expected, label, *, atol=1e-9, rtol=1e-9):
    demand(np.isfinite(actual).all() and np.isfinite(expected).all()
           and np.allclose(actual, expected, atol=atol, rtol=rtol), "arithmetic_mismatch", f"{label} differs from the saved-force reconstruction")


def integer(value, label):
    demand(type(value) is int, "invalid_integer", f"{label} must be an integer")
    return value


def index_list(value, count, label):
    demand(isinstance(value, list) and all(type(i) is int and 0 <= i < count for i in value)
           and len(value) == len(set(value)), "invalid_indices", f"{label} must have unique, in-range integer atom indices")
    return value


def finite_json_tree(value):
    """Reject nonfinite metadata too, including overflowing JSON exponents."""
    if isinstance(value, float):
        demand(math.isfinite(value), "nonfinite_json", "Record contains a nonfinite JSON number")
    elif isinstance(value, dict):
        for child in value.values():
            finite_json_tree(child)
    elif isinstance(value, list):
        for child in value:
            finite_json_tree(child)


def report_base():
    return {"schema_version": 1, "status": "unavailable", "numerical_reconstruction_verified": False,
            "scientific_model_validated": False, "electronic_state_verified": False,
            "transition_state_verified": False, "connectivity_verified": False,
            "findings": [], "limitations": [
                "Agreement checks saved-array consistency, not whether forces or the electronic method describe the molecule accurately.",
                "Recorded masses are used as supplied; isotope assignment and external geometry identity require separate evidence.",
                "No electronic event log is read; recorded call IDs are metadata, not proof that an electronic calculation completed.",
                "Step sensitivity, electronic-state identity, reaction connectivity and operating reliability are not established."]}


def reconstruct(record):
    """Validate the current standalone stationary schema and return diagnostics."""
    evidence = record["finite_difference_evidence"]
    demand(isinstance(evidence, dict), "invalid_evidence", "finite_difference_evidence must be an object")
    demand(evidence["units"] == {"length": "angstrom", "force": "eV/angstrom"}, "wrong_units", "Displacement evidence must use angstrom and eV/angstrom")
    demand(record["coordinate_order"] == "x,y,z for each atom in free_atom_indices", "coordinate_order", "Unsupported coordinate ordering")
    positions = evidence["reference_positions_angstrom"]
    demand(isinstance(positions, list) and len(positions) > 0, "invalid_shape", "Reference geometry must be nonempty")
    n = len(positions)
    reference = array(positions, (n, 3), "reference geometry")
    baseline = array(evidence["baseline_forces_ev_per_angstrom"], (n, 3), "baseline forces")
    free = index_list(record["free_atom_indices"], n, "free atoms")
    frozen = index_list(record["frozen_atom_indices"], n, "frozen atoms")
    demand(free and not set(free).intersection(frozen) and set(free + frozen) == set(range(n)),
           "coordinate_partition", "Free and frozen atoms must partition the complete reference geometry")
    d = len(free) * 3
    demand(d <= MAX_COORDINATES, "audit_size_limit", f"Audit supports at most {MAX_COORDINATES} free coordinates")
    demand(integer(record["free_coordinate_count"], "free coordinate count") == d, "coordinate_count", "Free-coordinate count disagrees with atom indices")
    demand(integer(record["force_requests"], "force requests") == 1 + 2 * d, "force_request_count", "Central stencil needs one baseline request and two per free coordinate")
    settings = record["settings"]
    step = number(settings["step_angstrom"], "step", positive=True)
    force_tolerance = number(settings["force_tolerance_ev_per_angstrom"], "force tolerance", positive=True)
    frequency_tolerance = number(settings["frequency_tolerance_cm1"], "frequency tolerance", positive=True)
    rows = evidence["displacements"]
    demand(isinstance(rows, list) and len(rows) == 2 * d, "stencil_count", "Exactly two displacement records per free coordinate are required")
    raw = np.empty((d, d))
    call_ids = [evidence.get("baseline_calculation_call_id")]
    representation_errors = []
    for column in range(d):
        atom, axis = free[column // 3], column % 3
        pair = []
        for side, sign in enumerate((1, -1)):
            row = rows[2 * column + side]
            demand(isinstance(row, dict), "invalid_stencil", "Every displacement must be an object")
            demand(integer(row["atom_index"], "displaced atom") == atom and integer(row["axis"], "displaced axis") == axis,
                   "stencil_order", "Displacements must be plus then minus for each ordered free coordinate")
            requested = number(row["requested_offset_angstrom"], "requested offset")
            actual = number(row["actual_offset_angstrom"], "actual offset")
            displaced = number(row["displaced_coordinate_angstrom"], "displaced coordinate")
            demand(requested == sign * step, "requested_offset", "Requested signed offsets must equal settings.step_angstrom")
            # These serialized binary64 values came from exactly these operations.
            # Avoid a loose absolute tolerance that would accept a zero stencil.
            demand(displaced == reference[atom, axis] + requested and actual == displaced - reference[atom, axis],
                   "displaced_coordinate", "Reference, requested/actual offset and displaced coordinate disagree")
            relative_error = abs(actual - requested) / step
            demand(sign * actual > 0 and math.isfinite(relative_error) and relative_error <= 1e-6,
                   "unresolved_displacement", "Actual displacement is missing, reversed or exceeds the current one-ppm representation guard")
            near(number(row["relative_step_representation_error"], "step representation error"), relative_error,
                 "step representation error", atol=1e-15, rtol=1e-8)
            representation_errors.append(relative_error)
            full_forces = array(row["forces_ev_per_angstrom"], (n, 3), "displaced full forces")
            pair.append(full_forces[free].ravel())
            call_ids.append(row.get("calculation_call_id"))
        # Match the recorded production stencil: denominator is REQUESTED 2h.
        # Actual offsets are audited above, not substituted into the estimator.
        raw[:, column] = -(pair[0] - pair[1]) / (2 * step)
    demand(np.isfinite(raw).all(), "nonfinite_reconstruction", "Reconstructed Hessian is nonfinite")
    demand(all(value is None or isinstance(value, str) for value in call_ids), "invalid_call_id", "Calculation IDs must be strings or explicit nulls")
    symmetric = (raw + raw.T) / 2
    asymmetry = raw - raw.T
    saved_hessian = array(record["hessian_ev_per_angstrom2"], (d, d), "saved Hessian")
    near(saved_hessian, symmetric, "symmetric Hessian")
    demand(record["hessian_symmetrized"] is True, "hessian_scope", "Current schema must label the saved Hessian as symmetrized")
    max_asymmetry = float(np.abs(asymmetry).max())
    raw_norm = float(np.linalg.norm(raw))
    relative_asymmetry = float(np.linalg.norm(asymmetry)) / raw_norm if raw_norm else 0.0
    near(number(record["hessian_asymmetry_max_abs_ev_per_angstrom2"], "maximum asymmetry"), max_asymmetry, "maximum raw asymmetry")
    near(number(record["hessian_asymmetry_relative_frobenius"], "relative asymmetry"), relative_asymmetry, "relative raw asymmetry")
    masses = array(record["free_masses_amu"], (len(free),), "free masses")
    demand((masses > 0).all(), "invalid_mass", "All free masses must be positive")
    repeated_mass = np.repeat(masses, 3)
    weighted = symmetric / np.sqrt(np.outer(repeated_mass, repeated_mass))
    eigenvalues = np.linalg.eigvalsh(weighted)
    saved_eigenvalues = array(record["mass_weighted_eigenvalues_ev_per_angstrom2_amu"], (d,), "eigenvalues")
    near(saved_eigenvalues, eigenvalues, "mass-weighted eigenvalues")
    modes = array(record["cartesian_modes_per_sqrt_amu"], (d, len(free), 3), "Cartesian modes")
    vectors = modes.reshape(d, d) * np.sqrt(repeated_mass)
    near(vectors @ vectors.T, np.eye(d), "mode mass orthonormality", atol=1e-8, rtol=0)
    residual = weighted @ vectors.T - vectors.T * saved_eigenvalues
    near(residual, np.zeros_like(residual), "mode eigenvector equation", atol=1e-8 * max(1.0, float(np.abs(weighted).max())), rtol=0)
    frequencies = array(record["frequencies_cm1"], (d,), "signed frequencies")
    # Use saved eigenvalues after independent eigenspectrum/eigenvector checks:
    # near-zero eigenvalue signs may differ by platform at machine precision.
    calculated_frequencies = np.sign(saved_eigenvalues) * np.sqrt(np.abs(saved_eigenvalues)) * WAVENUMBER_FACTOR
    near(frequencies, calculated_frequencies, "signed frequency conversion", atol=1e-6, rtol=1e-8)
    labels = ["resolved_negative" if f < -frequency_tolerance else "resolved_positive" if f > frequency_tolerance
              else "unresolved_near_zero" for f in frequencies]
    demand(record["mode_classifications"] == labels, "mode_classification", "Mode classifications disagree with signed frequencies and threshold")
    for key, label in (("negative_mode_count", "resolved_negative"), ("resolved_positive_mode_count", "resolved_positive"),
                       ("unresolved_near_zero_mode_count", "unresolved_near_zero")):
        demand(integer(record[key], key) == labels.count(label), "mode_count", f"{key} disagrees with the frequency threshold")
    demand(integer(record["raw_negative_eigenvalue_count"], "raw negative count") == int((saved_eigenvalues < 0).sum()),
           "mode_count", "Raw negative count disagrees with saved eigenvalues")
    gradient = array(record["free_gradient_ev_per_angstrom"], (len(free), 3), "free gradient")
    near(gradient, -baseline[free], "baseline free gradient")
    force_max = float(np.linalg.norm(baseline[free], axis=1).max())
    near(number(record["free_force_max_ev_per_angstrom"], "free force maximum"), force_max, "baseline free force maximum")
    near(number(record["free_force_cartesian_rms_ev_per_angstrom"], "free force RMS"), float(np.sqrt(np.mean(baseline[free] ** 2))), "baseline free force RMS")
    stationary = force_max <= force_tolerance
    negative = labels.count("resolved_negative")
    classification = "not_stationary_within_force_tolerance" if not stationary else (
        "stationary_point_with_no_resolved_negative_modes" if negative == 0 else
        "first_order_saddle_candidate" if negative == 1 else "higher_order_saddle_candidate")
    demand(record["stationary_within_force_tolerance"] is stationary and record["classification"] == classification,
           "stationarity_classification", "Stationarity/classification contradict baseline forces and resolved modes")
    demand(record["transition_state_verified"] is False, "unsupported_validation", "Local force differences alone do not establish a verified transition state")
    if "geometry_angstrom" in record:
        near(array(record["geometry_angstrom"], (n, 3), "attached geometry"), reference, "attached/reference geometry", atol=1e-10, rtol=0)
    resolution = record.get("displacement_resolution")
    if resolution is not None:
        demand(isinstance(resolution, dict) and resolution["checked"] is True, "displacement_resolution", "Resolution record must indicate completed checks")
        near(number(resolution["max_relative_step_representation_error_observed"], "observed representation error"), max(representation_errors), "maximum displacement representation error", atol=1e-15)
        near(number(resolution["max_relative_step_representation_error_allowed"], "allowed representation error"), 1e-6, "current schema representation guard", atol=0, rtol=0)
    present_ids = [value for value in call_ids if value is not None]
    return {"atoms": n, "free_coordinates": d, "displacements": len(rows),
            "hessian_asymmetry_max_abs_ev_per_angstrom2": max_asymmetry,
            "hessian_asymmetry_relative_frobenius": relative_asymmetry,
            "max_eigenvector_residual": float(np.abs(residual).max()),
            "stationary_within_recorded_force_tolerance": stationary,
            "resolved_negative_mode_count": negative, "unresolved_mode_count": labels.count("unresolved_near_zero"),
            "classification": classification,
            "call_ids": {"present": len(present_ids), "missing_or_null": len(call_ids) - len(present_ids),
                         "duplicated_values": [value for value, count in Counter(present_ids).items() if count > 1],
                         "external_event_linkage": "not_checked"},
            "geometry_binding": "attached_geometry_matches" if "geometry_angstrom" in record else "external_identity_not_checked"}


def verify_stationary(record):
    """Public API for a standalone stationary record; never mutates its input."""
    report = report_base()
    try:
        demand(isinstance(record, dict), "invalid_record", "Stationary record must be an object")
        finite_json_tree(record)
        if "finite_difference_evidence" not in record:
            report["findings"].append({"code": "raw_force_evidence_unavailable", "message": "No displaced-force evidence is saved; a full reconstruction was not performed."})
            return report
        details = reconstruct(record)
        report.update(status="passed", numerical_reconstruction_verified=True, details=details)
    except (EvidenceError, KeyError, TypeError, ValueError, OverflowError, RecursionError, np.linalg.LinAlgError) as error:
        report.update(status="failed", findings=[{"code": getattr(error, "code", "malformed_record"), "message": str(error)}])
    return report


def inspect_file(path):
    """Select a raw record, workflow stationary block, or S1 characterization."""
    report = report_base()
    source = {"path": str(Path(path).resolve())}
    try:
        with Path(path).open("rb") as stream:
            raw = stream.read(MAX_INPUT_BYTES + 1)
        demand(len(raw) <= MAX_INPUT_BYTES, "audit_size_limit", "Input exceeds 64 MiB read limit")
        source["sha256"] = hashlib.sha256(raw).hexdigest()
        def reject_constant(value):
            raise EvidenceError("nonfinite_json", f"Nonstandard nonfinite JSON value: {value}")
        data = json.loads(raw, parse_constant=reject_constant)
        demand(isinstance(data, dict), "invalid_record", "Result JSON must be an object")
        finite_json_tree(data)
        fields = [key for key in ("stationary", "transition_state_characterization") if key in data]
        demand(len(fields) <= 1 and not (fields and "finite_difference_evidence" in data),
               "ambiguous_record", "Result contains multiple characterization blocks")
        if fields:
            source["selected_field"] = fields[0]
            source["parent_status"] = data.get("status")
            selected = data[fields[0]]
            if selected is None:
                report["findings"] = [{"code": "characterization_unavailable", "message": "Selected characterization is null; no reconstruction performed."}]
            else:
                report = verify_stationary(selected)
        else:
            source["selected_field"] = "<root>"
            report = verify_stationary(data)
        # Wrapper completion is separate from reconstruction of its saved block.
        # A failed workflow can still contain mathematically consistent partial
        # evidence, but must never acquire a completed-workflow certificate here.
        source["parent_completion_assessed"] = False
    except (EvidenceError, OSError, ValueError, TypeError, OverflowError, RecursionError) as error:
        report.update(status="failed", findings=[{"code": getattr(error, "code", "unreadable_record"), "message": str(error)}])
    return {**report, "source": source, "audit_implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "numpy_version": np.__version__}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path, help="Explicit saved JSON result (read-only)")
    parser.add_argument("--output", type=Path, help="Create a NEW report file; never overwrite")
    args = parser.parse_args()
    report = inspect_file(args.result)
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(encoded)
    print(encoded, end="")
    return {"passed": 0, "failed": 1, "unavailable": 2}[report["status"]]


if __name__ == "__main__":
    sys.exit(main())
