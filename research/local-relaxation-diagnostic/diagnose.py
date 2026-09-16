"""Read-only local quadratic diagnostics; never relaxes a molecular structure."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MAX_COORDINATES = 600
MAX_BYTES = 64 * 1024 * 1024


class DiagnosticError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise DiagnosticError(message)


def positive(value, name):
    require(type(value) in (int, float) and math.isfinite(value) and value > 0,
            f"{name} must be a finite positive real number")
    return float(value)


def finite_array(value, name):
    def numbers(item):
        return all(numbers(x) for x in item) if isinstance(item, list) else (
            type(item) in (int, float) and math.isfinite(item))
    require(isinstance(value, list) and numbers(value), f"{name} must contain finite real numbers, without booleans/strings")
    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError, OverflowError) as error:
        raise DiagnosticError(f"{name} is not a rectangular numerical array") from error
    return result


def report_base():
    return {
        "schema_version": 1, "status": "unavailable", "findings": [], "estimate": None,
        "energy_error_bound_established": False, "actual_relaxation_performed": False,
        "electronic_state_verified": False, "transition_state_verified": False,
        "connectivity_verified": False, "scientific_accuracy_established": False,
        "limitations": [
            "This is a local quadratic stationary-point estimate, not an actual relaxation or a nonlinear energy-error bound.",
            "The Hessian and gradient must describe the same smooth electronic branch; saved input agreement and force reconstruction do not prove that.",
            "Force noise, finite-step curvature error, anharmonicity and branch changes can invalidate the predicted displacement and energy change.",
            "Caller thresholds are declared diagnostic choices, not experimentally justified accuracy or operating tolerances.",
            "Only the recorded free Cartesian space is used. No external-mode projection, eigenvalue clipping or pseudoinverse is performed.",
        ],
    }


def diagnose_quadratic(hessian_ev_per_angstrom2, gradient_ev_per_angstrom, *,
                       curvature_floor_ev_per_angstrom2, max_condition_number,
                       max_atom_shift_angstrom, expected_index):
    """Solve a full-rank Cartesian quadratic model, with explicit refusal limits.

    Gradient is (n_free_atoms, 3), Hessian is (3*n_free_atoms, 3*n_free_atoms).
    A result labelled estimated is still conditional on the local model.
    Malformed inputs raise DiagnosticError; scientific/numerical gates return
    refused or inconclusive with machine-readable findings.
    """
    floor = positive(curvature_floor_ev_per_angstrom2, "curvature floor")
    condition_limit = positive(max_condition_number, "maximum condition number")
    require(condition_limit >= 1, "maximum condition number must be at least one")
    shift_limit = positive(max_atom_shift_angstrom, "maximum atom shift")
    require(type(expected_index) is int and expected_index in (0, 1), "expected_index must be 0 (minimum) or 1 (first-order saddle)")
    gradient = finite_array(gradient_ev_per_angstrom, "gradient")
    require(gradient.ndim == 2 and gradient.shape[1] == 3 and 0 < gradient.size <= MAX_COORDINATES,
            f"gradient must have shape (n_free_atoms, 3), at most {MAX_COORDINATES} coordinates")
    hessian = finite_array(hessian_ev_per_angstrom2, "Hessian")
    dimension = gradient.size
    require(hessian.shape == (dimension, dimension), "Hessian shape must match all gradient coordinates")
    # The API takes an already symmetrized Hessian. Never conceal an asymmetric
    # matrix by silently changing it before the eigensolve.
    require(np.array_equal(hessian, hessian.T), "Hessian must be explicitly symmetric")
    report = report_base()
    report.update(status="refused", input_scope="unbound_numerical_arrays", settings={
        "curvature_floor_ev_per_angstrom2": floor, "max_condition_number": condition_limit,
        "max_atom_shift_angstrom": shift_limit, "expected_index": expected_index},
        coordinate_count=dimension, free_atom_count=gradient.shape[0],
        formulas={"gradient": "g = -F_free", "displacement": "delta = -inverse(H) g",
                  "energy_change": "E_quadratic(q+delta)-E_quadratic(q) = -0.5 g^T inverse(H) g"})
    try:
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            eigenvalues, vectors = np.linalg.eigh(hessian)
            require(np.isfinite(eigenvalues).all() and np.isfinite(vectors).all(), "Nonfinite eigensystem")
            magnitudes = np.abs(eigenvalues)
            smallest, largest = float(magnitudes.min()), float(magnitudes.max())
            condition = largest / smallest if smallest > 0 else None
            if condition is not None and not math.isfinite(condition):
                condition = None
            soft = int(np.count_nonzero(magnitudes <= floor))
            negative = int(np.count_nonzero(eigenvalues < -floor))
            report["spectrum"] = {
                "basis": "Euclidean free Cartesian coordinates; not mass-weighted modes",
                "eigenvalues_ev_per_angstrom2": eigenvalues.tolist(),
                "negative_count": negative, "raw_negative_count": int(np.count_nonzero(eigenvalues < 0)),
                "soft_count": soft, "absolute_condition_number": condition,
                "minimum_absolute_curvature_ev_per_angstrom2": smallest,
                "maximum_absolute_curvature_ev_per_angstrom2": largest,
            }
            if soft:
                report["findings"].append({"code": "soft_curvature", "message": "At least one curvature is at/below the declared absolute floor; no modes were discarded and no inverse was computed."})
            if condition is None or condition > condition_limit:
                report["findings"].append({"code": "ill_conditioned", "message": "The absolute spectral condition exceeds the declared limit or is undefined/nonfinite."})
            if report["findings"]:
                return report
            g = gradient.ravel()
            projected = vectors.T @ g
            mode_steps = -projected / eigenvalues
            delta = vectors @ mode_steps
            contributions = 0.5 * projected * mode_steps
            residual = g + hessian @ delta
            shifts = np.linalg.norm(delta.reshape(gradient.shape), axis=1)
            estimate = {
                "displacement_free_atoms_angstrom": delta.reshape(gradient.shape).tolist(),
                "signed_quadratic_energy_change_ev": float(np.sum(contributions)),
                "absolute_mode_contribution_sum_ev": float(np.sum(np.abs(contributions))),
                "positive_curvature_energy_change_ev": float(np.sum(contributions[eigenvalues > 0])),
                "negative_curvature_energy_change_ev": float(np.sum(contributions[eigenvalues < 0])),
                "max_atom_shift_angstrom": float(shifts.max()),
                "rms_atom_shift_angstrom": float(np.linalg.norm(delta) / math.sqrt(len(gradient))),
                "linearized_gradient_residual_norm_ev_per_angstrom": float(np.linalg.norm(residual)),
                "quadratic_identity_residual_ev": float(g @ delta + .5 * delta @ hessian @ delta - np.sum(contributions)),
                "per_cartesian_eigenmode": {
                    "gradient_projection_ev_per_angstrom": projected.tolist(),
                    "displacement_amplitude_angstrom": mode_steps.tolist(),
                    "signed_energy_contribution_ev": contributions.tolist(),
                    "interpretation": "Individual components depend on basis choice inside degenerate subspaces; totals and displacement do not.",
                },
            }
            # Do not emit Infinity/NaN if scalar conversion or an untrapped
            # backend operation overflowed.
            json.dumps(estimate, allow_nan=False)
            report.update(status="estimated", estimate=estimate,
                          quadratic_character="minimum" if negative == 0 else "first_order_saddle" if negative == 1 else "higher_order_saddle")
            if negative != expected_index:
                report["findings"].append({"code": "index_mismatch", "message": f"The Cartesian quadratic index is {negative}; caller requested {expected_index}."})
            if estimate["max_atom_shift_angstrom"] > shift_limit:
                report["findings"].append({"code": "displacement_limit_exceeded", "message": "The predicted maximum free-atom shift exceeds the caller limit; the unmodified estimate is retained only as a diagnostic."})
            if report["findings"]:
                report["status"] = "inconclusive"
            return report
    except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError) as error:
        report.update(status="refused", estimate=None)
        report["findings"].append({"code": "numerical_failure", "message": str(error)})
        return report


def load_dependency(name, relative_path):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, hashlib.sha256(path.read_bytes()).hexdigest()


def diagnose_saved_result(path, *, structure_path=None, max_relative_hessian_asymmetry, **settings):
    """Bind actual workflow/H1 records with N1 and reconstruct forces with E2.

    No quantum or production modules are imported by this adapter. Legacy files
    lacking full force evidence are unavailable, never assigned invented forces.
    """
    report = report_base()
    path = Path(path)
    report["source"] = {"result_path": str(path.resolve())}
    try:
        asymmetry_limit = positive(max_relative_hessian_asymmetry, "maximum relative Hessian asymmetry")
        n1, n1_hash = load_dependency("l1_n1_adapter", "research/mode-comparison/compare_modes.py")
        e2, e2_hash = load_dependency("l1_e2_verifier", "research/evidence-audit/verify_stationary.py")
        bound = n1.load_saved_result(path, structure_path)
        report.update(source=bound["source"], geometry_binding=bound["geometry_binding"],
                      source_notes=bound["notes"], quantum_settings=bound["quantum_settings"],
                      implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      dependencies_sha256={"N1_adapter": n1_hash, "E2_verifier": e2_hash})
        n1.checked_modes(bound)
        stationary = bound["stationary"]
        verification = e2.verify_stationary(stationary)
        report["saved_force_reconstruction"] = verification
        if verification["status"] != "passed":
            report["status"] = "unavailable" if verification["status"] == "unavailable" else "refused"
            report["findings"].append({"code": "unverified_force_evidence", "message": "Full same-reference force/Hessian reconstruction is required; see saved_force_reconstruction."})
            return report
        require(bound["geometry_binding"]["full_precision_force_reference_available"], "Exact saved force-reference geometry is required")
        relative_asymmetry = stationary["hessian_asymmetry_relative_frobenius"]
        report["hessian_asymmetry"] = {"relative_frobenius": relative_asymmetry,
            "maximum_allowed": asymmetry_limit,
            "maximum_absolute_ev_per_angstrom2": stationary["hessian_asymmetry_max_abs_ev_per_angstrom2"]}
        if relative_asymmetry > asymmetry_limit:
            report["status"] = "refused"
            report["findings"].append({"code": "hessian_asymmetry_limit_exceeded", "message": "The raw force-derived Hessian exceeds the caller asymmetry limit; symmetrization does not repair inconsistent derivatives."})
            return report
        numerical = diagnose_quadratic(stationary["hessian_ev_per_angstrom2"],
                                      stationary["free_gradient_ev_per_angstrom"], **settings)
        report.update(numerical)
        report.update(input_scope="hash_bound_workflow_or_H1_with_reconstructed_force_evidence",
                      free_atom_indices=stationary["free_atom_indices"],
                      frozen_atom_indices=stationary["frozen_atom_indices"],
                      reference_positions_angstrom=bound["positions_angstrom"],
                      recorded_force_tolerance_ev_per_angstrom=stationary["settings"]["force_tolerance_ev_per_angstrom"],
                      recorded_max_free_force_ev_per_angstrom=stationary["free_force_max_ev_per_angstrom"])
        report["settings"]["max_relative_hessian_asymmetry"] = asymmetry_limit
        # The returned displacement has free-atom rows only. Explicit zeroes
        # bind anchors while preserving the original complete atom order.
        if report["estimate"] is not None:
            full = np.zeros((len(bound["symbols"]), 3))
            full[stationary["free_atom_indices"]] = report["estimate"]["displacement_free_atoms_angstrom"]
            report["estimate"]["displacement_all_atoms_angstrom"] = full.tolist()
        return report
    except (OSError, KeyError, TypeError, ValueError, OverflowError, IndexError, np.linalg.LinAlgError) as error:
        report.update(status="refused", estimate=None)
        report["findings"].append({"code": "invalid_or_incompatible_evidence", "message": str(error)})
        return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("--structure", type=Path)
    parser.add_argument("--curvature-floor", type=float, required=True, help="Absolute Cartesian curvature floor, eV/angstrom^2")
    parser.add_argument("--max-condition", type=float, required=True)
    parser.add_argument("--max-atom-shift", type=float, required=True, help="Maximum per-free-atom displacement, angstrom")
    parser.add_argument("--max-relative-asymmetry", type=float, required=True)
    parser.add_argument("--expected-index", type=int, choices=(0, 1), required=True)
    parser.add_argument("--output", type=Path, help="Optional NEW report; never overwrite")
    args = parser.parse_args(argv)
    report = diagnose_saved_result(args.result, structure_path=args.structure,
        curvature_floor_ev_per_angstrom2=args.curvature_floor,
        max_condition_number=args.max_condition, max_atom_shift_angstrom=args.max_atom_shift,
        max_relative_hessian_asymmetry=args.max_relative_asymmetry, expected_index=args.expected_index)
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        try:
            with args.output.open("x", encoding="utf-8") as stream:
                stream.write(text)
        except OSError as error:
            print(str(error), file=sys.stderr)
            return 2
    else:
        print(text, end="")
    return 0 if report["status"] == "estimated" else 2


if __name__ == "__main__":
    sys.exit(main())
