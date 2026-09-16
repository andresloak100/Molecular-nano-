"""Offline comparison of saved total energies and raw analytical forces.

No calculator is instantiated. Numerical agreement is not a state or chemistry
certificate. The optional H1 adapter reads a validated native checkpoint.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

UNITS = {"length": "angstrom", "energy": "eV", "force": "eV/angstrom"}
MAX_BYTES = 64 * 1024 * 1024
MAX_STENCILS = 1800
MAX_RELATIVE_OFFSET_ERROR = 1e-6
SETTING_KEYS = {"charge", "spin", "xc", "basis", "dispersion", "grid_level",
                "conv_tol", "max_cycle", "threads", "memory_mb", "density_fit",
                "scf_initial_guess"}


class EvidenceError(ValueError):
    pass


class EvidenceUnavailable(EvidenceError):
    pass


def require(condition, message):
    if not condition:
        raise EvidenceError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def number(value, label):
    require(type(value) in (float, int), f"{label} must be a real JSON number")
    try:
        result = float(value)
    except OverflowError as error:
        raise EvidenceError(f"{label} overflows binary64") from error
    require(math.isfinite(result), f"{label} must be finite")
    return result


def text(value, label):
    require(type(value) is str and bool(value.strip()), f"{label} must be a nonempty string")
    return value


def sha(value, label):
    require(type(value) is str and len(value) == 64
            and all(c in "0123456789abcdef" for c in value), f"Invalid {label}")


def matrix(value, n, label):
    require(type(value) is list and len(value) == n, f"{label} must have {n} rows")
    for row in value:
        require(type(row) is list and len(row) == 3, f"{label} must have three columns")
        for component in row:
            number(component, label)
    return value


def _context(context):
    require(type(context) is dict, "Missing evaluation context")
    for key in ("atomic_numbers", "fixed_indices", "pbc", "units", "force_scope", "energy_scope", "quantum_settings"):
        if context.get(key) is None:
            raise EvidenceUnavailable(f"Required calculation binding {key} is unavailable")
    numbers = context["atomic_numbers"]
    require(type(numbers) is list and 0 < len(numbers) <= 600
            and all(type(z) is int and 1 <= z <= 118 for z in numbers), "Invalid atomic numbers")
    fixed = context["fixed_indices"]
    require(type(fixed) is list and all(type(i) is int and 0 <= i < len(numbers) for i in fixed)
            and len(set(fixed)) == len(fixed), "Invalid fixed indices")
    pbc = context["pbc"]
    require(type(pbc) is list and len(pbc) == 3 and all(x is False for x in pbc),
            "Only explicitly nonperiodic geometries are supported")
    require(context["units"] == UNITS, "Unsupported energy, force or length units")
    require(context["force_scope"] == "raw_unconstrained", "Constraint-masked forces are unsupported")
    require(context["energy_scope"] == "total_potential_energy", "Total potential energy is required")
    settings = context["quantum_settings"]
    require(type(settings) is dict, "Quantum settings must be an object")
    if not SETTING_KEYS <= settings.keys():
        raise EvidenceUnavailable("Explicit quantum settings are incomplete")
    for key in ("charge", "spin", "grid_level", "max_cycle", "threads", "memory_mb"):
        require(type(settings[key]) is int, f"{key} must be an integer")
    require(settings["spin"] >= 0 and 0 <= settings["grid_level"] <= 9, "Invalid spin/grid level")
    require(all(settings[key] > 0 for key in ("max_cycle", "threads", "memory_mb")), "Invalid solver controls")
    require(0 < number(settings["conv_tol"], "conv_tol") < 1, "Invalid conv_tol")
    for key in ("xc", "basis"):
        text(settings[key], key)
    require(settings["dispersion"] in (None, "d3bj", "d3zero"), "Invalid dispersion")
    require(type(settings["density_fit"]) is bool, "density_fit must be boolean")
    require(settings["scf_initial_guess"] in ("minao", "atom", "1e", "huckel"), "Invalid explicit starting guess")
    versions = context["software_versions"]
    require(type(versions) is dict and all(type(k) is str and (v is None or type(v) is str)
                                         for k, v in versions.items()), "Invalid software versions")
    require(context["resolved_numerics"] is None or type(context["resolved_numerics"]) is dict,
            "resolved_numerics must be an object or explicit null")
    return len(numbers)


def _evaluation(value, calls):
    require(type(value) is dict, "Evaluation must be an object")
    n = _context(value["context"])
    matrix(value["positions_angstrom"], n, "positions")
    if value["source_sha256"] is not None:
        sha(value["source_sha256"], "source_sha256")
    if value["source_record_id"] is not None:
        text(value["source_record_id"], "source_record_id")
    require(value["status"] in ("completed", "failed", "interrupted", "not_run"), "Invalid evaluation status")
    require(type(value["scf_converged"]) is bool and type(value["gradient_completed"]) is bool,
            "Convergence flags must be booleans")
    identifier = value["calculation_call_id"]
    if identifier is not None:
        text(identifier, "calculation_call_id")
    state = value["state_evidence_id"]
    if state is not None:
        text(state, "state_evidence_id")
    if value["energy_ev"] is not None:
        number(value["energy_ev"], "total energy")
    if "source_energy_conversion" in value:
        conversion = value["source_energy_conversion"]
        require(type(conversion) is dict and conversion.get("from") == "Hartree"
                and conversion.get("to") == "eV", "Malformed source-energy conversion")
        source, factor = conversion.get("total_energy_hartree"), conversion.get("hartree_eV")
        if source is not None:
            number(source, "original total Hartree energy")
        if factor is not None:
            require(27.0 < number(factor, "Hartree/eV conversion") < 28.0, "Implausible Hartree/eV conversion")
        if value["energy_ev"] is not None:
            require(source is not None and factor is not None
                    and float(source) * float(factor) == value["energy_ev"], "Converted energy contradicts its source binding")
    if value["forces_ev_per_angstrom"] is not None:
        matrix(value["forces_ev_per_angstrom"], n, "raw forces")
    if identifier is not None:
        # A resumed record may occur in several source snapshots. Its measured
        # payload must remain identical even when those source digests differ.
        measured = {k: v for k, v in value.items() if k not in ("source_sha256", "source_record_id")}
        binding = digest(measured)
        require(identifier not in calls or calls[identifier] == binding,
                "One calculation call ID is bound to contradictory evaluations")
        calls[identifier] = binding


def _missing(value, role, *, forces=False):
    if value is None:
        return [f"{role}: evaluation absent"]
    result = []
    if value["status"] != "completed" or not value["scf_converged"] or not value["gradient_completed"]:
        result.append(f"{role}: successful energy/gradient evidence absent ({value['status']})")
    if value["energy_ev"] is None:
        result.append(f"{role}: total energy absent")
    if forces and value["forces_ev_per_angstrom"] is None:
        result.append(f"{role}: raw analytical forces absent")
    if value["calculation_call_id"] is None:
        result.append(f"{role}: calculation call identity absent")
    if value["source_sha256"] is None or value["source_record_id"] is None:
        result.append(f"{role}: source artifact binding absent")
    return result


def _point_binding(point):
    if point is None:
        return None
    result = {key: point[key] for key in ("source_sha256", "source_record_id",
                                         "calculation_call_id", "state_evidence_id")}
    if "source_energy_conversion" in point:
        result["source_energy_conversion"] = point["source_energy_conversion"]
    return result


def _report():
    return {"schema_version": 1, "status": "unavailable", "comparisons": [],
            "unavailable_stencils": [], "step_sensitivity": [], "findings": [],
            "electronic_state_identity_verified": False, "electronic_branch_continuity_verified": False,
            "scientific_model_validated": False, "numerical_convergence_certified": False,
            "all_comparisons_within_declared_tolerance": None,
            "limitations": [
                "Saved source digests and call IDs are bindings, not physical authentication.",
                "Matching starting guesses, spin diagnostics or opaque state IDs do not prove branch continuity.",
                "Finite-step residuals and observed step sensitivity are not error bounds or chemical validation.",
                "Binary64 input-rounding sensitivity excludes solver noise, SCF error and truncation error.",
                "SCF conv_tol is not an energy-error bound; no tolerance is inferred from it."]}


def analyze_stencils(document, *, force_tolerance_ev_per_angstrom=None):
    """Pure numerical analysis; malformed evidence returns an invalid report.

    ``reference`` and every available displaced evaluation explicitly repeat the
    calculation context. Stencils may have unequal positive/negative offsets.
    A supplied tolerance compares observed residuals only; it closes no gate.
    """
    report = _report()
    try:
        canonical(document)  # Forbid nonfinite metadata and non-JSON values.
        require(type(document) is dict and type(document["schema_version"]) is int
                and document["schema_version"] == 1, "Unsupported schema")
        stencils = document["stencils"]
        require(type(stencils) is list and 0 < len(stencils) <= MAX_STENCILS, "Invalid stencil count")
        tolerance = force_tolerance_ev_per_angstrom
        if tolerance is not None:
            tolerance = number(tolerance, "force tolerance")
            require(tolerance >= 0, "Force tolerance must be nonnegative")
        report["force_tolerance_ev_per_angstrom"] = tolerance
        report["input_document_sha256"] = digest(document)
        calls, seen, groups = {}, set(), defaultdict(list)
        common = None
        for index, stencil in enumerate(stencils):
            reference, minus, plus = (stencil[key] for key in ("reference", "minus", "plus"))
            if reference is None:
                raise EvidenceUnavailable("Reference geometry/context is unavailable")
            for point in (reference, minus, plus):
                if point is not None:
                    _evaluation(point, calls)
            context, positions = reference["context"], reference["positions_angstrom"]
            binding = canonical({"context": context, "positions": positions})
            if common is None:
                common = binding
                if not context["resolved_numerics"]:
                    report["findings"].append("Resolved quadrature/auxiliary-basis policy was not recorded.")
                if not context["software_versions"] or any(v is None for v in context["software_versions"].values()):
                    report["findings"].append("Some software-version identity is unavailable.")
            require(common == binding, "Reference geometry, atom order, boundary or calculation context differs across stencils")
            atom, axis = stencil["atom_index"], stencil["axis"]
            require(type(atom) is int and 0 <= atom < len(positions)
                    and type(axis) is int and 0 <= axis <= 2, "Invalid displaced coordinate")
            require(atom not in context["fixed_indices"], "Displacements of fixed atoms are unsupported")
            requested_minus = number(stencil["requested_minus_offset_angstrom"], "requested minus offset")
            requested_plus = number(stencil["requested_plus_offset_angstrom"], "requested plus offset")
            require(requested_minus < 0 < requested_plus, "Displacement signs must bracket the reference")
            offsets = []
            for role, point, requested in (("minus", minus, requested_minus), ("plus", plus, requested_plus)):
                # The intended positions also define a pending slot; never use
                # the request alone as evidence that the calculation occurred.
                expected = deepcopy(positions)
                expected[atom][axis] = number(positions[atom][axis], "reference coordinate") + requested
                actual = expected[atom][axis] - positions[atom][axis]
                error = abs(actual - requested) / abs(requested)
                require(math.isfinite(actual) and actual * requested > 0 and error <= MAX_RELATIVE_OFFSET_ERROR,
                        "Requested displacement is unresolved or distorted in represented coordinates")
                if point is not None:
                    require(canonical(point["context"]) == canonical(context), "Evaluation calculation contexts differ")
                    require(point["positions_angstrom"] == expected,
                            "Displaced geometry differs from its exact single-coordinate request")
                offsets.append(actual)
            a, b = -offsets[0], offsets[1]
            key = (atom, axis, a, b)
            require(key not in seen, "Duplicate represented stencil")
            seen.add(key)
            identities = {role: _point_binding(point) for role, point in
                          (("reference", reference), ("minus", minus), ("plus", plus))}
            missing = _missing(reference, "reference", forces=True) + _missing(minus, "minus") + _missing(plus, "plus")
            if missing:
                report["unavailable_stencils"].append({"stencil_index": index, "atom_index": atom,
                    "axis": axis, "reasons": missing, "evidence": identities})
                continue
            e0, em, ep = (float(point["energy_ev"]) for point in (reference, minus, plus))
            total = a + b
            require(math.isfinite(total) and total > 0, "Represented stencil span overflows binary64")
            # Derivative of the quadratic interpolant at the BASELINE. A secant
            # divided by a+b would instead have a curvature bias if a != b.
            wm, wp = b / total, a / total
            derivative = math.fsum((wm * ((e0 - em) / a), wp * ((ep - e0) / b)))
            fd_force = -derivative
            analytic = float(reference["forces_ev_per_angstrom"][atom][axis])
            residual = fd_force - analytic
            cm, cp = -wm / a, wp / b
            c0 = -(cm + cp)
            rounding = math.fsum((abs(cm) * math.ulp(em), abs(c0) * math.ulp(e0), abs(cp) * math.ulp(ep))) / 2
            require(all(math.isfinite(x) for x in (fd_force, residual, rounding, cm, c0, cp)),
                    "Finite-difference arithmetic overflowed; use resolvable energies/coordinates")
            state_ids = [point["state_evidence_id"] for point in (reference, minus, plus)]
            branch_notes = ["No electronic branch continuity was verified."]
            if any(value is None for value in state_ids):
                branch_notes.append("One or more per-geometry electronic-state evidence IDs are absent.")
            provided = [value for value in state_ids if value is not None]
            if len(set(provided)) < len(provided):
                branch_notes.append("An opaque state evidence ID is reused across geometries; this is not a per-geometry state binding.")
            row = {"stencil_index": index, "atom_index": atom, "axis": axis,
                "requested_minus_offset_angstrom": requested_minus, "requested_plus_offset_angstrom": requested_plus,
                "actual_minus_distance_angstrom": a, "actual_plus_distance_angstrom": b,
                "represented_asymmetry_angstrom": b - a,
                "estimator": "three_point_nonuniform_at_reference",
                "derivative_energy_weights_per_angstrom": {"minus": cm, "reference": c0, "plus": cp},
                "analytic_force_ev_per_angstrom": analytic, "finite_difference_force_ev_per_angstrom": fd_force,
                "residual_ev_per_angstrom": residual, "absolute_residual_ev_per_angstrom": abs(residual),
                "within_declared_tolerance": None if tolerance is None else abs(residual) <= tolerance,
                "binary64_input_rounding_sensitivity_ev_per_angstrom": rounding,
                "input_rounding_exceeds_declared_tolerance": None if tolerance is None else rounding > tolerance,
                "numerical_resolution_findings": [
                    "Binary64 input-rounding sensitivity exceeds the declared tolerance; a small observed residual is not an informative precision check."
                ] if tolerance is not None and rounding > tolerance else [],
                "energies_ev": {"reference": e0, "minus": em, "plus": ep},
                "evidence": identities, "branch_findings": branch_notes}
            report["comparisons"].append(row)
            groups[(atom, axis)].append(row)
        for (atom, axis), rows in groups.items():
            fd = [r["finite_difference_force_ev_per_angstrom"] for r in rows]
            analytical = [r["analytic_force_ev_per_angstrom"] for r in rows]
            report["step_sensitivity"].append({"atom_index": atom, "axis": axis,
                "distinct_stencils": len(rows), "stencil_indices": [r["stencil_index"] for r in rows],
                "finite_difference_force_spread_ev_per_angstrom": None if len(rows) < 2 else max(fd) - min(fd),
                "analytic_force_spread_ev_per_angstrom": None if len(rows) < 2 else max(analytical) - min(analytical),
                "scope": "Observed sensitivity among available stencils; no convergence or error bound."})
        count, absent = len(report["comparisons"]), len(report["unavailable_stencils"])
        report["status"] = "partial" if count and absent else "comparison_produced" if count else "unavailable"
        if count and not absent and tolerance is not None:
            report["all_comparisons_within_declared_tolerance"] = all(r["within_declared_tolerance"] for r in report["comparisons"])
        canonical(report)
    except EvidenceUnavailable as error:
        report = _report()
        report["findings"] = [str(error)]
    except (EvidenceError, KeyError, TypeError, ValueError, OverflowError, ZeroDivisionError) as error:
        # Do not leave a partial success visible after a contradictory binding.
        report = _report()
        report.update(status="invalid", findings=[str(error)])
    return report


def _pairs_without_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def read_stencil_document(path):
    """Read bounded exact JSON bytes; strict parsing, no source mutations."""
    with Path(path).open("rb") as handle:
        raw = handle.read(MAX_BYTES + 1)
    require(len(raw) <= MAX_BYTES, "Input exceeds the 64 MiB limit")
    document = json.loads(raw, object_pairs_hook=_pairs_without_duplicates)
    canonical(document)
    return document


def _h1_module():
    path = Path(__file__).resolve().parents[1] / "characterization-resume" / "prototype.py"
    spec = importlib.util.spec_from_file_location("f1_h1_checkpoint_reader", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def from_h1_checkpoints(paths):
    """Convert validated H1 checkpoints without synthesizing missing energies.

    Native loader validates captured bytes under a shared lock. Total Hartree
    energy and its recorded conversion are required; SCF-only energy is ignored.
    This imports the producer's read-only inspector, never its acquisition API.
    """
    require(paths and len(paths) <= 32, "Supply one to 32 H1 checkpoint directories")
    h1 = _h1_module()
    stencils, sources, known_conversions = [], [], set()
    for path in paths:
        checkpoint = Path(path) / "force_checkpoint.json"
        require(checkpoint.stat().st_size <= MAX_BYTES, "Checkpoint exceeds the 64 MiB limit")
        payload, source_hash = h1.load_checkpoint(path)
        manifest = payload["manifest"]
        geometry = manifest["geometry"]
        records = {r["request"]["id"]: r for r in payload["records"]}
        base_context = {"atomic_numbers": geometry["atomic_numbers"], "fixed_indices": geometry["fixed_indices"],
            "pbc": geometry["pbc"], "quantum_settings": manifest["quantum_settings"],
            "software_versions": {"python": manifest["software"]["python"], **manifest["software"]["packages"]},
            "source_software_manifest": manifest["software"], "resolved_numerics": None,
            "energy_scope": "total_potential_energy", "force_scope": "raw_unconstrained", "units": UNITS,
            "checkpoint_geometry_context": {k: v for k, v in geometry.items() if k != "positions_angstrom"},
            "implementation_sha256": manifest["implementation_sha256"], "prototype_sha256": manifest["prototype_sha256"]}
        def convert(request):
            record = records.get(request["id"])
            if record is None:
                return None
            diagnostics = record["quantum_diagnostics"]
            positions = deepcopy(geometry["positions_angstrom"])
            if request["kind"] == "displacement":
                positions[request["atom_index"]][request["axis"]] = request["displaced_coordinate_angstrom"]
            energy = None
            if diagnostics.get("total_energy_hartree") is not None and diagnostics.get("hartree_eV") is not None:
                hartree = number(diagnostics["total_energy_hartree"], "total Hartree energy")
                conversion = number(diagnostics["hartree_eV"], "Hartree/eV conversion")
                require(27.0 < conversion < 28.0, "Implausible Hartree/eV conversion")
                known_conversions.add(conversion)
                require(len(known_conversions) == 1, "Recorded Hartree/eV conversions differ across evaluations")
                energy = number(hartree * conversion, "converted total energy")
            result = {"context": deepcopy(base_context), "positions_angstrom": positions,
                "source_sha256": source_hash, "source_record_id": request["id"],
                "calculation_call_id": record["calculation_call_id"], "status": record["status"],
                "scf_converged": diagnostics["scf_converged"], "gradient_completed": diagnostics["gradient_completed"],
                "energy_ev": energy, "forces_ev_per_angstrom": record["forces_ev_per_angstrom"],
                "state_evidence_id": diagnostics.get("state_evidence_id")}
            if diagnostics.get("versions") is not None:
                require(type(diagnostics["versions"]) is dict, "Malformed per-calculation software versions")
                result["context"]["software_versions"].update(diagnostics["versions"])
            result["context"]["recorded_response_controls"] = {
                key: diagnostics.get(key) for key in
                ("grid_response", "auxiliary_basis_response", "dispersion_three_body")}
            # Preserve absent conversion/energy explicitly; do not manufacture a
            # context contradiction merely because one evaluation lacks energy.
            result["source_energy_conversion"] = {"from": "Hartree", "to": "eV",
                "total_energy_hartree": diagnostics.get("total_energy_hartree"),
                "hartree_eV": diagnostics.get("hartree_eV")}
            return result
        requests = manifest["requests"]
        baseline = convert(requests[0])
        if baseline is None:
            baseline = {"context": deepcopy(base_context), "positions_angstrom": geometry["positions_angstrom"],
                "source_sha256": source_hash, "source_record_id": "baseline", "calculation_call_id": None,
                "status": "not_run", "scf_converged": False, "gradient_completed": False,
                "energy_ev": None, "forces_ev_per_angstrom": None, "state_evidence_id": None}
        for index in range(1, len(requests), 2):
            plus_request, minus_request = requests[index:index+2]
            stencils.append({"atom_index": plus_request["atom_index"], "axis": plus_request["axis"],
                "requested_minus_offset_angstrom": minus_request["requested_offset_angstrom"],
                "requested_plus_offset_angstrom": plus_request["requested_offset_angstrom"],
                "reference": deepcopy(baseline), "minus": convert(minus_request), "plus": convert(plus_request)})
        sources.append({"path": str(Path(path).resolve()), "checkpoint_sha256": source_hash,
            "manifest_sha256": payload["manifest_sha256"], "attempt_id": payload["attempt_id"],
            "status": payload["status"], "error": payload.get("error"),
            "failed_calculation": payload.get("failed_calculation"),
            "input_context": manifest["input_context"]})
    return {"schema_version": 1, "stencils": stencils, "source_checkpoints": sources}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input", type=Path, help="Normalized stencil JSON")
    inputs.add_argument("--checkpoints", nargs="+", type=Path, help="H1 checkpoint directories")
    parser.add_argument("--force-tolerance", type=float, default=None, help="Explicit numerical residual threshold in eV/angstrom")
    parser.add_argument("--output", type=Path, help="New report file; existing files are refused")
    args = parser.parse_args(argv)
    try:
        document = read_stencil_document(args.input) if args.input else from_h1_checkpoints(args.checkpoints)
        report = analyze_stencils(document, force_tolerance_ev_per_angstrom=args.force_tolerance)
    except (OSError, ValueError, KeyError, TypeError) as error:
        report = _report()
        report.update(status="invalid", findings=[str(error)])
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as output:
            output.write(encoded)
    else:
        print(encoded, end="")
    return 2 if report["status"] == "invalid" else 0 if report["status"] == "comparison_produced" else 1


if __name__ == "__main__":
    raise SystemExit(main())
