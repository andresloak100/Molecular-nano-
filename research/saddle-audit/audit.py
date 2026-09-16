"""Read S1 records without executing chemistry or importing the S1 implementation."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.units import Hartree
from mode_evidence import transfer_projection

KCAL_PER_EV = 23.060547830618307
REPO = Path(__file__).resolve().parents[2]


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and np.isfinite(value)


def audit_record(record):
    """Separate record consistency, geometric evidence and unverified chemistry.

    Missing work is a limitation, not fabricated failure or a completed result.
    No outcome of this audit certifies an electronic state or transition state.
    """
    findings, missing = [], []

    def issue(code, message):
        findings.append({"code": code, "message": message})

    pieces = ("methane", "ethynyl", "transition_state")
    statuses = {name: record.get(name, {}).get("status", "not_recorded") for name in pieces}
    complete = all(status == "converged" for status in statuses.values())
    for name, status in statuses.items():
        if status != "converged":
            missing.append(f"{name}: {status}")
    scans = {}
    for species in ("ethynyl", "transition_state"):
        scan = record.get(f"{species}_initial_guess_scan")
        if not scan:
            missing.append(f"{species}: initial-guess scan not recorded")
            scans[species] = {"status": "not_recorded"}
            continue
        energies = scan.get("energies_hartree", {})
        declared_attempts = scan.get("attempted_guesses")
        if declared_attempts is not None and set(declared_attempts) != set(energies):
            issue(f"{species}_scan_manifest_mismatch", "Attempted-guess list disagrees with recorded energy keys; all recorded keys are audited.")
        attempted = list(dict.fromkeys([*(declared_attempts or []), *energies]))
        valid = {name: energies[name] for name in attempted if finite(energies.get(name))}
        failed = [name for name in attempted if name not in valid]
        spread = (max(valid.values()) - min(valid.values())) * Hartree * KCAL_PER_EV if len(valid) >= 2 else None
        # 0.01 kcal/mol is S1's declared agreement threshold, not an accuracy bound.
        agreement = spread is not None and spread < 0.01
        full = bool(attempted) and len(valid) == len(attempted) and len(valid) >= 2
        scans[species] = {
            "attempted": attempted, "converged": list(valid), "missing_or_failed": failed,
            "spread_kcal_per_mol": spread, "agreement_over_converged_subset": agreement,
            "complete_recorded_scan": full,
            "status": "different_solutions" if spread is not None and not agreement else
                      "recorded_guesses_agree" if full else "incomplete",
        }
        if spread is not None and not agreement:
            issue(f"{species}_electronic_solution_ambiguity", f"Converged guesses differ by {spread:.9g} kcal/mol at this geometry.")
        if not full:
            missing.append(f"{species}: initial-guess scan incomplete")
        if scan.get("guess_independent") is True and not (full and agreement):
            issue(f"{species}_unsupported_guess_independence", "Reported independence is unsupported by the recorded converged energies.")
        if scan.get("complete_four_guess_check") is True and not (full and len(valid) == 4):
            issue(f"{species}_incorrect_scan_completeness", "The complete-four-guess flag disagrees with the recorded energies.")

    barrier = record.get("barrier", {})
    recomputed_delta = None
    if complete:
        values = [record[name].get("energy_ev") for name in pieces]
        if all(finite(x) for x in values):
            recomputed_delta = (values[2] - values[0] - values[1]) * KCAL_PER_EV
            reported = barrier.get("value_kcal_per_mol")
            if finite(reported) and not np.isclose(recomputed_delta, reported, atol=1e-7, rtol=1e-8):
                issue("energy_arithmetic_mismatch", "Reported energy difference does not match the recorded component energies.")
        else:
            issue("missing_component_energy", "A converged component lacks a finite electronic energy.")
        if all(isinstance(record[name].get("chemical_symbols"), list) for name in pieces):
            reactants = Counter(record["methane"]["chemical_symbols"]) + Counter(record["ethynyl"]["chemical_symbols"])
            if reactants != Counter(record["transition_state"]["chemical_symbols"]):
                issue("atom_balance_mismatch", "Transition geometry and separated reactants have different atom counts.")
        else:
            missing.append("component atom accounting not recorded")

    char = record.get("transition_state_characterization") or {}
    geometric_candidate = False
    curvature = {"status": "not_recorded"}
    transfer = {"assessed": False, "reason": "characterization not recorded"}
    if not char:
        missing.append("transition geometry: completed Hessian characterization not recorded")
    else:
        try:
            tolerance = char["settings"]["force_tolerance_ev_per_angstrom"]
            freq_tolerance = char["settings"]["frequency_tolerance_cm1"]
            free = char["free_atom_indices"]
            masses = np.asarray(char["free_masses_amu"], dtype=float)
            frequencies = np.asarray(char["frequencies_cm1"], dtype=float)
            modes = np.asarray(char["cartesian_modes_per_sqrt_amu"], dtype=float)
            gradient = np.asarray(char["free_gradient_ev_per_angstrom"], dtype=float)
            d = len(free) * 3
            if not d or len(set(free)) != len(free) or any(type(x) is not int or x < 0 for x in free):
                raise ValueError("invalid free atom indices")
            ts = record.get("transition_state", {})
            symbols = ts["chemical_symbols"]
            optimized_geometry = np.asarray(ts["final_geometry_angstrom"], dtype=float)
            characterized_geometry = np.asarray(char["geometry_angstrom"], dtype=float)
            if (not symbols or optimized_geometry.shape != (len(symbols), 3)
                    or characterized_geometry.shape != optimized_geometry.shape
                    or not np.isfinite(optimized_geometry).all()
                    or not np.isfinite(characterized_geometry).all()):
                raise ValueError("missing or invalid transition/characterization geometry")
            if not np.allclose(optimized_geometry, characterized_geometry, atol=1e-10, rtol=0):
                raise ValueError("characterization geometry differs from optimized transition geometry")
            if sorted(free) != list(range(len(symbols))) or char.get("frozen_atom_indices", []):
                raise ValueError("S1 unconstrained characterization must include every transition atom")
            if masses.shape != (len(free),) or modes.shape != (d, len(free), 3) or frequencies.shape != (d,) or gradient.shape != (len(free), 3):
                raise ValueError("inconsistent masses, modes, frequencies or gradient dimensions")
            if not all(np.isfinite(x).all() for x in (masses, frequencies, modes, gradient)) or np.any(masses <= 0):
                raise ValueError("nonfinite values or invalid masses")
            if not all(finite(x) and x > 0 for x in (tolerance, freq_tolerance)):
                raise ValueError("invalid force or frequency tolerance")
            force = float(np.linalg.norm(gradient, axis=1).max())
            stationary = force <= tolerance
            negative = int(np.count_nonzero(frequencies < -freq_tolerance))
            if char.get("stationary_within_force_tolerance") is not stationary:
                issue("stationarity_flag_mismatch", "Stationarity flag disagrees with the saved gradient and force tolerance.")
            if char.get("negative_mode_count") != negative:
                issue("negative_mode_count_mismatch", "Mode count disagrees with signed frequencies and resolution tolerance.")
            weighted = modes.reshape(d, d) * np.repeat(np.sqrt(masses), 3)
            orthonormal = bool(np.allclose(weighted @ weighted.T, np.eye(d), atol=1e-6, rtol=1e-6))
            if not orthonormal:
                issue("mode_mass_normalization_mismatch", "Saved Cartesian modes are not orthonormal in the recorded mass metric.")
            geometric_candidate = bool(complete and stationary and negative == 1 and orthonormal)
            curvature = {"status": "inspected", "free_force_max_ev_per_angstrom": force,
                         "stationary": stationary, "resolved_negative_modes": negative,
                         "mass_orthonormal_modes": orthonormal}
        except (KeyError, TypeError, ValueError) as error:
            issue("invalid_characterization", str(error))
            curvature = {"status": "invalid_record"}
        if char.get("finite_difference_step_convergence_checked") is not True:
            missing.append("finite-difference step convergence not established in characterization")
        if char.get("transition_state_verified") is True:
            issue("unsupported_characterization_verification", "S1 local Hessian characterization alone cannot verify transition-state connectivity.")
        if record.get("reaction_indices") and curvature["status"] == "inspected":
            try:
                transfer = transfer_projection(char, record["reaction_indices"], record["transition_state"]["chemical_symbols"])
                claimed = record.get("transition_mode_analysis", {}).get("single_transfer_like_mode")
                if claimed is True and transfer.get("single_negative_mode_passes_overlap_heuristic") is not True:
                    issue("transfer_heuristic_disagreement", "S1 labels a transfer-like mode, but it does not pass the independent declared mass-metric overlap threshold.")
            except (KeyError, TypeError, ValueError) as error:
                transfer = {"assessed": False, "reason": str(error)}
                issue("invalid_transfer_evidence", str(error))
        else:
            transfer = {"assessed": False, "reason": "valid geometry-bound characterization and reaction indices required"}

    # A geometric saddle candidate remains distinct from the intended chemical transition.
    if barrier.get("transition_state_verified") is True or "verified" in barrier.get("status", "").replace("unverified", ""):
        issue("unsupported_transition_state_verification", "This S1 schema contains no verified endpoint-connectivity evidence.")
    if barrier.get("status") == "barrier_on_saddle_candidate" and not geometric_candidate:
        issue("unsupported_saddle_candidate", "The claimed saddle candidate lacks a consistent stationary, single-negative-mode record.")
    if finite(barrier.get("value_kcal_per_mol")) and not geometric_candidate:
        issue("unconfirmed_energy_difference_named_barrier", "An electronic energy difference is present without established local saddle evidence; display it as an unconfirmed energy difference.")

    missing.extend(["IRC or equivalent endpoint connectivity not established by this audit",
                    "electronic-state identity not established by agreement among guesses",
                    "zero-point, thermal and tunnelling corrections are outside this audit"])
    return {
        "component_status": statuses,
        "geometric_first_order_saddle_candidate": geometric_candidate,
        "curvature": curvature,
        "independent_transfer_evidence": transfer,
        "electronic_guess_evidence": scans,
        "recomputed_electronic_energy_difference_kcal_per_mol": recomputed_delta,
        "transition_state_verified": False,
        "findings": findings,
        "missing_evidence_or_limits": missing,
        "audit_interpretation": "Read-only record consistency and evidence review; no chemistry was computed or certified.",
    }


def inspect_file(path):
    raw = path.read_bytes()
    source = {"path": str(path.resolve()), "sha256": hashlib.sha256(raw).hexdigest()}
    try:
        record = json.loads(raw, parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"nonfinite JSON: {value}")))
        if not isinstance(record, dict):
            raise ValueError("record must be a JSON object")
        result = audit_record(record)
    except (ValueError, TypeError, KeyError, AttributeError) as error:
        # An actively written file can be temporarily incomplete; do not alter it.
        result = {"findings": [{"code": "unreadable_record", "message": str(error)}],
                  "transition_state_verified": False}
    return {"source": source, **result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records", nargs="*", type=Path)
    parser.add_argument("--output", type=Path, help="Create a NEW audit JSON file; existing files are never overwritten.")
    args = parser.parse_args()
    paths = args.records or sorted((REPO / "research/reference-saddle").rglob("saddle_study.json"))
    report = {"schema_version": 1, "audited_utc": datetime.now(timezone.utc).isoformat(),
              "quantum_jobs_started": 0, "source_files_modified": False,
              "audit_implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "mode_evidence_implementation_sha256": hashlib.sha256(Path(__file__).with_name("mode_evidence.py").read_bytes()).hexdigest(),
              "records": [inspect_file(path) for path in paths]}
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(encoded)
    else:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
