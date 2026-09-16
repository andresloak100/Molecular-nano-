"""Read-only consistency audit of the E1 archives; never imports a quantum solver.

NumPy is used only for linear algebra on archived 6-by-6 H2 Hessians. No model
is accepted as chemically accurate by passing these checks. Paired CC evidence
is deliberately left to the independently owned Q1 auditor.
"""
from __future__ import annotations

import argparse
from collections import Counter
import contextlib
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_ROOT = REPO / "data" / "validation"
# SI exact definitions; atomic mass constant is CODATA 2022, as used by the
# SciPy version that generated the archived frequencies. No core import.
EV_J = 1.602176634e-19
AMU_KG = 1.66053906892e-27
LIGHT_M_S = 299792458
EV_PER_KCAL_MOL = 4184 / (6.02214076e23 * EV_J)
# Historical benchmark export constant: retained for reproduction, alongside
# the independently calculated current-SI value. This differs by 0.255 ppm.
ARCHIVED_BENCHMARK_EV_PER_KCAL_MOL = 0.0433641153087705
WAVENUMBER = math.sqrt(EV_J / (1e-20 * AMU_KG)) / (2 * math.pi * LIGHT_M_S * 100)
SPECIES = {
    "methane": ("CH4", Counter(C=1, H=4), 0),
    "ethynyl_radical": ("CCH", Counter(C=2, H=1), 1),
    "methane_ethynyl_ts": ("CH3-H-CCH", Counter(C=3, H=5), 1),
    "acetylene": ("HCCH", Counter(C=2, H=2), 0),
    "methyl_radical": ("CH3", Counter(C=1, H=3), 1),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, label, atol=1e-9):
    require(not isinstance(actual, bool) and isinstance(actual, (int, float)),
            f"{label}: expected a number")
    require(math.isfinite(actual) and math.isclose(actual, expected, rel_tol=0, abs_tol=atol),
            f"{label}: saved {actual!r}, recomputed {expected!r}")


def finite_tree(value):
    if isinstance(value, float):
        require(math.isfinite(value), "non-finite JSON number")
    elif isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)


class Archive:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.files = {}

    def raw(self, relative):
        path = (self.root / relative).resolve()
        require(path.is_relative_to(self.root), f"path escapes archive: {relative}")
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        key = str(path.relative_to(self.root))
        require(key not in self.files or self.files[key] == digest, f"input changed during audit: {key}")
        self.files[key] = digest
        return data

    def json(self, relative):
        result = json.loads(self.raw(relative))
        finite_tree(result)
        return result

    def digest(self, relative, expected):
        actual = hashlib.sha256(self.raw(relative)).hexdigest()
        require(actual == expected, f"SHA-256 mismatch: {relative}")
        return actual

    def xyz(self, relative):
        lines = self.raw(relative).decode().splitlines()
        count = int(lines[0])
        rows = [row.split() for row in lines[2:] if row.strip()]
        require(len(rows) == count, f"single XYZ frame/count mismatch: {relative}")
        positions = np.array([[float(v) for v in row[1:4]] for row in rows])
        require(positions.shape == (count, 3) and np.isfinite(positions).all(),
                f"invalid XYZ coordinates: {relative}")
        return [row[0] for row in rows], positions


def audit_reference(archive, directory):
    base = Path(directory)
    manifest = archive.json(base / "input_manifest.json")
    benchmark = archive.json(base / "benchmark.json")
    require(benchmark["status"] == "completed", "reference benchmark is incomplete")
    require(benchmark["input_hashes"] == manifest["file_sha256"], "manifest/benchmark input hashes disagree")
    for name, digest in manifest["file_sha256"].items():
        archive.digest(base / "inputs" / name, digest)
    provenance = archive.json(base / "inputs/provenance.json")
    require(provenance == manifest["provenance"], "embedded provenance differs from saved input")
    require(provenance["license"] == "CC BY-NC 4.0", "separate source license missing")
    conversion = provenance["conversion"]
    require(conversion["source_length_unit"] == "bohr" and conversion["output_length_unit"] == "angstrom",
            "reference length units changed")
    close(conversion["angstrom_per_bohr"], 0.529177210544, "recorded bohr conversion", 1e-14)
    source = archive.raw(base / "inputs/temelso_2006_si.txt").decode()
    archive.digest(base / "inputs/temelso_2006_si.txt", provenance["source_sha256"])
    require(set(benchmark["species"]) == set(SPECIES), "reference species missing or unexpected")
    energies = {}
    for name, (source_label, elements, spin) in SPECIES.items():
        row = archive.json(base / "species" / f"{name}.json")
        require(row == benchmark["species"][name], f"{name}: species and summary disagree")
        require(row["status"] == "completed" and row["scf_converged"] is True,
                f"{name}: incomplete electronic calculation")
        require(row["quantum_diagnostics"]["gradient_completed"] is True, f"{name}: gradient incomplete")
        archive.digest(base / row["file"], row["file_sha256"])
        symbols, positions = archive.xyz(base / row["file"])
        require(Counter(symbols) == elements and row["atoms"] == len(symbols), f"{name}: atom accounting")
        # Extract source coordinates independently; no published energy is a target.
        blocks = re.findall(r"(?m)^\s*" + re.escape(source_label)
                            + r"\s+-\s*UCCSD\(T\)/cc-pVDZ[^\n]*\n((?:[CH]\s+[^\n]+\n)+)", source)
        require(len(blocks) == 1, f"{name}: expected one source geometry")
        rows = [line.split() for line in blocks[0].strip().splitlines()]
        require([row[0] for row in rows] == symbols, f"{name}: source atom order")
        converted = np.array([[float(v) for v in row[1:]] for row in rows]) * conversion["angstrom_per_bohr"]
        require(np.allclose(positions, converted, rtol=0, atol=5e-11), f"{name}: source coordinate conversion")
        settings = dict(benchmark["quantum_settings"], charge=0, spin=spin)
        require(row["quantum_settings"] == settings and row["charge"] == 0 and row["spin_2s"] == spin,
                f"{name}: inconsistent state/settings")
        diagnostics = row["quantum_diagnostics"]
        require(diagnostics["scf_converged"] is True, f"{name}: diagnostic SCF convergence contradicts completion")
        require(diagnostics["settings"] == settings, f"{name}: diagnostics settings mismatch")
        close(row["energy_ev"], row["energy_hartree"] * diagnostics["hartree_eV"], f"{name}: energy units", 1e-8)
        close(row["energy_hartree"], diagnostics["total_energy_hartree"], f"{name}: diagnostic energy")
        forces = np.array(row["forces_ev_per_angstrom"])
        require(forces.shape == (len(symbols), 3), f"{name}: force shape")
        close(row["max_force_ev_per_angstrom"], float(np.linalg.norm(forces, axis=1).max()), f"{name}: max force")
        close(row["expected_s2"], spin / 2 * (spin / 2 + 1), f"{name}: nominal spin square")
        close(row["spin_contamination"], row["s2"] - row["expected_s2"], f"{name}: spin deviation")
        energies[name] = row["energy_ev"]
    reactants = energies["methane"] + energies["ethynyl_radical"]
    relative = energies["methane_ethynyl_ts"] - reactants
    reaction = energies["acetylene"] + energies["methyl_radical"] - reactants
    for label, value in (("reference_geometry_energy", relative), ("reaction_energy", reaction)):
        close(benchmark["computed"][label + "_ev"], value, label)
        close(benchmark["computed"][label + "_kcal_per_mol"], value / ARCHIVED_BENCHMARK_EV_PER_KCAL_MOL, label + " archived units")
    for flag in ("method_validated", "saddle_verified", "tool_validated"):
        require(benchmark[flag] is False, f"unsupported {flag} flag")
    return {"species": len(energies), "source_coordinates_and_units_checked": True,
            "fixed_geometry_relative_energy_kcal_per_mol_archived": relative / ARCHIVED_BENCHMARK_EV_PER_KCAL_MOL,
            "fixed_geometry_relative_energy_kcal_per_mol_current_si": relative / EV_PER_KCAL_MOL,
            "unit_conversion_difference_kcal_per_mol": relative / EV_PER_KCAL_MOL - relative / ARCHIVED_BENCHMARK_EV_PER_KCAL_MOL,
            "unit_note": "Historical kcal/mol conversion differs by about 0.255 ppm from current exact-SI definitions; original records preserved."}


def audit_h2(archive):
    base = Path("h2-integration")
    relax = archive.json(base / "relax/result.json")
    comparison = archive.json(base / "step-comparison.json")
    require(relax["status"] == "completed" and relax["geometry_converged"] is True, "H2 relaxation incomplete")
    require(relax["structure"]["quantum_diagnostics"]["scf_converged"] is True
            and relax["structure"]["quantum_diagnostics"]["gradient_completed"] is True, "H2 relaxation electronic evidence incomplete")
    relaxed_forces = np.array(relax["structure"]["forces_ev_per_angstrom"])
    require(relaxed_forces.shape == (2, 3), "H2 relaxed force dimensions")
    relaxed_force = float(np.linalg.norm(relaxed_forces, axis=1).max())
    close(relax["structure"]["free_force_max_ev_per_angstrom"], relaxed_force, "H2 relaxed force maximum")
    require(relaxed_force <= relax["optimization"]["fmax_ev_per_angstrom"], "H2 relaxation not within claimed force tolerance")
    symbols, geometry = archive.xyz(base / "relax/structure.extxyz")
    require(symbols == ["H", "H"], "H2 atom identity")
    records = []
    max_residuals = []
    for directory in ("modes", "modes-half-step"):
        record = archive.json(base / directory / "result.json")
        require(record["status"] == "completed", f"{directory}: incomplete characterization")
        require(record["quantum_settings"] == relax["quantum_settings"], "H2 method changed between steps")
        diagnostics = record["quantum_diagnostics"]
        require(diagnostics["scf_converged"] is True and diagnostics["gradient_completed"] is True,
                "H2 base electronic evidence incomplete")
        require(diagnostics["settings"] == record["quantum_settings"], "H2 diagnostic settings mismatch")
        for key, file in (("design_sha256", "inputs/design.json"), ("initial_sha256", "inputs/h2.xyz"),
                          ("final_sha256", "inputs/h2.xyz"), ("structure_sha256", f"{directory}/input.extxyz")):
            archive.digest(base / file, record["input_hashes"][key])
            if key != "structure_sha256":
                require(record["input_hashes"][key] == relax["input_hashes"][key], "H2 design changed")
        mode_symbols, mode_geometry = archive.xyz(base / directory / "input.extxyz")
        require(mode_symbols == symbols and np.array_equal(mode_geometry, geometry), "H2 mode geometry differs from relaxation")
        stationary = record["stationary"]
        require(stationary["free_atom_indices"] == [0, 1] and stationary["frozen_atom_indices"] == [], "H2 coordinate subspace changed")
        require(stationary["free_coordinate_count"] == 6 and stationary["force_requests"] == 13
                and record["expected_force_evaluations"] == 13, "H2 displacement evaluation count")
        masses = np.array(stationary["free_masses_amu"])
        require(masses.shape == (2,) and (masses > 0).all(), "H2 masses invalid")
        require(np.array_equal(masses, [1.008, 1.008]), "H2 isotope masses changed")
        hessian = np.array(stationary["hessian_ev_per_angstrom2"])
        modes = np.array(stationary["cartesian_modes_per_sqrt_amu"])
        saved_values = np.array(stationary["mass_weighted_eigenvalues_ev_per_angstrom2_amu"])
        require(hessian.shape == (6, 6) and modes.shape == (6, 2, 3) and saved_values.shape == (6,), "H2 mode array dimensions")
        require(np.allclose(hessian, hessian.T, rtol=0, atol=1e-10), "saved symmetrized Hessian is asymmetric")
        mass = np.repeat(masses, 3)
        weighted = hessian / np.sqrt(np.outer(mass, mass))
        values = np.linalg.eigvalsh(weighted)
        require(np.allclose(saved_values, values, rtol=0, atol=1e-10), "Hessian eigenvalues disagree")
        vectors = modes.reshape(6, 6) * np.sqrt(mass)
        require(np.allclose(vectors @ vectors.T, np.eye(6), rtol=0, atol=1e-8), "mass-normalized modes are not orthonormal")
        residual = float(np.max(np.abs(weighted @ vectors.T - vectors.T * saved_values)))
        require(residual < 1e-8, "mode eigenvector residual too large")
        max_residuals.append(residual)
        # Check saved eigenvalue-to-frequency arithmetic. Near-zero eigensolver
        # signs may vary by platform, so never demand matching numerical noise.
        frequencies = np.array(stationary["frequencies_cm1"])
        calculated = np.sign(saved_values) * np.sqrt(np.abs(saved_values)) * WAVENUMBER
        require(frequencies.shape == (6,) and np.allclose(frequencies, calculated, rtol=0, atol=1e-7), "frequency unit conversion")
        threshold = stationary["settings"]["frequency_tolerance_cm1"]
        close(threshold, comparison["frequency_threshold_cm1"], "frequency threshold")
        labels = ["resolved_negative" if f < -threshold else "resolved_positive" if f > threshold
                  else "unresolved_near_zero" for f in frequencies]
        require(labels == stationary["mode_classifications"], "frequency classifications disagree")
        for key, label in (("negative_mode_count", "resolved_negative"), ("resolved_positive_mode_count", "resolved_positive"),
                           ("unresolved_near_zero_mode_count", "unresolved_near_zero")):
            require(stationary[key] == labels.count(label), f"incorrect {key}")
        require(stationary["raw_negative_eigenvalue_count"] == int((saved_values < 0).sum()), "raw negative count")
        gradient = np.array(stationary["free_gradient_ev_per_angstrom"])
        require(gradient.shape == (2, 3), "H2 gradient shape")
        force = float(np.linalg.norm(gradient, axis=1).max())
        close(stationary["free_force_max_ev_per_angstrom"], force, "H2 maximum force")
        close(stationary["free_force_cartesian_rms_ev_per_angstrom"], float(np.sqrt(np.mean(gradient**2))), "H2 RMS force")
        require(stationary["stationary_within_force_tolerance"] is (force <= stationary["settings"]["force_tolerance_ev_per_angstrom"]), "stationarity flag disagrees with force")
        require(stationary["stationary_within_force_tolerance"] is True
                and stationary["negative_mode_count"] == 0
                and stationary["classification"] == "stationary_point_with_no_resolved_negative_modes",
                "H2 stationary classification disagrees with mode/force evidence")
        require(stationary["transition_state_verified"] is False and record["validation"]["design_validated"] is False
                and record["validation"]["transition_state_validated"] is False,
                "H2 integration evidence promoted to validated design/transition state")
        # Link the recorded base and displaced computations to their completion events.
        events = [json.loads(line) for line in archive.raw(base / directory / "electronic.jsonl").splitlines() if line.strip()]
        finite_tree(events)
        done = [event for event in events if event["event"] == "calculation_completed"]
        identifiers = {event["call_id"] for event in done}
        require(len(done) == len(identifiers) == 13, "H2 gradient completion events missing or duplicated")
        require(all(event["scf_converged"] is True and event["gradient_completed"] is True for event in done),
                "H2 completion flags indicate failure")
        require({event["call_id"] for event in events} == identifiers, "H2 unfinished event calls")
        for identifier in identifiers:
            group = [event for event in events if event["call_id"] == identifier]
            scf = [event for event in group if event["event"] == "scf_completed"]
            started = [event for event in group if event["event"] == "gradient_started"]
            require(len(scf) == len(started) == 1 and scf[0]["converged"] is True,
                    "H2 per-call SCF/gradient evidence missing")
            completed = next(event for event in group if event["event"] == "calculation_completed")
            require(group.index(scf[0]) < group.index(started[0]) < group.index(completed), "H2 per-call event order")
        require(diagnostics["call_id"] in identifiers, "H2 base call absent in events")
        base_event = next(event for event in done if event["call_id"] == diagnostics["call_id"])
        close(base_event["energy_eV"], diagnostics["total_energy_hartree"] * diagnostics["hartree_eV"], "H2 base-event energy", 1e-8)
        records.append(stationary)
    for i, record in enumerate(records):
        close(comparison["steps_angstrom"][i], record["settings"]["step_angstrom"], "saved displacement step")
        close(comparison["stretch_frequencies_cm1"][i], max(record["frequencies_cm1"]), "saved stretch")
        require(comparison["negative_mode_counts"][i] == record["negative_mode_count"], "step comparison negative counts")
        require(comparison["unresolved_mode_counts"][i] == record["unresolved_near_zero_mode_count"], "step comparison unresolved counts")
    delta = max(records[1]["frequencies_cm1"]) - max(records[0]["frequencies_cm1"])
    close(comparison["stretch_frequency_change_cm1"], delta, "H2 stretch change")
    return {"stretch_change_cm1": delta, "max_mode_equation_residual": max(max_residuals),
            "unresolved_modes_per_step": [r["unresolved_near_zero_mode_count"] for r in records],
            "raw_displaced_forces_available": False,
            "limitation": "Raw displaced force arrays and the unsymmetrized Hessian are not archived; their derivation/asymmetry cannot be re-audited from these files."}


def audit_direct_df(archive):
    for directory in ("h-abstraction-direct-initial", "h-abstraction-df-initial"):
        base = Path(directory)
        metadata = archive.json(base / "interpretation.json")
        result = archive.json(base / "result.json")
        for name, digest in metadata["file_sha256"].items():
            archive.digest(base / name, digest)
        for name in ("initial", "final"):
            archive.digest(base / f"input-{name}.extxyz", result["input_hashes"][f"{name}_sha256"])
            symbols, _ = archive.xyz(base / f"input-{name}.extxyz")
            require(Counter(symbols) == Counter(C=22, H=31), "53-atom candidate accounting")
        require(result["design"]["fixed_indices"] == [7, 8, 9, 35, 36, 37], "candidate anchors changed")
        require(result["validation"]["design_validated"] is False, "single point promoted to validated design")
        reproduction = metadata.get("reproduction_input")
        if reproduction:
            archive.digest(base / reproduction["file"], reproduction["sha256"])
            design = archive.json(base / reproduction["file"])
            original = result["design"]
            expected = dict(original, length_unit="angstrom", initial="input-initial.extxyz", final="input-final.extxyz")
            require(design == expected, "reproduction design contains undocumented changes")
    archive.json("h-abstraction-df-initial/comparison.json")
    # Reuse the existing owner's independent arithmetic implementation. Disable
    # Python -O, which would remove that implementation's assertion checks.
    require(__debug__, "run without Python -O; the reused DF checker uses assertions")
    spec = importlib.util.spec_from_file_location("saved_df_checker", DEFAULT_ROOT / "verify_density_fitting_comparison.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ROOT = archive.root
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        module.main()
    result = json.loads(output.getvalue())
    result["limitation"] = "Original direct-run input JSON bytes and event/thread instrumentation are not archived; no reconstruction of its design hash was attempted."
    return result


def audit(root=DEFAULT_ROOT):
    archive = Archive(root)
    checks = []
    for name, function in (("pbe0_svp", lambda: audit_reference(archive, "pbe0-svp")),
                           ("pbe0_tzvp", lambda: audit_reference(archive, "pbe0-tzvp")),
                           ("h2_characterization", lambda: audit_h2(archive)),
                           ("direct_df", lambda: audit_direct_df(archive))):
        try:
            checks.append({"check": name, "status": "passed", "details": function()})
        except (ValueError, KeyError, TypeError, IndexError, OSError, AssertionError) as error:
            checks.append({"check": name, "status": "failed", "error": f"{type(error).__name__}: {error}"})
    return {"schema_version": 1, "audited_utc": datetime.now(timezone.utc).isoformat(),
            "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "reused_df_checker_sha256": hashlib.sha256((DEFAULT_ROOT / "verify_density_fitting_comparison.py").read_bytes()).hexdigest(),
            "numpy_version": np.__version__,
            "status": "passed" if all(c["status"] == "passed" for c in checks) else "failed",
            "scope": "Archived evidence consistency only; no quantum calculations or scientific validation",
            "scientific_model_validated": False, "checks": checks,
            "limitations": ["H2 raw displaced forces are unavailable; only saved Hessian/eigenpair consistency is checked.",
                            "Portable input manifests omit source_directory; the original manifest hash is retained, not reconstructed.",
                            "Checksums detect disagreement with saved digests, not authenticity or correctness of a quantum calculation.",
                            "Paired CC archives belong to the separate Q1 audit; this report does not cover them."],
            "audited_file_sha256": dict(sorted(archive.files.items()))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output", type=Path, help="New audit report file; existing files are never overwritten")
    args = parser.parse_args()
    report = audit(args.validation_root)
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(encoded)
    print(encoded, end="")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
