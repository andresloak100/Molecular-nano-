#!/usr/bin/env python3
"""Read-only consistency audit for the two archived paired reference runs.

Standard library only: does not import production calculation code, execute
quantum calculations, change evidence, or establish physical validation.
Run with no arguments from any directory, or pass one or more archive roots.
JSON goes to stdout; any consistency failure gives exit status 1. Historical
schema-1 omissions are warnings, never silently filled with modern defaults.
"""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys


# Pinned archive conventions, deliberately independent of nanodesign/ASE.
# These check reproducibility of archived conversions, not metrology accuracy.
HARTREE_TO_EV = 27.211386024367243
EV_PER_KCAL_PER_MOL = 0.0433641153087705
ENERGY_ABS_TOL = 1e-8
SPECIES = {
    "methane": ({"C": 1, "H": 4}, 0),
    "ethynyl_radical": ({"C": 2, "H": 1}, 1),
    "methane_ethynyl_ts": ({"C": 3, "H": 5}, 1),
    "acetylene": ({"C": 2, "H": 2}, 0),
    "methyl_radical": ({"C": 1, "H": 3}, 1),
}
ATOMIC_NUMBERS = {"C": 6, "H": 1}
UNITS = {"energy": "eV", "secondary_energy": "kcal/mol", "length": "angstrom"}
COMMON_FALSE_FLAGS = (
    "method_accuracy_validated", "saddle_verified", "tool_validated",
    "published_barrier_used_as_acceptance_target",
)
STATE_FALSE_FLAGS = (
    "electronic_state_identity_verified", "ground_state_verified",
    "initial_guess_scan_performed",
)


class EvidenceError(ValueError):
    pass


class Audit:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.failures = []
        self.warnings = []
        self.checks = 0

    def check(self, condition, location, message):
        self.checks += 1
        if not condition:
            self.failures.append({"location": location, "message": message})

    def equal(self, actual, expected, location):
        self.check(actual == expected and type(actual) is type(expected), location,
                   f"Expected {expected!r}; found {actual!r}")

    def number(self, value, location):
        if type(value) not in (int, float) or not math.isfinite(value):
            raise EvidenceError(f"{location}: missing or non-finite number (booleans excluded)")
        return value

    def near(self, actual, expected, location, tolerance=ENERGY_ABS_TOL):
        actual = self.number(actual, location)
        expected = self.number(expected, location + " expected")
        self.check(math.isclose(actual, expected, abs_tol=tolerance, rel_tol=0),
                   location, f"Expected {expected:.16g}; found {actual:.16g}")

    def finite_tree(self, value, location):
        if isinstance(value, dict):
            for key, child in value.items():
                self.finite_tree(child, f"{location}.{key}")
        elif isinstance(value, list):
            for index, child in enumerate(value):
                self.finite_tree(child, f"{location}[{index}]")
        elif isinstance(value, float):
            self.check(math.isfinite(value), location, "Non-finite archived number")

    def path(self, relative):
        if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
            raise EvidenceError(f"Invalid relative evidence path: {relative!r}")
        resolved = (self.root / relative).resolve()
        if not resolved.is_relative_to(self.root):
            raise EvidenceError(f"Evidence path escapes archive: {relative!r}")
        return resolved

    def read_json(self, relative):
        value = json.loads(self.path(relative).read_text())
        if not isinstance(value, dict):
            raise EvidenceError(f"{relative}: expected JSON object")
        self.finite_tree(value, relative)
        return value

    def false_flags(self, record, flags, location, historical=False):
        for flag in flags:
            if historical and flag not in record:
                self.warnings.append({"location": f"{location}.{flag}",
                                      "message": "Not recorded in historical schema; remains unknown"})
            else:
                self.equal(record.get(flag), False, f"{location}.{flag}")

    def read_xyz(self, relative):
        rows = self.path(relative).read_text().splitlines()
        count = int(rows[0])
        atom_rows = [line.split() for line in rows[2:] if line.strip()]
        if count <= 0 or len(atom_rows) != count or any(len(row) != 4 for row in atom_rows):
            raise EvidenceError(f"{relative}: XYZ count/row shape mismatch")
        symbols, positions = [], []
        for index, row in enumerate(atom_rows):
            if row[0] not in ATOMIC_NUMBERS:
                raise EvidenceError(f"{relative}: unsupported element {row[0]!r}")
            symbols.append(row[0])
            positions.append([self.number(float(x), f"{relative} atom {index}") for x in row[1:]])
        return symbols, positions

    def hash(self, relative):
        return hashlib.sha256(self.path(relative).read_bytes()).hexdigest()


def relative_energies(energies):
    reactants = energies["methane"] + energies["ethynyl_radical"]
    ts = energies["methane_ethynyl_ts"] - reactants
    reaction = energies["acetylene"] + energies["methyl_radical"] - reactants
    return {
        "nominal_ts_relative_energy_ev": ts,
        "nominal_ts_relative_energy_kcal_per_mol": ts / EV_PER_KCAL_PER_MOL,
        "reaction_energy_ev": reaction,
        "reaction_energy_kcal_per_mol": reaction / EV_PER_KCAL_PER_MOL,
    }


def audit_species(audit, name, row, comparison, manifest, benchmark, schema):
    location = f"species.{name}"
    composition, spin = SPECIES[name]
    geometry_file = row["geometry_file"]
    symbols, positions = audit.read_xyz(geometry_file)
    audit.equal(dict(Counter(symbols)), composition, location + ".composition")
    xyz_hash = audit.hash(geometry_file)
    filename = Path(geometry_file).name
    audit.equal(manifest["file_sha256"].get(filename), xyz_hash, location + ".manifest_hash")
    for key in ("dft_xyz_sha256", "cc_xyz_sha256"):
        audit.equal(row.get(key), xyz_hash, location + "." + key)
    for key, expected in (("charge", 0), ("spin_2s", spin), ("basis", comparison["basis"])):
        audit.equal(row.get(key), expected, location + "." + key)
    electrons = sum(ATOMIC_NUMBERS[symbol] for symbol in symbols)
    audit.check(electrons >= spin and (electrons - spin) % 2 == 0,
                location + ".spin_accounting", "Electron count and nominal spin are incompatible")
    dft = audit.read_json(row["dft_record"])
    cc = audit.read_json(row["cc_record"])
    audit.equal(benchmark["species"].get(name), dft, location + ".benchmark_species_copy")
    audit.equal(dft.get("file"), "inputs/" + filename, location + ".dft.file")
    audit.equal(dft.get("file_sha256"), xyz_hash, location + ".dft.file_sha256")
    audit.equal(dft.get("atoms"), len(symbols), location + ".dft.atoms")
    formula = "".join(symbol + (str(composition[symbol]) if composition[symbol] > 1 else "")
                      for symbol in ("C", "H"))
    audit.equal(dft.get("formula"), formula, location + ".dft.formula")
    for key in ("geometry_file", "dft_xyz_sha256", "cc_xyz_sha256"):
        audit.equal(cc.get(key), row.get(key), location + ".cc." + key)
    cc_geometry = cc["quantum_diagnostics"]["geometry"]
    audit.equal(cc_geometry.get("symbols"), symbols, location + ".cc.geometry.symbols")
    audit.equal(cc_geometry.get("positions_angstrom"), positions, location + ".cc.geometry.positions")
    structures = [s for s in manifest["provenance"]["structures"] if s.get("file") == filename]
    audit.equal(len(structures), 1, location + ".provenance_entry_count")
    if len(structures) == 1:
        for key, expected in (("atoms", len(symbols)), ("charge", 0), ("multiplicity", spin + 1)):
            audit.equal(structures[0].get(key), expected, location + ".provenance." + key)
    for method, record, prefix in (("dft", dft, ""), ("cc", cc, "hf_")):
        loc = location + "." + method
        diag = record["quantum_diagnostics"]
        settings = diag["settings"]
        for key, expected in (("species", name), ("status", "completed"), ("charge", 0), ("spin_2s", spin)):
            audit.equal(record.get(key), expected, loc + "." + key)
        for key, expected in (("charge", 0), ("spin", spin), ("basis", comparison["basis"])):
            audit.equal(settings.get(key), expected, loc + ".settings." + key)
        template = dict(comparison["dft_settings" if method == "dft" else "cc_settings"])
        template.update(charge=0, spin=spin)
        audit.equal(settings, template, loc + ".settings_vs_comparison")
        if method == "dft":
            audit.equal(record.get("quantum_settings"), settings, loc + ".quantum_settings")
        for key, expected in (("electron_count", electrons), ("alpha_electrons", (electrons + spin)//2),
                              ("beta_electrons", (electrons - spin)//2)):
            audit.equal(diag.get(key), expected, loc + "." + key)
        audit.equal(diag.get("scf_converged"), True, loc + ".scf_converged")
        energy = audit.number(record.get("energy_ev"), loc + ".energy_ev")
        hartree = audit.number(record.get("energy_hartree"), loc + ".energy_hartree")
        audit.near(row.get(method + "_energy_ev"), energy, loc + ".summary_energy")
        audit.near(energy, hartree * HARTREE_TO_EV, loc + ".hartree_conversion")
        audit.near(diag.get("total_energy_hartree"), hartree, loc + ".diagnostic_energy")
        conversion_key = "hartree_eV" if method == "dft" else "hartree_to_ev"
        audit.near(diag.get(conversion_key), HARTREE_TO_EV, loc + ".conversion_constant")
        expected_s2 = (spin / 2) * (spin / 2 + 1)
        audit.near(diag.get(prefix + "expected_s2"), expected_s2, loc + ".expected_s2")
        s2 = audit.number(diag.get(prefix + "s2"), loc + ".s2")
        audit.near(diag.get(prefix + "s2_deviation"), s2 - expected_s2, loc + ".s2_deviation")
        for suffix, expected in (("s2", s2), ("expected_s2", expected_s2), ("s2_deviation", s2 - expected_s2)):
            audit.near(row.get(method + "_" + prefix + suffix), expected, loc + ".summary_" + suffix)
        if method == "cc":
            audit.equal(record.get("basis"), comparison["basis"], loc + ".basis")
            audit.equal(diag.get("reference"), "UHF" if spin else "RHF", loc + ".reference")
            audit.equal(diag.get("variant"), "UCCSD(T)" if spin else "RCCSD(T)", loc + ".variant")
            audit.equal(row.get("cc_reference"), diag.get("reference"), loc + ".summary_reference")
            audit.equal(row.get("cc_variant"), diag.get("variant"), loc + ".summary_variant")
            for flag in ("ccsd_converged", "triples_completed"):
                audit.equal(diag.get(flag), True, loc + "." + flag)
            audit.false_flags(diag, ("geometry_optimized", "forces_computed", "chemical_accuracy_validated",
                                     "cc_wavefunction_s2_evaluated"), loc)
            audit.false_flags(diag, STATE_FALSE_FLAGS, loc, historical=schema == 1)
            audit.near(diag.get("total_energy_ev"), energy, loc + ".diagnostic_energy_ev")
            hf = audit.number(diag.get("hf_energy_hartree"), loc + ".hf_energy")
            corr = audit.number(diag.get("ccsd_correlation_energy_hartree"), loc + ".ccsd_correlation")
            triples = audit.number(diag.get("triples_correction_hartree"), loc + ".triples")
            audit.near(diag.get("ccsd_total_energy_hartree"), hf + corr, loc + ".ccsd_sum")
            audit.near(diag.get("total_correlation_energy_hartree"), corr + triples, loc + ".correlation_sum")
            audit.near(hartree, hf + corr + triples, loc + ".total_sum")
            if schema == 2:
                guess = comparison["cc_settings"]["scf_initial_guess"]
                audit.check(isinstance(guess, str) and bool(guess), loc + ".guess", "Missing guess")
                for source, key in ((row, "cc_scf_initial_guess"), (record, "scf_initial_guess"), (diag, "scf_initial_guess")):
                    audit.equal(source.get(key), guess, loc + "." + key)
                audit.false_flags(record, ("electronic_state_identity_verified",), loc + ".record")
        else:
            audit.equal(record.get("scf_converged"), True, loc + ".record.scf_converged")
            audit.near(hartree, audit.number(diag.get("dft_energy_hartree"), loc + ".dft_energy") +
                       audit.number(diag.get("dispersion_energy_hartree"), loc + ".dispersion"), loc + ".total_sum")
            forces = record["forces_ev_per_angstrom"]
            if len(forces) != len(symbols) or any(not isinstance(v, list) or len(v) != 3 for v in forces):
                raise EvidenceError(loc + ": invalid forces shape")
            maximum = max(math.sqrt(sum(audit.number(x, loc + ".force") ** 2 for x in v)) for v in forces)
            audit.near(record.get("max_force_ev_per_angstrom"), maximum, loc + ".max_force")
            events = [json.loads(line) for line in audit.path("dft/" + record["event_log"]).read_text().splitlines() if line.strip()]
            audit.finite_tree(events, loc + ".event_log")
            completions = [e for e in events if e.get("event") == "calculation_completed" and e.get("call_id") == diag["call_id"]]
            audit.equal(len(completions), 1, loc + ".completion_event_count")
            if len(completions) == 1:
                audit.near(completions[0].get("energy_eV"), energy, loc + ".completion_energy")
                audit.equal(completions[0].get("scf_converged"), True, loc + ".completion_convergence")
    audit.near(row.get("dft_minus_cc_energy_ev"), dft["energy_ev"] - cc["energy_ev"], location + ".difference")
    if schema == 2:
        audit.equal(row.get("nominal_charge_spin_inputs_matched"), True, location + ".nominal_inputs_matched")
        audit.false_flags(row, ("electronic_state_identity_verified",), location)
    return {"atoms": len(symbols), "composition": dict(Counter(symbols)), "charge": 0,
            "spin_2s": spin, "electrons": electrons, "xyz_sha256": xyz_hash,
            "dft_energy_ev": dft["energy_ev"], "cc_energy_ev": cc["energy_ev"],
            "cc_initial_guess_recorded": cc["quantum_diagnostics"]["settings"].get("scf_initial_guess")}


def audit_archive(root):
    audit = Audit(root)
    report = {"archive": str(audit.root), "species": {}, "recomputed": {}}
    try:
        comparison = audit.read_json("method_comparison.json")
        schema = comparison.get("schema_version")
        if type(schema) is not int or schema not in (1, 2):
            raise EvidenceError(f"Unsupported comparison schema: {schema!r}")
        report["schema_version"] = schema
        audit.equal(comparison.get("status"), "completed", "comparison.status")
        audit.equal(comparison.get("units"), UNITS, "comparison.units")
        audit.equal(comparison.get("geometry_basis_pairing_verified"), True, "comparison.geometry_basis_pairing_verified")
        audit.false_flags(comparison, COMMON_FALSE_FLAGS, "comparison")
        audit.false_flags(comparison, STATE_FALSE_FLAGS, "comparison", historical=schema == 1)
        if schema == 1 and "scf_initial_guess" not in comparison["cc_settings"]:
            audit.warnings.append({"location": "comparison.cc_settings.scf_initial_guess",
                                   "message": "Not recorded; archive nickname/default behavior cannot establish the actual guess"})
        benchmark = audit.read_json(comparison["dft_evidence"])
        manifest = audit.read_json(comparison["input_manifest"])
        audit.equal(benchmark.get("status"), "completed", "benchmark.status")
        audit.equal(benchmark.get("units"), UNITS, "benchmark.units")
        audit.false_flags(benchmark, ("method_validated", "saddle_verified", "tool_validated"), "benchmark")
        audit.equal(manifest.get("archived_inputs"), "inputs", "manifest.archived_inputs")
        audit.equal(benchmark.get("input_hashes"), manifest["file_sha256"], "benchmark.input_hashes")
        audit.equal(audit.read_json("dft/inputs/provenance.json"), manifest.get("provenance"), "manifest.provenance_copy")
        for filename, expected_hash in manifest["file_sha256"].items():
            audit.equal(audit.hash("dft/inputs/" + filename), expected_hash, "manifest.hash." + filename)
        audit.equal(set(comparison["species"]), set(SPECIES), "comparison.species_names")
        audit.equal(set(benchmark["species"]), set(SPECIES), "benchmark.species_names")
        for name in SPECIES:
            try:
                report["species"][name] = audit_species(audit, name, comparison["species"][name], comparison, manifest, benchmark, schema)
            except (OSError, ValueError, KeyError, TypeError, IndexError, OverflowError, AttributeError) as exc:
                audit.failures.append({"location": "species." + name, "message": str(exc)})
        audit.equal(manifest.get("reaction_conservation"),
                    {"elements": {"C": 3, "H": 5}, "charge": 0, "electrons": 23, "spin_2s": 1},
                    "manifest.reaction_conservation")
        if len(report["species"]) == len(SPECIES):
            for method, field in (("dft", "dft_energy_ev"), ("ccsd_t", "cc_energy_ev")):
                values = relative_energies({name: item[field] for name, item in report["species"].items()})
                report["recomputed"][method] = values
            report["recomputed"]["dft_minus_ccsd_t"] = {
                key: report["recomputed"]["dft"][key] - report["recomputed"]["ccsd_t"][key]
                for key in report["recomputed"]["dft"]}
            for method, values in report["recomputed"].items():
                for key, expected in values.items():
                    audit.near(comparison["computed"][method].get(key), expected, f"computed.{method}.{key}")
            for key, expected in report["recomputed"]["dft"].items():
                benchmark_key = key.replace("nominal_ts_relative", "reference_geometry")
                audit.near(benchmark["computed"].get(benchmark_key), expected, "benchmark.computed." + benchmark_key)
    except (OSError, ValueError, KeyError, TypeError, IndexError, OverflowError, AttributeError) as exc:
        audit.failures.append({"location": "archive", "message": str(exc)})
    report.update(consistency_passed=not audit.failures, physical_validation_established=False,
                  checks=audit.checks, failures=audit.failures, warnings=audit.warnings)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path, nargs="*")
    args = parser.parse_args(argv)
    validation = Path(__file__).resolve().parents[2] / "data" / "validation"
    archives = args.archives or [validation / "paired-ccpvdz", validation / "paired-ccpvdz-atom"]
    reports = [audit_archive(root) for root in archives]
    passed = all(report["consistency_passed"] for report in reports)
    output = {
        "audit": "independent_paired_archive_consistency", "consistency_passed": passed,
        "physical_validation_established": False,
        "scope": "Archived bytes, record linkage, nominal species accounting, and arithmetic only; not a quantum recomputation or authentication of original calculations.",
        "limitations": ["Hash agreement checks internal consistency; coordinated alteration of all records is not detectable without an external trusted digest.",
                        "Nominal charge/spin and convergence do not establish electronic-state identity, ground states, chemical accuracy, a saddle, reaction connectivity, or a molecular machine.",
                        "Source-literature correctness, source-to-XYZ conversion, and historical initial guesses are not independently established by this audit."],
        "conversion_convention": {"hartree_to_ev": HARTREE_TO_EV,
                                  "ev_per_kcal_per_mol": EV_PER_KCAL_PER_MOL,
                                  "energy_absolute_tolerance": ENERGY_ABS_TOL,
                                  "scope": "Pinned archived conversion convention; independent arithmetic, no production-code import."},
        "archives": reports,
    }
    print(json.dumps(output, indent=2, allow_nan=False))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
