"""Bounded corruption tests: all writes target temporary archive copies."""

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("audit_paired_evidence.py")
SPEC = importlib.util.spec_from_file_location("audit_paired_evidence", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
VALIDATION = SCRIPT.resolve().parents[2] / "data" / "validation"


class PairedEvidenceTests(unittest.TestCase):
    def copied_archive(self, name="paired-ccpvdz-atom"):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        destination = Path(temporary.name) / name
        shutil.copytree(VALIDATION / name, destination)
        return destination

    def change(self, root, filename, keys, value):
        path = root / filename
        document = json.loads(path.read_text())
        target = document
        for key in keys[:-1]:
            target = target[key]
        target[keys[-1]] = value
        path.write_text(json.dumps(document))

    def assert_failure_at(self, root, fragment):
        result = MODULE.audit_archive(root)
        self.assertFalse(result["consistency_passed"])
        self.assertFalse(result["physical_validation_established"])
        failures = json.dumps(result["failures"])
        self.assertIn(fragment, failures, failures)

    def test_both_archives_and_historical_unknowns(self):
        old = MODULE.audit_archive(VALIDATION / "paired-ccpvdz")
        new = MODULE.audit_archive(VALIDATION / "paired-ccpvdz-atom")
        self.assertTrue(old["consistency_passed"], old["failures"])
        self.assertTrue(new["consistency_passed"], new["failures"])
        self.assertFalse(new["physical_validation_established"])
        self.assertTrue(old["warnings"])
        self.assertEqual(new["warnings"], [])
        self.assertIsNone(old["species"]["ethynyl_radical"]["cc_initial_guess_recorded"])
        self.assertEqual(new["species"]["ethynyl_radical"]["cc_initial_guess_recorded"], "atom")
        self.assertAlmostEqual(new["recomputed"]["ccsd_t"]["nominal_ts_relative_energy_kcal_per_mol"], 2.398615614424693)

    def test_changed_coordinate_bytes_and_missing_record(self):
        root = self.copied_archive()
        path = root / "dft/inputs/methane.xyz"
        path.write_text(path.read_text() + "\n")
        self.assert_failure_at(root, "manifest.hash.methane.xyz")
        root = self.copied_archive()
        (root / "ccsd_t/ethynyl_radical.json").unlink()
        self.assert_failure_at(root, "ethynyl_radical")

    def test_numeric_geometry_accounting_and_arithmetic_corruption(self):
        mutations = [
            ("method_comparison.json", ["computed", "ccsd_t", "reaction_energy_ev"], 20.0, "computed.ccsd_t.reaction_energy_ev"),
            ("method_comparison.json", ["computed", "dft", "reaction_energy_kcal_per_mol"], -1.1622528631892237, "computed.dft.reaction_energy_kcal_per_mol"),
            ("ccsd_t/ethynyl_radical.json", ["energy_ev"], True, "energy_ev"),
            ("ccsd_t/ethynyl_radical.json", ["energy_ev"], float("nan"), "Non-finite"),
            ("ccsd_t/ethynyl_radical.json", ["quantum_diagnostics", "electron_count"], 12, "electron_count"),
            ("ccsd_t/ethynyl_radical.json", ["quantum_diagnostics", "geometry", "positions_angstrom", 0, 2], 0.99, "geometry.positions"),
            ("ccsd_t/ethynyl_radical.json", ["quantum_diagnostics", "hartree_to_ev"], 1.0, "conversion_constant"),
            ("dft/species/ethynyl_radical.json", ["energy_ev"], -2100.0, "benchmark_species_copy"),
        ]
        for filename, keys, value, expected_failure in mutations:
            with self.subTest(filename=filename, keys=keys):
                root = self.copied_archive()
                self.change(root, filename, keys, value)
                self.assert_failure_at(root, expected_failure)

    def test_scientific_claim_and_missing_current_schema_flag(self):
        for flag, value in (("saddle_verified", True), ("electronic_state_identity_verified", None)):
            with self.subTest(flag=flag):
                root = self.copied_archive()
                self.change(root, "method_comparison.json", [flag], value)
                self.assert_failure_at(root, "comparison." + flag)
        root = self.copied_archive()
        path = root / "method_comparison.json"
        record = json.loads(path.read_text())
        del record["ground_state_verified"]
        path.write_text(json.dumps(record))
        self.assert_failure_at(root, "comparison.ground_state_verified")

    def test_cli_json_and_failure_exit_without_quantum_imports(self):
        good = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True, check=False)
        self.assertEqual(good.returncode, 0, good.stderr)
        self.assertTrue(json.loads(good.stdout)["consistency_passed"])
        root = self.copied_archive()
        (root / "method_comparison.json").write_text("{ invalid JSON")
        bad = subprocess.run([sys.executable, str(SCRIPT), str(root)], capture_output=True, text=True, check=False)
        self.assertEqual(bad.returncode, 1, bad.stderr)
        self.assertFalse(json.loads(bad.stdout)["consistency_passed"])


if __name__ == "__main__":
    unittest.main()
