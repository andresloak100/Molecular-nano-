"""Read-only workbench integration and path-boundary checks; no quantum jobs."""

from collections import Counter
import hashlib
from html.parser import HTMLParser
import json
import math
import os
from pathlib import Path
import tempfile
import threading
import unittest
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from workbench.server import (
    EvidenceError, PROJECT_ROOT, SafeRoot, WorkbenchStore, decode_json,
    completed_result_errors, make_server, parse_xyz,
)


class RecordedEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.store = WorkbenchStore()

    def test_actual_candidate_coordinates_identity_and_default_distance(self):
        structure = self.store.structure(self.store.snapshot.catalog["default_structure_id"])
        symbols, positions = parse_xyz((PROJECT_ROOT / "data/validation/h-abstraction-direct-initial/structure.extxyz").read_bytes())
        self.assertEqual(structure["symbols"], symbols)
        self.assertEqual(structure["positions"], positions)
        self.assertEqual(Counter(symbols), {"C": 22, "H": 31})
        self.assertEqual(structure["fixed_indices"], [7, 8, 9, 35, 36, 37])
        self.assertEqual(structure["hydrogen_transfer"], {"donor": 0, "hydrogen": 10, "acceptor": 26})
        self.assertAlmostEqual(math.dist(positions[0], positions[26]), 3.6, places=12)
        self.assertIn("not computed electronic bond orders", structure["bond_source"])

    def test_nine_pending_poses_have_distinct_exact_geometry_and_endpoints(self):
        campaign = self.store.snapshot.catalog["campaigns"][0]
        rows = campaign["entries"]
        self.assertEqual(len(rows), 9)
        self.assertEqual({row["status"] for row in rows}, {"pending"})
        self.assertTrue(campaign["plan_hash_matches_ledger"])
        coordinate_sets = set()
        for row in rows:
            initial = self.store.structure(row["initial_structure_id"])
            final = self.store.structure(row["final_structure_id"])
            coordinate_sets.add(tuple(tuple(point) for point in initial["positions"]))
            expected = math.hypot(row["pose"]["separation_angstrom"], row["pose"]["lateral_offset_angstrom"])
            self.assertAlmostEqual(math.dist(initial["positions"][0], initial["positions"][26]), expected, places=8)
            self.assertNotEqual(initial["positions"][10], final["positions"][10])
            self.assertIsNone(initial["recorded_result"])
        self.assertEqual(len(coordinate_sets), 9)

    def test_energy_force_spin_status_and_reference_values_are_recorded_not_synthetic(self):
        catalog = self.store.snapshot.catalog
        for result in catalog["calculations"]:
            raw, _ = self.store.evidence(result["id"])
            record = json.loads(raw)
            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["structure"], record["structure"])
            self.assertEqual(result["quantum_settings"], record["quantum_settings"])
            self.assertEqual(result["elapsed_seconds"], record["elapsed_seconds"])
            self.assertFalse(result["validation"]["design_validated"])
            self.assertEqual(result["path_status"], "not computed")
            self.assertEqual(result["vibrations_status"], "not computed")
            self.assertIn("timing_caveat", result["interpretation"])
        for result in catalog["references"]:
            raw, _ = self.store.evidence(result["id"])
            record = json.loads(raw)
            self.assertEqual(result["computed"], record["computed"])
            self.assertEqual(result["species"], record["species"])
            self.assertIsNot(result.get("electronic_state_identity_verified"), True)
            self.assertFalse(result["saddle_verified"])
        minao, atom = catalog["references"]
        self.assertNotEqual(minao["computed"]["ccsd_t"], atom["computed"]["ccsd_t"])

    def test_retrospective_annotation_is_separate_and_bound_to_exact_original(self):
        record = self.store.snapshot.catalog["references"][0]
        raw, _ = self.store.evidence(record["id"])
        annotation = record["retrospective_annotation"]
        self.assertEqual(record["retrospective_annotation_status"], "verified source binding")
        self.assertEqual(annotation["source_sha256"], hashlib.sha256(raw).hexdigest())
        self.assertEqual(record["interpretation"], json.loads(raw)["interpretation"])
        self.assertEqual(annotation["inferred_cc_initial_guess"], "minao")
        self.assertIsNone(record["cc_settings"].get("scf_initial_guess"))
        self.assertIn("inferred minao", record["label"])


class BoundaryAndImportTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.boundary = SafeRoot(self.root)

    def test_forbidden_paths_and_symlinks_are_rejected(self):
        for path in ("../secret", "/etc/passwd", "a/../../secret", "a\\secret", "", "."):
            with self.subTest(path=path), self.assertRaises(EvidenceError):
                self.boundary.read_bytes(path)
        (self.root / "real.json").write_text("{}")
        (self.root / "linked.json").symlink_to(self.root / "real.json")
        (self.root / "linked-directory").symlink_to(self.root, target_is_directory=True)
        for path in ("linked.json", "linked-directory/real.json"):
            with self.subTest(path=path), self.assertRaises(EvidenceError):
                self.boundary.read_bytes(path)
        with self.assertRaises(EvidenceError):
            self.boundary.directory("linked-directory")

    def test_nonregular_file_is_rejected_without_blocking(self):
        os.mkfifo(self.root / "pipe.json")
        with self.assertRaises(EvidenceError):
            self.boundary.read_bytes("pipe.json")

    def test_extended_xyz_property_order_and_coordinate_precision(self):
        symbols, positions = parse_xyz(b'2\nProperties=pos:R:3:species:S:1:forces:R:3\n0.123456789123456 0 0 C 0 0 0\n1.223456789123456 0 0 H 0 0 0\n')
        self.assertEqual(symbols, ["C", "H"])
        self.assertEqual(positions[0][0], 0.123456789123456)
        self.assertEqual(positions[1][0], 1.223456789123456)

    def test_invalid_xyz_and_nonfinite_json_are_rejected(self):
        for value in (b"", b"2\n\nC 0 0 0\n", b"1\n\nH nan 0 0\n", b"1\nProperties=species:S:1\nH\n"):
            with self.subTest(value=value), self.assertRaises(EvidenceError):
                parse_xyz(value)
        for value in (b"", b"[]", b'{"x": NaN}', b'{"x": 1e999}'):
            with self.subTest(value=value), self.assertRaises(EvidenceError):
                decode_json(value)

    def test_empty_missing_and_failed_records_remain_explicit(self):
        store = WorkbenchStore(self.root, campaigns=[])
        self.assertIsNone(store.snapshot.catalog["default_structure_id"])
        self.assertEqual({row["status"] for row in store.snapshot.catalog["calculations"]}, {"missing"})
        direct = self.root / "data/validation/h-abstraction-direct-initial"
        direct.mkdir(parents=True)
        (direct / "result.json").write_text("")
        self.assertEqual(store.catalog()["calculations"][0]["status"], "unreadable")
        failed = {"status": "failed", "error": {"message": "SCF did not converge"}}
        (direct / "result.json").write_text(json.dumps(failed))
        result = store.catalog()["calculations"][0]
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["error"], failed["error"])
        self.assertNotIn("structure", result)

    def write_campaign(self, *, design_path="design/design.json", attempts=None, quantum=None):
        campaign = self.root / "selected"
        campaign.mkdir()
        quantum = {"basis": "sto-3g", "spin": 1} if quantum is None else quantum
        (campaign / "design").mkdir()
        (campaign / "design/design.json").write_text(json.dumps({
            "schema_version": 1, "length_unit": "angstrom", "quantum": quantum,
            "initial": "initial.xyz", "final": "final.xyz", "fixed_indices": [0],
            "metadata": {"description": "<script>not executable</script>"},
        }))
        for name in ("initial.xyz", "final.xyz"):
            (campaign / "design" / name).write_text("2\n\nC 0 0 0\nH 1.1 0 0\n")
        hashes = {f"{key}_sha256": hashlib.sha256((campaign / "design" / name).read_bytes()).hexdigest()
                  for key, name in (("design", "design.json"), ("initial", "initial.xyz"), ("final", "final.xyz"))}
        plan = {"schema_version": 1, "quantum_settings": quantum,
                "controls": {"stage": "singlepoint", "state": "initial", "fmax": .03, "steps": 200, "images": 7},
                "designs": [{"id": "one", "design": design_path, "pose": {}, "input_hashes": hashes}]}
        encoded = json.dumps(plan).encode()
        (campaign / "plan.json").write_bytes(encoded)
        ledger = {"schema_version": 1, "plan_sha256": hashlib.sha256(encoded).hexdigest(), "attempts": {"one": [] if attempts is None else attempts}}
        (campaign / "campaign.json").write_text(json.dumps(ledger))
        return campaign

    def result_record(self, campaign, status="completed"):
        plan = json.loads((campaign / "plan.json").read_text())
        return {"status": status, "stage": "singlepoint", "state": "initial",
                "quantum_settings": plan["quantum_settings"], "input_hashes": plan["designs"][0]["input_hashes"],
                "optimization": {"fmax_ev_per_angstrom": .03, "max_steps_per_stage": 200, "images": 7},
                "structure": {"energy_ev": -1.0, "free_force_max_ev_per_angstrom": .1,
                    "forces_ev_per_angstrom": [[0, 0, 2], [0, 0, .1]],
                    "quantum_diagnostics": {"scf_converged": True, "gradient_completed": True}}}

    def campaign_row(self):
        return WorkbenchStore(PROJECT_ROOT, campaign_root=self.root, campaigns=["selected"]).snapshot.catalog["campaigns"][0]["entries"][0]

    def write_attempt(self, campaign, record):
        attempt = campaign / "runs/one/attempt-1"
        attempt.mkdir(parents=True, exist_ok=True)
        encoded = json.dumps(record).encode()
        (attempt / "result.json").write_bytes(encoded)
        ledger = json.loads((campaign / "campaign.json").read_text())
        ledger["attempts"]["one"] = [{"status": record["status"], "output": "runs/one/attempt-1",
                                        "result_sha256": hashlib.sha256(encoded).hexdigest()}]
        (campaign / "campaign.json").write_text(json.dumps(ledger))

    def test_selected_campaign_boundary_and_attempt_failure_visibility(self):
        campaign = self.write_campaign(attempts=[{"status": "failed", "output": "runs/one/attempt-1"}])
        attempt = campaign / "runs/one/attempt-1"
        attempt.mkdir(parents=True)
        result = self.result_record(campaign, status="failed")
        result["error"] = {"message": "numerical failure"}
        (attempt / "result.json").write_text(json.dumps(result))
        store = WorkbenchStore(PROJECT_ROOT, campaign_root=self.root, campaigns=["selected"])
        row = store.snapshot.catalog["campaigns"][0]["entries"][0]
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["calculation_status"], "failed")
        self.assertEqual(row["recorded_result"]["error"]["message"], "numerical failure")
        with self.assertRaises(EvidenceError):
            WorkbenchStore(PROJECT_ROOT, campaign_root=self.root, campaigns=["../outside"])

    def test_imported_design_traversal_does_not_escape_selected_campaign(self):
        self.write_campaign(design_path="../outside.json")
        (self.root / "outside.json").write_text('{"secret": true}')
        store = WorkbenchStore(PROJECT_ROOT, campaign_root=self.root, campaigns=["selected"])
        row = store.snapshot.catalog["campaigns"][0]["entries"][0]
        self.assertEqual(row["status"], "invalid evidence")
        self.assertIn("traversal", row["load_error"])
        self.assertNotIn(b'{"secret": true}', [value[0] for value in store.snapshot.raw.values()])

    def test_imported_coordinate_units_and_schema_must_be_explicit(self):
        campaign = self.write_campaign()
        path = campaign / "design/design.json"
        original = json.loads(path.read_text())
        for change in ({"length_unit": "bohr"}, {"length_unit": None}, {"schema_version": 2}):
            with self.subTest(change=change):
                path.write_text(json.dumps({**original, **change}))
                row = self.campaign_row()
                self.assertEqual(row["status"], "invalid evidence")
                self.assertIn("units are never inferred", row["load_error"])
                self.assertIsNone(row["initial_structure_id"])

    def test_tampered_plan_or_input_bytes_invalidate_pending_status(self):
        campaign = self.write_campaign()
        self.assertEqual(self.campaign_row()["status"], "pending")
        (campaign / "design/initial.xyz").write_text("2\n\nC 0 0 0\nH 1.2 0 0\n")
        row = self.campaign_row()
        self.assertEqual(row["recorded_status"], "pending")
        self.assertEqual(row["status"], "invalid evidence")
        self.assertFalse(row["input_hashes_match"])
        (campaign / "plan.json").write_text((campaign / "plan.json").read_text() + "\n")
        self.assertTrue(any("Plan hash" in error for error in self.campaign_row()["integrity_errors"]))

    def test_completed_ledger_without_result_is_invalid_not_completed(self):
        self.write_campaign(attempts=[{"status": "completed", "output": "runs/missing", "result_sha256": "0" * 64}])
        row = self.campaign_row()
        self.assertEqual(row["recorded_status"], "completed")
        self.assertEqual(row["status"], "invalid evidence")
        self.assertEqual(row["calculation_status"], "missing")

    def test_completed_attempt_requires_matching_digest_and_provenance(self):
        campaign = self.write_campaign()
        record = self.result_record(campaign)
        self.write_attempt(campaign, record)
        self.assertEqual(self.campaign_row()["status"], "completed")
        path = campaign / "runs/one/attempt-1/result.json"
        path.write_text(path.read_text() + "\n")
        row = self.campaign_row()
        self.assertEqual(row["status"], "invalid evidence")
        self.assertFalse(row["result_hash_matches"])
        for field, value in (("input_hashes", {}), ("stage", "relax"), ("quantum_settings", {"basis": "wrong"}), ("optimization", {})):
            with self.subTest(field=field):
                self.write_attempt(campaign, {**record, field: value})
                row = self.campaign_row()
                self.assertEqual(row["status"], "invalid evidence")
                self.assertTrue(row["result_hash_matches"])

    def test_completed_attempt_cannot_hide_missing_or_inconsistent_force_evidence(self):
        campaign = self.write_campaign()
        record = self.result_record(campaign)
        record["structure"]["free_force_max_ev_per_angstrom"] = 0
        self.write_attempt(campaign, record)
        row = self.campaign_row()
        self.assertEqual(row["status"], "invalid evidence")
        self.assertTrue(any("residual" in error for error in row["integrity_errors"]))

    def test_selected_root_replacement_cannot_follow_a_new_directory(self):
        (self.root / "allowed/selected").mkdir(parents=True)
        (self.root / "outside").mkdir()
        (self.root / "allowed/selected/value.json").write_text('{"inside": true}')
        (self.root / "outside/value.json").write_text('{"outside": true}')
        selected = SafeRoot(self.root / "allowed").directory("selected")
        (self.root / "allowed/selected").rename(self.root / "allowed/original")
        (self.root / "allowed/selected").symlink_to(self.root / "outside", target_is_directory=True)
        with self.assertRaises((EvidenceError, OSError)):
            selected.read_bytes("value.json")

    def test_old_snapshot_ids_are_rejected_instead_of_returning_new_bytes(self):
        store = WorkbenchStore()
        first_id = store.snapshot.catalog["default_structure_id"]
        first_source = store.structure(first_id)["source_href"].split("id=", 1)[1]
        self.assertTrue(store.evidence(first_source)[0])
        newer = store.catalog()
        self.assertNotEqual(first_id, newer["default_structure_id"])
        with self.assertRaises(KeyError):
            store.structure(first_id)
        with self.assertRaises(KeyError):
            store.evidence(first_source)

    def test_hash_matched_nonnumeric_or_nonstring_quantum_controls_are_invalid(self):
        campaign = self.write_campaign(quantum={"basis": 7, "spin": "banana"})
        self.write_attempt(campaign, self.result_record(campaign))
        row = self.campaign_row()
        self.assertEqual(row["status"], "invalid evidence")
        self.assertFalse(row["status_verified"])
        self.assertTrue(any("quantum" in value for value in row["integrity_errors"]))

    def test_null_path_records_are_errors_not_uncaught_exceptions(self):
        controls = {"stage": "path", "images": 7, "fmax": .03}
        summary = {"energy_ev": -1.0, "free_force_max_ev_per_angstrom": 0,
                   "forces_ev_per_angstrom": [[0, 0, 0]],
                   "quantum_diagnostics": {"scf_converged": True, "gradient_completed": True},
                   "geometry_converged": True, "topology_screen": {"preserved": True}}
        for endpoint in (None, {**summary, "topology_screen": None}):
            result = {"endpoints_converged": True, "neb_converged": True,
                      "endpoints": [endpoint, summary], "images": [summary] * 7,
                      "relative_energies_ev": [0] * 7}
            self.assertTrue(completed_result_errors(result, controls, 1, []))

    def test_annotation_cannot_infer_a_guess_for_changed_source_bytes(self):
        archive = self.root / "data/validation/paired-ccpvdz"
        archive.mkdir(parents=True)
        original = b'{"status":"completed","interpretation":"Original wording","cc_settings":{}}'
        (archive / "method_comparison.json").write_bytes(original + b"\n")
        (archive / "interpretation.json").write_text(json.dumps({
            "schema_version": 1, "kind": "retrospective_interpretation", "source": "method_comparison.json",
            "source_sha256": hashlib.sha256(original).hexdigest(), "inferred_cc_initial_guess": "minao",
        }))
        reference = WorkbenchStore(self.root, campaigns=[]).snapshot.catalog["references"][0]
        self.assertEqual(reference["retrospective_annotation_status"], "invalid source binding")
        self.assertTrue(reference["annotation_errors"])
        self.assertNotIn("inferred", reference["label"])
        self.assertEqual(reference["interpretation"], "Original wording")


class HTTPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = make_server(WorkbenchStore(), port=0)
        cls.worker = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.worker.start()
        cls.base = f"http://127.0.0.1:{cls.server.server_port}"

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.worker.join(timeout=2)

    def test_loopback_catalog_structure_and_source_routes(self):
        self.assertEqual(self.server.server_address[0], "127.0.0.1")
        with urlopen(self.base + "/api/catalog") as response:
            catalog = json.load(response)
            self.assertEqual(response.headers["X-Content-Type-Options"], "nosniff")
            self.assertIn("script-src 'self'", response.headers["Content-Security-Policy"])
        with urlopen(self.base + "/api/structure?id=" + catalog["default_structure_id"]) as response:
            structure = json.load(response)
        self.assertEqual(structure["atom_count"], 53)
        with urlopen(self.base + structure["source_href"]) as response:
            self.assertEqual(hashlib.sha256(response.read()).hexdigest(), structure["source_sha256"])

    def test_page_references_only_locally_served_assets(self):
        class Assets(HTMLParser):
            def __init__(self):
                super().__init__()
                self.paths = []

            def handle_starttag(self, tag, attrs):
                values = dict(attrs)
                if tag == "script" and "src" in values:
                    self.paths.append(values["src"])
                if tag == "link" and values.get("rel") == "stylesheet":
                    self.paths.append(values["href"])

        parser = Assets()
        with urlopen(self.base + "/") as response:
            parser.feed(response.read().decode())
        self.assertEqual(set(parser.paths), {"/app.css", "/app.js"})
        for path in parser.paths:
            with urlopen(self.base + path) as response:
                self.assertEqual(response.status, 200)
                self.assertGreater(len(response.read()), 100)

    def test_arbitrary_files_unregistered_ids_and_writes_are_unavailable(self):
        for endpoint in ("/../README.md", "/api/evidence?id=../../README.md", "/api/structure?id=/etc/passwd"):
            with self.subTest(endpoint=endpoint), self.assertRaises(HTTPError) as caught:
                urlopen(self.base + endpoint)
            self.assertEqual(caught.exception.code, 404)
            caught.exception.close()
        for method in ("POST", "PUT", "DELETE", "PATCH"):
            with self.subTest(method=method), self.assertRaises(HTTPError) as caught:
                urlopen(Request(self.base + "/api/catalog", data=b"{}", method=method))
            self.assertEqual(caught.exception.code, 405)
            caught.exception.close()

    def test_foreign_host_header_is_rejected(self):
        with self.assertRaises(HTTPError) as caught:
            urlopen(Request(self.base + "/api/catalog", headers={"Host": "untrusted.example"}))
        self.assertEqual(caught.exception.code, 403)
        caught.exception.close()


if __name__ == "__main__":
    unittest.main()
