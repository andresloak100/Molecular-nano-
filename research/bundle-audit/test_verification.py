"""Independent B2 transport-integrity probes; all evidence is synthetic."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from nanodesign import bundle
from fixtures import file_bytes, make_source


class VerificationBoundaryTests(unittest.TestCase):
    def test_roundtrip_stays_portable_and_read_only_after_original_removed(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            source, exported = base / "source", base / "exported"
            expected = make_source(source)
            manifest = bundle.export_bundle(source, exported)
            self.assertEqual(file_bytes(source), expected)
            self.assertEqual(file_bytes(exported / "payload"), expected)
            self.assertEqual(manifest["scientific_validation"], "not_assessed")
            renamed = base / "copied-for-review"
            exported.rename(renamed)
            source.rename(base / "original-no-longer-at-selected-path")
            before = file_bytes(renamed)
            report = bundle.verify_bundle(renamed)
            self.assertTrue(report["integrity_verified"], report)
            self.assertEqual(report["scientific_validation"], "not_assessed")
            self.assertEqual(file_bytes(renamed), before)

    def test_output_inside_source_via_symlink_alias_does_not_modify_source(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            source = base / "source"
            expected = make_source(source)
            alias = base / "alias"
            alias.symlink_to(source, target_is_directory=True)
            with self.assertRaises(bundle.BundleError):
                bundle.export_bundle(source, alias / "exported")
            self.assertEqual(file_bytes(source), expected)
            self.assertFalse((source / "exported").exists())

    def test_payload_changed_after_its_verified_read_is_detected(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            source, exported = base / "source", base / "exported"
            make_source(source)
            bundle.export_bundle(source, exported)
            original_read = bundle._read
            changed = False

            def read_then_modify(root, relative, *args, **kwargs):
                nonlocal changed
                raw = original_read(root, relative, *args, **kwargs)
                if relative == "payload/opaque.bin" and not changed:
                    changed = True
                    # Same length: a size-only final comparison would miss it.
                    (exported / relative).write_bytes(bytes(reversed(raw)))
                return raw

            with patch.object(bundle, "_read", side_effect=read_then_modify):
                report = bundle.verify_bundle(exported)
            self.assertTrue(changed)
            self.assertFalse(report["integrity_verified"], report)
            self.assertTrue(report["issues"])

    def test_manifest_replacement_between_initial_read_and_inventory_is_detected(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            source, exported = base / "source", base / "exported"
            make_source(source)
            bundle.export_bundle(source, exported)
            original_read = bundle._read
            changed = False

            def read_then_replace(root, relative, *args, **kwargs):
                nonlocal changed
                raw = original_read(root, relative, *args, **kwargs)
                if relative == "manifest.json" and not changed:
                    changed = True
                    record = json.loads(raw)
                    record["source_name"] = "a different captured source"
                    (exported / relative).write_text(json.dumps(record))
                return raw

            with patch.object(bundle, "_read", side_effect=read_then_replace):
                report = bundle.verify_bundle(exported)
            self.assertTrue(changed)
            self.assertFalse(report["integrity_verified"], report)
            self.assertTrue(report["issues"])

    def test_duplicate_json_key_cannot_replace_recorded_inventory(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            source, exported = base / "source", base / "exported"
            make_source(source)
            bundle.export_bundle(source, exported)
            path = exported / "manifest.json"
            original = path.read_text()
            path.write_text('{"files": [],' + original.lstrip()[1:])
            before = file_bytes(exported)
            report = bundle.verify_bundle(exported)
            self.assertFalse(report["integrity_verified"], report)
            self.assertTrue(any("Duplicate JSON key" in issue["message"] for issue in report["issues"]))
            self.assertEqual(file_bytes(exported), before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
