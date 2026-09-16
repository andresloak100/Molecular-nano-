"""Independent B2 export mutation/publication checks with tiny synthetic files.

No quantum dependency, existing evidence, or production test helper is used.
Post-publication cleanup failures are checked as incomplete bundles requiring
inspection; publication may already have occurred when cleanup raises.
"""
from __future__ import annotations

import errno
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import nanodesign.bundle as bundle


def setup_tree(base):
    source = base / "source"
    source.mkdir()
    (source / "record.txt").write_bytes(b"original\n")
    return source, base / "bundle"


class ExportMutationTests(unittest.TestCase):
    def test_same_size_source_edit_after_first_capture_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = setup_tree(Path(directory))
            real_read = bundle._read
            changed = False

            def read_then_change(root, relative, expected, limit, **kwargs):
                nonlocal changed
                raw = real_read(root, relative, expected, limit, **kwargs)
                if not changed and root.path == source:
                    changed = True
                    (source / relative).write_bytes(b"modified\n")
                return raw

            with patch.object(bundle, "_read", side_effect=read_then_change):
                with self.assertRaises(bundle.BundleError):
                    bundle.export_bundle(source, output)
            self.assertTrue(changed)
            self.assertFalse((output / "manifest.json").exists())
            self.assertFalse(bundle.verify_bundle(output)["integrity_verified"])

    def test_source_root_replacement_before_capture_is_rejected(self):
        for replacement in ("directory", "symlink"):
            with self.subTest(replacement=replacement), tempfile.TemporaryDirectory() as directory:
                base = Path(directory)
                source, output = setup_tree(base)
                outside = base / "outside"
                outside.mkdir()
                (outside / "record.txt").write_bytes(b"outside secret\n")
                real_read = bundle._read
                changed = False

                def replace_then_read(root, relative, expected, limit, **kwargs):
                    nonlocal changed
                    if not changed and root.path == source:
                        changed = True
                        source.rename(base / "original-source")
                        if replacement == "symlink":
                            source.symlink_to(outside, target_is_directory=True)
                        else:
                            source.mkdir()
                            (source / "record.txt").write_bytes(b"replacement\n")
                    return real_read(root, relative, expected, limit, **kwargs)

                with patch.object(bundle, "_read", side_effect=replace_then_read):
                    with self.assertRaises(bundle.BundleError):
                        bundle.export_bundle(source, output)
                self.assertTrue(changed)
                self.assertFalse((output / "manifest.json").exists())
                self.assertFalse((output / "payload" / "record.txt").exists())
                self.assertEqual((outside / "record.txt").read_bytes(), b"outside secret\n")

    def test_partial_payload_write_cannot_publish_manifest(self):
        self.check_interrupted_write("payload/record.txt")

    def test_partial_pending_manifest_write_cannot_publish_manifest(self):
        self.check_interrupted_write(".manifest.pending")

    def check_interrupted_write(self, target):
        with tempfile.TemporaryDirectory() as directory:
            source, output = setup_tree(Path(directory))
            real_write = bundle._write
            interrupted = False

            def partial_write(root, relative, raw):
                nonlocal interrupted
                if relative == target:
                    interrupted = True
                    real_write(root, relative, raw[:3])
                    raise OSError(errno.EIO, "B2 injected interrupted write")
                return real_write(root, relative, raw)

            with patch.object(bundle, "_write", side_effect=partial_write):
                with self.assertRaises(bundle.BundleError):
                    bundle.export_bundle(source, output)
            self.assertTrue(interrupted)
            self.assertFalse((output / "manifest.json").exists())
            self.assertFalse(bundle.verify_bundle(output)["integrity_verified"])
            self.assertEqual((source / "record.txt").read_bytes(), b"original\n")

    def test_link_publication_failure_cannot_publish_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = setup_tree(Path(directory))
            with patch.object(bundle.os, "link", side_effect=OSError(errno.EIO, "B2 injected link failure")):
                with self.assertRaises(bundle.BundleError):
                    bundle.export_bundle(source, output)
            self.assertFalse((output / "manifest.json").exists())
            self.assertFalse(bundle.verify_bundle(output)["integrity_verified"])
            self.assertEqual((source / "record.txt").read_bytes(), b"original\n")

    def test_post_publication_cleanup_error_is_not_silently_verified(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = setup_tree(Path(directory))
            real_unlink = bundle.os.unlink

            def fail_pending_cleanup(path, *args, **kwargs):
                if path == ".manifest.pending":
                    raise OSError(errno.EIO, "B2 injected pending cleanup failure")
                return real_unlink(path, *args, **kwargs)

            with patch.object(bundle.os, "unlink", side_effect=fail_pending_cleanup):
                with self.assertRaises(bundle.BundleError):
                    bundle.export_bundle(source, output)
            # The verifier must reject the leftover pending file even if the
            # exporter has already linked a syntactically valid final manifest.
            self.assertFalse(bundle.verify_bundle(output)["integrity_verified"])
            self.assertEqual((source / "record.txt").read_bytes(), b"original\n")
            if (output / "manifest.json").exists():
                self.assertEqual((output / "manifest.json").read_bytes(),
                                 (output / ".manifest.pending").read_bytes())


if __name__ == "__main__":
    source_path = Path(bundle.__file__)
    print("bundle_sha256_before=" + hashlib.sha256(source_path.read_bytes()).hexdigest(), flush=True)
    result = unittest.main(verbosity=2, exit=False).result
    print("bundle_sha256_after=" + hashlib.sha256(source_path.read_bytes()).hexdigest(), flush=True)
    raise SystemExit(not result.wasSuccessful())
