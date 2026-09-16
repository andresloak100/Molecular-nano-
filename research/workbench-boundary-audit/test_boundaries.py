"""Independent V2 desired-behavior checks using only synthetic temporary files.

No solver imports, repository-evidence edits, or permanent server processes.
A failure identifies a reproducible workbench contract gap, not a chemistry result.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from workbench.server import EvidenceError, SafeRoot, WorkbenchStore


def write_json(path, value):
    content = json.dumps(value, sort_keys=True).encode()
    path.write_bytes(content)
    return hashlib.sha256(content).hexdigest()


def write_campaign(root, separation):
    root.mkdir(parents=True, exist_ok=True)
    design = {"schema_version": 1, "length_unit": "angstrom",
              "quantum": {"basis": "sto-3g", "spin": 0},
              "initial": "initial.xyz", "final": "final.xyz",
              "fixed_indices": [], "metadata": {}}
    design_hash = write_json(root / "design.json", design)
    xyz = f"2\nsynthetic V2 fixture\nH 0 0 0\nH {separation} 0 0\n".encode()
    for name in ("initial.xyz", "final.xyz"):
        (root / name).write_bytes(xyz)
    xyz_hash = hashlib.sha256(xyz).hexdigest()
    plan = {"schema_version": 1, "quantum_settings": design["quantum"],
            "controls": {"stage": "singlepoint", "state": "initial", "fmax": .03, "steps": 2, "images": 7},
            "designs": [{"id": "probe", "design": "design.json", "pose": {"separation_angstrom": separation},
                         "input_hashes": {"design_sha256": design_hash, "initial_sha256": xyz_hash, "final_sha256": xyz_hash}}]}
    plan_hash = write_json(root / "plan.json", plan)
    write_json(root / "campaign.json", {"schema_version": 1, "plan_sha256": plan_hash, "attempts": {"probe": []}})


class BoundaryRegressionTests(unittest.TestCase):
    def test_selected_directory_replacement_cannot_escape_original_boundary(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            allowed, outside = base / "allowed", base / "outside"
            selected = allowed / "selected"
            selected.mkdir(parents=True)
            outside.mkdir()
            (selected / "probe.json").write_bytes(b'{"source":"inside"}')
            (outside / "probe.json").write_bytes(b'{"source":"outside"}')
            boundary = SafeRoot(allowed).directory("selected")
            self.assertEqual(boundary.read_bytes("probe.json"), b'{"source":"inside"}')
            selected.rename(allowed / "original-selected")
            selected.symlink_to(outside, target_is_directory=True)
            try:
                content = boundary.read_bytes("probe.json")
            except (EvidenceError, OSError):
                return  # Rejecting a replaced root is also safe.
            self.assertEqual(content, b'{"source":"inside"}',
                             "A previously selected campaign boundary followed its replacement symlink outside the configured root")

    def test_catalog_link_remains_consistent_after_another_reader_refreshes(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            campaign = project / "selected"
            write_campaign(campaign, .74)
            store = WorkbenchStore(project, campaign_root=project, campaigns=["selected"])
            first_catalog = store.snapshot.catalog
            first_row = first_catalog["campaigns"][0]["entries"][0]
            self.assertEqual(first_row["status"], "pending")
            first_id = first_row["initial_structure_id"]
            first_structure = store.structure(first_id)
            first_source = first_structure["source_href"].split("id=", 1)[1]
            first_source_bytes = store.evidence(first_source)[0]
            write_campaign(campaign, 1.48)
            second_catalog = store.catalog()  # Another browser window refreshes.
            self.assertEqual(second_catalog["campaigns"][0]["entries"][0]["status"], "pending")
            try:
                fetched = store.structure(first_id)
                source_bytes = store.evidence(first_source)[0]
            except (KeyError, EvidenceError):
                return  # Explicit stale-link rejection is better than mixing versions.
            self.assertEqual(fetched["source_sha256"], first_structure["source_sha256"],
                             "An ID in the first catalog silently returns coordinates from a later catalog")
            self.assertEqual(source_bytes, first_source_bytes,
                             "A source link silently changes bytes after another reader refreshes")


if __name__ == "__main__":
    unittest.main(verbosity=2)
