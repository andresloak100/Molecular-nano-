"""V2 desired-behavior checks for imported evidence, using synthetic files only.

The malformed-quantum and null-path tests intentionally assert the desired safe
behavior and failed before the recorded V1 repairs. No solver is
imported or run, and no saved scientific evidence is changed.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from workbench.server import WorkbenchStore


VALID_QUANTUM = {
    "charge": 0, "spin": 0, "xc": "pbe0", "basis": "sto-3g",
    "dispersion": None, "grid_level": 3, "conv_tol": 1e-9,
    "max_cycle": 150, "threads": 1, "memory_mb": 2000,
    "density_fit": False, "scf_initial_guess": "minao",
}


def write_json(path, value):
    content = json.dumps(value, sort_keys=True).encode()
    path.write_bytes(content)
    return hashlib.sha256(content).hexdigest()


def write_completed_campaign(project, *, quantum=None, stage="singlepoint", damage=None):
    """Write consistent digests so each test isolates its named schema defect."""
    root = project / "selected"
    root.mkdir()
    quantum = deepcopy(VALID_QUANTUM if quantum is None else quantum)
    design = {
        "schema_version": 1, "length_unit": "angstrom", "quantum": quantum,
        "initial": "initial.xyz", "final": "final.xyz",
        "fixed_indices": [], "metadata": {},
    }
    hashes = {"design_sha256": write_json(root / "design.json", design)}
    xyz = b"2\nsynthetic V2 evidence fixture\nH 0 0 0\nH 0.74 0 0\n"
    for state in ("initial", "final"):
        (root / f"{state}.xyz").write_bytes(xyz)
        hashes[f"{state}_sha256"] = hashlib.sha256(xyz).hexdigest()
    controls = {"stage": stage, "state": "initial", "fmax": .03, "steps": 2, "images": 7}
    plan = {
        "schema_version": 1, "quantum_settings": quantum, "controls": controls,
        "designs": [{"id": "probe", "design": "design.json", "pose": {}, "input_hashes": hashes}],
    }
    plan_hash = write_json(root / "plan.json", plan)
    summary = {
        "energy_ev": -1., "free_force_max_ev_per_angstrom": .01,
        "forces_ev_per_angstrom": [[0, 0, .01], [0, 0, .01]],
        "quantum_diagnostics": {"scf_converged": True, "gradient_completed": True, "settings": quantum},
        "geometry_converged": True, "endpoint_identity_ok": True,
        "topology_screen": {"preserved": True},
    }
    result = {
        "status": "completed", "stage": stage, "state": "initial",
        "quantum_settings": quantum, "input_hashes": hashes,
        "optimization": {"fmax_ev_per_angstrom": .03, "max_steps_per_stage": 2, "images": 7},
        "structure": summary,
    }
    if stage == "path":
        result.update(
            endpoints_converged=True, neb_converged=True,
            neb_force_max_ev_per_angstrom=.01,
            endpoints=[deepcopy(summary), deepcopy(summary)],
            images=[deepcopy(summary) for _ in range(7)],
            relative_energies_ev=[0.] * 7,
        )
    if damage:
        damage(result)
    output = root / "runs" / "probe"
    output.mkdir(parents=True)
    result_hash = write_json(output / "result.json", result)
    write_json(root / "campaign.json", {
        "schema_version": 1, "plan_sha256": plan_hash,
        "attempts": {"probe": [{"status": "completed", "output": "runs/probe", "result_sha256": result_hash}]},
    })


def imported_row(project):
    store = WorkbenchStore(project, campaign_root=project, campaigns=["selected"])
    return store.snapshot.catalog["campaigns"][0]["entries"][0]


class ImportedEvidenceRegressionTests(unittest.TestCase):
    def test_well_formed_fixture_is_completed(self):
        for stage in ("singlepoint", "path"):
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as temp:
                project = Path(temp)
                write_completed_campaign(project, stage=stage)
                row = imported_row(project)
                self.assertEqual(row["status"], "completed")
                self.assertTrue(row["status_verified"])
                self.assertEqual(row["integrity_errors"], [])

    def test_matching_hashes_do_not_validate_malformed_quantum_settings(self):
        # Before the repair, matching invalid fields still produced "completed".
        for bad_fields in ({"basis": 7}, {"spin": "banana"}):
            with self.subTest(bad_fields=bad_fields), tempfile.TemporaryDirectory() as temp:
                project = Path(temp)
                write_completed_campaign(project, quantum={**VALID_QUANTUM, **bad_fields})
                row = imported_row(project)
                self.assertEqual(row["status"], "invalid evidence")
                self.assertFalse(row["status_verified"])
                self.assertTrue(row["integrity_errors"])

    def test_null_path_endpoint_is_invalid_evidence_without_aborting_catalog(self):
        # Before the repair, this raised an uncaught AttributeError.
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            write_completed_campaign(project, stage="path",
                                     damage=lambda result: result["endpoints"].__setitem__(0, None))
            row = imported_row(project)
            self.assertEqual(row["status"], "invalid evidence")
            self.assertFalse(row["status_verified"])
            self.assertTrue(row["integrity_errors"])

    def test_null_path_topology_is_invalid_evidence_without_aborting_catalog(self):
        # Before the repair, this raised an uncaught AttributeError.
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            write_completed_campaign(project, stage="path",
                                     damage=lambda result: result["endpoints"][0].__setitem__("topology_screen", None))
            row = imported_row(project)
            self.assertEqual(row["status"], "invalid evidence")
            self.assertFalse(row["status_verified"])
            self.assertTrue(row["integrity_errors"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
