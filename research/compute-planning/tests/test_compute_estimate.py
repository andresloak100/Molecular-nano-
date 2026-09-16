"""Software-only planning checks; synthetic fixtures are not chemical evidence."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


DIRECTORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("compute_estimate", DIRECTORY / "estimate.py")
planner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(planner)


def fixtures(root):
    """Deliberately synthetic two-atom records, confined to a temporary root."""
    coordinates = b"2\nsynthetic software fixture\nH 0 0 0\nH 0 0 0.74\n"
    for name, relative in planner.ARCHIVES.items():
        folder = root / relative
        folder.mkdir(parents=True)
        (folder / "input-initial.extxyz").write_bytes(coordinates)
        record = {
            "status": "completed", "stage": "singlepoint", "state": "initial",
            "elapsed_seconds": 20 if name == "direct" else 10,
            "quantum_settings": {"density_fit": name == "density_fitting", "threads": 2},
            "structure": {"quantum_diagnostics": {
                "scf_converged": True, "gradient_completed": True,
                "basis_functions": 2, "elapsed_seconds": 9,
            }},
            "design": {"fixed_indices": []},
            "input_hashes": {"initial_sha256": hashlib.sha256(coordinates).hexdigest()},
        }
        if name == "density_fitting":
            record["structure"]["quantum_diagnostics"]["effective_pyscf_threads"] = 1
        save_fixture(folder, record)


def save_fixture(folder, record):
    raw = json.dumps(record).encode()
    (folder / "result.json").write_bytes(raw)
    (folder / "interpretation.json").write_text(json.dumps({
        "file_sha256": {"result.json": hashlib.sha256(raw).hexdigest()},
        "timing_caveat": "Synthetic fixture; not a measured scientific result.",
    }))


class PlanningTests(unittest.TestCase):
    def test_default_path_accounts_for_stage_starts_and_cache(self):
        result = planner.path_counts()
        self.assertEqual(result["iteration_only_evaluations"], 2400)
        self.assertEqual(result["stage_start_geometry_requests"], 12)
        self.assertEqual(result["ordinary_to_climbing_cache_credit"], 5)
        self.assertEqual(result["geometry_evaluations_with_cross_stage_cache"], 2407)
        self.assertEqual(result["geometry_evaluations_without_cross_stage_cache_credit"], 2412)

    def test_already_converged_stages_still_need_initial_geometries(self):
        result = planner.path_counts(endpoint_steps=0, neb_steps=0, climbing_steps=0)
        self.assertEqual(result["geometry_evaluations_with_cross_stage_cache"], 7)
        separate = planner.path_counts(endpoint_steps=3, neb_steps=2, climbing_steps=4)
        self.assertEqual(separate["geometry_evaluations_with_cross_stage_cache"], 43)

    def test_default_hessian_guard_and_request_calculation_distinction(self):
        result = planner.hessian_counts(53, [7, 8, 9, 35, 36, 37])
        self.assertEqual(result["free_cartesian_coordinates"], 141)
        self.assertEqual(result["full_stencil_force_requests"], 283)
        self.assertEqual(result["workflow_force_requests_including_guard"], 284)
        self.assertEqual(result["cold_start_quantum_evaluations_if_all_guards_pass"], 283)
        self.assertFalse(result["coordinate_guard_allows_run"])
        self.assertTrue(planner.hessian_counts(53, [7, 8, 9, 35, 36, 37], 141)["coordinate_guard_allows_run"])

    def test_invalid_scenarios_are_not_coerced(self):
        for kwargs in ({"images": 4}, {"images": True}, {"endpoint_steps": -1},
                       {"neb_steps": 1.5}, {"climbing_steps": False}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                planner.path_counts(**kwargs)
        for fixed in ([0, 0], [True], [2], [0, 1]):
            with self.subTest(fixed=fixed), self.assertRaises(ValueError):
                planner.hessian_counts(2, fixed)

    def test_duration_rejects_invalid_or_overflowing_input(self):
        self.assertEqual(planner.duration(3, 1200)["hours"], 1)
        for seconds in (None, True, 0, -1, float("nan"), float("inf")):
            with self.subTest(seconds=seconds), self.assertRaises(ValueError):
                planner.duration(3, seconds)
        with self.assertRaises(ValueError):
            planner.duration(10 ** 1000, 20)

    def test_missing_evidence_does_not_fall_back_to_estimated_timings(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(FileNotFoundError):
                planner.estimate(temporary)

    def test_fixture_report_preserves_missing_thread_data_and_source_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixtures(root)
            before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
            result = planner.estimate(root, steps=0)
            self.assertIsNone(result["measurements"]["direct"]["effective_pyscf_threads_as_recorded"])
            self.assertEqual(result["measurements"]["density_fitting"]["effective_pyscf_threads_as_recorded"], 1)
            self.assertEqual(result["projections"]["direct"]["path_with_cross_stage_cache"]["seconds"], 140)
            self.assertEqual(result["hessian"]["free_cartesian_coordinates"], 6)
            self.assertEqual(result["quantum_jobs_launched"], 0)
            self.assertFalse(result["runtime_guaranteed"])
            self.assertEqual(before, {p: p.read_bytes() for p in root.rglob("*") if p.is_file()})
            json.dumps(result, allow_nan=False)

    def test_changed_archive_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixtures(root)
            path = root / planner.ARCHIVES["direct"] / "result.json"
            path.write_bytes(path.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "result hash"):
                planner.load_measurements(root)

    def test_failed_or_nonfinite_measurement_is_rejected_even_with_valid_hash(self):
        for field, value in (("status", "failed"), ("elapsed_seconds", float("nan"))):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                fixtures(root)
                folder = root / planner.ARCHIVES["direct"]
                record = json.loads((folder / "result.json").read_text())
                record[field] = value
                save_fixture(folder, record)
                with self.assertRaises(ValueError):
                    planner.load_measurements(root)

    def test_real_archives_supply_counts_and_timings_without_quantum(self):
        result = planner.estimate()
        self.assertEqual(result["measurements"]["direct"]["elements"], {"C": 22, "H": 31})
        self.assertEqual(result["measurements"]["direct"]["basis_functions"], 463)
        self.assertAlmostEqual(result["projections"]["density_fitting"]["full_hessian_if_guards_pass"]["hours"], 55.14813890606863)

    def test_count_matches_live_workflow_with_controlled_updates_and_ase_cache(self):
        # This integration check uses an analytic counting fixture, NEVER a
        # scientific potential or quantum backend. Prescribed updates isolate
        # cost accounting from physical convergence; real ASE cache/NEB/IDPP
        # behavior and production stage transitions are exercised.
        sys.path.insert(0, str(DIRECTORY.parents[1]))
        try:
            import numpy as np
            from ase import Atoms
            from ase.calculators.calculator import Calculator, all_changes
            from nanodesign import workflow
        except ImportError:
            self.skipTest("optional integration check requires the project's installed ASE/NumPy")
        evaluations = []

        class CountingCalculator(Calculator):
            implemented_properties = ["energy", "forces"]
            diagnostics = {}

            def __init__(self, *args, **kwargs):
                super().__init__()

            def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
                super().calculate(atoms, properties, system_changes)
                evaluations.append(self.atoms.positions.copy())
                self.results = {"energy": float((self.atoms.positions ** 2).sum()),
                                "forces": -2 * self.atoms.positions}

        class PrescribedFIRE:
            def __init__(self, atoms, **kwargs):
                self.atoms = atoms

            def run(self, *, fmax, steps):
                self.atoms.get_forces()
                for _ in range(steps):
                    positions = self.atoms.get_positions()
                    positions[:, 0] += 0.003
                    self.atoms.set_positions(positions)
                    self.atoms.get_forces()
                return True  # software fixture only, not a convergence claim

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, distance in (("initial", 0.74), ("final", 1.1)):
                (root / f"{name}.xyz").write_text(f"2\nsynthetic fixture\nH 0 0 0\nH 0 0 {distance}\n")
            design = {"schema_version": 1, "length_unit": "angstrom", "initial": "initial.xyz",
                      "final": "final.xyz", "fixed_indices": [0], "quantum": {"spin": 0}}
            (root / "design.json").write_text(json.dumps(design))
            with patch.object(workflow, "PySCFCalculator", CountingCalculator), patch.object(workflow, "FIRE", PrescribedFIRE):
                workflow.run(root / "design.json", root / "out", stage="path", steps=2, images=7)
            self.assertEqual(len(evaluations), planner.path_counts(endpoint_steps=2, neb_steps=2, climbing_steps=2)["geometry_evaluations_with_cross_stage_cache"])


if __name__ == "__main__":
    unittest.main()
