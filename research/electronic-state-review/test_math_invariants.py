"""Independent E4 AO-metric cases with analytic answers; no solver is run."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import math
from pathlib import Path
import unittest

import numpy as np

import nanodesign.electronic_state as electronic_state
from e4_fixtures import build_snapshot
from e4_reference_cases import metric_axes, restricted_case, unrestricted_case


class MathInvariantTests(unittest.TestCase):
    def compare(self, left, right):
        result = electronic_state.compare_snapshots(left, right)
        self.assertEqual(result["comparison_status"], "compared", result.get("reasons"))
        return result["metrics"]

    def assert_channel(self, channel, singular, distance_squared):
        np.testing.assert_allclose(channel["singular_values"], singular, rtol=0, atol=1e-12)
        self.assertAlmostEqual(channel["density_distance_squared"], distance_squared, places=12)
        expected_angles = [math.acos(value) for value in singular]
        # arccos is ill-conditioned at 1: O(epsilon) singular-value roundoff
        # becomes O(sqrt(epsilon)) angular roundoff for identical subspaces.
        np.testing.assert_allclose(channel["principal_angles_radians"], expected_angles,
                                   rtol=0, atol=5e-8)

    def test_nonorthogonal_metric_preserves_spin_electron_counts(self):
        snapshot = build_snapshot(unrestricted_case())
        diagnostics = electronic_state.validate_snapshot(snapshot)
        self.assertAlmostEqual(diagnostics["alpha"]["density_electron_trace"], 2., places=12)
        self.assertAlmostEqual(diagnostics["beta"]["density_electron_trace"], 1., places=12)
        alpha = np.asarray(snapshot["channels"]["alpha"]["occupied_coefficients"])
        beta = np.asarray(snapshot["channels"]["beta"]["occupied_coefficients"])
        # Euclidean traces intentionally disagree: omitting the AO metric must
        # not accidentally pass because the fixture uses an identity overlap.
        self.assertAlmostEqual(float(np.trace(alpha @ alpha.T)), 2.5, places=12)
        self.assertAlmostEqual(float(np.trace(beta @ beta.T)), .25, places=12)
        metrics = self.compare(snapshot, deepcopy(snapshot))
        self.assert_channel(metrics["alpha"], [1., 1.], 0.)
        self.assert_channel(metrics["beta"], [1.], 0.)

    def test_occupied_rotation_and_sign_do_not_change_the_determinant(self):
        left = build_snapshot(unrestricted_case())
        right = build_snapshot(unrestricted_case(gauge_rotation=True))
        metrics = self.compare(left, right)
        self.assert_channel(metrics["alpha"], [1., 1.], 0.)
        self.assert_channel(metrics["beta"], [1.], 0.)

    def test_orbital_permutation_moves_occupations_without_changing_the_state(self):
        case = unrestricted_case()
        permuted = deepcopy(case)
        for spin, order in enumerate(([1, 0, 2], [2, 1, 0])):
            permuted["mo_coeff"][spin] = case["mo_coeff"][spin][:, order]
            permuted["mo_occ"][spin] = case["mo_occ"][spin][order]
        permuted["mo_coeff"][0, :, 0] *= -1
        permuted["mo_coeff"][1, :, 0] *= -1
        metrics = self.compare(build_snapshot(case), build_snapshot(permuted))
        self.assert_channel(metrics["alpha"], [1., 1.], 0.)
        self.assert_channel(metrics["beta"], [1.], 0.)

    def test_equal_total_density_does_not_hide_changed_spin_channels(self):
        left = build_snapshot(unrestricted_case())
        right = build_snapshot(unrestricted_case(hidden_spin_change=True))

        def total_density(snapshot):
            coefficients = [np.asarray(snapshot["channels"][name]["occupied_coefficients"])
                            for name in ("alpha", "beta")]
            return sum(c @ c.T for c in coefficients)

        np.testing.assert_allclose(total_density(left), total_density(right), rtol=0, atol=1e-15)
        metrics = self.compare(left, right)
        self.assert_channel(metrics["alpha"], [1., 0.], 2.)
        self.assert_channel(metrics["beta"], [0.], 2.)

    def test_partial_subspace_change_has_analytic_principal_angle_and_distance(self):
        original = unrestricted_case()
        changed = deepcopy(original)
        _, axes = metric_axes()
        # Keep alpha u1, rotate occupied u2 toward virtual u3 by atan(4/3).
        changed["mo_coeff"][0, :, 1] = .6 * axes[:, 1] + .8 * axes[:, 2]
        changed["mo_coeff"][0, :, 2] = -.8 * axes[:, 1] + .6 * axes[:, 2]
        metrics = self.compare(build_snapshot(original), build_snapshot(changed))
        self.assert_channel(metrics["alpha"], [1., .6], 32 / 25)
        self.assert_channel(metrics["beta"], [1.], 0.)

    def test_restricted_double_occupation_normalizes_to_one_electron_per_spin(self):
        snapshot = build_snapshot(restricted_case(), reference="RKS")
        self.assertEqual(snapshot["source_occupations"]["spatial"], [2., 0., 0.])
        diagnostics = electronic_state.validate_snapshot(snapshot)
        for name in ("alpha", "beta"):
            self.assertEqual(snapshot["channels"][name]["n_electrons"], 1)
            self.assertAlmostEqual(diagnostics[name]["density_electron_trace"], 1., places=12)
        self.assertEqual(snapshot["channels"]["alpha"], snapshot["channels"]["beta"])
        metrics = self.compare(snapshot, deepcopy(snapshot))
        self.assert_channel(metrics["alpha"], [1.], 0.)
        self.assert_channel(metrics["beta"], [1.], 0.)


if __name__ == "__main__":
    source_path = Path(electronic_state.__file__)
    print("electronic_state_sha256_before=" + hashlib.sha256(source_path.read_bytes()).hexdigest(), flush=True)
    result = unittest.main(verbosity=2, exit=False).result
    print("electronic_state_sha256_after=" + hashlib.sha256(source_path.read_bytes()).hexdigest(), flush=True)
    raise SystemExit(not result.wasSuccessful())
