"""E4 exact synthetic AO-metric fixtures; no molecular solver or state claim."""
from __future__ import annotations

import numpy as np


def metric_axes():
    """Columns are orthonormal in S, not in the Euclidean array metric."""
    overlap = np.diag([2.0, 0.5, 4.0])
    axes = np.diag([1.0 / np.sqrt(2.0), np.sqrt(2.0), 0.5])
    return overlap, axes


def occupied_rotation():
    """Exact 3/5,4/5 rotation of the first two orbital columns."""
    return np.array([[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]])


def unrestricted_case(*, hidden_spin_change=False, gauge_rotation=False):
    """Three electrons: two alpha and one beta in synthetic AO coordinates.

    Default alpha spans axes 1/2 and beta occupies axis 3. The hidden-spin
    variant has alpha axes 1/3 and beta axis 2: equal total density, changed
    spin channels. Expected squared metric distances are 2 per spin; occupied
    overlap singular values are alpha [1,0], beta [0].
    """
    overlap, axes = metric_axes()
    alpha = axes.copy()
    beta = axes.copy()
    alpha_occ = np.array([1.0, 1.0, 0.0])
    beta_occ = np.array([0.0, 0.0, 1.0])
    if hidden_spin_change:
        alpha_occ = np.array([1.0, 0.0, 1.0])
        beta_occ = np.array([0.0, 1.0, 0.0])
    if gauge_rotation:
        if hidden_spin_change:
            raise ValueError("Gauge fixture rotates only the default occupied span.")
        alpha = axes @ occupied_rotation()
        beta[:, 2] *= -1.0
    return {
        "overlap": overlap,
        "mo_coeff": np.array([alpha, beta]),
        "mo_occ": np.array([alpha_occ, beta_occ]),
        "nelec": (2, 1),
    }


def restricted_case():
    """Two electrons in one spatial orbital; each spin density is u1 u1.T."""
    overlap, axes = metric_axes()
    return {
        "overlap": overlap,
        "mo_coeff": axes,
        "mo_occ": np.array([2.0, 0.0, 0.0]),
        "nelec": (1, 1),
    }
