"""Explicit synthetic snapshot fixtures, independent of the E3 producer.

The AO metric is an exact algebra example, not Gaussian integrals for the
placeholder atom/basis metadata. No fixture represents a quantum calculation.
"""
from copy import deepcopy

import numpy as np

from e4_reference_cases import restricted_case, unrestricted_case


def build_snapshot(case=None, *, reference="UKS", call_id="e4-synthetic-call"):
    if case is None:
        case = restricted_case() if reference == "RKS" else unrestricted_case()
    alpha_count, beta_count = case["nelec"]
    total = alpha_count + beta_count
    symbol = {1: "H", 2: "He", 3: "Li"}[total]
    settings = {
        "charge": 0, "spin": alpha_count - beta_count, "xc": "pbe0",
        "basis": "e4-synthetic-three-s", "dispersion": None, "grid_level": 3,
        "conv_tol": 1e-9, "max_cycle": 150, "threads": 1,
        "memory_mb": 128, "density_fit": False, "scf_initial_guess": "minao",
    }
    coeff, occ = case["mo_coeff"], case["mo_occ"]
    if reference == "RKS":
        source = {"spatial": occ.tolist()}
        spin_coeff = [coeff, coeff]
        spin_occ = [occ / 2, occ / 2]
    else:
        source = {"alpha": occ[0].tolist(), "beta": occ[1].tolist()}
        spin_coeff, spin_occ = coeff, occ
    channels = {}
    for name, count, c, o in zip(("alpha", "beta"), case["nelec"], spin_coeff, spin_occ):
        indices = np.flatnonzero(o == 1)
        channels[name] = {
            "n_electrons": count, "occupied_indices": indices.tolist(),
            "occupied_coefficients": c[:, indices].tolist(),
        }
    return {
        "schema_version": 1, "kind": "mean_field_electronic_snapshot",
        "reference": reference, "phase": "converged_scf_only",
        "producer": {"engine": "PySCF", "version": "E4-synthetic-no-solver",
                     "class": "E4.synthetic." + reference},
        "call_id": call_id, "settings": settings, "scf_converged": True,
        "scf_energy_hartree": -1.0,
        "geometry": {"symbols": [symbol], "atomic_numbers": [total],
                     "positions_bohr": [[0.0, 0.0, 0.0]], "charge": 0,
                     "spin": alpha_count - beta_count, "frame": "solver_cartesian",
                     "bohr_angstrom": 0.52917721092, "coordinate_frame_id": "e4-synthetic-frame",
                     "pbc": [False, False, False]},
        "basis": {"n_ao": 3, "representation": "spherical",
                  "normalize_gto": True, "ecp": False, "pseudo": False,
                  "expanded_basis": {symbol: [[0, [exponent, 1.0]] for exponent in (1.0, 0.5, 0.25)]},
                  "ao_labels": [[0, symbol, f"{i + 1}s", ""] for i in range(3)],
                  "shells": [{"atom_index": 0, "angular_momentum": 0, "kappa": 0,
                              "exponents": [exponent], "contraction_coefficients": [[1.0]],
                              "libcint_coefficients": [[1.0]]}
                             for exponent in (1.0, 0.5, 0.25)]},
        "ao_overlap": case["overlap"].tolist(),
        "source_occupations": source, "channels": deepcopy(channels),
        "electronic_state_identity_verified": False, "ground_state_verified": False,
        "stability_assessed_by_capture": False,
    }
