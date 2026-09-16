"""Independent transfer-coordinate projection; no electronic calculations.

The score measures a direction in a chosen mass metric. It is a declared
screening heuristic, never proof of endpoint connectivity or state identity.
"""
from __future__ import annotations

import numpy as np


def transfer_projection(characterization, reaction_indices, symbols, *, minimum_overlap=0.25):
    """Return squared normalized overlap with grad(q), q = r(H,A)-r(H,D).

    Cartesian mode v and coordinate gradient g use the mass metric M:
      overlap = (g dot v)^2 / ((v.T M v) * (g.T inv(M) g)).
    This treats isotope masses explicitly and is invariant to mode sign,
    mode amplitude, translation and simultaneous proper rotations.
    """
    if isinstance(minimum_overlap, bool) or not np.isfinite(minimum_overlap) or not 0 < minimum_overlap <= 1:
        raise ValueError("minimum_overlap must be finite and in (0, 1]")
    positions = np.asarray(characterization["geometry_angstrom"], dtype=float)
    masses = np.asarray(characterization["free_masses_amu"], dtype=float)
    free = characterization["free_atom_indices"]
    modes = np.asarray(characterization["cartesian_modes_per_sqrt_amu"], dtype=float)
    frequencies = np.asarray(characterization["frequencies_cm1"], dtype=float)
    tolerance = characterization["settings"]["frequency_tolerance_cm1"]
    if positions.shape != (len(symbols), 3) or not np.isfinite(positions).all():
        raise ValueError("geometry must contain one finite Cartesian position per symbol")
    if (not free or any(type(i) is not int or i < 0 or i >= len(symbols) for i in free)
            or len(set(free)) != len(free)):
        raise ValueError("invalid free atom indices")
    if masses.shape != (len(free),) or not np.isfinite(masses).all() or np.any(masses <= 0):
        raise ValueError("invalid free atom masses")
    if (frequencies.ndim != 1 or modes.shape != (len(frequencies), len(free), 3)
            or not np.isfinite(modes).all() or not np.isfinite(frequencies).all()):
        raise ValueError("invalid mode/frequency arrays")
    if isinstance(tolerance, bool) or not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("invalid frequency resolution tolerance")
    a, h, d = [reaction_indices[key] for key in ("acceptor_carbon", "transferring_hydrogen", "donor_carbon")]
    if (any(type(i) is not int or i < 0 or i >= len(symbols) for i in (a,h,d))
            or len({a,h,d}) != 3):
        raise ValueError("reaction indices must identify three distinct atoms")
    if [symbols[i] for i in (a,h,d)] != ["C", "H", "C"]:
        raise ValueError("reaction indices must identify carbon, hydrogen, carbon")
    if any(i not in free for i in (a,h,d)):
        return {"assessed": False, "reason": "reaction atom is outside the free-mode space",
                "establishes_connectivity": False}
    ha, hd = positions[h] - positions[a], positions[h] - positions[d]
    la, ld = np.linalg.norm(ha), np.linalg.norm(hd)
    if min(la, ld) <= 1e-12:
        raise ValueError("reaction coordinate is undefined at coincident reaction atoms")
    ua, ud = ha/la, hd/ld
    full_gradient = np.zeros_like(positions)
    full_gradient[h], full_gradient[a], full_gradient[d] = ua-ud, -ua, ud
    gradient = full_gradient[free]
    inverse_mass_norm = float(np.sum(gradient**2 / masses[:,None]))
    if not np.isfinite(inverse_mass_norm) or inverse_mass_norm <= 0:
        raise ValueError("reaction gradient has no finite norm in the mass metric")
    rows = []
    local = {atom: row for row, atom in enumerate(free)}
    for index, frequency in enumerate(frequencies):
        if frequency >= -tolerance:
            continue
        mode = modes[index]
        norm = float(np.sum(masses[:,None] * mode**2))
        if not np.isfinite(norm) or norm <= 0:
            raise ValueError("resolved negative mode has zero or invalid mass norm")
        derivative = float(np.sum(gradient * mode))
        score = derivative**2 / (norm * inverse_mass_norm)
        if not np.isfinite(score) or score > 1 + 1e-10:
            raise ValueError("transfer projection violates normalized overlap bounds")
        score = min(1., max(0., score))
        h_share = float(masses[local[h]] * np.sum(mode[local[h]]**2) / norm)
        to_acceptor = float(np.dot(mode[local[h]] - mode[local[a]], ua))
        to_donor = float(np.dot(mode[local[h]] - mode[local[d]], ud))
        rows.append({"mode_index": index, "frequency_cm1": float(frequency),
                     "squared_mass_metric_overlap": score,
                     "hydrogen_mass_weighted_displacement_share": h_share,
                     "h_acceptor_distance_derivative": to_acceptor,
                     "h_donor_distance_derivative": to_donor,
                     "exceeds_declared_overlap_threshold": score >= minimum_overlap,
                     "dq_per_unit_mode": derivative})
    return {"assessed": True, "coordinate": "q = distance(H,acceptor)-distance(H,donor)",
            "metric": "squared normalized overlap using recorded free-atom masses",
            "minimum_overlap": float(minimum_overlap),
            "threshold_is_declared_heuristic": True,
            "negative_modes": rows,
            "single_negative_mode_passes_overlap_heuristic": len(rows) == 1 and rows[0]["exceeds_declared_overlap_threshold"],
            "establishes_connectivity": False, "establishes_electronic_state": False,
            "limitations": ["A chosen local reaction coordinate does not establish reactant/product connectivity.",
                            "Coordinate overlap can reflect carbon motion; hydrogen participation and bond-distance derivatives are separate diagnostics.",
                            "The threshold is a screening choice, not calibrated chemical accuracy or a transfer probability."]}
