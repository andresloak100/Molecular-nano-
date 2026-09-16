"""Exact polynomial oracles for F2 software tests; not molecular evidence.

Energies are evaluated as rational numbers at the actual supplied binary-float
coordinates, then rounded once into the saved float. Large offsets deliberately
allow tests to expose information lost during that final energy rounding.
"""
from copy import deepcopy
from fractions import Fraction
import hashlib


REFERENCE = [[0.3125, 0.0, 0.0], [1.5, 0.25, -0.125]]


def rational(value):
    return Fraction(value)


def polynomial_record(positions, *, reference=None, offset=-100.0, slope=1.25,
                      curvature=4.0, cubic=0.0, quartic=0.0, anchor_slope=3.0):
    """Return an energy and raw forces from an independently known function.

    Atom 0 x is the derivative coordinate; atom 1 has a nonzero raw x load.
    The other coordinates have harmonic energies about the reference geometry.
    Marking atom 1 as fixed belongs to the evidence adapter, not this potential.
    """
    reference = REFERENCE if reference is None else reference
    delta = [[rational(x) - rational(x0) for x, x0 in zip(row, ref)]
             for row, ref in zip(positions, reference)]
    x = delta[0][0]
    g, k, c, q = map(rational, (slope, curvature, cubic, quartic))
    energy = rational(offset) + g*x + k*x*x/2 + c*x**3 + q*x**4
    forces = [[Fraction(0) for _ in range(3)] for _ in reference]
    forces[0][0] = -(g + k*x + 3*c*x*x + 4*q*x**3)
    for atom, row in enumerate(delta):
        for axis, value in enumerate(row):
            if (atom, axis) == (0, 0):
                continue
            energy += value*value
            forces[atom][axis] = -2*value
    energy += rational(anchor_slope)*delta[1][0]
    forces[1][0] -= rational(anchor_slope)
    return {"energy_ev": float(energy),
            "forces_ev_per_angstrom": [[float(x) for x in row] for row in forces],
            "exact_energy_fraction": str(energy)}


def displaced(offset, *, reference=None, atom=0, axis=0):
    positions = deepcopy(REFERENCE if reference is None else reference)
    positions[atom][axis] += offset
    return positions


def exact_three_point_force(*, a, b, slope=1.25, cubic=0.0, quartic=0.0):
    """Symbolic result at x0 for this polynomial, independent of F1 code.

    Quadratic terms cancel for any positive a,b. The cubic/quartic residuals
    expose truncation dependence; they are not chemical uncertainty bounds.
    """
    a, b = rational(a), rational(b)
    return float(-rational(slope) - rational(cubic)*a*b - rational(quartic)*a*b*(b-a))


def stencil_document(pairs=((0.125, 0.125),), *, state_evidence_id=None, **potential):
    """Encode independent polynomial values in F1's proposed evidence schema.

    Hashes and IDs below identify explicit synthetic fixture records. They are
    not asserted to authenticate a solver execution or an electronic state.
    """
    from nanodesign.quantum import QuantumSettings
    context = {
        "atomic_numbers": [1, 1], "fixed_indices": [1], "pbc": [False, False, False],
        "quantum_settings": QuantumSettings(spin=0, basis="sto-3g", dispersion=None, threads=1).to_dict(),
        "software_versions": {"fixture": "1"},
        "resolved_numerics": None,
        "energy_scope": "total_potential_energy",
        "units": {"length": "angstrom", "energy": "eV", "force": "eV/angstrom"},
        "force_scope": "raw_unconstrained",
    }

    def record(positions, identity):
        values = polynomial_record(positions, **potential)
        return {
            "context": deepcopy(context), "positions_angstrom": deepcopy(positions),
            "source_sha256": hashlib.sha256(f"F2 synthetic fixture {identity}".encode()).hexdigest(),
            "source_record_id": identity, "calculation_call_id": identity,
            "status": "completed", "scf_converged": True, "gradient_completed": True,
            "energy_ev": values["energy_ev"],
            "forces_ev_per_angstrom": values["forces_ev_per_angstrom"],
            "state_evidence_id": state_evidence_id,
        }

    baseline = record(REFERENCE, "baseline")
    stencils = []
    for index, (a, b) in enumerate(pairs):
        stencils.append({
            "atom_index": 0, "axis": 0,
            "requested_minus_offset_angstrom": -a,
            "requested_plus_offset_angstrom": b,
            "reference": deepcopy(baseline),
            "minus": record(displaced(-a), f"minus-{index}"),
            "plus": record(displaced(b), f"plus-{index}"),
        })
    return {"schema_version": 1, "stencils": stencils}
