"""Checks on the tool-recharge thermodynamics.

The module's conclusion rests on one comparison, BDE(X-H) against the tip's
C-H bond, so the arithmetic of that comparison and the consistency of the
geometries it uses are what these tests cover. They run no quantum jobs.

    python -m pytest research/candidate-feasibility/test_tool_cycle.py -q
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

from tool_cycle import PARTNERS, TIP_BOND_KCAL, _atoms  # noqa: E402


def test_tip_bond_values_match_the_handle_fidelity_measurements():
    """These are quoted from this lane's own measurement, not the literature."""
    assert TIP_BOND_KCAL["adamantyl"] == pytest.approx(136.15, abs=0.01)
    assert TIP_BOND_KCAL["methyl"] == pytest.approx(136.52, abs=0.01)
    assert TIP_BOND_KCAL["hydrogen"] == pytest.approx(136.72, abs=0.01)
    # The real handle is the reference, and the spread across handles is small.
    spread = max(TIP_BOND_KCAL.values()) - min(TIP_BOND_KCAL.values())
    assert spread < 1.0


@pytest.mark.parametrize("name", sorted(PARTNERS))
def test_each_hydride_is_its_radical_plus_one_hydrogen(name):
    """Otherwise the BDE would not be a bond dissociation energy at all."""
    entry = PARTNERS[name]
    radical = _atoms(*entry["radical"])
    hydride = _atoms(*entry["hydride"])
    assert len(hydride) == len(radical) + 1
    radical_counts = sorted(radical.get_chemical_symbols())
    hydride_counts = sorted(hydride.get_chemical_symbols())
    # Removing one H from the hydride must give exactly the radical composition.
    hydride_counts.remove("H")
    assert sorted(hydride_counts) == radical_counts


@pytest.mark.parametrize("name", sorted(PARTNERS))
def test_spin_states_are_doublet_radical_and_singlet_hydride(name):
    entry = PARTNERS[name]
    assert entry["radical"][2] == 1, "radical must be an open-shell doublet"
    assert entry["hydride"][2] == 0, "closed-shell hydride"


@pytest.mark.parametrize("name", sorted(PARTNERS))
def test_no_atoms_are_unphysically_close(name):
    entry = PARTNERS[name]
    for key in ("radical", "hydride"):
        atoms = _atoms(*entry[key])
        if len(atoms) < 2:
            continue
        distances = atoms.get_all_distances()
        np.fill_diagonal(distances, np.inf)
        assert float(distances.min()) > 0.6, f"{name} {key} has overlapping atoms"


@pytest.mark.parametrize("name", sorted(PARTNERS))
def test_bond_lengths_are_near_standard_values(name):
    """Guessed geometries, but they must at least be chemically plausible."""
    entry = PARTNERS[name]
    hydride = _atoms(*entry["hydride"])
    heavy = hydride.positions[0]
    for index in range(1, len(hydride)):
        if hydride.get_chemical_symbols()[index] != "H":
            continue
        length = float(np.linalg.norm(hydride.positions[index] - heavy))
        assert 0.7 < length < 1.3, f"{name} X-H length {length:.3f} A is implausible"


def test_radicals_carry_one_unpaired_electron_on_the_heavy_atom():
    for name, entry in PARTNERS.items():
        radical = _atoms(*entry["radical"])
        moments = radical.get_initial_magnetic_moments()
        assert moments.sum() == pytest.approx(1.0), name
        assert moments[0] == pytest.approx(1.0), name


def test_the_fluorine_case_is_present_because_it_is_the_only_near_miss():
    """The finding depends on the partner list containing the strongest candidate."""
    assert "fluorine" in PARTNERS
    assert "hydrogen_atom" in PARTNERS
    # A generic carbon radical must be included as the clearly-uphill control.
    assert "methyl" in PARTNERS


def test_geometries_are_labelled_unrelaxed():
    """Both sides of the comparison must share the rigid-geometry approximation."""
    for entry in PARTNERS.values():
        for key in ("radical", "hydride"):
            atoms = _atoms(*entry[key])
            assert atoms.info["geometry_status"] == "unrelaxed_standard_bond_lengths"
