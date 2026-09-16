"""Checks that the reduced models really are the candidate with a smaller handle.

The whole reduced-model argument depends on the substrate, the tip and the
pose being untouched. If truncation quietly moved the tip or changed the
target C-H, the cost numbers would describe a different system than the one
they are offered as a substitute for. These are structural checks only; they
establish nothing about chemistry.

    python -m pytest research/candidate-feasibility/test_reduced_models.py -q
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

from nanodesign.candidates import make_h_abstraction  # noqa: E402
from reduced_models import HANDLES, TIP_APEX, TIP_DISTAL, reduced_candidate, tool_fragment  # noqa: E402

EXPECTED_ATOMS = {"hydrogen": 29, "methyl": 32, "adamantyl": 53}


@pytest.mark.parametrize("handle", HANDLES)
def test_substrate_and_tip_are_bit_identical_to_the_candidate(handle):
    """Truncation must not move a single substrate or tip atom."""
    reference, _, _ = make_h_abstraction(3.6, 0.0)
    reduced, _ = reduced_candidate(handle)
    shared = list(range(26)) + [TIP_APEX, TIP_DISTAL]
    assert [reduced.symbols[i] for i in shared] == [reference.symbols[i] for i in shared]
    np.testing.assert_allclose(reduced.positions[shared], reference.positions[shared], atol=0.0)


@pytest.mark.parametrize("handle,count", EXPECTED_ATOMS.items())
def test_atom_counts(handle, count):
    reduced, info = reduced_candidate(handle)
    assert len(reduced) == count
    assert info["handle"] == handle


@pytest.mark.parametrize("handle", HANDLES)
def test_pose_is_preserved(handle):
    """The donor-to-apex distance is the pose; it must survive truncation."""
    reduced, _ = reduced_candidate(handle)
    distance = float(np.linalg.norm(reduced.positions[TIP_APEX] - reduced.positions[0]))
    assert distance == pytest.approx(3.6, abs=1e-9)


@pytest.mark.parametrize("handle", HANDLES)
def test_no_atoms_are_too_close(handle):
    reduced, _ = reduced_candidate(handle)
    distances = reduced.get_all_distances()
    np.fill_diagonal(distances, np.inf)
    assert float(distances.min()) > 0.7


@pytest.mark.parametrize("handle", HANDLES)
def test_reduced_models_are_neutral_doublets_with_the_spin_on_the_apex(handle):
    reduced, info = reduced_candidate(handle)
    assert info["total_charge"] == 0
    assert info["spin"] == 1
    moments = reduced.get_initial_magnetic_moments()
    assert moments[TIP_APEX] == pytest.approx(1.0)
    assert moments.sum() == pytest.approx(1.0)


@pytest.mark.parametrize("handle", HANDLES)
def test_anchor_scheme_is_declared_as_different_not_scaled(handle):
    """The report leans on this distinction, so the metadata must carry it."""
    _, info = reduced_candidate(handle)
    assert "anchor_scheme" in info
    if handle != "adamantyl":
        assert "different mechanical boundary condition" in info["anchor_scheme"]
        assert len(info["fixed_indices"]) == 4  # three substrate carbons + one handle atom
    else:
        assert len(info["fixed_indices"]) == 6


@pytest.mark.parametrize("handle", HANDLES)
def test_tool_fragments_have_consistent_composition(handle):
    """R-CC-H must be exactly R-CC* plus one hydrogen."""
    radical, radical_info = tool_fragment(handle, hydrogenated=False)
    hydrogenated, hydrogenated_info = tool_fragment(handle, hydrogenated=True)
    assert len(hydrogenated) == len(radical) + 1
    assert radical_info["spin"] == 1 and hydrogenated_info["spin"] == 0
    assert hydrogenated.get_chemical_symbols()[-1] == "H"
    # The shared atoms must be identical, so the affinity is not contaminated
    # by a geometry change between the two species.
    np.testing.assert_allclose(hydrogenated.positions[: len(radical)], radical.positions, atol=0.0)


@pytest.mark.parametrize("handle", HANDLES)
def test_tool_fragment_tip_matches_the_candidate_tip(handle):
    reference, _, _ = make_h_abstraction(3.6, 0.0)
    radical, _ = tool_fragment(handle, hydrogenated=False)
    np.testing.assert_allclose(radical.positions[0], reference.positions[TIP_APEX], atol=0.0)
    np.testing.assert_allclose(radical.positions[1], reference.positions[TIP_DISTAL], atol=0.0)


def test_added_hydrogen_sits_on_the_apex_away_from_the_handle():
    for handle in HANDLES:
        radical, _ = tool_fragment(handle, hydrogenated=False)
        hydrogenated, info = tool_fragment(handle, hydrogenated=True)
        apex = hydrogenated.positions[0]
        added = hydrogenated.positions[info["added_hydrogen_index"]]
        assert float(np.linalg.norm(added - apex)) == pytest.approx(1.06, abs=1e-9)
        # Anti-parallel to the tool axis: the H arrives from the substrate side.
        axis = radical.positions[1] - radical.positions[0]
        axis /= np.linalg.norm(axis)
        assert float(np.dot(added - apex, axis)) < 0.0


def test_methyl_hydrogens_are_tetrahedral_and_equidistant():
    fragment, _ = tool_fragment("methyl", hydrogenated=False)
    carbon = fragment.positions[2]
    hydrogens = fragment.positions[3:6]
    lengths = [float(np.linalg.norm(h - carbon)) for h in hydrogens]
    assert lengths == pytest.approx([1.09, 1.09, 1.09], abs=1e-9)
    axis = fragment.positions[1] - carbon
    axis /= np.linalg.norm(axis)
    for hydrogen in hydrogens:
        direction = (hydrogen - carbon) / np.linalg.norm(hydrogen - carbon)
        angle = np.degrees(np.arccos(np.clip(np.dot(direction, axis), -1.0, 1.0)))
        assert angle == pytest.approx(109.4712206, abs=1e-6)


def test_unknown_handle_is_rejected():
    with pytest.raises(ValueError, match="handle must be one of"):
        reduced_candidate("buckyball")
    with pytest.raises(ValueError, match="handle must be one of"):
        tool_fragment("buckyball", hydrogenated=False)
