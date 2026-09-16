"""Check structural bookkeeping, not chemical or mechanosynthetic feasibility."""

from collections import Counter

import numpy as np
import pytest

from nanodesign.candidates import make_h_abstraction


def test_composition_atom_conservation_and_radical_bookkeeping():
    reactant, product, meta = make_h_abstraction()
    assert len(reactant) == 53
    assert Counter(reactant.get_chemical_symbols()) == {"C": 22, "H": 31}
    assert np.array_equal(reactant.numbers, product.numbers)
    assert meta["charge"] == 0
    assert meta["multiplicity"] == 2
    assert int(reactant.numbers.sum()) % 2 == 1
    assert np.isclose(reactant.get_initial_magnetic_moments().sum(), 1)
    assert np.isclose(product.get_initial_magnetic_moments().sum(), 1)
    moved = np.where(np.linalg.norm(reactant.positions - product.positions, axis=1) > 1e-12)[0]
    assert moved.tolist() == [meta["transferred_hydrogen"]]
    assert reactant[meta["transferred_hydrogen"]].symbol == "H"
    assert set(meta["substrate_indices"]).isdisjoint(meta["tool_indices"])
    assert sorted(meta["substrate_indices"] + meta["tool_indices"]) == list(range(53))


@pytest.mark.parametrize("endpoint_name", ["reactant", "product"])
def test_explicit_bond_valences_and_endpoint_lengths(endpoint_name):
    reactant, product, meta = make_h_abstraction()
    endpoint = reactant if endpoint_name == "reactant" else product
    valences = np.zeros(len(endpoint), dtype=int)
    bonds = meta[f"{endpoint_name}_bonds"]
    pairs = {frozenset((i, j)) for i, j, _ in bonds}
    assert len(pairs) == len(bonds)
    for i, j, order in bonds:
        valences[i] += order
        valences[j] += order
        pair = {i, j}
        if pair == {meta["tip_apex"], meta["tip_distal"]}:
            expected = 1.21
            assert order == 3
        elif pair == {meta["tip_distal"], meta["handle_bridgehead"]}:
            expected = 1.46
        elif "H" in {endpoint[i].symbol, endpoint[j].symbol}:
            expected = 1.06 if endpoint_name == "product" and meta["tip_apex"] in pair else 1.09
        else:
            expected = 1.54
        assert endpoint.get_distance(i, j) == pytest.approx(expected, abs=1e-10)
    radical = meta[f"{endpoint_name}_radical_index"]
    for i, atom in enumerate(endpoint):
        expected_valence = 1 if atom.symbol == "H" else 3 if i == radical else 4
        assert valences[i] == expected_valence


def test_adamantane_is_subdivided_tetrahedron_with_tetrahedral_bond_angles():
    reactant, _, meta = make_h_abstraction()
    cage_bonds = [(i, j) for i, j, _ in meta["reactant_bonds"] if i < 10 and j < 10]
    assert len(cage_bonds) == 12
    neighbors = {i: [] for i in range(10)}
    for i, j in cage_bonds:
        neighbors[i].append(j)
        neighbors[j].append(i)
    assert sorted(map(len, neighbors.values())) == [2] * 6 + [3] * 4
    assert {frozenset(neighbors[i]) for i in range(4, 10)} == {
        frozenset((a, b)) for a in range(4) for b in range(a + 1, 4)
    }
    for center, adjacent in neighbors.items():
        for a in range(len(adjacent)):
            for b in range(a + 1, len(adjacent)):
                assert reactant.get_angle(adjacent[a], center, adjacent[b]) == pytest.approx(
                    np.degrees(np.arccos(-1 / 3)), abs=1e-9
                )


def test_distal_carbon_anchors_match_in_both_endpoints_without_gross_overlap():
    reactant, product, meta = make_h_abstraction()
    fixed = meta["fixed_indices"]
    assert len(fixed) == 6
    assert all(reactant[i].symbol == "C" for i in fixed)
    assert not set(fixed) & {meta["target_carbon"], meta["tip_apex"], meta["transferred_hydrogen"]}
    assert np.array_equal(reactant.positions[fixed], product.positions[fixed])
    for atoms in (reactant, product):
        assert set(atoms.constraints[0].get_indices()) == set(fixed)
        distances = atoms.get_all_distances()
        np.fill_diagonal(distances, np.inf)
        assert distances.min() >= 0.7
        assert not atoms.pbc.any()


def test_separation_and_lateral_offset_translate_tool_and_preserve_target():
    initial_r, initial_p, meta = make_h_abstraction()
    shifted_r, shifted_p, shifted_meta = make_h_abstraction(4.2, 0.35)
    displacement = np.array([0.35, 0, 0.6])
    tool = meta["tool_indices"]
    target = meta["substrate_indices"]
    assert np.allclose(shifted_r.positions[tool] - initial_r.positions[tool], displacement)
    assert np.allclose(shifted_p.positions[tool] - initial_p.positions[tool], displacement)
    assert np.allclose(shifted_r.positions[target], initial_r.positions[target])
    unchanged = [i for i in target if i != meta["transferred_hydrogen"]]
    assert np.allclose(shifted_p.positions[unchanged], initial_p.positions[unchanged])
    assert np.allclose(
        shifted_p.positions[meta["transferred_hydrogen"]] - initial_p.positions[meta["transferred_hydrogen"]],
        displacement,
    )
    assert shifted_r.positions[meta["target_carbon"]] == pytest.approx([0, 0, 0])
    assert shifted_r.positions[meta["tip_apex"]] == pytest.approx([0.35, 0, 4.2])
    assert shifted_meta["target_apex_distance_angstrom"] == pytest.approx(np.hypot(4.2, 0.35))


@pytest.mark.parametrize("separation,lateral", [(0, 0), (-1, 0), (np.nan, 0), (np.inf, 0), (3.6, np.nan), (3.6, np.inf)])
def test_invalid_dimensions_fail(separation, lateral):
    with pytest.raises(ValueError):
        make_h_abstraction(separation, lateral)


def test_colliding_starting_guess_fails():
    with pytest.raises(ValueError, match="closer than 0.7"):
        make_h_abstraction(1.1)
