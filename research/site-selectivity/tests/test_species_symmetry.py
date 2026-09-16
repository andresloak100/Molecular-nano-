"""Why stage 1 needs exactly two adamantyl radicals, and not a site survey.

There are two different degeneracies in this lane and they must not be confused,
because one is rigorous and the other is an artifact of an input geometry.

1. **Census shells.** The 15 competing hydrogens fall into four groups by how far
   the tool apex must move to reach them.  Those groups are a property of the
   *tool's placement relative to the cage* in an unrelaxed candidate, and
   relaxation or a different tool orientation could lift them by an amount the
   census cannot bound.  Treating them as a symmetry reduction would be a
   screening shortcut.

2. **Species identity in stage 1.** Isolated adamantane has T_d symmetry. Its 4
   bridgehead hydrogens form a single orbit and its 12 methylene hydrogens form
   a single orbit, so there are exactly *two* distinct C-H environments and
   therefore exactly two distinct adamantyl radicals.  Deleting H14 or H16 or
   H15 does not produce three similar molecules to be sampled; it produces the
   *same molecule* in three orientations.  This is rigorous, survives
   relaxation, and is the reason stage 1 computes two radicals rather than
   surveying twelve sites.

This file verifies (2) numerically from the actual coordinates, rather than
asserting the point group.  The test is a congruence check on a
rotation/translation-invariant fingerprint: two structures with identical sorted
interatomic distance spectra to machine precision are congruent, hence the same
molecule, hence exactly equal in energy for any method.
"""

from __future__ import annotations

import numpy as np
import pytest

from species import _radical, _target_cage, build_species


def fingerprint(atoms):
    """Sorted interatomic distances, split by element pair.

    Invariant to translation, rotation and reflection, and to atom ordering.
    """
    symbols = atoms.get_chemical_symbols()
    positions = np.asarray(atoms.positions, dtype=float)
    groups: dict[str, list[float]] = {}
    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            key = "".join(sorted((symbols[i], symbols[j])))
            groups.setdefault(key, []).append(
                float(np.linalg.norm(positions[i] - positions[j]))
            )
    return {key: np.sort(np.array(value)) for key, value in sorted(groups.items())}


def assert_congruent(first, second, tolerance=1e-10):
    left, right = fingerprint(first), fingerprint(second)
    assert set(left) == set(right)
    for key in left:
        assert left[key].shape == right[key].shape
        assert np.abs(left[key] - right[key]).max() < tolerance, key


def cage_environments():
    """Group the cage's hydrogens by the substitution of their carbon."""
    cage = _target_cage()
    symbols = cage.get_chemical_symbols()
    positions = np.asarray(cage.positions, dtype=float)
    distances = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=-1)
    bridgehead, methylene = [], []
    for hydrogen, symbol in enumerate(symbols):
        if symbol != "H":
            continue
        carbon = int(
            np.argmin([distances[hydrogen, k] if symbols[k] == "C" else np.inf
                       for k in range(len(symbols))])
        )
        carbon_neighbours = sum(
            1 for k in range(len(symbols))
            if symbols[k] == "C" and k != carbon and distances[carbon, k] < 1.8
        )
        (bridgehead if carbon_neighbours == 3 else methylene).append((hydrogen, carbon))
    return bridgehead, methylene


def test_adamantane_has_four_bridgehead_and_twelve_methylene_hydrogens():
    bridgehead, methylene = cage_environments()
    assert len(bridgehead) == 4
    assert len(methylene) == 12


def test_every_methylene_radical_is_the_same_molecule():
    """All 12 methylene hydrogens are one symmetry orbit, so one species.

    If this holds, computing a second, third or sixth methylene representative
    cannot produce a different energy, and a "survey" of the nearest shell would
    be confirming arithmetic rather than measuring chemistry.
    """
    _bridgehead, methylene = cage_environments()
    reference = _radical(*methylene[0])
    for hydrogen, carbon in methylene[1:]:
        assert_congruent(reference, _radical(hydrogen, carbon))


def test_every_bridgehead_radical_is_the_same_molecule():
    bridgehead, _methylene = cage_environments()
    reference = _radical(*bridgehead[0])
    for hydrogen, carbon in bridgehead[1:]:
        assert_congruent(reference, _radical(hydrogen, carbon))


def test_the_two_radical_classes_are_genuinely_different_molecules():
    """The reduction must not be so aggressive that it merges the two sites."""
    bridgehead, methylene = cage_environments()
    with pytest.raises(AssertionError):
        assert_congruent(_radical(*bridgehead[0]), _radical(*methylene[0]))


def test_stage_one_species_set_is_exactly_the_distinct_chemistry():
    """Six species, and none of them is a redundant site representative."""
    species = build_species()
    assert set(species) == {
        "adamantane", "adamantyl_bridgehead", "adamantyl_methylene",
        "hydrogen_atom", "ethynyl", "acetylene",
    }
    # Exactly two of the six are adamantyl radicals: one per C-H environment.
    radicals = [name for name in species if name.startswith("adamantyl_")]
    assert len(radicals) == 2
    assert {species[name].site_type for name in radicals} == {
        "bridgehead_tertiary", "methylene_secondary"
    }
    formulas = {name: species[name].formula for name in species}
    assert formulas["adamantyl_bridgehead"] == formulas["adamantyl_methylene"] == "C10H15"
    assert formulas["adamantane"] == "C10H16"


def test_the_two_radicals_are_isomers_so_differences_cancel_cleanly():
    """Identical composition is what makes the site difference well conditioned."""
    species = build_species()
    left = species["adamantyl_bridgehead"].atoms.get_chemical_symbols()
    right = species["adamantyl_methylene"].atoms.get_chemical_symbols()
    assert sorted(left) == sorted(right)
