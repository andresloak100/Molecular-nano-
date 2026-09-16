"""The six species needed for both site-selectivity thermochemical cycles.

The adamantane cage is sliced out of the 53-atom candidate itself, via the
public ``nanodesign.candidates.make_h_abstraction`` entry point, so the isolated
molecule studied here is constructed identically to the target cage in the
design being screened.  Nothing in this module optimizes a geometry or performs
an electronic-structure calculation; the coordinates are the same unoptimized
ideal-lattice starting guesses the candidate uses.

In isolated adamantane all four bridgehead hydrogens are equivalent and all
twelve methylene hydrogens are equivalent under the molecular point group, so
exactly two distinct adamantyl radicals exist and one calculation covers each
class.  The site census in ``site_census.py`` confirms that symmetry numerically
from the actual coordinates.

``spin`` in this repository follows PySCF: it is N_alpha - N_beta, i.e. 2S, not
the multiplicity.  Every doublet radical therefore carries ``spin=1``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase import Atoms

from nanodesign.candidates import make_h_abstraction

# Ethynyl and acetylene reference geometries, Angstrom.  These are starting
# guesses only; every species in this study is relaxed before its energy is
# used, so the guess affects cost and not the reported result.
_ETHYNYL_CC = 1.21
_ETHYNYL_CH = 1.06
_ACETYLENE_CC = 1.20
_ACETYLENE_CH = 1.06


@dataclass(frozen=True)
class Species:
    """One species in a thermochemical cycle, with its explicit electronic state."""

    name: str
    atoms: Atoms
    spin: int
    charge: int
    description: str
    site_type: str | None = None

    @property
    def formula(self) -> str:
        return self.atoms.get_chemical_formula()


def _target_cage() -> Atoms:
    """Isolated adamantane, exactly as the candidate's target cage is built."""
    reactant, _product, metadata = make_h_abstraction()
    substrate = list(metadata["substrate_indices"])
    cage = reactant[substrate]
    cage.set_constraint()
    cage.set_initial_magnetic_moments(np.zeros(len(cage)))
    cage.info.clear()
    if cage.get_chemical_formula() != "C10H16":
        raise ValueError(f"Expected adamantane C10H16, built {cage.get_chemical_formula()}")
    return cage


def _radical(hydrogen_index: int, carbon_index: int, cage: Atoms | None = None) -> Atoms:
    """Adamantyl radical made by deleting one cage hydrogen.

    ``cage`` optionally supplies an already relaxed adamantane geometry, which
    only changes the starting point of the radical's own relaxation and so
    affects cost rather than the converged result.
    """
    cage = _target_cage() if cage is None else cage.copy()
    symbols = cage.get_chemical_symbols()
    if symbols[hydrogen_index] != "H":
        raise ValueError(f"Index {hydrogen_index} is {symbols[hydrogen_index]}, not H")
    keep = [i for i in range(len(cage)) if i != hydrogen_index]
    radical = cage[keep]
    if radical.get_chemical_formula() != "C10H15":
        raise ValueError(f"Expected C10H15, built {radical.get_chemical_formula()}")
    # The unpaired electron sits on the carbon that lost its hydrogen.  Index
    # shifts by one for atoms after the deleted hydrogen; cage carbons all
    # precede the hydrogens, so the carbon index is unchanged, but assert it.
    if radical.get_chemical_symbols()[carbon_index] != "C":
        raise ValueError(f"Radical centre index {carbon_index} is not carbon after deletion")
    moments = np.zeros(len(radical))
    moments[carbon_index] = 1.0
    radical.set_initial_magnetic_moments(moments)
    return radical


def build_species(relaxed_cage: Atoms | None = None) -> dict[str, Species]:
    """Return every species keyed by name, with its state fixed explicitly.

    ``relaxed_cage`` optionally seeds the two radicals from an already relaxed
    adamantane instead of the ideal-lattice guess.  This is a cost optimization
    only: each radical is still relaxed independently to the same force
    threshold, so the converged energies do not depend on the seed.
    """
    cage = _target_cage()

    # Bridgehead (tertiary) C-H: carbon 0, hydrogen 10 -- the site the 53-atom
    # candidate targets.  Methylene (secondary) C-H: carbon 4, hydrogen 14 --
    # the nearest alternative site found by the geometric census.
    bridgehead = _radical(hydrogen_index=10, carbon_index=0, cage=relaxed_cage)
    methylene = _radical(hydrogen_index=14, carbon_index=4, cage=relaxed_cage)

    hydrogen = Atoms("H", positions=[[0.0, 0.0, 0.0]], pbc=False)
    hydrogen.set_initial_magnetic_moments([1.0])

    ethynyl = Atoms(
        "CCH",
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, _ETHYNYL_CC], [0.0, 0.0, _ETHYNYL_CC + _ETHYNYL_CH]],
        pbc=False,
    )
    ethynyl.set_initial_magnetic_moments([1.0, 0.0, 0.0])

    acetylene = Atoms(
        "HCCH",
        positions=[
            [0.0, 0.0, -_ACETYLENE_CH],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, _ACETYLENE_CC],
            [0.0, 0.0, _ACETYLENE_CC + _ACETYLENE_CH],
        ],
        pbc=False,
    )

    species = [
        Species(
            name="adamantane",
            atoms=cage,
            spin=0,
            charge=0,
            description="Isolated adamantane C10H16, the candidate's target cage, closed-shell singlet",
        ),
        Species(
            name="adamantyl_bridgehead",
            atoms=bridgehead,
            spin=1,
            charge=0,
            description="1-adamantyl radical: bridgehead tertiary C-H broken at cage carbon 0",
            site_type="bridgehead_tertiary",
        ),
        Species(
            name="adamantyl_methylene",
            atoms=methylene,
            spin=1,
            charge=0,
            description="2-adamantyl radical: methylene secondary C-H broken at cage carbon 4",
            site_type="methylene_secondary",
        ),
        Species(
            name="hydrogen_atom",
            atoms=hydrogen,
            spin=1,
            charge=0,
            description="Isolated hydrogen atom, doublet; the free-atom dissociation partner",
        ),
        Species(
            name="ethynyl",
            atoms=ethynyl,
            spin=1,
            charge=0,
            description="Ethynyl radical C2H, doublet; the abstracting tip species",
        ),
        Species(
            name="acetylene",
            atoms=acetylene,
            spin=0,
            charge=0,
            description="Acetylene C2H2, closed-shell singlet; the abstraction product",
        ),
    ]
    return {item.name: item for item in species}


# Species that contain an unpaired electron and therefore require the
# initial-guess scan described in the shared coordination notes.
OPEN_SHELL = ("adamantyl_bridgehead", "adamantyl_methylene", "hydrogen_atom", "ethynyl")
