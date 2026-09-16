"""Independent software fixtures for interruption/reuse accounting, not chemistry."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
import numpy as np


def anchored_structure():
    """Six free Cartesian coordinates, so one full stencil requires 13 records."""
    atoms = Atoms("H3", positions=[[0, 0, 0], [1.2, 0.1, 0], [0.2, 1.3, 0.4]],
                  masses=[1.0, 2.0, 3.0])
    atoms.set_constraint(FixAtoms(indices=[1]))
    return atoms


class CountingForceFixture(Calculator):
    """Linear force field with controllable failure at one attempted evaluation.

    The caller owns settings/provenance configuration. This fixture supplies no
    electronic-state evidence. Each instance records actual calculator cache
    misses separately from public force requests and previously saved records.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, reference, *, fail_at=None, failure=KeyboardInterrupt, label="fixture", settings=None):
        super().__init__()
        self.reference = np.asarray(reference, dtype=float).copy()
        self.fail_at = fail_at
        self.failure = failure
        self.label = label
        self.settings = settings
        self.calls = []
        self.diagnostics = {"synthetic_fixture": True}

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {}
        self.calls.append(atoms.positions.copy())
        self.diagnostics = {"synthetic_fixture": True, "call_id": f"{self.label}-{len(self.calls)}",
                            "settings": asdict(self.settings), "scf_converged": False,
                            "gradient_completed": False}
        if self.fail_at is not None and len(self.calls) == self.fail_at:
            self.diagnostics["deliberate_failure"] = True
            raise self.failure(f"synthetic interruption at evaluation {self.fail_at}")
        delta = atoms.positions - self.reference
        forces = -2 * delta
        forces[1] += [3, 0, -2]  # deliberately nonzero load on the fixed atom
        self.results = {"energy": float((delta * delta).sum()), "forces": forces}
        self.diagnostics.update(scf_converged=True, gradient_completed=True)


def snapshot_tree(directory):
    """Read exact fixture bytes to check that resume never changes its source."""
    root = Path(directory)
    return {path.relative_to(root).as_posix(): path.read_bytes()
            for path in sorted(root.rglob("*")) if path.is_file()}
