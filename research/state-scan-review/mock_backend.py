"""Synthetic solver for S3 tests; never imports or invokes quantum chemistry."""

from dataclasses import asdict
from uuid import uuid4

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree
import numpy as np


def install_probe(monkeypatch, module):
    """Replace the scan's calculator with an observable, controllable ASE fake."""
    control = {"calls": [], "instances": [], "faults": {}, "energies": {}}

    class ProbeCalculator(Calculator):
        implemented_properties = ["energy", "forces"]

        def __init__(self, settings, *, event_log=None, **kwargs):
            super().__init__(**kwargs)
            self.settings = settings
            self.event_log = event_log
            self.diagnostics = {}
            control["instances"].append(self)

        def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            guess = self.settings.scf_initial_guess
            control["calls"].append({
                "calculator": self,
                "settings": asdict(self.settings),
                "positions": atoms.positions.copy(),
                "numbers": atoms.numbers.copy(),
            })
            fault = control["faults"].get(guess)
            self.diagnostics = {
                "call_id": str(uuid4()), "settings": asdict(self.settings),
                "scf_initial_guess": guess, "scf_converged": False,
                "gradient_completed": False, "initial_guess_scan_performed": False,
                "ground_state_verified": False, "electronic_state_identity_verified": False,
            }
            if fault == "failure":
                raise RuntimeError("S3 deliberately failed synthetic calculation")
            if fault == "interrupt":
                raise KeyboardInterrupt("S3 deliberately interrupted synthetic calculation")
            energy = control["energies"].get(guess, -1.0)
            forces = np.arange(len(atoms) * 3, dtype=float).reshape(-1, 3) / 10 + 0.1
            self.diagnostics.update(
                scf_converged=True, gradient_completed=True,
                total_energy_hartree=energy / Hartree,
                hartree_eV=float(Hartree),
                s2=self.settings.spin / 2 * (self.settings.spin / 2 + 1),
                expected_s2=self.settings.spin / 2 * (self.settings.spin / 2 + 1),
                s2_deviation=0.0,
            )
            if fault == "nonfinite_energy":
                energy = float("nan")
            elif fault == "malformed_forces":
                forces = np.zeros((len(atoms), 2))
            elif fault == "missing_gradient":
                del self.diagnostics["gradient_completed"]
            elif fault == "wrong_settings":
                self.diagnostics["settings"]["basis"] = "unexpected-basis"
            elif fault == "mutate_atoms":
                atoms.positions[0, 0] += 0.5
            self.results = {"energy": energy, "forces": forces}

    monkeypatch.setattr(module, "PySCFCalculator", ProbeCalculator)
    return control
