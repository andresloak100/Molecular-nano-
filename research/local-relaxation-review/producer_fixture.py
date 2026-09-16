"""Independent synthetic H1 acquisition; no electronic-structure calculation."""
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.calculators.calculator import Calculator, all_changes
from nanodesign.quantum import QuantumSettings

ROOT = Path(__file__).resolve().parents[2]

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

class HarmonicAcquisition(Calculator):
    implemented_properties = ['forces']
    def __init__(self, settings, reference):
        super().__init__()
        self.settings = settings
        self.reference = np.array(reference, copy=True)
        self.calls = 0
        self.diagnostics = {}
        self.H = np.diag([2., 3., 4., 5., 6., 7.])
        self.H[0, 3] = self.H[3, 0] = .4
        self.g = np.array([.012, -.009, .006, -.015, .008, -.004])
    def calculate(self, atoms=None, properties=('forces',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls += 1
        delta = (atoms.positions[[0, 2]] - self.reference[[0, 2]]).ravel()
        forces = np.zeros((3, 3))
        forces[[0, 2]] = -(self.g + self.H @ delta).reshape(2, 3)
        forces[1] = [100., -200., 300.]  # fixed support load must not enter Newton solve
        self.results = {'forces': forces}
        self.diagnostics = {'settings': asdict(self.settings),
            'scf_initial_guess': self.settings.scf_initial_guess,
            'scf_converged': True, 'gradient_completed': True,
            'call_id': f'l2-synthetic-{self.calls}',
            'versions': {'fixture': 'analytic polynomial, not quantum evidence'}}

def produce(tmp_path, monkeypatch, masses=(1.5, 20., 3.5), hessian=None):
    h1 = load_module('l2_actual_h1', ROOT / 'research/characterization-resume/prototype.py')
    class Acquisition(HarmonicAcquisition):
        def __init__(self, settings, reference):
            super().__init__(settings, reference)
            if hessian is not None:
                self.H = np.array(hessian, copy=True)
    monkeypatch.setattr(h1, 'PySCFCalculator', Acquisition)
    atoms = Atoms('HHH', positions=[[.123456789123, 0., 0.], [2., 0., 0.], [4.234567891234, .3, 0.]])
    atoms.set_masses(masses)
    atoms.set_constraint(FixAtoms(indices=[1]))
    settings = QuantumSettings(spin=1, dispersion=None, threads=1)
    atoms.calc = Acquisition(settings, atoms.positions)
    source = tmp_path / 'source.xyz'
    from ase.io import write
    write(source, atoms, format='extxyz')
    context = {'input_hashes': {'structure_sha256': hashlib.sha256(source.read_bytes()).hexdigest()},
               'design': {'fixed_indices': [1], 'quantum': asdict(settings)},
               'synthetic_fixture': True}
    output = tmp_path / 'produced'
    result = h1.run_characterization_checkpoint(atoms, output, settings=settings, input_context=context)
    assert result['status'] == 'completed', result
    assert atoms.calc.calls == 13
    assert json.loads((output / 'result.json').read_text()) == result
    return output / 'result.json', atoms.calc.H.copy(), atoms.calc.g.copy()
