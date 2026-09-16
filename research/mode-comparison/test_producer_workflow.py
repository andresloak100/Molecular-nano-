"""Actual workflow adapter integration using a synthetic calculator, never DFT."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write
import numpy as np
import pytest

HERE=Path(__file__).parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(HERE))
from compare_modes import ComparisonError, compare_bound, load_saved_result
from nanodesign import workflow
from nanodesign.quantum import QuantumSettings


class SyntheticProducerCalculator(Calculator):
    """Known diagonal curvature, with exact fixture reference geometry."""
    implemented_properties=["energy","forces"]
    reference=np.array([[0.,0.,2.],[.741234567891234,0.,0.]])

    def __init__(self,settings,event_log=None):
        super().__init__()
        self.settings=settings
        self.calls=0
        self.diagnostics={}

    def calculate(self,atoms=None,properties=("energy","forces"),system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        self.calls+=1
        displacement=atoms.positions-self.reference
        self.results={"energy":float(np.sum(displacement**2)/2),"forces":-displacement}
        self.diagnostics={"settings":asdict(self.settings),"scf_initial_guess":self.settings.scf_initial_guess,
            "scf_converged":True,"gradient_completed":True,"call_id":f"synthetic-{self.calls}",
            "versions":{"test":"synthetic adapter regression only"}}


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def actual_workflow_pair(tmp_path,monkeypatch):
    monkeypatch.setattr(workflow,"PySCFCalculator",SyntheticProducerCalculator)
    atoms=Atoms("HeH",positions=SyntheticProducerCalculator.reference)
    write(tmp_path/"initial.xyz",atoms);write(tmp_path/"final.xyz",atoms)
    settings=QuantumSettings(spin=1,dispersion=None,threads=1)
    design=tmp_path/"design.json"
    design.write_text(json.dumps({"schema_version":1,"length_unit":"angstrom","initial":"initial.xyz",
        "final":"final.xyz","fixed_indices":[0],"quantum":asdict(settings)}))
    other=atoms.copy();other.positions[1,0]+=.2
    source=tmp_path/"two-frames.traj"
    write(source,[other,atoms])
    source_hash=checksum(source)
    outputs=[]
    for i,step in enumerate((.002,.001)):
        out=tmp_path/f"step-{i}"
        record=workflow.run_characterization(design,source,out,image=1,step=step,fmax=.01)
        assert record["status"]=="completed"
        assert record["units"]["force"]=="eV/angstrom"
        assert record["input_hashes"]["structure_sha256"]==source_hash
        assert record["input_hashes"]["input_snapshot_sha256"]==checksum(out/"input.extxyz")
        assert record["input_hashes"]["input_snapshot_sha256"]!=source_hash
        outputs.append(out/"result.json")
    assert checksum(source)==source_hash
    return outputs


def test_actual_workflow_selected_binary_frame_loads_without_rewriting(actual_workflow_pair):
    paths=actual_workflow_pair
    before=[path.read_bytes() for path in paths]
    a,b=map(load_saved_result,paths)
    out=compare_bound(a,b)
    assert out["status"]=="compared",out["findings"]
    assert a["source"]["snapshot_hash_binding"]=="input_snapshot_sha256"
    assert a["geometry_binding"]["snapshot_relation"]=="reference_rounded_to_eight_decimal_angstrom"
    assert a["positions_angstrom"]==SyntheticProducerCalculator.reference.tolist()
    assert np.max(np.abs(out["spectral_rank_frequency_changes_cm1"])) < 1e-5
    assert [path.read_bytes() for path in paths]==before


def test_old_source_hash_cannot_silently_bind_reserialized_snapshot(actual_workflow_pair):
    path=actual_workflow_pair[0]
    data=json.loads(path.read_text())
    del data["input_hashes"]["input_snapshot_sha256"]
    path.write_text(json.dumps(data))
    with pytest.raises(ComparisonError,match="distinct input_snapshot_sha256"):
        load_saved_result(path)


def test_valid_source_hash_cannot_override_wrong_explicit_snapshot(actual_workflow_pair):
    path=actual_workflow_pair[0]
    data=json.loads(path.read_text())
    data["input_hashes"]["structure_sha256"]=checksum(path.with_name("input.extxyz"))
    data["input_hashes"]["input_snapshot_sha256"]=None
    path.write_text(json.dumps(data))
    with pytest.raises(ComparisonError,match="Explicit input snapshot hash"):
        load_saved_result(path)
