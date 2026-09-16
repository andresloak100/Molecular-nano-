"""Synthetic linear-algebra fixtures plus read-only historical H2 comparison."""
from copy import deepcopy
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

HERE=Path(__file__).parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from compare_modes import ComparisonError, compare_bound, load_saved_result


def bound(step=.005):
    masses=np.array([1.,12.])
    return {"symbols":["H","C"],"positions_angstrom":[[0.,0.,0.],[1.,0.,0.]],
            "masses_amu":masses.tolist(),"quantum_settings":{
                "charge":0,"spin":1,"xc":"pbe0","basis":"def2-svp","dispersion":"d3bj",
                "grid_level":3,"conv_tol":1e-9,"max_cycle":150,"density_fit":False,"scf_initial_guess":"minao"},
            "software_versions":{"fixture":"synthetic, not a quantum result"},"notes":[],
            "stationary":{"method":"synthetic regression fixture", "free_atom_indices":[0,1],
                "frozen_atom_indices":[],"coordinate_order":"x,y,z for each atom in free_atom_indices",
                "external_modes_removed":False, "free_masses_amu":masses.tolist(),
                "cartesian_modes_per_sqrt_amu":(np.eye(6)/np.repeat(np.sqrt(masses),3)).reshape(6,2,3).tolist(),
                "frequencies_cm1":[100.,105.,250.,400.,600.,800.],
                "settings":{"step_angstrom":step,"force_tolerance_ev_per_angstrom":.03,
                            "frequency_tolerance_cm1":20.,"max_free_coordinates":120}}}


def rotate_rows(data, i, j, angle):
    modes=np.asarray(data["stationary"]["cartesian_modes_per_sqrt_amu"])
    old=modes.copy()
    modes[i]=np.cos(angle)*old[i]+np.sin(angle)*old[j]
    modes[j]=-np.sin(angle)*old[i]+np.cos(angle)*old[j]
    data["stationary"]["cartesian_modes_per_sqrt_amu"]=modes.tolist()


def test_known_rotated_degenerate_subspace():
    a,b=bound(),bound(.0025)
    rotate_rows(b,0,1,np.pi/4)
    out=compare_bound(a,b)
    assert out["status"]=="compared"
    assert np.array(out["squared_mass_metric_overlap_matrix"])[:2,:2]==pytest.approx(np.full((2,2),.5))
    assert out["groups"][0]["principal_cosines_squared"]==pytest.approx([1.,1.])
    assert out["groups"][0]["comparison_kind"]=="subspace_only"
    assert not out["groups"][0]["individual_identity_assigned"]
    assert not out["step_convergence_certified"]


def test_known_subspace_leakage():
    a,b=bound(),bound(.0025)
    rotate_rows(b,1,2,np.pi/6)
    out=compare_bound(a,b)
    assert out["groups"][0]["principal_cosines_squared"]==pytest.approx([1.,.75])
    assert out["groups"][0]["maximum_principal_angle_degrees"]==pytest.approx(30.)


def test_sign_invariance():
    a,b=bound(),bound(.0025)
    modes=np.array(b["stationary"]["cartesian_modes_per_sqrt_amu"])
    modes[[0,2,5]]*=-1
    b["stationary"]["cartesian_modes_per_sqrt_amu"]=modes.tolist()
    out=compare_bound(a,b)
    assert np.array(out["squared_mass_metric_overlap_matrix"])==pytest.approx(np.eye(6))


def test_rank_swap_is_exposed_not_identity_assigned():
    a,b=bound(),bound(.0025)
    b["stationary"]["frequencies_cm1"][2:4]=[400.,250.]
    out=compare_bound(a,b)
    assert out["squared_mass_metric_overlap_matrix"][2][3]==pytest.approx(1.)
    assert out["groups"][1]["principal_cosines_squared"]==pytest.approx([0.])
    assert not out["frequency_changes_are_identity_assignments"]


def test_connected_group_span_and_full_space_warning():
    a,b=bound(),bound(.0025)
    a["stationary"]["frequencies_cm1"]=[100,105,110,115,120,125]
    b["stationary"]["frequencies_cm1"]=[100,105,110,115,120,125]
    rotate_rows(b,0,5,.89)
    out=compare_bound(a,b,degeneracy_gap_cm1=6.)
    group=out["groups"][0]
    assert group["dimension"]==6
    assert group["left_frequency_range_cm1"]==[100,125]
    assert group["full_coordinate_space_group"]
    assert not group["directional_overlap_informative"]
    assert group["principal_cosines_squared"]==pytest.approx([1.]*6)


def test_resolved_sign_and_raw_sign_are_separate():
    a,b=bound(),bound(.0025)
    a["stationary"]["frequencies_cm1"][:2]=[-21.,-1.]
    b["stationary"]["frequencies_cm1"][:2]=[-19.,1.]
    out=compare_bound(a,b)
    assert out["raw_negative_frequency_counts"]==[2,1]
    assert out["resolved_negative_mode_counts"]==[1,0]
    assert out["changed_sign_spectral_ranks"]==[0]


@pytest.mark.parametrize("key,value",[("symbols",["C","H"]),("positions_angstrom",[[0,0,0],[1.01,0,0]]),
                                     ("masses_amu",[2.,12.]),("quantum_settings",{})])
def test_changed_binding_is_not_comparable(key,value):
    a,b=bound(),bound(.0025);b[key]=value
    out=compare_bound(a,b)
    assert out["status"]=="not_comparable"
    assert out["findings"]


@pytest.mark.parametrize("key,value",[("charge",1),("spin",3),("basis","other"),("xc","other"),
    ("scf_initial_guess","atom"),("grid_level",4),("dispersion",None),("density_fit",True),("conv_tol",1e-8)])
def test_method_and_state_input_changes_rejected(key,value):
    a,b=bound(),bound(.0025);b["quantum_settings"][key]=value
    assert compare_bound(a,b)["status"]=="not_comparable"


@pytest.mark.parametrize("key",["spin","charge","scf_initial_guess"])
def test_equal_missing_state_metadata_not_binding(key):
    a,b=bound(),bound(.0025)
    del a["quantum_settings"][key];del b["quantum_settings"][key]
    assert compare_bound(a,b)["status"]=="not_comparable"


def test_invalid_equal_metadata_rejected():
    a,b=bound(),bound(.0025)
    a["quantum_settings"]["spin"]=b["quantum_settings"]["spin"]=None
    assert compare_bound(a,b)["status"]=="not_comparable"


@pytest.mark.parametrize("key,value",[("grid_level",999),("conv_tol",4.),("scf_initial_guess","sap")])
def test_invalid_equal_electronic_controls_rejected(key,value):
    a,b=bound(),bound(.0025)
    a["quantum_settings"][key]=b["quantum_settings"][key]=value
    assert compare_bound(a,b)["status"]=="not_comparable"


def test_known_software_change_is_not_step_only_comparison():
    a,b=bound(),bound(.0025)
    b["software_versions"]={"fixture":"different version"}
    assert compare_bound(a,b)["status"]=="not_comparable"


def test_controls_beyond_step_cannot_change():
    a,b=bound(),bound(.0025)
    b["stationary"]["settings"]["frequency_tolerance_cm1"]=30.
    assert compare_bound(a,b)["status"]=="not_comparable"


def test_free_atom_order_cannot_change():
    a,b=bound(),bound(.0025)
    b["stationary"]["free_atom_indices"]=[1,0]
    b["stationary"]["free_masses_amu"]=[12.,1.]
    b["stationary"]["cartesian_modes_per_sqrt_amu"]=np.array(b["stationary"]["cartesian_modes_per_sqrt_amu"])[:,::-1].tolist()
    assert compare_bound(a,b)["status"]=="not_comparable"


def test_nonorthonormal_modes_rejected():
    a,b=bound(),bound(.0025)
    b["stationary"]["cartesian_modes_per_sqrt_amu"][0][0][0]=2.
    assert compare_bound(a,b)["status"]=="not_comparable"


@pytest.mark.parametrize("value", [True,"100",float("nan"),float("inf")])
def test_non_numeric_or_nonfinite_frequencies_rejected(value):
    a,b=bound(),bound(.0025)
    b["stationary"]["frequencies_cm1"][0]=value
    assert compare_bound(a,b)["status"]=="not_comparable"


def test_identical_steps_and_missing_versions_are_explicit():
    a,b=bound(),bound()
    b.pop("software_versions")
    out=compare_bound(a,b)
    assert out["status"]=="compared" and not out["distinct_steps"]
    assert out["legacy_or_missing_metadata"]


def test_inputs_remain_unchanged():
    a,b=bound(),bound(.0025);before=deepcopy((a,b))
    compare_bound(a,b)
    assert (a,b)==before


@pytest.fixture
def archive(tmp_path):
    source=ROOT/"data/validation/h2-integration"
    for name in ("modes","modes-half-step"):
        (tmp_path/name).mkdir()
        for filename in ("result.json","input.extxyz"):
            shutil.copyfile(source/name/filename,tmp_path/name/filename)
    return tmp_path


def test_existing_h2_comparison_without_new_calculations(archive):
    a=load_saved_result(archive/"modes/result.json")
    b=load_saved_result(archive/"modes-half-step/result.json")
    out=compare_bound(a,b)
    assert out["status"]=="compared"
    assert out["spectral_rank_frequency_changes_cm1"][-1]==pytest.approx(-.08958370070286037)
    assert out["resolved_negative_mode_counts"]==[0,0]
    assert out["legacy_or_missing_metadata"]
    assert "legacy absent" in str(out["legacy_or_missing_metadata"])


def edit_record(path, mutate):
    data=json.loads(path.read_text());mutate(data);path.write_text(json.dumps(data))


def test_structure_hash_binding(archive):
    path=archive/"modes/input.extxyz"
    path.write_text(path.read_text()+"\n")
    with pytest.raises(ComparisonError,match="hash"):
        load_saved_result(archive/"modes/result.json")


def test_snapshot_rounding_uses_exact_recorded_reference_for_pair_binding(archive):
    left_path=archive/"modes/result.json"; right_path=archive/"modes-half-step/result.json"
    a=load_saved_result(left_path)
    reference=np.array(a["positions_angstrom"])
    reference[1,2]+=1e-9
    for path in (left_path,right_path):
        edit_record(path,lambda d:d["stationary"].update(finite_difference_evidence={"reference_positions_angstrom":reference.tolist()}))
    a,b=load_saved_result(left_path),load_saved_result(right_path)
    assert a["geometry_binding"]["snapshot_relation"]=="reference_rounded_to_eight_decimal_angstrom"
    assert compare_bound(a,b)["status"]=="compared"
    reference[1,2]+=1e-9
    edit_record(right_path,lambda d:d["stationary"].update(finite_difference_evidence={"reference_positions_angstrom":reference.tolist()}))
    assert compare_bound(a,load_saved_result(right_path))["status"]=="not_comparable"


def test_snapshot_rounding_does_not_hide_larger_reference_change(archive):
    path=archive/"modes/result.json"
    reference=load_saved_result(path)["positions_angstrom"]
    reference[1][2]+=.001
    edit_record(path,lambda d:d["stationary"].update(finite_difference_evidence={"reference_positions_angstrom":reference}))
    with pytest.raises(ComparisonError,match="serialization"):
        load_saved_result(path)


@pytest.mark.parametrize("mutation",[
    lambda d:d.update(status="failed"),
    lambda d:d["quantum_settings"].update(scf_initial_guess=None),
    lambda d:d["quantum_diagnostics"]["settings"].update(spin=2),
    lambda d:d["design"].update(fixed_indices=[0]),
    lambda d:d.update(units={}),
    lambda d:d.update(schema_version=999),
    lambda d:d["quantum_diagnostics"].update(scf_initial_guess="atom"),
    lambda d:d["units"].update(force="Hartree/bohr"),
])
def test_loader_rejects_inconsistent_records(archive,mutation):
    path=archive/"modes/result.json";edit_record(path,mutation)
    with pytest.raises(ComparisonError):load_saved_result(path)


def test_cli_creates_report_and_refuses_overwrite(archive,tmp_path):
    output=tmp_path/"comparison.json"
    cmd=[sys.executable,str(HERE/"compare_modes.py"),str(archive/"modes/result.json"),str(archive/"modes-half-step/result.json"),"--output",str(output)]
    run=subprocess.run(cmd,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    before=output.read_bytes()
    assert json.loads(before)["quantum_jobs_started"]==0
    run=subprocess.run(cmd,capture_output=True,text=True)
    assert run.returncode!=0 and output.read_bytes()==before
