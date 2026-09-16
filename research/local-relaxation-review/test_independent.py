"""Independent L2 analytic and actual-producer tests, all synthetic."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).parent))
from producer_fixture import ROOT, load_module, produce
l1 = load_module('l2_reviewed_diagnostic', ROOT / 'research/local-relaxation-diagnostic/diagnose.py')
OPTIONS = dict(curvature_floor_ev_per_angstrom2=1e-9, max_condition_number=1e6,
               max_atom_shift_angstrom=10., expected_index=0)

def run(H, g, **kw):
    return l1.diagnose_quadratic(np.asarray(H).tolist(), np.asarray(g).reshape(-1, 3).tolist(), **(OPTIONS | kw))

def estimate_check(result, H, g):
    assert result['status'] == 'estimated', result
    expected = np.linalg.solve(H, -g)
    e = result['estimate']
    np.testing.assert_allclose(np.array(e['displacement_free_atoms_angstrom']).ravel(), expected, atol=1e-13)
    assert e['signed_quadratic_energy_change_ev'] == pytest.approx(g @ expected + .5 * expected @ H @ expected, abs=1e-13)
    assert not result['energy_error_bound_established']
    assert not result['actual_relaxation_performed']
    assert not result['electronic_state_verified']
    return e

def test_coupled_cartesian_solve_sign_and_immutability():
    H = [[2., .5, 0.], [.5, 3., .2], [0., .2, 4.]]
    g = [[.2, -.3, .1]]
    before = deepcopy((H, g))
    a = l1.diagnose_quadratic(H, g, **OPTIONS)
    ea = estimate_check(a, np.array(H), np.array(g).ravel())
    b = run(H, -np.array(g))
    np.testing.assert_allclose(b['estimate']['displacement_free_atoms_angstrom'], -np.array(ea['displacement_free_atoms_angstrom']))
    assert b['estimate']['signed_quadratic_energy_change_ev'] == pytest.approx(ea['signed_quadratic_energy_change_ev'])
    assert (H, g) == before

@pytest.mark.parametrize('g,expected', [([1., 0., 0.], .5), ([0., 1., 0.], -.5), ([1., 1., 0.], 0.)])
def test_saddle_signed_change_and_cancellation(g, expected):
    result = run(np.diag([-1., 1., 2.]), g, expected_index=1)
    e = estimate_check(result, np.diag([-1., 1., 2.]), np.array(g))
    assert e['signed_quadratic_energy_change_ev'] == pytest.approx(expected)
    assert e['positive_curvature_energy_change_ev'] <= 0
    assert e['negative_curvature_energy_change_ev'] >= 0
    if expected == 0:
        assert e['absolute_mode_contribution_sum_ev'] == pytest.approx(1.)
        assert e['max_atom_shift_angstrom'] == pytest.approx(2**.5)

@pytest.mark.parametrize('diagonal,options,code', [
    ([0., 1., 2.], {}, 'soft_curvature'),
    ([1e-9, 1., 2.], {}, 'soft_curvature'),
    ([1e-7, 1., 2.], {}, 'ill_conditioned'),
    ([-1e-10, 1., 2.], {}, 'soft_curvature'),
])
def test_unforced_soft_direction_still_refused(diagonal, options, code):
    r = run(np.diag(diagonal), [0., .01, 0.], **options)
    assert r['status'] == 'refused' and r['estimate'] is None
    assert code in [f['code'] for f in r['findings']]

def test_small_force_not_energy_or_displacement_bound():
    r = run(np.diag([1e-4, 1., 1.]), [.029, 0., 0.], max_atom_shift_angstrom=.1)
    assert r['status'] == 'inconclusive'
    assert r['estimate']['max_atom_shift_angstrom'] == pytest.approx(290.)
    assert r['estimate']['signed_quadratic_energy_change_ev'] == pytest.approx(-4.205)
    assert not r['energy_error_bound_established']

def test_index_mismatch_retains_signed_diagnostic():
    r = run(np.diag([-1., -2., 3.]), [.1, .2, .3], expected_index=1)
    assert r['status'] == 'inconclusive' and r['estimate'] is not None
    assert r['quadratic_character'] == 'higher_order_saddle'

@pytest.mark.parametrize('key,value', [('curvature_floor_ev_per_angstrom2',True),
    ('max_condition_number',.9), ('max_atom_shift_angstrom',float('inf')),
    ('expected_index',True), ('expected_index',1.0)])
def test_invalid_controls(key,value):
    with pytest.raises(l1.DiagnosticError):
        run(np.eye(3), [0., 0., 0.], **{key:value})

def test_asymmetric_input_rejected():
    with pytest.raises(l1.DiagnosticError):
        run([[1., .1, 0.], [0., 2., 0.], [0., 0., 3.]], [0., 0., 0.])

def test_rotation_invariance():
    Q, _ = np.linalg.qr(np.array([[1., 2., 3.], [2., -1., 1.], [1., 0., -2.]]))
    H = np.diag([1., 1., 3.]); g = np.array([.1, .2, .3])
    rotated = Q @ H @ Q.T; rotated = (rotated + rotated.T) / 2
    a = run(H,g)['estimate']; b = run(rotated,Q@g)['estimate']
    np.testing.assert_allclose(np.array(b['displacement_free_atoms_angstrom']).ravel(), Q @ np.array(a['displacement_free_atoms_angstrom']).ravel())
    assert b['signed_quadratic_energy_change_ev'] == pytest.approx(a['signed_quadratic_energy_change_ev'])

@pytest.fixture
def actual(tmp_path,monkeypatch):
    return produce(tmp_path,monkeypatch)

def bound(path):
    return l1.diagnose_saved_result(path, max_relative_hessian_asymmetry=1e-6, **OPTIONS)

def test_actual_nonzero_producer_anchor_mapping_exact_reference_and_immutability(actual):
    path,H,g=actual
    before={p:p.read_bytes() for p in path.parent.rglob('*') if p.is_file()}
    r=bound(path); e=estimate_check(r,H,g)
    assert r['free_atom_indices']==[0,2] and r['frozen_atom_indices']==[1]
    np.testing.assert_array_equal(e['displacement_all_atoms_angstrom'][1],[0.,0.,0.])
    np.testing.assert_allclose(np.array(e['displacement_all_atoms_angstrom'])[[0,2]], e['displacement_free_atoms_angstrom'])
    assert r['geometry_binding']['full_precision_force_reference_available']
    assert r['reference_positions_angstrom'][0][0] == .123456789123
    assert r['saved_force_reconstruction']['status']=='passed'
    assert all(p.read_bytes()==data for p,data in before.items())

def test_changed_masses_leave_cartesian_estimate_unchanged(tmp_path,monkeypatch):
    estimates=[]
    for i,masses in enumerate([(1.5,20.,3.5),(15.,2.,.35)]):
        folder=tmp_path/str(i);folder.mkdir()
        path,H,g=produce(folder,monkeypatch,masses)
        estimates.append(estimate_check(bound(path),H,g))
    np.testing.assert_allclose(estimates[0]['displacement_free_atoms_angstrom'],estimates[1]['displacement_free_atoms_angstrom'],atol=1e-14)

@pytest.mark.parametrize('tamper',['snapshot','force','reference','settings','legacy'])
def test_actual_binding_and_reconstruction_refuse_tampering(actual,tamper):
    path,H,g=actual
    data=json.loads(path.read_text())
    if tamper=='snapshot':
        snapshot=path.with_name('input.extxyz');snapshot.write_bytes(snapshot.read_bytes()+b'\n')
    elif tamper=='force':
        data['stationary']['finite_difference_evidence']['baseline_forces_ev_per_angstrom'][0][0]+=.1
    elif tamper=='reference':
        data['stationary']['finite_difference_evidence']['reference_positions_angstrom'][0][0]+=.1
    elif tamper=='settings':
        data['quantum_settings']['spin']=3
    else:
        del data['stationary']['finite_difference_evidence']
    if tamper!='snapshot': path.write_text(json.dumps(data))
    before=path.read_bytes(); r=bound(path)
    assert r['status'] in ('refused','unavailable') and r['estimate'] is None, r
    assert path.read_bytes()==before

def test_actual_raw_asymmetry_gate_precedes_symmetric_estimate(tmp_path,monkeypatch):
    H=np.diag([2.,3.,4.,5.,6.,7.]);H[0,3]=.2;H[3,0]=-.2
    path,_,_=produce(tmp_path,monkeypatch,hessian=H)
    r=bound(path)
    assert r['saved_force_reconstruction']['status']=='passed'
    assert r['status']=='refused' and r['estimate'] is None
    assert 'hessian_asymmetry_limit_exceeded' in [f['code'] for f in r['findings']]
