"""Cheap synthetic record tests. Fixtures are NOT molecular calculation results."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).parent))
from audit import audit_record, inspect_file


def characterization(force=0.0, negative=1):
    # Synthetic two-free-atom modes with deliberately unequal masses.
    masses = np.array([1., 12.])
    modes = (np.eye(6) / np.repeat(np.sqrt(masses), 3)).reshape(6, 2, 3)
    return {
        "settings": {"force_tolerance_ev_per_angstrom": .03, "frequency_tolerance_cm1": 20.},
        "free_atom_indices": [0, 1], "free_masses_amu": masses.tolist(),
        "geometry_angstrom": [[0, 0, 0], [1, 0, 0]],
        "frequencies_cm1": [-100.] * negative + [100.] * (6-negative),
        "cartesian_modes_per_sqrt_amu": modes.tolist(),
        "free_gradient_ev_per_angstrom": [[force, 0, 0], [0, 0, 0]],
        "stationary_within_force_tolerance": bool(force <= .03),
        "negative_mode_count": negative,
    }


def record():
    # Deliberately artificial atom counts/energies, not chemical species models.
    return {"methane": {"status": "converged", "energy_ev": -2., "chemical_symbols": ["H"]},
            "ethynyl": {"status": "converged", "energy_ev": -3., "chemical_symbols": ["C"]},
            "transition_state": {"status": "converged", "energy_ev": -4.,
                                 "chemical_symbols": ["H", "C"],
                                 "final_geometry_angstrom": [[0, 0, 0], [1, 0, 0]]},
            "transition_state_characterization": characterization(),
            "barrier": {"value_kcal_per_mol": 23.060547830618307}}


def codes(result):
    return {finding["code"] for finding in result["findings"]}


def test_pending_preserved_without_fake_values():
    out = audit_record({})
    assert out["component_status"]["transition_state"] == "not_recorded"
    assert out["recomputed_electronic_energy_difference_kcal_per_mol"] is None
    assert not out["geometric_first_order_saddle_candidate"]


@pytest.mark.parametrize("force,negative,expected", [(0,1,True),(.05,1,False),(0,0,False),(0,2,False)])
def test_geometric_candidate_requires_stationarity_and_one_mode(force, negative, expected):
    data = record()
    data["transition_state_characterization"] = characterization(force, negative)
    out = audit_record(data)
    assert out["geometric_first_order_saddle_candidate"] is expected
    assert out["transition_state_verified"] is False


def test_false_stationarity_flag_does_not_override_gradient():
    data = record()
    data["transition_state_characterization"] = characterization(.1)
    data["transition_state_characterization"]["stationary_within_force_tolerance"] = True
    out = audit_record(data)
    assert not out["geometric_first_order_saddle_candidate"]
    assert "stationarity_flag_mismatch" in codes(out)


def test_incomplete_guesses_do_not_report_complete_agreement():
    out = audit_record({"ethynyl_initial_guess_scan": {"energies_hartree": {"minao": -1., "atom": -1., "1e": None, "huckel": None}, "guess_independent": True}})
    scan = out["electronic_guess_evidence"]["ethynyl"]
    assert scan["agreement_over_converged_subset"]
    assert not scan["complete_recorded_scan"]
    assert "ethynyl_unsupported_guess_independence" in codes(out)


def test_one_guess_does_not_establish_agreement():
    scan = audit_record({"ethynyl_initial_guess_scan": {"energies_hartree": {"minao": -1.}}})["electronic_guess_evidence"]["ethynyl"]
    assert scan["spread_kcal_per_mol"] is None
    assert not scan["agreement_over_converged_subset"]


def test_distinct_solutions_remain_visible():
    result = audit_record({"ethynyl_initial_guess_scan": {"energies_hartree": {"minao": -1., "atom": -1.01}}})
    assert "ethynyl_electronic_solution_ambiguity" in codes(result)


def test_mismatched_guess_manifest_does_not_hide_recorded_energy():
    result = audit_record({"ethynyl_initial_guess_scan": {
        "attempted_guesses": ["minao", "atom"],
        "energies_hartree": {"minao": -1., "atom": -1., "huckel": -.5}}})
    assert {"ethynyl_electronic_solution_ambiguity", "ethynyl_scan_manifest_mismatch"} <= codes(result)


@pytest.mark.parametrize("field,value", [("free_atom_indices", [98, 99]),
                                         ("geometry_angstrom", [[0, 0, 0], [500, 0, 0]])])
def test_characterization_must_belong_to_optimized_geometry(field, value):
    data = record()
    data["transition_state_characterization"][field] = value
    result = audit_record(data)
    assert "invalid_characterization" in codes(result)
    assert not result["geometric_first_order_saddle_candidate"]


@pytest.mark.parametrize("field,value", [("free_masses_amu",[1,-12]), ("frequencies_cm1",[-100]), ("free_atom_indices",[0,0])])
def test_malformed_characterization_fails_closed(field, value):
    data = record()
    data["transition_state_characterization"][field] = value
    out = audit_record(data)
    assert "invalid_characterization" in codes(out)
    assert not out["geometric_first_order_saddle_candidate"]


def test_cartesian_normalization_is_not_mass_normalization():
    data = record()
    data["transition_state_characterization"]["cartesian_modes_per_sqrt_amu"] = np.eye(6).reshape(6,2,3).tolist()
    out = audit_record(data)
    assert "mode_mass_normalization_mismatch" in codes(out)
    assert not out["geometric_first_order_saddle_candidate"]


def test_failed_optimizer_cannot_produce_candidate():
    data = record()
    data["transition_state"]["status"] = "failed"
    out = audit_record(data)
    assert not out["geometric_first_order_saddle_candidate"]
    assert out["recomputed_electronic_energy_difference_kcal_per_mol"] is None


def test_barrier_arithmetic_and_false_verification():
    data = record()
    data["barrier"] = {"value_kcal_per_mol": 999., "transition_state_verified": True}
    result = audit_record(data)
    assert {"energy_arithmetic_mismatch", "unsupported_transition_state_verification"} <= codes(result)


def test_interrupted_json_preserved(tmp_path):
    path = tmp_path / "partial.json"
    path.write_text('{"transition_state":')
    raw = path.read_bytes()
    result = inspect_file(path)
    assert "unreadable_record" in codes(result)
    assert path.read_bytes() == raw


def test_nan_json_rejected(tmp_path):
    path = tmp_path / "nan.json"
    path.write_text('{"energy": NaN}')
    assert "unreadable_record" in codes(inspect_file(path))


def test_valid_file_is_hashed_without_modification(tmp_path):
    path = tmp_path / "valid.json"
    path.write_text(json.dumps(record()))
    raw = path.read_bytes()
    result = inspect_file(path)
    assert len(result["source"]["sha256"]) == 64
    assert path.read_bytes() == raw


def test_saved_transfer_claim_gets_independent_projection_check():
    data = record()
    positions = [[1.,0,0],[0.,0,0],[-1.,0,0]]
    data["transition_state"].update(chemical_symbols=["C","H","C"], final_geometry_angstrom=positions)
    data["ethynyl"]["chemical_symbols"] = ["C","C"]
    masses = np.array([12.,1.,12.])
    modes = np.eye(9)[[4,0,1,2,3,5,6,7,8]] / np.repeat(np.sqrt(masses),3)
    modes[0,3] = 1e-14
    char = characterization()
    char.update(free_atom_indices=[0,1,2], free_masses_amu=masses.tolist(),
                geometry_angstrom=positions, frequencies_cm1=[-100.]+[100.]*8,
                free_gradient_ev_per_angstrom=[[0,0,0]]*3,
                cartesian_modes_per_sqrt_amu=modes.reshape(9,3,3).tolist())
    data["transition_state_characterization"] = char
    data["reaction_indices"] = {"acceptor_carbon":0, "transferring_hydrogen":1, "donor_carbon":2}
    data["transition_mode_analysis"] = {"single_transfer_like_mode":True}
    result = audit_record(data)
    assert result["geometric_first_order_saddle_candidate"]
    assert "transfer_heuristic_disagreement" in codes(result)
    assert not result["independent_transfer_evidence"]["single_negative_mode_passes_overlap_heuristic"]
    assert not result["transition_state_verified"]
