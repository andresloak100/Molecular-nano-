"""Tiny actual-PySCF capture roundtrips intended for an uncongested CI runner.

Each marked case performs one SCF and one analytical gradient. H2 and H with
STO-3G have two and one AO functions respectively. Serialization, comparison,
and cache checks add no SCF or gradient evaluations. These are integration
checks of saved numerical evidence, not tests of chemical design accuracy.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from uuid import UUID

from ase import Atoms
from ase.units import Hartree
import numpy as np
import pytest

from nanodesign.electronic_state import (
    compare_snapshots, read_snapshot, save_snapshot, validate_snapshot,
)
from nanodesign.quantum import PySCFCalculator, QuantumSettings


@pytest.mark.quantum
@pytest.mark.parametrize("density_fit", [False, True], ids=["direct", "density-fit"])
@pytest.mark.parametrize("reference", ["RKS", "UKS"], ids=["h2-rks", "h-uks"])
def test_real_capture_roundtrip_preserves_solver_evidence_and_accepted_lifecycle(
    tmp_path, reference, density_fit,
):
    # Importing this test module does not import or run the electronic solver.
    from pyscf import lib

    # Real paths also work on systems whose temporary-directory alias is a
    # symlink, since the snapshot writer deliberately rejects symlink parents.
    directory = tmp_path.resolve()
    state_directory = directory / "states"
    event_log = directory / "electronic.jsonl"
    if reference == "RKS":
        atoms = Atoms("H2", positions=[[-0.18, 0.22, -0.40], [-0.18, 0.22, 0.39]])
        spin, expected_ao, expected_alpha, expected_beta = 0, 2, 1, 1
    else:
        # This also exercises a genuine zero-electron beta channel.
        atoms = Atoms("H", positions=[[0.27, -0.43, 0.61]])
        spin, expected_ao, expected_alpha, expected_beta = 1, 1, 1, 0
    input_positions = atoms.positions.copy()
    settings = QuantumSettings(
        charge=0, spin=spin, xc="pbe0", basis="sto-3g", dispersion=None,
        grid_level=1, conv_tol=1e-10, max_cycle=50, threads=1,
        memory_mb=256, density_fit=density_fit, scf_initial_guess="minao",
    )
    previous_threads = lib.num_threads()
    atoms.calc = calculator = PySCFCalculator(
        settings, event_log=event_log, electronic_state_directory=state_directory,
    )
    forces = atoms.get_forces()
    energy_ev = atoms.get_potential_energy()  # ASE cache, no second solve.
    diagnostics = calculator.diagnostics
    assert lib.num_threads() == previous_threads
    assert np.isfinite(energy_ev)
    assert forces.shape == (len(atoms), 3) and np.isfinite(forces).all()
    np.testing.assert_array_equal(atoms.positions, input_positions)
    assert diagnostics["calculation_status"] == "completed"
    assert diagnostics["scf_converged"] is diagnostics["gradient_completed"] is True
    assert diagnostics["reference"] == reference
    assert diagnostics["auxiliary_basis_response"] is density_fit
    assert diagnostics["requested_pyscf_threads"] == 1
    assert diagnostics["basis_functions"] == expected_ao
    assert diagnostics["electron_count"] == expected_alpha + expected_beta

    saved = diagnostics["electronic_state_snapshot"]
    UUID(saved["call_id"])
    assert saved["call_id"] == diagnostics["call_id"]
    assert saved["phase"] == "converged_scf_only"
    assert saved["whole_evaluation_accepted"] is True
    snapshot_path = Path(saved["path"])
    assert snapshot_path == state_directory / f"{diagnostics['call_id']}.json"
    raw = snapshot_path.read_bytes()
    assert saved["sha256"] == hashlib.sha256(raw).hexdigest()
    assert saved["size_bytes"] == len(raw)
    snapshot = read_snapshot(snapshot_path)
    assert snapshot["call_id"] == diagnostics["call_id"]
    assert snapshot["reference"] == reference
    assert snapshot["phase"] == "converged_scf_only"
    assert snapshot["settings"] == settings.to_dict()
    assert snapshot["scf_converged"] is True
    assert snapshot["scf_energy_hartree"] == pytest.approx(diagnostics["dft_energy_hartree"], rel=0, abs=1e-12)
    assert snapshot["scf_energy_hartree"] * Hartree == pytest.approx(energy_ev, rel=0, abs=1e-10)
    assert snapshot["producer"]["engine"] == "PySCF"
    assert snapshot["electronic_state_identity_verified"] is False
    assert snapshot["ground_state_verified"] is False
    assert snapshot["stability_assessed_by_capture"] is False

    geometry = snapshot["geometry"]
    assert geometry["symbols"] == atoms.get_chemical_symbols()
    assert geometry["atomic_numbers"] == atoms.numbers.tolist()
    assert geometry["charge"] == 0 and geometry["spin"] == spin
    assert geometry["pbc"] == [False, False, False]
    assert geometry["bohr_angstrom"] == float(lib.param.BOHR)
    # Use the solver's recorded conversion, not a separately rounded ASE Bohr.
    np.testing.assert_allclose(np.asarray(geometry["positions_bohr"]) * geometry["bohr_angstrom"],
                               input_positions, rtol=0, atol=2e-14)
    assert snapshot["basis"]["n_ao"] == expected_ao
    assert len(snapshot["basis"]["ao_labels"]) == expected_ao
    assert snapshot["basis"]["shells"] and snapshot["basis"]["expanded_basis"]
    assert np.asarray(snapshot["ao_overlap"]).shape == (expected_ao, expected_ao)
    arithmetic = validate_snapshot(snapshot)
    for channel, expected in (("alpha", expected_alpha), ("beta", expected_beta)):
        evidence = snapshot["channels"][channel]
        assert evidence["n_electrons"] == expected
        assert len(evidence["occupied_indices"]) == expected
        assert np.asarray(evidence["occupied_coefficients"]).shape == (expected_ao, expected)
        assert arithmetic[channel]["density_electron_trace"] == pytest.approx(expected, rel=0, abs=1e-8)
        assert arithmetic[channel]["metric_orthonormality_max_abs_residual"] <= 1e-8

    # Exercise the real serializer again with already-read evidence only.
    roundtrip_path = directory / "roundtrip.json"
    roundtrip_reference = save_snapshot(roundtrip_path, snapshot)
    assert roundtrip_path.read_bytes() == raw
    assert roundtrip_reference["sha256"] == saved["sha256"]
    compared = compare_snapshots(snapshot, read_snapshot(roundtrip_path))
    assert compared["comparison_status"] == "compared" and compared["reasons"] == []
    assert compared["call_ids"] == {"left": diagnostics["call_id"], "right": diagnostics["call_id"]}
    assert compared["electronic_state_identity_verified"] is False
    assert compared["ground_state_verified"] is False
    for channel, expected in (("alpha", expected_alpha), ("beta", expected_beta)):
        metric = compared["metrics"][channel]
        assert metric["n_occupied"] == expected
        np.testing.assert_allclose(metric["singular_values"], np.ones(expected), rtol=0, atol=1e-8)
        np.testing.assert_allclose(metric["principal_angles_radians"], np.zeros(expected), rtol=0, atol=2e-6)
        assert metric["density_distance_squared"] == pytest.approx(0.0, abs=1e-12)

    events_before = event_log.read_bytes()
    events = [json.loads(line) for line in events_before.splitlines()]
    assert {record["call_id"] for record in events} == {diagnostics["call_id"]}
    assert sum(record["event"] == "scf_completed" for record in events) == 1
    assert sum(record["event"] == "calculation_completed" for record in events) == 1
    assert events[-1]["event"] == "calculation_completed"
    assert events[-1]["scf_converged"] is events[-1]["gradient_completed"] is True
    assert list(state_directory.iterdir()) == [snapshot_path]
    np.testing.assert_array_equal(atoms.get_forces(), forces)
    assert atoms.get_potential_energy() == energy_ev
    assert event_log.read_bytes() == events_before
    assert snapshot_path.read_bytes() == raw
    assert list(state_directory.iterdir()) == [snapshot_path]
