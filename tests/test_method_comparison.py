"""Paired-method arithmetic and durable identity/convergence rejection checks."""

from dataclasses import asdict, replace
import json
from pathlib import Path

from ase.units import Hartree
import pytest

from nanodesign import benchmark, method_comparison
from nanodesign.highlevel import CoupledClusterCalculationError


DFT_ENERGIES = {"methane": -10.0, "ethynyl_radical": -20.0, "methane_ethynyl_ts": -30.5,
                "acetylene": -21.0, "methyl_radical": -11.0}
CC_ENERGIES = {"CH4": -10.2, "C2H": -20.3, "C3H5": -30.0, "C2H2": -21.5, "CH3": -11.4}


@pytest.fixture
def paired_backends(monkeypatch):
    """Deterministic bookkeeping fixtures, never claimed chemical results."""
    calls = []
    def fake_dft(settings, directory, output):
        output = Path(output)
        output.mkdir()
        _, manifest = benchmark._snapshot_and_validate(directory or benchmark.REFERENCE_DIRECTORY, output)
        species = {}
        for name, spec in benchmark.SPECIES.items():
            state = replace(settings, charge=spec["charge"], spin=spec["spin"])
            species[name] = {"status": "completed", "file_sha256": manifest["file_sha256"][spec["file"]],
                             "charge": spec["charge"], "spin_2s": spec["spin"], "quantum_settings": asdict(state),
                             "energy_ev": DFT_ENERGIES[name], "energy_hartree": DFT_ENERGIES[name] / Hartree,
                             "scf_converged": True, "quantum_diagnostics": {"gradient_completed": True},
                             "s2": .75 if spec["spin"] else 0,
                             "expected_s2": .75 if spec["spin"] else 0, "spin_contamination": 0}
        return {"status": "completed", "input_hashes": manifest["file_sha256"], "species": species}
    def fake_cc(atoms, settings):
        calls.append((atoms.copy(), settings))
        energy = CC_ENERGIES[atoms.get_chemical_formula()]
        return {"status": "completed", "settings": asdict(settings), "total_energy_ev": energy,
                "total_energy_hartree": energy / Hartree, "scf_converged": True, "ccsd_converged": True,
                "triples_completed": True, "geometry": {"symbols": atoms.get_chemical_symbols(), "positions_angstrom": atoms.positions.tolist()},
                "hf_s2": .75 if settings.spin else 0, "hf_expected_s2": .75 if settings.spin else 0,
                "hf_s2_deviation": 0, "reference": "UHF" if settings.spin else "RHF",
                "variant": "UCCSD(T)" if settings.spin else "RCCSD(T)"}
    monkeypatch.setattr(method_comparison, "run_benchmark", fake_dft)
    monkeypatch.setattr(method_comparison, "coupled_cluster_energy", fake_cc)
    return fake_dft, fake_cc, calls


def test_paired_arithmetic_geometry_states_and_immutable_output(tmp_path, paired_backends):
    _, _, calls = paired_backends
    out = tmp_path / "paired"
    result = method_comparison.run_method_comparison(out)
    assert result["status"] == "completed"
    assert result["geometry_basis_pairing_verified"] is True
    assert result["electronic_state_identity_verified"] is False
    assert result["scientific_status"] == "electronic_state_identity_unverified"
    values = result["computed"]
    assert values["dft"]["nominal_ts_relative_energy_ev"] == pytest.approx(-.5)
    assert values["ccsd_t"]["nominal_ts_relative_energy_ev"] == pytest.approx(.5)
    assert values["dft_minus_ccsd_t"]["nominal_ts_relative_energy_ev"] == pytest.approx(-1)
    assert values["dft_minus_ccsd_t"]["reaction_energy_ev"] == pytest.approx(.4)
    assert values["dft_minus_ccsd_t"]["reaction_energy_kcal_per_mol"] == pytest.approx(.4 / benchmark.EV_PER_KCAL_PER_MOL)
    assert result["method_accuracy_validated"] is result["saddle_verified"] is result["tool_validated"] is False
    assert result["published_barrier_used_as_acceptance_target"] is False
    assert len(calls) == 5
    for name, record in result["species"].items():
        assert record["dft_xyz_sha256"] == record["cc_xyz_sha256"] == benchmark.sha256(out / record["geometry_file"])
        assert record["basis"] == "cc-pvdz"
        assert record["spin_2s"] == benchmark.SPECIES[name]["spin"]
        assert record["nominal_charge_spin_inputs_matched"] is True
        assert record["electronic_state_identity_verified"] is False
        cc_record = json.loads((out / record["cc_record"]).read_text())
        assert cc_record["status"] == "completed"
        assert cc_record["electronic_state_identity_verified"] is False
    assert (out / "dft/inputs/temelso_2006_si.txt").read_bytes() == (benchmark.REFERENCE_DIRECTORY / "temelso_2006_si.txt").read_bytes()
    assert json.loads((out / "method_comparison.json").read_text()) == result
    with pytest.raises(FileExistsError):
        method_comparison.run_method_comparison(out)


def test_large_basis_rejected_before_dft_or_output(tmp_path, paired_backends):
    out = tmp_path / "never-started"
    with pytest.raises(ValueError, match="only cc-pvdz"):
        method_comparison.run_method_comparison(out, basis="cc-pvtz")
    assert not out.exists()


def test_explicit_cc_guess_applies_to_every_species_and_is_archived(tmp_path, paired_backends):
    _, _, calls = paired_backends
    out = tmp_path / "explicit-guess"
    result = method_comparison.run_method_comparison(out, cc_initial_guess="atom")
    assert result["cc_settings"]["scf_initial_guess"] == "atom"
    assert result["ground_state_verified"] is result["initial_guess_scan_performed"] is False
    assert len(calls) == 5
    assert all(settings.scf_initial_guess == "atom" for _, settings in calls)
    for record in result["species"].values():
        assert record["cc_scf_initial_guess"] == "atom"
        cc_record = json.loads((out / record["cc_record"]).read_text())
        assert cc_record["scf_initial_guess"] == "atom"
        assert cc_record["quantum_diagnostics"]["settings"]["scf_initial_guess"] == "atom"


def test_invalid_cc_guess_rejected_before_dft_or_output(tmp_path, paired_backends):
    out = tmp_path / "bad-guess"
    with pytest.raises(ValueError, match="scf_initial_guess"):
        method_comparison.run_method_comparison(out, cc_initial_guess="typo")
    assert not out.exists()


@pytest.mark.parametrize("failure,status", [(RuntimeError, "failed"), (KeyboardInterrupt, "interrupted")])
def test_cc_failure_keeps_completed_species_and_failure_record(tmp_path, paired_backends, monkeypatch, failure, status):
    _, original_cc, _ = paired_backends
    out = tmp_path / "partial"
    def fail_second(atoms, settings):
        if atoms.get_chemical_formula() == "C2H":
            saved = json.loads((out / "method_comparison.json").read_text())
            assert saved["status"] == "running"
            assert set(saved["species"]) == {"methane"}
            assert json.loads((out / "ccsd_t/ethynyl_radical.json").read_text())["status"] == "running"
            raise failure("deliberate partial calculation")
        return original_cc(atoms, settings)
    monkeypatch.setattr(method_comparison, "coupled_cluster_energy", fail_second)
    with pytest.raises(failure, match="deliberate"):
        method_comparison.run_method_comparison(out)
    saved = json.loads((out / "method_comparison.json").read_text())
    assert saved["status"] == status
    assert saved["geometry_basis_pairing_verified"] is False
    assert "computed" not in saved
    assert set(saved["species"]) == {"methane"}
    failed = json.loads((out / "ccsd_t/ethynyl_radical.json").read_text())
    assert failed["status"] == status
    assert "energy_ev" not in failed


def test_quantum_failure_diagnostics_preserved(tmp_path, paired_backends, monkeypatch):
    def fail_cc(atoms, settings):
        raise CoupledClusterCalculationError("CCSD not converged", {"scf_converged": True, "ccsd_converged": False})
    monkeypatch.setattr(method_comparison, "coupled_cluster_energy", fail_cc)
    out = tmp_path / "failure-evidence"
    with pytest.raises(CoupledClusterCalculationError):
        method_comparison.run_method_comparison(out)
    record = json.loads((out / "ccsd_t/methane.json").read_text())
    assert record["quantum_diagnostics"]["ccsd_converged"] is False


@pytest.mark.parametrize("fault", ["settings", "geometry", "convergence", "units", "nonfinite", "nonfinite_spin", "wrong_spin", "reference", "extra_nonfinite"])
def test_bad_cc_result_cannot_produce_comparison(tmp_path, paired_backends, monkeypatch, fault):
    _, original_cc, _ = paired_backends
    def invalid_cc(atoms, settings):
        record = original_cc(atoms, settings)
        if fault == "settings":
            record["settings"]["basis"] = "sto-3g"
        elif fault == "geometry":
            record["geometry"]["positions_angstrom"][0][0] += .001
        elif fault == "convergence":
            record["ccsd_converged"] = False
        elif fault == "units":
            record["total_energy_ev"] = record["total_energy_hartree"]
        elif fault == "nonfinite_spin":
            record["hf_s2"] = float("nan")
        elif fault == "wrong_spin":
            record["hf_expected_s2"] = 7
        elif fault == "reference":
            record["reference"] = "ROHF"
        elif fault == "extra_nonfinite":
            record["additional_diagnostic"] = float("inf")
        else:
            record["total_energy_ev"] = float("nan")
        return record
    monkeypatch.setattr(method_comparison, "coupled_cluster_energy", invalid_cc)
    out = tmp_path / fault
    with pytest.raises(method_comparison.MethodComparisonError):
        method_comparison.run_method_comparison(out)
    saved = json.loads((out / "method_comparison.json").read_text())
    assert saved["status"] == "failed" and "computed" not in saved


@pytest.mark.parametrize("fault", ["xyz", "state", "basis", "dft_failed", "nonfinite_spin", "units", "convergence"])
def test_dft_pairing_mismatch_rejected_before_cc(tmp_path, paired_backends, monkeypatch, fault):
    original_dft, _, calls = paired_backends
    def invalid_dft(settings, directory, output):
        result = original_dft(settings, directory, output)
        if fault == "xyz":
            path = output / "inputs/methane.xyz"
            path.write_bytes(path.read_bytes() + b"\n")
        elif fault == "state":
            result["species"]["methane"]["spin_2s"] = 2
        elif fault == "basis":
            result["species"]["methane"]["quantum_settings"]["basis"] = "sto-3g"
        elif fault == "nonfinite_spin":
            result["species"]["methane"]["s2"] = float("nan")
        elif fault == "units":
            result["species"]["methane"]["energy_ev"] = result["species"]["methane"]["energy_hartree"]
        elif fault == "convergence":
            result["species"]["methane"]["scf_converged"] = False
        else:
            raise RuntimeError("DFT campaign failed")
        return result
    monkeypatch.setattr(method_comparison, "run_benchmark", invalid_dft)
    out = tmp_path / fault
    with pytest.raises(RuntimeError):
        method_comparison.run_method_comparison(out)
    assert calls == []
    assert json.loads((out / "method_comparison.json").read_text())["status"] == "failed"


def test_archived_geometry_change_during_cc_invalidates_comparison(tmp_path, paired_backends, monkeypatch):
    _, original_cc, _ = paired_backends
    out = tmp_path / "late-tamper"
    def tamper_after_last(atoms, settings):
        if atoms.get_chemical_formula() == "CH3":
            path = out / "dft/inputs/methane.xyz"
            path.write_bytes(path.read_bytes() + b"\n")
        return original_cc(atoms, settings)
    monkeypatch.setattr(method_comparison, "coupled_cluster_energy", tamper_after_last)
    with pytest.raises(method_comparison.MethodComparisonError, match="hash mismatch"):
        method_comparison.run_method_comparison(out)
    saved = json.loads((out / "method_comparison.json").read_text())
    assert saved["status"] == "failed" and "computed" not in saved
    assert len(saved["species"]) == 5


def test_backend_geometry_mutation_cannot_be_certified_as_paired(tmp_path, paired_backends, monkeypatch):
    _, original_cc, _ = paired_backends
    def mutated_cc(atoms, settings):
        atoms.positions[0, 0] += .1
        return original_cc(atoms, settings)
    monkeypatch.setattr(method_comparison, "coupled_cluster_energy", mutated_cc)
    out = tmp_path / "mutated-backend"
    with pytest.raises(method_comparison.MethodComparisonError, match="geometry mismatch"):
        method_comparison.run_method_comparison(out)
    assert json.loads((out / "method_comparison.json").read_text())["status"] == "failed"


def test_overflow_after_unit_conversion_is_rejected_before_saving_results():
    values = dict.fromkeys(DFT_ENERGIES, 0.0)
    values["methane_ethynyl_ts"] = 1e307
    with pytest.raises(method_comparison.MethodComparisonError, match="unit conversion"):
        method_comparison._energy_differences(values)


def test_equal_determinant_spin_squares_do_not_certify_electronic_state(tmp_path, paired_backends):
    result = method_comparison.run_method_comparison(tmp_path / "same-spin")
    for entry in result["species"].values():
        assert entry["dft_s2"] == entry["cc_hf_s2"] == entry["dft_expected_s2"]
        assert entry["electronic_state_identity_verified"] is False
    assert "does not verify the same electronic state" in result["interpretation"]
