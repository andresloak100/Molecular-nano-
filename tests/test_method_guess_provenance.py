"""DFT guess provenance checks using archived inputs and mocked calculations."""

import copy
from dataclasses import asdict, replace
import json
from pathlib import Path
import shutil

from ase.units import Hartree
import pytest

from nanodesign import method_comparison
from nanodesign.benchmark import SPECIES


ARCHIVES = Path(__file__).resolve().parents[1] / "data" / "validation"
SOURCES = ("benchmark", "species", "diagnostic_settings", "diagnostic")


def guess_sources(dft, name="methane"):
    record = dft["species"][name]
    diagnostics = record["quantum_diagnostics"]
    return {
        "benchmark": dft["quantum_settings"],
        "species": record["quantum_settings"],
        "diagnostic_settings": diagnostics["settings"],
        "diagnostic": diagnostics,
    }


@pytest.fixture
def archived_dft():
    root = ARCHIVES / "paired-ccpvdz-atom"
    return root, json.loads((root / "dft" / "benchmark.json").read_text())


@pytest.mark.parametrize("archive", ["paired-ccpvdz", "paired-ccpvdz-atom"])
def test_actual_legacy_evidence_accepts_every_species_without_mutation(archive):
    root = ARCHIVES / archive
    files = [root / "dft" / "benchmark.json"]
    files += [root / "dft" / "inputs" / spec["file"] for spec in SPECIES.values()]
    original_bytes = {path: path.read_bytes() for path in files}
    dft = json.loads(original_bytes[files[0]])
    original = copy.deepcopy(dft)
    for name, spec in SPECIES.items():
        assert all("scf_initial_guess" not in source for source in guess_sources(dft, name).values())
        atoms, digest = method_comparison._paired_input(dft, root, name, "cc-pvdz")
        assert atoms.get_chemical_formula() == spec["formula"]
        assert digest == dft["input_hashes"][spec["file"]]
    assert dft == original
    assert {path: path.read_bytes() for path in files} == original_bytes


@pytest.mark.parametrize("mask", range(16))
def test_missing_and_explicit_minao_are_compatible_without_mutation(archived_dft, mask):
    root, dft = archived_dft
    # Benchmark spin is intentionally overridden for methane, but not ethynyl.
    for name in ("methane", "ethynyl_radical"):
        for index, source in enumerate(guess_sources(dft, name).values()):
            if mask & (1 << index):
                source["scf_initial_guess"] = "minao"
        before = copy.deepcopy(dft)
        method_comparison._paired_input(dft, root, name, "cc-pvdz")
        assert dft == before


@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("value", ["atom", "1e", "huckel", None, False, 1, [], {}, "", "invented", "MINAO", " minao "])
def test_each_conflicting_or_invalid_guess_is_rejected(archived_dft, source, value):
    root, dft = archived_dft
    guess_sources(dft)[source]["scf_initial_guess"] = value
    with pytest.raises(method_comparison.MethodComparisonError, match="scf_initial_guess"):
        method_comparison._paired_input(dft, root, "methane", "cc-pvdz")


def test_agreeing_metadata_cannot_override_requested_dft_guess(archived_dft):
    root, dft = archived_dft
    for source in guess_sources(dft).values():
        source["scf_initial_guess"] = "atom"
    with pytest.raises(method_comparison.MethodComparisonError, match="expected 'minao'"):
        method_comparison._paired_input(dft, root, "methane", "cc-pvdz")


@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("value", [None, []])
def test_explicit_malformed_containers_are_not_legacy_absence(archived_dft, source, value):
    root, dft = archived_dft
    record = dft["species"]["methane"]
    if source == "benchmark":
        dft["quantum_settings"] = value
    elif source == "species":
        record["quantum_settings"] = value
    elif source == "diagnostic_settings":
        record["quantum_diagnostics"]["settings"] = value
    else:
        record["quantum_diagnostics"] = value
    with pytest.raises(method_comparison.MethodComparisonError, match="Invalid DFT"):
        method_comparison._paired_input(dft, root, "methane", "cc-pvdz")


def test_absent_optional_settings_containers_keep_existing_fixture_compatibility(archived_dft):
    root, dft = archived_dft
    del dft["quantum_settings"]
    del dft["species"]["methane"]["quantum_diagnostics"]["settings"]
    method_comparison._paired_input(dft, root, "methane", "cc-pvdz")


def test_unrelated_missing_method_field_is_not_filled_from_current_defaults(archived_dft):
    root, dft = archived_dft
    del dft["species"]["methane"]["quantum_settings"]["xc"]
    with pytest.raises(method_comparison.MethodComparisonError, match="Unexpected DFT method"):
        method_comparison._paired_input(dft, root, "methane", "cc-pvdz")


@pytest.fixture
def mock_pairing(monkeypatch, archived_dft):
    """Only bookkeeping runs; synthetic CC energies are not chemical evidence."""
    root, legacy = archived_dft
    control = {"calls": [], "mutate": lambda dft: None, "legacy": False}

    def fake_dft(settings, directory, output):
        inputs = Path(output) / "inputs"
        inputs.mkdir(parents=True)
        for spec in SPECIES.values():
            shutil.copyfile(root / "dft" / "inputs" / spec["file"], inputs / spec["file"])
        dft = copy.deepcopy(legacy)
        if not control["legacy"]:
            dft["quantum_settings"] = asdict(settings)
            for name, spec in SPECIES.items():
                record = dft["species"][name]
                record["quantum_settings"] = asdict(replace(settings, charge=spec["charge"], spin=spec["spin"]))
                record["quantum_diagnostics"]["settings"] = copy.deepcopy(record["quantum_settings"])
                record["quantum_diagnostics"]["scf_initial_guess"] = settings.scf_initial_guess
        control["mutate"](dft)
        control["dft"] = dft
        control["original"] = copy.deepcopy(dft)
        return dft

    def fake_cc(atoms, settings):
        control["calls"].append(settings)
        if control.get("late_conflict") and len(control["calls"]) == len(SPECIES):
            guess_sources(control["dft"])["diagnostic"]["scf_initial_guess"] = "atom"
        return {
            "status": "completed", "settings": asdict(settings),
            "geometry": {"symbols": atoms.get_chemical_symbols(), "positions_angstrom": atoms.positions.tolist()},
            "total_energy_ev": -1.0, "total_energy_hartree": -1.0 / Hartree,
            "scf_converged": True, "ccsd_converged": True, "triples_completed": True,
            "hf_s2": .75 if settings.spin else 0, "hf_expected_s2": .75 if settings.spin else 0,
            "hf_s2_deviation": 0, "reference": "UHF" if settings.spin else "RHF",
            "variant": "UCCSD(T)" if settings.spin else "RCCSD(T)",
        }

    monkeypatch.setattr(method_comparison, "run_benchmark", fake_dft)
    monkeypatch.setattr(method_comparison, "coupled_cluster_energy", fake_cc)
    return control


@pytest.mark.parametrize("source", SOURCES)
def test_conflicting_first_species_preserves_failure_and_never_starts_cc(tmp_path, mock_pairing, source):
    def conflict(dft):
        guess_sources(dft)[source]["scf_initial_guess"] = "atom"
    mock_pairing["mutate"] = conflict
    output = tmp_path / "conflict"
    with pytest.raises(method_comparison.MethodComparisonError, match="scf_initial_guess"):
        method_comparison.run_method_comparison(output)
    assert mock_pairing["calls"] == []
    summary = json.loads((output / "method_comparison.json").read_text())
    failed = json.loads((output / "ccsd_t" / "methane.json").read_text())
    assert summary["status"] == failed["status"] == "failed"
    assert summary["geometry_basis_pairing_verified"] is False
    assert summary["electronic_state_identity_verified"] is False
    assert "computed" not in summary and "energy_ev" not in failed


@pytest.mark.parametrize("legacy", [False, True])
def test_cc_atom_does_not_change_dft_minao_or_certify_states(tmp_path, mock_pairing, legacy):
    mock_pairing["legacy"] = legacy
    result = method_comparison.run_method_comparison(tmp_path / "paired", cc_initial_guess="atom")
    assert result["status"] == "completed"
    assert result["dft_settings"]["scf_initial_guess"] == "minao"
    assert result["cc_settings"]["scf_initial_guess"] == "atom"
    assert len(mock_pairing["calls"]) == len(SPECIES)
    assert all(settings.scf_initial_guess == "atom" for settings in mock_pairing["calls"])
    assert mock_pairing["dft"] == mock_pairing["original"]
    for flag in ("electronic_state_identity_verified", "ground_state_verified", "initial_guess_scan_performed", "tool_validated"):
        assert result[flag] is False
    assert all(record["electronic_state_identity_verified"] is False for record in result["species"].values())


def test_late_provenance_conflict_blocks_completed_comparison(tmp_path, mock_pairing):
    mock_pairing["late_conflict"] = True
    output = tmp_path / "late-conflict"
    with pytest.raises(method_comparison.MethodComparisonError, match="scf_initial_guess"):
        method_comparison.run_method_comparison(output)
    summary = json.loads((output / "method_comparison.json").read_text())
    assert len(mock_pairing["calls"]) == len(summary["species"]) == len(SPECIES)
    assert summary["status"] == "failed" and "computed" not in summary
    assert summary["geometry_basis_pairing_verified"] is False
