"""Comparison arithmetic, published identity integrity, and durable failure evidence."""

from dataclasses import replace
import json
from pathlib import Path
import shutil

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree
import numpy as np
import pytest

from nanodesign import benchmark
from nanodesign.quantum import QuantumSettings


class KnownEnergies(Calculator):
    """Analytic bookkeeping fixture only; never a claimed chemical result."""

    implemented_properties = ["energy", "forces"]
    # Deliberately negative TS-relative energy: the comparison must not clamp it.
    energies = {"CH4": -10.0, "C2H": -20.0, "C3H5": -30.5, "C2H2": -21.0, "CH3": -11.0}
    states = {"CH4": 0, "C2H": 1, "C3H5": 1, "C2H2": 0, "CH3": 1}
    observed = []

    def __init__(self, settings, *, event_log):
        super().__init__()
        self.settings = settings
        self.event_log = Path(event_log)
        self.diagnostics = {}

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        formula = atoms.get_chemical_formula()
        assert self.settings.spin == self.states[formula]
        assert self.settings.charge == 0
        self.observed.append(self.settings)
        energy = self.energies[formula]
        self.results = {"energy": energy, "forces": np.full((len(atoms), 3), 0.02)}
        self.diagnostics = {"scf_converged": True, "gradient_completed": True, "total_energy_hartree": energy / Hartree,
                            "s2": 0.75 if self.settings.spin else 0.0, "expected_s2": 0.75 if self.settings.spin else 0.0,
                            "s2_deviation": 0.0, "stability_checked": False, "test_fixture": True}
        self.event_log.write_text(json.dumps({"event": "test_fixture", "energy_ev": energy}) + "\n")


@pytest.fixture
def reference(tmp_path):
    path = tmp_path / "reference"
    shutil.copytree(benchmark.REFERENCE_DIRECTORY, path)
    return path


@pytest.fixture
def calculator(monkeypatch):
    KnownEnergies.observed = []
    monkeypatch.setattr(benchmark, "PySCFCalculator", KnownEnergies)
    return KnownEnergies


def test_fixed_geometry_arithmetic_state_policy_and_evidence(tmp_path, reference, calculator):
    output = tmp_path / "comparison"
    base = QuantumSettings(charge=7, spin=4, grid_level=5, basis="cc-pvtz", dispersion=None)
    result = benchmark.run_benchmark(base, reference, output)
    assert result["status"] == "completed"
    energy = result["computed"]
    assert energy["reference_geometry_energy_ev"] == pytest.approx(-0.5)
    assert energy["reference_geometry_energy_kcal_per_mol"] == pytest.approx(-0.5 / benchmark.EV_PER_KCAL_PER_MOL)
    assert energy["reaction_energy_ev"] == pytest.approx(-2.0)
    assert energy["reaction_energy_kcal_per_mol"] == pytest.approx(-2.0 / benchmark.EV_PER_KCAL_PER_MOL)
    assert "barrier" not in " ".join(energy)
    assert result["comparison"]["signed_discrepancy_kcal_per_mol"] == pytest.approx(-0.5 / benchmark.EV_PER_KCAL_PER_MOL - 2.2)
    assert result["comparison"]["absolute_discrepancy_kcal_per_mol"] == pytest.approx(0.5 / benchmark.EV_PER_KCAL_PER_MOL + 2.2)
    assert result["published_reference"]["value_kcal_per_mol"] == 2.2
    assert result["method_validated"] is result["saddle_verified"] is result["tool_validated"] is False
    assert "method, basis and geometry" in result["comparison"]["interpretation"]
    assert len(calculator.observed) == 5
    for settings in calculator.observed:
        assert replace(settings, charge=7, spin=4) == base
    assert json.loads((output / "benchmark.json").read_text()) == result
    manifest = json.loads((output / "input_manifest.json").read_text())
    assert manifest["reaction_conservation"] == {"elements": {"C": 3, "H": 5}, "charge": 0, "electrons": 23, "spin_2s": 1}
    assert manifest["data_license"].startswith("CC BY-NC")
    for filename, digest in manifest["file_sha256"].items():
        assert benchmark.sha256(output / "inputs" / filename) == digest
        assert (output / "inputs" / filename).read_bytes() == (reference / filename).read_bytes()
    for name, record in result["species"].items():
        assert json.loads((output / "species" / f"{name}.json").read_text()) == record
        assert (output / record["event_log"]).exists()
        assert record["max_force_ev_per_angstrom"] == pytest.approx(np.sqrt(3) * 0.02)
    with pytest.raises(FileExistsError):
        benchmark.run_benchmark(base, reference, output)


def test_persistent_output_is_required():
    with pytest.raises(ValueError, match="output is required"):
        benchmark.run_benchmark()


@pytest.mark.parametrize("error_type,status", [(RuntimeError, "failed"), (KeyboardInterrupt, "interrupted")])
def test_later_failure_preserves_completed_species(tmp_path, reference, monkeypatch, error_type, status):
    output = tmp_path / "failed-run"

    class StopsAtRadical(KnownEnergies):
        def calculate(self, atoms=None, properties=None, system_changes=all_changes):
            if atoms.get_chemical_formula() == "C2H":
                # A later failure must not erase the first expensive result.
                saved = json.loads((output / "benchmark.json").read_text())
                assert saved["status"] == "running"
                assert saved["active_species"] == "ethynyl_radical"
                assert set(saved["species"]) == {"methane"}
                assert (output / "species" / "methane.json").is_file()
                self.diagnostics = {"scf_converged": False, "error": "deliberate fixture failure"}
                raise error_type("deliberate fixture failure")
            return super().calculate(atoms, properties, system_changes)

    monkeypatch.setattr(benchmark, "PySCFCalculator", StopsAtRadical)
    with pytest.raises(error_type, match="deliberate fixture failure"):
        benchmark.run_benchmark(directory=reference, output=output)
    saved = json.loads((output / "benchmark.json").read_text())
    assert saved["status"] == status
    assert "computed" not in saved and "comparison" not in saved
    assert set(saved["species"]) == {"methane"}
    failed = json.loads((output / "species" / "ethynyl_radical.json").read_text())
    assert failed["status"] == status
    assert failed["error_type"] == error_type.__name__
    assert "energy_ev" not in failed


@pytest.mark.parametrize("fault", ["unconverged", "gradient_incomplete", "wrong_energy_unit", "nonfinite_energy", "nonfinite_force"])
def test_invalid_electronic_result_cannot_produce_comparison(tmp_path, reference, monkeypatch, fault):
    class InvalidResult(KnownEnergies):
        def calculate(self, *args, **kwargs):
            super().calculate(*args, **kwargs)
            if fault == "unconverged":
                self.diagnostics["scf_converged"] = False
            elif fault == "gradient_incomplete":
                self.diagnostics["gradient_completed"] = False
            elif fault == "wrong_energy_unit":
                self.diagnostics["total_energy_hartree"] = self.results["energy"]
            elif fault == "nonfinite_energy":
                self.results["energy"] = float("nan")
            else:
                self.results["forces"][0, 0] = float("inf")

    monkeypatch.setattr(benchmark, "PySCFCalculator", InvalidResult)
    output = tmp_path / "invalid"
    with pytest.raises(benchmark.BenchmarkError):
        benchmark.run_benchmark(directory=reference, output=output)
    result = json.loads((output / "benchmark.json").read_text())
    assert result["status"] == "failed"
    assert "computed" not in result
    assert json.loads((output / "species" / "methane.json").read_text())["status"] == "failed"


@pytest.mark.parametrize("fault", ["source", "charge", "multiplicity", "bool_charge", "float_count", "duplicate_metadata", "formula", "coordinates", "units", "extra_frame", "missing_file"])
def test_reference_tampering_rejected_before_quantum_work(tmp_path, reference, calculator, fault):
    metadata_file = reference / "provenance.json"
    metadata = json.loads(metadata_file.read_text())
    methane = next(entry for entry in metadata["structures"] if entry["file"] == "methane.xyz")
    xyz = reference / "methane.xyz"
    if fault == "source":
        with (reference / "temelso_2006_si.txt").open("a") as stream:
            stream.write("\nmodified\n")
    elif fault in ("charge", "multiplicity"):
        methane[fault] += 1
    elif fault == "bool_charge":
        methane["charge"] = False
    elif fault == "float_count":
        methane["atoms"] = 5.0
    elif fault == "duplicate_metadata":
        metadata["structures"].append(methane)
    elif fault == "formula":
        xyz.write_text(xyz.read_text().replace("C ", "N "))
    elif fault == "coordinates":
        # Moving a carbon preserves formula but must break source validation.
        lines = xyz.read_text().splitlines()
        row = lines[3].split()
        row[1] = str(float(row[1]) + 0.001)
        lines[3] = " ".join(row)
        xyz.write_text("\n".join(lines) + "\n")
    elif fault == "units":
        metadata["conversion"]["output_length_unit"] = "bohr"
    elif fault == "extra_frame":
        xyz.write_text(xyz.read_text() * 2)
    else:
        xyz.unlink()
    metadata_file.write_text(json.dumps(metadata))
    output = tmp_path / "tampered"
    with pytest.raises(benchmark.BenchmarkError):
        benchmark.run_benchmark(directory=reference, output=output)
    assert calculator.observed == []
    saved = json.loads((output / "benchmark.json").read_text())
    assert saved["status"] == "failed" and "computed" not in saved


def test_energies_are_computed_from_archived_geometry(tmp_path, reference, monkeypatch):
    class MutatesOriginal(KnownEnergies):
        def calculate(self, *args, **kwargs):
            (reference / "ethynyl_radical.xyz").write_text("corrupted after snapshot\n")
            return super().calculate(*args, **kwargs)

    monkeypatch.setattr(benchmark, "PySCFCalculator", MutatesOriginal)
    result = benchmark.run_benchmark(directory=reference, output=tmp_path / "archived")
    assert result["status"] == "completed"
    assert (tmp_path / "archived" / "inputs" / "ethynyl_radical.xyz").read_text() != (reference / "ethynyl_radical.xyz").read_text()


def test_unverified_source_energies_are_never_targets(tmp_path, reference, calculator):
    metadata_file = reference / "provenance.json"
    metadata = json.loads(metadata_file.read_text())
    for entry in metadata["structures"]:
        entry["source_energy_hartree"] = 123456.0
    metadata_file.write_text(json.dumps(metadata))
    result = benchmark.run_benchmark(directory=reference, output=tmp_path / "no-source-energy-target")
    assert result["computed"]["reference_geometry_energy_ev"] == pytest.approx(-0.5)
    assert result["published_reference"]["value_kcal_per_mol"] == 2.2
    assert "unverified and unused" in result["source_energy_caveat"]


def test_reaction_must_remain_expected_neutral_doublet(tmp_path, reference, calculator, monkeypatch):
    # If a future edit changes a species state, matching metadata alone must
    # not let an inconsistent combined reaction through to costly calculations.
    methane = dict(benchmark.SPECIES["methane"], spin=2)
    monkeypatch.setitem(benchmark.SPECIES, "methane", methane)
    metadata_file = reference / "provenance.json"
    metadata = json.loads(metadata_file.read_text())
    next(entry for entry in metadata["structures"] if entry["file"] == "methane.xyz")["multiplicity"] = 3
    metadata_file.write_text(json.dumps(metadata))
    with pytest.raises(benchmark.BenchmarkError, match="conserve"):
        benchmark.run_benchmark(directory=reference, output=tmp_path / "wrong-total-state")
    assert calculator.observed == []
