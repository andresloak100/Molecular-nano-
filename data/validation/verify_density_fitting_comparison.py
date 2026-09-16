"""Independent, standard-library audit of archived DF/direct comparison arithmetic.

Run with Python 3.11+; no nanodesign, NumPy or quantum-chemistry imports required.
This checks saved evidence and arithmetic, not the validity of the DFT model.
"""
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text())


def close(actual, expected, *, absolute=1e-11):
    assert math.isclose(actual, expected, rel_tol=0, abs_tol=absolute), (actual, expected)


def main():
    direct_dir = ROOT / "h-abstraction-direct-initial"
    fitted_dir = ROOT / "h-abstraction-df-initial"
    direct, fitted = [load(directory / "result.json") for directory in (direct_dir, fitted_dir)]
    comparison = load(fitted_dir / "comparison.json")
    for directory in (direct_dir, fitted_dir):
        metadata = load(directory / "interpretation.json")
        for filename, digest in metadata["file_sha256"].items():
            assert hashlib.sha256((directory / filename).read_bytes()).hexdigest() == digest
    assert hashlib.sha256((direct_dir / "result.json").read_bytes()).hexdigest() == comparison["direct_result_sha256"]
    assert direct["status"] == fitted["status"] == "completed"
    assert dict(direct["quantum_settings"], density_fit=True) == fitted["quantum_settings"]
    for name in ("initial_sha256", "final_sha256"):
        assert direct["input_hashes"][name] == fitted["input_hashes"][name]
    assert direct["design"]["fixed_indices"] == fitted["design"]["fixed_indices"]
    ds, fs = direct["structure"], fitted["structure"]
    de = ds["quantum_diagnostics"]
    fe = fs["quantum_diagnostics"]
    assert de["hartree_eV"] == fe["hartree_eV"]
    for structure in (ds, fs):
        diagnostics = structure["quantum_diagnostics"]
        assert diagnostics["scf_converged"] and diagnostics["gradient_completed"]
        close(diagnostics["total_energy_hartree"] * diagnostics["hartree_eV"], structure["energy_ev"], absolute=1e-9)
    energy_difference = float(Decimal(str(fs["energy_ev"])) - Decimal(str(ds["energy_ev"])))
    close(comparison["energy_difference_ev"], energy_difference)
    close(comparison["energy_difference_hartree"] * de["hartree_eV"], energy_difference)
    # Exact SI definitions, calculated independently from the application.
    ev_per_kcal_per_mol = Decimal(4184) / (Decimal("6.02214076e23") * Decimal("1.602176634e-19"))
    kcal_difference = float(Decimal(str(energy_difference)) / ev_per_kcal_per_mol)
    close(comparison["energy_difference_kcal_per_mol"], kcal_difference)
    close(comparison["unit_conversion"]["ev_per_kcal_per_mol"], float(ev_per_kcal_per_mol))
    force_rows = []
    for fitted_row, direct_row, saved in zip(fs["forces_ev_per_angstrom"], ds["forces_ev_per_angstrom"], comparison["force_differences_ev_per_angstrom"], strict=True):
        assert len(fitted_row) == len(direct_row) == len(saved) == 3
        row = [f - d for f, d in zip(fitted_row, direct_row, strict=True)]
        for actual, expected in zip(saved, row, strict=True):
            close(actual, expected, absolute=1e-15)
        force_rows.append(row)
    assert len(force_rows) == 53
    norms = [math.sqrt(sum(value * value for value in row)) for row in force_rows]
    flattened = [value for row in force_rows for value in row]
    max_vector = max(norms)
    rms_component = math.sqrt(sum(value * value for value in flattened) / len(flattened))
    close(comparison["max_atomic_force_vector_difference_ev_per_angstrom"], max_vector, absolute=1e-15)
    close(comparison["max_force_component_difference_ev_per_angstrom"], max(map(abs, flattened)), absolute=1e-15)
    close(comparison["rms_force_component_difference_ev_per_angstrom"], rms_component, absolute=1e-15)
    fixed = set(direct["design"]["fixed_indices"])
    close(comparison["max_free_atom_force_vector_difference_ev_per_angstrom"], max(norm for atom, norm in enumerate(norms) if atom not in fixed), absolute=1e-15)
    close(comparison["observed_walltime_ratio_direct_over_density_fitting"], direct["elapsed_seconds"] / fitted["elapsed_seconds"])
    print(json.dumps({"arithmetic_and_saved_hashes_verified": True, "atoms": len(force_rows),
                      "energy_difference_ev": energy_difference, "energy_difference_kcal_per_mol": kcal_difference,
                      "max_force_vector_difference_ev_per_angstrom": max_vector,
                      "rms_force_component_difference_ev_per_angstrom": rms_component,
                      "scientific_model_validated": False}, indent=2))


if __name__ == "__main__":
    main()
