"""Synthetic campaign fixtures for software checks; never scientific evidence.

All artifacts are created in a caller-supplied temporary directory. No quantum
solver runs. The numbers only satisfy the campaign reader's document contract.
"""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

from ase import Atoms
from ase.io import write

from nanodesign.campaign import create_campaign
from nanodesign.design import load_design, sha256
from nanodesign.quantum import QuantumSettings


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def fixture(root, *, saved_status="completed", input_guess=None, result_guess=None,
            drop_setting=None, extra_setting=None):
    """None means absent legacy key; other values are literal fixture contents."""
    root = Path(root)
    source = root / "source"
    source.mkdir(parents=True)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, .74]])
    write(source / "initial.xyz", atoms)
    write(source / "final.xyz", atoms)
    quantum = asdict(QuantumSettings(spin=0, threads=1))
    quantum.pop("scf_initial_guess", None)
    if input_guess is not None:
        quantum["scf_initial_guess"] = input_guess
    save(source / "design.json", {
        "schema_version": 1, "length_unit": "angstrom", "initial": "initial.xyz",
        "final": "final.xyz", "fixed_indices": [], "quantum": quantum,
        "metadata": {"fixture_only": True},
    })
    campaign_root = root / "campaign"
    create_campaign([source / "design.json"], campaign_root)
    # Construct a legacy plan document in the fixture only. Production archives
    # are never modified or rehashed by this review.
    plan_path = campaign_root / "plan.json"
    plan = json.loads(plan_path.read_text())
    plan["quantum_settings"] = dict(quantum)
    save(plan_path, plan)
    state_path = campaign_root / "campaign.json"
    state = json.loads(state_path.read_text())
    state["plan_sha256"] = sha256(plan_path)
    if saved_status is not None:
        entry = plan["designs"][0]
        output = "runs/pose-0001/attempt-0001"
        result_path = campaign_root / output / "result.json"
        result_path.parent.mkdir(parents=True)
        saved_settings = dict(quantum)
        saved_settings.pop("scf_initial_guess", None)
        if result_guess is not None:
            saved_settings["scf_initial_guess"] = result_guess
        if drop_setting:
            saved_settings.pop(drop_setting)
        if extra_setting:
            saved_settings.update(extra_setting)
        controls = plan["controls"]
        save(result_path, {
            "fixture_only": True, "status": saved_status,
            "stage": controls["stage"], "state": controls["state"],
            "input_hashes": load_design(campaign_root / entry["design"])[4],
            "quantum_settings": saved_settings,
            "optimization": {"fmax_ev_per_angstrom": controls["fmax"],
                             "max_steps_per_stage": controls["steps"], "images": controls["images"]},
            "structure": {"energy_ev": -1.0, "free_force_max_ev_per_angstrom": 0.0,
                "forces_ev_per_angstrom": [[0, 0, 0], [0, 0, 0]],
                "quantum_diagnostics": {"scf_converged": saved_status == "completed",
                                        "gradient_completed": saved_status == "completed"}},
            "validation": {"design_validated": False},
        })
        state["attempts"]["pose-0001"] = [{
            "status": saved_status, "output": output, "result_sha256": sha256(result_path),
        }]
    save(state_path, state)
    return campaign_root


def file_hashes(root):
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in Path(root).rglob("*") if path.is_file()}
