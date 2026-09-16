"""Compare saved characterization steps without running electronic calculations."""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
from ase.io import read

MAX_BYTES = 64 * 1024 * 1024
MAX_COORDINATES = 600
REQUIRED_QUANTUM = {"charge", "spin", "xc", "basis", "dispersion", "grid_level",
                    "conv_tol", "max_cycle", "density_fit", "scf_initial_guess"}


class ComparisonError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ComparisonError(message)


def positive(value, name):
    require(type(value) in (int, float) and math.isfinite(value) and value > 0,
            f"{name} must be finite and positive")
    return float(value)


def read_bytes(path):
    require(path.stat().st_size <= MAX_BYTES, f"Input exceeds {MAX_BYTES} byte limit")
    return path.read_bytes()


def strict_json(raw):
    def reject(value):
        raise ComparisonError(f"Nonfinite JSON value: {value}")
    return json.loads(raw, parse_constant=reject)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def finite_array(value, shape, name):
    def numbers(item):
        return all(numbers(x) for x in item) if isinstance(item,list) else type(item) in (int,float) and math.isfinite(item)
    require(isinstance(value,list) and numbers(value), f"{name} must contain only finite real numbers, not strings/booleans")
    result=np.asarray(value,dtype=float)
    require(result.shape==shape, f"{name} has incorrect shape")
    return result


def canonical_quantum(settings, notes, scope):
    require(isinstance(settings, dict), f"Missing {scope} quantum settings")
    value = deepcopy(settings)
    if "scf_initial_guess" not in value:
        value["scf_initial_guess"] = "minao"
        notes.append(f"{scope}: legacy absent scf_initial_guess interpreted as documented minao default; original bytes unchanged")
    require(REQUIRED_QUANTUM <= value.keys(), f"Incomplete {scope} method/state-input metadata")
    require(type(value["charge"]) is int and type(value["spin"]) is int and value["spin"] >= 0,
            f"Invalid {scope} charge/spin input")
    require(all(isinstance(value[k],str) and value[k].strip() for k in ("xc","basis")), f"Invalid {scope} method/basis")
    require(value["dispersion"] in (None,"d3bj","d3zero") and type(value["density_fit"]) is bool,
            f"Invalid {scope} dispersion/density-fitting input")
    require(type(value["grid_level"]) is int and 0 <= value["grid_level"] <= 9
            and type(value["max_cycle"]) is int and value["max_cycle"] > 0, f"Invalid {scope} electronic controls")
    require(positive(value["conv_tol"], f"{scope} electronic convergence tolerance") < 1, f"Invalid {scope} convergence tolerance")
    for key in ("threads","memory_mb"):
        if key in value:
            require(type(value[key]) is int and value[key]>0, f"Invalid {scope} {key}")
    require(value["scf_initial_guess"] in ("minao", "atom", "1e", "huckel"),
            f"Invalid {scope} initial guess")
    return value


def load_saved_result(path, structure_path=None):
    """Load workflow result.json plus an explicitly hash-bound XYZ snapshot.

    H1's successful workflow-shaped results use this same adapter. A standalone
    Hessian without method/geometry provenance is intentionally insufficient.
    """
    path = Path(path)
    raw = read_bytes(path)
    record = strict_json(raw)
    require(isinstance(record, dict), "Result must be a JSON object")
    require(type(record.get("schema_version")) is int and record["schema_version"]==1, "Only workflow schema_version 1 is supported")
    require(record.get("stage") == "characterize" and record.get("status") == "completed",
            "A completed workflow characterize result is required")
    units=record.get("units")
    required_units={"length":"angstrom", "energy":"eV", "frequency":"cm^-1"}
    require(isinstance(units,dict) and all(units.get(k)==v for k,v in required_units.items())
            and set(units) <= set(required_units)|{"force"}
            and ("force" not in units or units["force"]=="eV/angstrom"), "Unsupported or missing result units")
    structure_path = Path(structure_path) if structure_path else path.with_name("input.extxyz")
    structure_bytes = read_bytes(structure_path)
    hashes=record.get("input_hashes",{})
    require(isinstance(hashes,dict),"Input hashes must be an object")
    if "input_snapshot_sha256" in hashes:
        expected=hashes["input_snapshot_sha256"]
        hash_basis="input_snapshot_sha256"
        require(expected==digest(structure_bytes),"Explicit input snapshot hash is invalid or mismatched; source hash cannot override it")
    else:
        expected=hashes.get("structure_sha256")
        hash_basis="legacy_structure_sha256_equals_snapshot_bytes"
        require(expected==digest(structure_bytes),"Legacy structure/source hash does not bind this input snapshot; a distinct input_snapshot_sha256 is required")
    frames = read(io.StringIO(structure_bytes.decode("utf-8")), format="extxyz", index=":")
    require(len(frames) == 1, "A one-frame characterization input snapshot is required")
    atoms = frames[0]
    require(not atoms.pbc.any(), "Periodic characterization comparisons are unsupported")
    notes = []
    if hash_basis.startswith("legacy"):
        notes.append("Legacy source hash exactly matches snapshot bytes; no separate input_snapshot_sha256 was recorded")
    quantum = canonical_quantum(record.get("quantum_settings"), notes, "result")
    design_quantum = canonical_quantum(record.get("design", {}).get("quantum"), notes, "design")
    require(quantum == design_quantum, "Design and result quantum settings disagree")
    diagnostics = record.get("quantum_diagnostics")
    if diagnostics is not None:
        require(isinstance(diagnostics, dict), "Quantum diagnostics must be an object")
        actual = canonical_quantum(diagnostics.get("settings"), notes, "diagnostics")
        require(quantum == actual, "Diagnostics and result quantum settings disagree")
        if "scf_initial_guess" in diagnostics:
            require(diagnostics["scf_initial_guess"]==quantum["scf_initial_guess"], "Top-level diagnostic initial guess contradicts normalized settings")
    else:
        notes.append("Baseline quantum diagnostics absent; settings bind declared input, not independently linked electronic events")
    stationary = record.get("stationary")
    require(isinstance(stationary, dict), "Missing stationary characterization")
    require(record.get("design",{}).get("fixed_indices") == stationary.get("frozen_atom_indices"),
            "Design fixed atoms differ from characterization frozen atoms")
    reference = stationary.get("finite_difference_evidence", {}).get("reference_positions_angstrom")
    positions = atoms.positions
    geometry_binding = {"primary_source":"hash_bound_structure_snapshot", "full_precision_force_reference_available":False}
    if reference is not None:
        positions=finite_array(reference,(len(atoms),3),"recorded force-reference positions")
        exact=np.array_equal(positions,atoms.positions)
        rounded=np.array([[float(f"{coordinate:.8f}") for coordinate in row] for row in positions])
        require(exact or np.array_equal(rounded,atoms.positions), "Recorded force reference differs from snapshot beyond known eight-decimal serialization")
        geometry_binding={"primary_source":"recorded_full_precision_force_reference",
                          "full_precision_force_reference_available":True,
                          "snapshot_relation":"exact" if exact else "reference_rounded_to_eight_decimal_angstrom",
                          "maximum_snapshot_rounding_difference_angstrom":float(np.abs(positions-atoms.positions).max())}
        if not exact:
            notes.append("Snapshot matches force-reference coordinates rounded to eight decimal angstrom; full-precision recorded reference is used for exact cross-record binding")
    else:
        notes.append("Historical displaced-force evidence absent; raw Hessian reconstruction unavailable in this comparison")
        notes.append("Full-precision force-reference geometry unavailable; cross-record geometry binding is limited to saved snapshot precision")
    if "geometry_angstrom" in stationary:
        require(np.array_equal(np.asarray(stationary["geometry_angstrom"]), positions), "Characterization geometry differs from bound reference")
    return {"stationary": stationary, "symbols": atoms.get_chemical_symbols(),
            "positions_angstrom": positions.tolist(), "masses_amu": atoms.get_masses().tolist(),
            "geometry_binding":geometry_binding,
            "quantum_settings": quantum, "software_versions": (diagnostics or {}).get("versions"),
            "notes": notes,
            "source": {"result_path": str(path.resolve()), "result_sha256": digest(raw),
                       "input_snapshot_path": str(structure_path.resolve()), "input_snapshot_sha256": digest(structure_bytes),
                       "snapshot_hash_binding":hash_basis,"source_structure_sha256":hashes.get("structure_sha256")}}


def checked_modes(bound):
    char = bound["stationary"]
    n = len(bound["symbols"])
    positions = finite_array(bound["positions_angstrom"],(n,3),"bound positions")
    all_masses = finite_array(bound["masses_amu"],(n,),"bound masses")
    require(n > 0 and positions.shape == (n,3) and np.isfinite(positions).all(), "Invalid bound geometry")
    require(all_masses.shape == (n,) and np.isfinite(all_masses).all() and np.all(all_masses>0), "Invalid bound isotope masses")
    free, frozen = char["free_atom_indices"], char["frozen_atom_indices"]
    require(isinstance(free,list) and isinstance(frozen,list) and bool(free), "Free/frozen indices required")
    require(all(type(i) is int for i in free+frozen) and sorted(free+frozen)==list(range(n)), "Free/frozen atoms must partition geometry")
    d = len(free)*3
    require(d <= MAX_COORDINATES, f"At most {MAX_COORDINATES} free coordinates supported")
    require(char["coordinate_order"] == "x,y,z for each atom in free_atom_indices", "Unsupported coordinate ordering")
    masses = finite_array(char["free_masses_amu"],(len(free),),"free masses")
    require(np.array_equal(masses,all_masses[free]), "Mode masses differ from bound snapshot masses")
    modes = finite_array(char["cartesian_modes_per_sqrt_amu"],(d,len(free),3),"modes")
    frequencies = finite_array(char["frequencies_cm1"],(d,),"frequencies")
    require(modes.shape==(d,len(free),3) and frequencies.shape==(d,), "Incomplete mode/frequency arrays")
    require(np.isfinite(modes).all() and np.isfinite(frequencies).all(), "Nonfinite mode/frequency arrays")
    weighted = modes.reshape(d,d)*np.repeat(np.sqrt(masses),3)
    require(np.allclose(weighted@weighted.T,np.eye(d),atol=1e-7,rtol=0), "Modes are not orthonormal in the saved mass metric")
    positive(char["settings"]["step_angstrom"], "displacement step")
    positive(char["settings"]["frequency_tolerance_cm1"], "frequency tolerance")
    require(char.get("external_modes_removed") is False, "External-mode projection metadata missing or unsupported")
    return frequencies, weighted


def compare_bound(left, right, *, degeneracy_gap_cm1=10.):
    """Compare explicitly bound records; reject changes beyond step/metadata.

    Grouping merges adjacent spectral ranks when either spectrum's gap is below
    the declared tolerance. A group may span more than the pairwise tolerance.
    No individual identity is assigned inside a multi-mode group.
    """
    result = {"schema_version":1, "status":"not_comparable", "findings":[],
              "sources":[left.get("source"),right.get("source")],
              "legacy_or_missing_metadata":left.get("notes",[])+right.get("notes",[]),
              "metadata_notes_by_source":{"left":left.get("notes",[]),"right":right.get("notes",[])},
              "step_convergence_certified":False, "electronic_state_verified":False,
              "transition_state_verified":False, "connectivity_verified":False,
              "limitations":["Two saved steps measure observed sensitivity, not an asymptotic convergence or chemical-accuracy bound.",
                             "Matching calculation inputs does not prove equal electronic solutions at displaced geometries.",
                             "Mode overlap and spectral rank do not establish physical reaction identity or endpoint connectivity.",
                             "This comparator does not reconstruct forces/Hessians; independent E2 verification is separate."]}
    try:
        gap = positive(degeneracy_gap_cm1,"degeneracy gap")
        f1,w1 = checked_modes(left); f2,w2 = checked_modes(right)
        for name, bound in (("left",left),("right",right)):
            require(isinstance(bound.get("quantum_settings"),dict) and REQUIRED_QUANTUM <= bound["quantum_settings"].keys(),
                    f"Incomplete explicitly bound {name} method/state-input settings")
            canonical_quantum(bound["quantum_settings"],[],name)
        require(left["symbols"]==right["symbols"], "Atom order/elements differ")
        require(np.array_equal(left["positions_angstrom"],right["positions_angstrom"]), "Geometries or coordinate frames differ")
        result["geometry_binding_details"]=[left.get("geometry_binding"),right.get("geometry_binding")]
        require(np.array_equal(left["masses_amu"],right["masses_amu"]), "Atomic/isotope masses differ")
        require(left["quantum_settings"]==right["quantum_settings"], "Quantum method/state-input settings differ")
        a,b = left["stationary"],right["stationary"]
        for key in ("free_atom_indices","frozen_atom_indices","coordinate_order","method","external_modes_removed"):
            require(a.get(key)==b.get(key) and key in a and key in b, f"Characterization {key} differs or is missing")
        ca,cb = deepcopy(a["settings"]),deepcopy(b["settings"])
        step_a,step_b = ca.pop("step_angstrom"),cb.pop("step_angstrom")
        require(ca==cb, "Characterization controls beyond displacement step differ")
        va,vb = left.get("software_versions"),right.get("software_versions")
        result["software_versions"]=[va,vb]
        require(va is None or vb is None or va==vb, "Known software versions differ; step-only comparison is confounded")
        result["recorded_software_versions_match"] = va is not None and va==vb
        if va!=vb or va is None:
            result["legacy_or_missing_metadata"].append("Software version metadata differs or is missing; observed differences may include software effects")
        order_a,order_b = np.argsort(f1,kind="stable"),np.argsort(f2,kind="stable")
        f1,f2,w1,w2 = f1[order_a],f2[order_b],w1[order_a],w2[order_b]
        overlap=np.clip((w1@w2.T)**2,0.,1.)
        groups=[[0]]
        for i in range(1,len(f1)):
            if min(f1[i]-f1[i-1],f2[i]-f2[i-1]) <= gap:
                groups[-1].append(i)
            else:
                groups.append([i])
        grouped=[]
        for ranks in groups:
            cosines=np.clip(np.linalg.svd(w1[ranks]@w2[ranks].T,compute_uv=False),0.,1.)
            grouped.append({"spectral_ranks":ranks, "left_mode_indices":order_a[ranks].tolist(),
                            "right_mode_indices":order_b[ranks].tolist(),"dimension":len(ranks),
                            "comparison_kind":"individual_direction" if len(ranks)==1 else "subspace_only",
                            "individual_identity_assigned":False,
                            "full_coordinate_space_group":len(ranks)==len(f1),
                            "directional_overlap_informative":len(ranks)<len(f1),
                            "left_frequency_range_cm1":[float(f1[ranks[0]]),float(f1[ranks[-1]])],
                            "right_frequency_range_cm1":[float(f2[ranks[0]]),float(f2[ranks[-1]])],
                            "principal_cosines_squared":(cosines**2).tolist(),
                            "maximum_principal_angle_degrees":float(np.degrees(np.arccos(min(cosines))))})
        tolerance=ca["frequency_tolerance_cm1"]
        def labels(freq):
            return ["negative" if f < -tolerance else "positive" if f > tolerance else "unresolved" for f in freq]
        sign_a,sign_b=labels(f1),labels(f2)
        result.update(status="compared", binding={"geometry_and_atom_order":"identical", "masses_and_free_space":"identical",
                       "quantum_method_and_state_inputs":"identical_normalized_settings_with_any_legacy_defaults_reported", "only_characterization_control_change":"step_angstrom" if step_a!=step_b else "none"},
                       steps_angstrom=[step_a,step_b], distinct_steps=step_a!=step_b,
                       degeneracy_grouping={"gap_cm1":gap,"rule":"connected adjacent spectral ranks; merge if the gap in either spectrum is <= declared threshold; transitive groups may span more than threshold"},
                       sorted_mode_indices=[order_a.tolist(),order_b.tolist()],
                       sorted_frequencies_cm1=[f1.tolist(),f2.tolist()],
                       spectral_rank_frequency_changes_cm1=(f2-f1).tolist(),
                       frequency_changes_are_identity_assignments=False,
                       resolved_sign_labels=[sign_a,sign_b],
                       raw_negative_frequency_counts=[int(np.count_nonzero(f1<0)),int(np.count_nonzero(f2<0))],
                       resolved_negative_mode_counts=[sign_a.count("negative"),sign_b.count("negative")],
                       changed_sign_spectral_ranks=[i for i,(x,y) in enumerate(zip(sign_a,sign_b)) if x!=y],
                       squared_mass_metric_overlap_matrix=overlap.tolist(), groups=grouped)
    except (KeyError,TypeError,ValueError,IndexError) as error:
        result["findings"].append(str(error))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left",type=Path);parser.add_argument("right",type=Path)
    parser.add_argument("--left-structure",type=Path);parser.add_argument("--right-structure",type=Path)
    parser.add_argument("--degeneracy-gap-cm1",type=float,default=10.)
    parser.add_argument("--output",type=Path,help="New JSON file only; existing files are never overwritten")
    args=parser.parse_args()
    try:
        left=load_saved_result(args.left,args.left_structure)
        right=load_saved_result(args.right,args.right_structure)
        result=compare_bound(left,right,degeneracy_gap_cm1=args.degeneracy_gap_cm1)
    except (OSError,KeyError,TypeError,ValueError) as error:
        result={"schema_version":1,"status":"not_comparable","findings":[str(error)],
                "step_convergence_certified":False,"electronic_state_verified":False,"transition_state_verified":False}
    result["comparison_implementation_sha256"]=digest(Path(__file__).read_bytes())
    result["quantum_jobs_started"]=0
    encoded=json.dumps(result,indent=2,allow_nan=False)+"\n"
    if args.output:
        with args.output.open("x") as stream:stream.write(encoded)
    else:print(encoded,end="")
    return 0 if result["status"]=="compared" else 2


if __name__=="__main__":
    raise SystemExit(main())
