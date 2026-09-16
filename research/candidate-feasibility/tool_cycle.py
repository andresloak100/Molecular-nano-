"""Can the tool be recharged? A tip that fires once is not a machine.

Rung 12 of this lane's completeness ladder, and the only one marked "not even
posed". Every result in this repository concerns a single half-operation on a
fresh tool: the ethynyl tip abstracts one hydrogen and becomes
acetylene-terminated. Mechanosynthesis needs *reusable* tools, so removing that
hydrogen again is a second reaction with its own barrier, its own selectivity
problem and its own error rate. It may well be the harder half, and nothing
models it.

The thermodynamic requirement is stateable immediately and is decided by a
number this lane already measured.

THE TRAP. Recharging by abstraction means

    R-C#C-H  +  X*  ->  R-C#C*  +  X-H

which is downhill only if the partner forms a STRONGER bond to hydrogen than
the tip does:

    BDE(X-H)  >  BDE(tip C-H)

This lane measured the tip's bond at 136.15 kcal/mol for the real adamantyl
handle (bare electronic, no zero-point). That is an extremely strong bond -
terminal alkyne C-H bonds are among the strongest in chemistry - and it is
strong for exactly the reason the tool works at all. The abstraction is driven
by forming it.

So the design has a structural tension, and it is the third face of the same
reactivity-selectivity trade-off already seen twice in this project:

    a hungry tip abstracts well          (good)
    a hungry tip discriminates poorly    (the selectivity problem)
    a hungry tip is hard to recharge     (this file)

All three follow from one number. You cannot soften one without softening the
others.

WHAT THIS MODULE DOES. It computes candidate partner X-H bond strengths at the
same level of theory as the tip's, so the comparison is internally consistent
rather than a mix of in-house and literature values. Mixing is what produced a
retracted conclusion in this project earlier today.

WHAT IT DOES NOT DO. Thermodynamics is not kinetics: a downhill recharge still
needs an accessible barrier, and an uphill one is not strictly impossible if
driven mechanically or photochemically. This bounds the chemistry, it does not
settle the engineering. Nor does it consider handing the hydrogen to a separate
purpose-built acceptor tool, which is a real architecture and a different
calculation.

Run from the repository root:

    python research/candidate-feasibility/tool_cycle.py
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

from nanodesign.quantum import PySCFCalculator, QuantumSettings  # noqa: E402
from handle_fidelity import guess_scan, species_energy  # noqa: E402
from timing import load_snapshot  # noqa: E402

EVIDENCE = LANE / "evidence"
EV_KCAL = 23.060548

# Candidate hydrogen acceptors, smallest first. Each entry is the radical X*
# and its hydride X-H, with geometries as simple internal-coordinate guesses at
# roughly standard bond lengths. These are unrelaxed, exactly like the tip
# fragments they are compared against, so the comparison is like-for-like: both
# sides carry the same rigid-geometry approximation.
PARTNERS = {
    "fluorine": {
        "radical": ("F", [[0.0, 0.0, 0.0]], 1),
        "hydride": ("FH", [[0.0, 0.0, 0.0], [0.0, 0.0, 0.917]], 0),
        "why": "The only common element whose X-H bond rivals a terminal alkyne C-H.",
    },
    "hydroxyl": {
        "radical": ("OH", [[0.0, 0.0, 0.0], [0.0, 0.0, 0.970]], 1),
        "hydride": ("OH2", [[0.0, 0.0, 0.0], [0.0, 0.0, 0.958], [0.0, 0.927, -0.240]], 0),
        "why": "Strong but probably not strong enough; also oxidises a diamond surface.",
    },
    "hydrogen_atom": {
        "radical": ("H", [[0.0, 0.0, 0.0]], 1),
        "hydride": ("HH", [[0.0, 0.0, 0.0], [0.0, 0.0, 0.741]], 0),
        "why": "Forms H2. The reference case for 'just pull it off with atomic hydrogen'.",
    },
    "methyl": {
        "radical": ("CH3", [[0.0, 0.0, 0.0], [1.079, 0.0, 0.0],
                            [-0.539, 0.934, 0.0], [-0.539, -0.934, 0.0]], 1),
        "hydride": ("CH4", [[0.0, 0.0, 0.0], [0.629, 0.629, 0.629], [-0.629, -0.629, 0.629],
                            [0.629, -0.629, -0.629], [-0.629, 0.629, -0.629]], 0),
        "why": "A generic carbon radical: the case where recharge is clearly uphill.",
    },
}

# Measured by this lane, same method, rigid candidate-pose geometry.
TIP_BOND_KCAL = {
    "adamantyl": 136.15,
    "methyl": 136.52,
    "hydrogen": 136.72,
}


def _atoms(symbols: str, positions, spin: int) -> Atoms:
    atoms = Atoms(symbols, positions=np.array(positions, dtype=float))
    moments = np.zeros(len(atoms))
    if spin:
        moments[0] = 1.0
    atoms.set_initial_magnetic_moments(moments)
    atoms.info["geometry_status"] = "unrelaxed_standard_bond_lengths"
    return atoms


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--xc", default="pbe0")
    parser.add_argument("--tip", default="adamantyl", choices=sorted(TIP_BOND_KCAL))
    parser.add_argument("--out", default="tool-cycle.json")
    args = parser.parse_args()

    EVIDENCE.mkdir(parents=True, exist_ok=True)
    out = EVIDENCE / args.out
    if out.exists():
        raise SystemExit(f"{out} exists; evidence is never overwritten. Choose --out.")

    base = QuantumSettings(
        xc=args.xc, basis=args.basis, dispersion="d3bj",
        density_fit=True, threads=1,
    )
    tip_bond = TIP_BOND_KCAL[args.tip]

    report: dict = {
        "schema_version": 1,
        "question": "Can the ethynyl tip be recharged after it abstracts a hydrogen?",
        "requirement": "BDE(X-H) > BDE(tip C-H) for abstraction recharge to be downhill",
        "tip_handle": args.tip,
        "tip_ch_bond_kcal_per_mol": tip_bond,
        "tip_bond_provenance": (
            "Measured by this lane in handle-fidelity, same method and same "
            "rigid candidate-pose geometry as the partners below."
        ),
        "method": f"{args.xc.upper()}-D3(BJ)/{args.basis}, density fitting",
        "geometry_status": "all species unrelaxed; both sides carry the same approximation",
        "excludes": [
            "zero-point energy, which differs between X-H and C-H and is not small",
            "kinetics: a downhill recharge still needs an accessible barrier",
            "mechanical or photochemical driving, which can pay for an uphill step",
            "handing the hydrogen to a separate purpose-built acceptor tool",
        ],
        "host_at_start": load_snapshot(),
        "species": {},
        "partners": {},
        "status": "running",
    }
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    def save() -> None:
        out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    print(f"tip C-H bond ({args.tip} handle): {tip_bond:.2f} kcal/mol")
    print(f"recharge is downhill only for partners above that\n")

    hydrogen = _atoms("H", [[0.0, 0.0, 0.0]], 1)
    report["species"]["H_atom"] = species_energy(hydrogen, replace(base, spin=1), "cycle_H_atom")
    save()
    h_energy = report["species"]["H_atom"]["energy_ev"]

    for name, entry in PARTNERS.items():
        print(f"\npartner: {name}")
        radical_symbols, radical_positions, radical_spin = entry["radical"]
        hydride_symbols, hydride_positions, hydride_spin = entry["hydride"]
        radical = _atoms(radical_symbols, radical_positions, radical_spin)
        hydride = _atoms(hydride_symbols, hydride_positions, hydride_spin)

        radical_record = species_energy(radical, replace(base, spin=radical_spin), f"cycle_{name}_radical")
        report["species"][f"{name}_radical"] = radical_record
        save()
        hydride_record = species_energy(hydride, replace(base, spin=hydride_spin), f"cycle_{name}_hydride")
        report["species"][f"{name}_hydride"] = hydride_record
        save()

        # BDE(X-H) = E(X*) + E(H) - E(X-H), positive for a bound hydride.
        bde_ev = radical_record["energy_ev"] + h_energy - hydride_record["energy_ev"]
        bde = bde_ev * EV_KCAL
        recharge = tip_bond - bde  # positive = uphill
        report["partners"][name] = {
            "bde_x_h_kcal_per_mol": bde,
            "recharge_reaction_energy_kcal_per_mol": recharge,
            "recharge_is_downhill": recharge < 0,
            "why_considered": entry["why"],
        }
        verdict = "DOWNHILL" if recharge < 0 else "uphill"
        print(f"    BDE(X-H) = {bde:7.2f} kcal/mol   recharge {recharge:+7.2f}  {verdict}")
        save()

    downhill = [n for n, v in report["partners"].items() if v["recharge_is_downhill"]]
    report["summary"] = {
        "partners_that_work": downhill,
        "n_partners_tested": len(report["partners"]),
        "finding": (
            "Recharging by abstraction needs a partner forming a stronger bond "
            "to hydrogen than a terminal alkyne C-H, which is a very short list."
            if len(downhill) <= 1 else
            "Several partners are thermodynamically capable of recharging the tip."
        ),
        "structural_point": (
            "The tip's bond strength is what drives the abstraction and is also "
            "what makes recharge hard. Reactivity, selectivity and "
            "recyclability are set by one number and cannot be tuned "
            "independently."
        ),
        "not_a_verdict_on_the_design": (
            "Thermodynamics only, no zero-point, no barriers, no mechanical "
            "driving, and no separate acceptor-tool architecture. An uphill "
            "recharge is a constraint on the architecture, not a proof of "
            "impossibility."
        ),
    }
    report["status"] = "completed"
    report["host_at_end"] = load_snapshot()
    save()

    print(f"\npartners able to recharge the tip: {downhill or 'none of those tested'}")
    print(f"wrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
