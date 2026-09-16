"""Does removing spin contamination move the +2.40 kcal/mol reference barrier?

The in-house reference uses a UHF reference (S^2 = 1.21 at the TS). The
published paper used an ROHF-based reference, which has no spin contamination
in its orbitals. The reliability diagnostics (see diagnose_reference.py) flag
the UHF open-shell species as strongly correlated/contaminated. This is the
decisive test of what that actually costs: recompute the bare electronic
barrier

    barrier = E(TS) - E(CH4) - E(C2H)

with CCSD(T) on BOTH an ROHF and a UHF reference, at cc-pVDZ on the archived
geometries, and compare. If the two agree, the UHF reference was fine despite
the diagnostic and +2.40 stands. If they diverge, the reference — and every
"DFT is wrong by 5 kcal/mol" statement resting on it — moves by the difference.

The closed-shell CH4 is identical under both references; only the open-shell
C2H and the TS differ. So the whole ROHF-vs-UHF barrier difference is carried
by the two contaminated species, which is exactly where the diagnostic said to
look.

Run:  python3 research/reference-reliability/rohf_vs_uhf_barrier.py
Compute is small (cc-pVDZ, <=8 atoms) and runs off-host, serial, one process.
"""

from __future__ import annotations

import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import Hartree
from pyscf import gto, scf, cc

ROOT = Path(__file__).resolve().parents[2]
REFDIR = ROOT / "data/reference"
HERE = Path(__file__).resolve().parent
KCAL = 627.5094740631  # Hartree -> kcal/mol

# species needed for the barrier, with (charge, spin=2S)
BARRIER = {
    "methane.xyz": (0, 0),
    "ethynyl_radical.xyz": (0, 1),
    "methane_ethynyl_ts.xyz": (0, 1),
}


def ccsd_t_energy(fname, charge, spin, reference):
    """Return CCSD(T) total energy (Hartree) and diagnostics for one species.

    reference in {"UHF", "ROHF"}; closed-shell (spin==0) always uses RHF and
    both branches coincide there.
    """
    atoms = read(REFDIR / fname)
    mol = gto.M(atom=[(a.symbol, tuple(a.position)) for a in atoms],
                unit="Angstrom", basis="cc-pvdz", charge=charge, spin=spin, verbose=0)
    if spin == 0:
        mf = scf.RHF(mol)
    elif reference == "ROHF":
        mf = scf.ROHF(mol)
    else:
        mf = scf.UHF(mol)
    mf.init_guess = "atom"
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    mf.kernel()
    ss = float(mf.spin_square()[0]) if spin else 0.0
    mycc = cc.CCSD(mf)
    mycc.conv_tol = 1e-9
    mycc.max_cycle = 150
    ecorr, t1, _ = mycc.kernel()
    et = mycc.ccsd_t()
    e_tot = float(mf.e_tot + ecorr + et)
    # T1 diagnostic under this reference
    if isinstance(t1, (tuple, list)):
        nocc = sum(int(x.shape[0]) for x in t1)
        norm2 = sum(float(np.sum(np.asarray(x) ** 2)) for x in t1)
    else:
        nocc = 2 * int(t1.shape[0])
        norm2 = float(np.sum(np.asarray(t1) ** 2))
    t1_diag = float(np.sqrt(norm2 / nocc)) if nocc else 0.0
    return {
        "reference": "RHF" if spin == 0 else reference,
        "scf_converged": bool(mf.converged),
        "ccsd_converged": bool(mycc.converged),
        "hf_s_squared": round(ss, 4),
        "ccsd_t_total_hartree": e_tot,
        "t1_diagnostic": round(t1_diag, 5),
    }


def barrier(records):
    e = {k: records[k]["ccsd_t_total_hartree"] for k in BARRIER}
    return (e["methane_ethynyl_ts.xyz"] - e["methane.xyz"] - e["ethynyl_radical.xyz"]) * KCAL


def main():
    started = time.perf_counter()
    out = {"generated_utc": datetime.now(timezone.utc).isoformat(),
           "compute_host": {"system": os.uname().sysname, "machine": os.uname().machine,
                            "loadavg": os.getloadavg()},
           "barrier_definition": "E(TS) - E(CH4) - E(C2H), CCSD(T)/cc-pVDZ, atom guess",
           "references": {}}
    for reference in ("UHF", "ROHF"):
        recs = {f: ccsd_t_energy(f, c, s, reference) for f, (c, s) in BARRIER.items()}
        b = barrier(recs)
        out["references"][reference] = {"species": recs, "barrier_kcal_per_mol": round(b, 4)}
        print(f"{reference}: barrier = {b:.4f} kcal/mol  "
              f"(TS T1={recs['methane_ethynyl_ts.xyz']['t1_diagnostic']}, "
              f"S^2={recs['methane_ethynyl_ts.xyz']['hf_s_squared']})")

    u = out["references"]["UHF"]["barrier_kcal_per_mol"]
    r = out["references"]["ROHF"]["barrier_kcal_per_mol"]
    out["comparison"] = {
        "uhf_barrier_kcal_per_mol": u,
        "rohf_barrier_kcal_per_mol": r,
        "rohf_minus_uhf_kcal_per_mol": round(r - u, 4),
        "published_table4_kcal_per_mol": 2.4,
        "reading": (
            "If |ROHF - UHF| is small (<~0.5 kcal/mol) the UHF reference was "
            "adequate despite the T1/D1 flags and +2.40 stands. If it is "
            "comparable to or larger than the ~5 kcal/mol DFT discrepancy, the "
            "reference itself is method-sensitive and the DFT comparison must "
            "carry that as an explicit reference uncertainty."
        ),
        "limits": [
            "ROHF-based CCSD here uses ROHF orbitals in PySCF's CCSD; it is close to but not bit-identical to the paper's RCCSD(T) variant.",
            "Neither reference resolves genuine multireference character; only a multireference method (CASPT2/NEVPT2) does. Agreement between ROHF and UHF bounds spin-contamination sensitivity, not static-correlation error.",
            "cc-pVDZ, fixed supplied collinear TS geometry with three reported imaginary modes.",
        ],
    }
    out["elapsed_seconds"] = round(time.perf_counter() - started, 1)
    (HERE / "rohf-vs-uhf-barrier.json").write_text(json.dumps(out, indent=1))
    print(f"\nROHF - UHF = {r - u:+.4f} kcal/mol   (published Table 4: +2.4)")
    print("written:", (HERE / "rohf-vs-uhf-barrier.json").relative_to(ROOT))


if __name__ == "__main__":
    main()
