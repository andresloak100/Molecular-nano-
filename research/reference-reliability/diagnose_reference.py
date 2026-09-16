"""Is the in-house CCSD(T) reference itself trustworthy at the transition structure?

The whole "DFT gets this barrier wrong by ~5 kcal/mol" conclusion is measured
against an in-house CCSD(T)/cc-pVDZ number (+2.40 kcal/mol) treated as a gold
standard. But CCSD(T) is a *single-reference* method: it is reliable only when
the Hartree-Fock determinant dominates the wavefunction. At the collinear
transition structure the UHF reference has S^2 = 1.21 versus 0.75 for a clean
doublet, and multireference character has been suspected but never quantified.

This script quantifies it with the standard coupled-cluster reliability
diagnostics, which are byproducts of the CCSD amplitudes the reference already
computes but never saved:

  T1 diagnostic (Lee & Taylor 1989): ||t1|| / sqrt(n_correlated_electrons).
      Closed-shell rule of thumb: > 0.02 means single-reference CCSD(T) is
      questionable. Open-shell references run higher intrinsically, so the
      honest reading here is RELATIVE: the transition structure against the
      well-behaved closed-shell species computed the same way.

  D1 diagnostic (Janssen & Nielsen 1998): largest singular value of the t1
      matrix. More sensitive to localized near-degeneracy. Closed-shell rule
      of thumb: > 0.05 questionable.

  Largest |t2| amplitude: a single dominant doubles amplitude is a direct
      fingerprint of a strongly correlated pair the single reference misses.

It recomputes the five cc-pVDZ species on the archived geometries at the
lowest SCF solution (guess 'atom', which the forensics lane established
reproduces the published energies), so the diagnostics attach to the exact
determinant the +2.40 reference uses. No archived record is modified; this
writes only into its own lane. Compute runs off-host (cloud container).

Run:  python3 research/reference-reliability/diagnose_reference.py
"""

from __future__ import annotations

import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from ase.io import read
from pyscf import gto, scf, cc

ROOT = Path(__file__).resolve().parents[2]
REFDIR = ROOT / "data/reference"
HERE = Path(__file__).resolve().parent

# Standard literature thresholds, cited not invented. They are closed-shell
# rules of thumb; for open-shell species the relative reading is what matters.
T1_CLOSED_SHELL_THRESHOLD = 0.02   # Lee & Taylor, Int. J. Quantum Chem. 1989
D1_CLOSED_SHELL_THRESHOLD = 0.05   # Janssen & Nielsen, CPL 1998

# geometry file -> (charge, spin=2S, role)
SPECIES = {
    "methane.xyz": (0, 0, "reactant (closed shell, control)"),
    "acetylene.xyz": (0, 0, "product (closed shell, control)"),
    "methyl_radical.xyz": (0, 1, "product radical (doublet)"),
    "ethynyl_radical.xyz": (0, 1, "reactant radical (doublet)"),
    "methane_ethynyl_ts.xyz": (0, 1, "TRANSITION STRUCTURE (doublet)"),
}


def frobenius_norm2(t1):
    if isinstance(t1, (tuple, list)):
        return float(sum(np.sum(np.asarray(x) ** 2) for x in t1))
    return float(np.sum(np.asarray(t1) ** 2))


def largest_singular_value(t1):
    def sv(x):
        x = np.asarray(x)
        return float(np.linalg.svd(x.reshape(x.shape[0], -1), compute_uv=False)[0]) if x.size else 0.0
    if isinstance(t1, (tuple, list)):
        return max(sv(x) for x in t1)
    return sv(t1)


def largest_amplitude(t):
    if isinstance(t, (tuple, list)):
        return max((float(np.max(np.abs(x))) if np.asarray(x).size else 0.0) for x in t)
    return float(np.max(np.abs(t))) if np.asarray(t).size else 0.0


def run_species(fname, charge, spin, guess="atom"):
    atoms = read(REFDIR / fname)
    mol = gto.M(
        atom=[(a.symbol, tuple(a.position)) for a in atoms],
        unit="Angstrom", basis="cc-pvdz", charge=charge, spin=spin, verbose=0,
    )
    unrestricted = spin != 0
    mf = scf.UHF(mol) if unrestricted else scf.RHF(mol)
    mf.init_guess = guess
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    mf.kernel()
    ss, _ = mf.spin_square() if unrestricted else (0.0, 1.0)
    solver = cc.UCCSD(mf) if unrestricted else cc.RCCSD(mf)
    solver.conv_tol = 1e-9
    solver.max_cycle = 150
    _, t1, t2 = solver.kernel()

    if unrestricted:
        nocc = int(solver.nocc[0]) + int(solver.nocc[1])
    else:
        nocc = 2 * int(solver.nocc)
    t1_diag = float(np.sqrt(frobenius_norm2(t1) / nocc))
    d1_diag = largest_singular_value(t1)
    return {
        "file": fname,
        "charge": charge,
        "spin_2S": spin,
        "reference": "UHF/UCCSD(T)" if unrestricted else "RHF/RCCSD(T)",
        "scf_initial_guess": guess,
        "scf_converged": bool(mf.converged),
        "ccsd_converged": bool(solver.converged),
        "hf_s_squared": round(float(ss), 4),
        "hf_s_squared_ideal": round(spin / 2 * (spin / 2 + 1), 4),
        "n_correlated_electrons": nocc,
        "t1_diagnostic": round(t1_diag, 5),
        "d1_diagnostic": round(d1_diag, 5),
        "largest_t2_amplitude": round(largest_amplitude(t2), 4),
        "t1_exceeds_closed_shell_threshold": t1_diag > T1_CLOSED_SHELL_THRESHOLD,
        "d1_exceeds_closed_shell_threshold": d1_diag > D1_CLOSED_SHELL_THRESHOLD,
    }


def main():
    started = time.perf_counter()
    results = {}
    for fname, (charge, spin, role) in SPECIES.items():
        rec = run_species(fname, charge, spin)
        rec["role"] = role
        results[fname] = rec
        print(f"{fname:28s} T1={rec['t1_diagnostic']:.4f}  D1={rec['d1_diagnostic']:.4f}"
              f"  maxT2={rec['largest_t2_amplitude']:.3f}  S^2={rec['hf_s_squared']}  ({role})")

    ts = results["methane_ethynyl_ts.xyz"]
    controls = [results["methane.xyz"], results["acetylene.xyz"]]
    control_max_t1 = max(c["t1_diagnostic"] for c in controls)

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "compute_host": {"system": os.uname().sysname, "machine": os.uname().machine,
                          "loadavg": os.getloadavg()},
        "elapsed_seconds": round(time.perf_counter() - started, 1),
        "question": (
            "Is single-reference CCSD(T) trustworthy at the transition structure "
            "that the in-house +2.40 kcal/mol reference is evaluated on?"
        ),
        "thresholds": {
            "t1_closed_shell": T1_CLOSED_SHELL_THRESHOLD,
            "d1_closed_shell": D1_CLOSED_SHELL_THRESHOLD,
            "note": "Closed-shell rules of thumb (Lee-Taylor 1989; Janssen-Nielsen 1998). "
                    "Open-shell references sit higher intrinsically; read the TS relative to the closed-shell controls.",
        },
        "species": results,
        "reading": {
            "ts_t1_diagnostic": ts["t1_diagnostic"],
            "ts_d1_diagnostic": ts["d1_diagnostic"],
            "ts_largest_t2": ts["largest_t2_amplitude"],
            "closed_shell_control_max_t1": round(control_max_t1, 5),
            "ts_over_control_ratio": round(ts["t1_diagnostic"] / control_max_t1, 1) if control_max_t1 else None,
        },
        "limits": [
            "T1/D1 flag single-reference inadequacy; they do not by themselves prove a multireference wavefunction or give a corrected barrier.",
            "Diagnostics are basis- and geometry-specific (cc-pVDZ, the supplied collinear TS with three reported imaginary modes).",
            "A confirmed diagnostic makes the +2.40 reference a lower-confidence anchor, not a wrong number; establishing the true barrier needs a multireference method (CASPT2/NEVPT2/MRCI) or ROHF-based RCCSD(T) as the paper used.",
            "This does not change any DFT result; it changes the confidence attached to the yardstick DFT is measured against.",
        ],
    }
    (HERE / "reference-reliability.json").write_text(json.dumps(report, indent=1))
    print("\nTS T1 =", ts["t1_diagnostic"], "vs closed-shell control max", round(control_max_t1, 4),
          f"({report['reading']['ts_over_control_ratio']}x)")
    print("written:", (HERE / "reference-reliability.json").relative_to(ROOT))


if __name__ == "__main__":
    main()
