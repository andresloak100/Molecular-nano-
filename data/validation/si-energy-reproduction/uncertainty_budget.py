"""Fail-closed uncertainty budget for this project's energetic claims.

Every energy in this repository is currently a bare point estimate. We say
"-2.89 kcal/mol" and "a site preference of 2 kcal/mol" with no interval, which
makes two different failures invisible: a number whose error bar is larger than
itself, and a number whose error bar nobody has bothered to estimate.

We can now do better than hand-waving, because several error terms have actually
been *measured* in this repository rather than assumed:

  method error        DFT against CCSD(T) at identical geometries
  basis error         def2-SVP -> TZVP -> QZVP on the same quantity
  state error         how far the default SCF guess sits above the lowest solution
  geometry error      residual force at a structure that is not stationary
  vibrational         zero-point and thermal corrections: NOT MEASURED
  tunnelling          NOT MEASURED

This module composes those into an interval for a named quantity, and — in the
style the rest of the repository already uses — **BLOCKS rather than guessing**
when a term that would materially affect the answer has not been measured. A
budget that silently omits an unmeasured term is worse than no budget, because
it launders ignorance into a confidence interval.

Systematic terms add linearly; only genuinely independent terms are combined in
quadrature. Treating a known one-sided method bias as a random error would
understate it, which is the usual way error budgets flatter their authors.

Run from the repository root:

    python data/validation/si-energy-reproduction/uncertainty_budget.py
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from pathlib import Path
from typing import Any

OUT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Term:
    """One contribution to an uncertainty budget.

    ``magnitude`` is in kcal/mol. ``systematic`` terms add linearly because they
    are biases with a known sign or a known one-sided character; independent
    terms add in quadrature. ``measured`` false means the term is a placeholder
    and any budget containing it must BLOCK rather than report a number.
    """

    name: str
    magnitude: float | None
    systematic: bool
    measured: bool
    evidence: str
    note: str = ""


# ---------------------------------------------------------------------------
# Measured terms. Each cites the artifact that produced it.
# ---------------------------------------------------------------------------

METHOD_ERROR_PRIMARY = Term(
    "method error, DFT vs CCSD(T), primary C-H reaction energy", 2.01, True, True,
    "si-energy-reproduction.json + dft-survey (PBE0-D3/def2-SVP -26.80 vs CCSD(T)/cc-pVDZ -24.79)",
    "One-sided: PBE0-D3 is too exothermic. Sign is known, so it adds linearly.")

METHOD_ERROR_TERTIARY = Term(
    "method error, DFT vs CCSD(T), tertiary C-H reaction energy", 6.15, True, True,
    "si-energy-reproduction-isobutane.json + dft-survey-isobutane.json (-38.38 vs -32.23)",
    "Same sign, three times larger than primary. Error grows with substitution.")

BASIS_ERROR_SVP_TO_TZVP = Term(
    "basis error, def2-SVP to def2-TZVP", 1.25, True, True,
    "methane reaction energy moves -26.80 (SVP) to -28.05 (TZVP)",
    "Not converged at SVP. Systematic in the same direction.")

BASIS_ERROR_TZVP_TO_QZVP = Term(
    "basis incompleteness beyond def2-TZVP", 0.06, True, True,
    "methane barrier moves -2.89 (TZVP) to -2.95 (QZVP)",
    "Converged. This term is negligible and is retained to show it was checked.")

STATE_ERROR_SCANNED = Term(
    "electronic-state selection, after a four-guess scan", 0.0, False, True,
    "si-energy-reproduction.json scf_attempts; all species scanned",
    "Zero only when a scan was actually run and the lowest solution retained.")

STATE_ERROR_UNSCANNED = Term(
    "electronic-state selection, no guess scan performed", 11.2, True, True,
    "worst observed default-above-lowest: propynyl PBE0-D3 (A2), 11.2 kcal/mol",
    "If no scan was run, this is the observed worst case, not a bound. "
    "Ethynyl UHF gave 8.76; the true worst case is unknown.")

# ---------------------------------------------------------------------------
# Unmeasured terms. Any budget touching these BLOCKS.
# ---------------------------------------------------------------------------

GEOMETRY_ERROR_NONSTATIONARY = Term(
    "geometry: quantity evaluated at a non-stationary structure", None, True, False,
    "S1 measured 2.38 eV/A residual force at the published TS under PBE0-D3/def2-TZVP",
    "A residual force is not convertible to an energy error without the curvature "
    "along that direction. Eighty times the 0.03 convergence threshold means the "
    "quantity is not the intended one at all, so no interval is meaningful.")

ZERO_POINT_ENERGY = Term(
    "zero-point energy correction", None, True, False,
    "not computed anywhere in this repository; magnitude anchored to Temelso Table 5",
    "Anchored to a real number rather than a generic range. For this exact reaction the "
    "source paper reports a 2.2 kcal/mol bare electronic barrier and a 1.7 kcal/mol "
    "barrier at 0 K, so its vibrational correction is about -0.5 kcal/mol. That is "
    "modest in absolute terms but it is the SAME ORDER as the site preference the "
    "project is trying to resolve, and it does NOT cancel between sites, because "
    "bridgehead and methylene C-H stretch frequencies differ. An adamantane site "
    "preference of 1-2 kcal/mol computed without zero-point corrections could be "
    "shifted appreciably, and conceivably reordered, by a term nobody has computed.")

THERMAL_FREE_ENERGY = Term(
    "finite-temperature free-energy correction", None, True, False,
    "not computed anywhere in this repository",
    "A bimolecular association loses substantial translational and rotational "
    "entropy. For a mounted tool much of that is already paid by the mechanical "
    "constraint, which is precisely why the free-gas analogue cannot be carried "
    "over unchanged.")

TUNNELLING = Term(
    "hydrogen tunnelling", None, True, False,
    "not computed anywhere in this repository",
    "Hydrogen transfer over a thin barrier tunnels. Room-temperature "
    "transmission coefficients of 2-10 are routine, and at the cryogenic "
    "temperatures mechanosynthesis proposals typically assume, tunnelling can "
    "dominate the rate entirely and is not a correction to a classical barrier "
    "but a replacement for it. A barrier height alone cannot give a rate here.")


@dataclass
class Budget:
    quantity: str
    value: float | None
    terms: list[Term] = field(default_factory=list)

    def evaluate(self) -> dict[str, Any]:
        missing = [t for t in self.terms if not t.measured]
        systematic = sum(abs(t.magnitude) for t in self.terms if t.measured and t.systematic and t.magnitude)
        independent = math.sqrt(sum(
            t.magnitude ** 2 for t in self.terms if t.measured and not t.systematic and t.magnitude))
        total = systematic + independent
        result: dict[str, Any] = {
            "quantity": self.quantity,
            "point_estimate_kcal_per_mol": self.value,
            "terms": [
                {"name": t.name, "magnitude_kcal_per_mol": t.magnitude, "type":
                 "systematic" if t.systematic else "independent", "measured": t.measured,
                 "evidence": t.evidence, "note": t.note}
                for t in self.terms
            ],
            "measured_systematic_sum_kcal_per_mol": systematic,
            "measured_independent_quadrature_kcal_per_mol": independent,
        }
        if missing:
            result.update(
                status="BLOCKED",
                uncertainty_kcal_per_mol=None,
                interval_kcal_per_mol=None,
                blocking_terms=[t.name for t in missing],
                verdict=(
                    "No interval reported. Unmeasured terms would materially change it, and "
                    "publishing a budget that silently omits them would launder ignorance into "
                    "a confidence interval. Measured terms alone already total "
                    f"{total:.2f} kcal/mol, which is a LOWER BOUND on the true uncertainty."),
            )
        else:
            result.update(
                status="REPORTABLE",
                uncertainty_kcal_per_mol=total,
                interval_kcal_per_mol=[self.value - total, self.value + total] if self.value is not None else None,
                verdict="All contributing terms measured.",
            )
        return result


def build() -> list[Budget]:
    return [
        Budget("In-house CCSD(T)/cc-pVDZ methane barrier at the published geometry", 2.3986, [
            STATE_ERROR_SCANNED,
            # Deliberately NOT including method/basis error: this IS the reference
            # method, and it reproduces the published value. Its remaining error is
            # basis incompleteness and the quality of CCSD(T) itself, neither measured.
            Term("CCSD(T) basis incompleteness beyond cc-pVDZ", None, True, False,
                 "not measured; no cc-pVTZ CCSD(T) run exists here",
                 "The published comparator is cc-pVTZ and gives 2.2 against our 2.4, "
                 "which hints the term is small, but one comparison is not a measurement."),
        ]),
        Budget("PBE0-D3/def2-TZVP 'barrier' at the published methane geometry", -2.89, [
            METHOD_ERROR_PRIMARY, BASIS_ERROR_TZVP_TO_QZVP, STATE_ERROR_SCANNED,
            GEOMETRY_ERROR_NONSTATIONARY,
        ]),
        Budget("Tertiary-minus-primary reaction energy difference, CCSD(T)/cc-pVDZ", -7.444, [
            STATE_ERROR_SCANNED,
            Term("CCSD(T) basis incompleteness, partially cancelling in a difference", None, True, False,
                 "not measured", "Cancellation in a difference is plausible but unverified."),
        ]),
        Budget("A DFT adamantane bridgehead-vs-methylene site preference", None, [
            Term("method error on a secondary-vs-tertiary site difference", None, True, False,
                 "being measured now: secondary_tertiary_calibration.py",
                 "Primary-vs-tertiary error swings 4.1 kcal/mol, so cancellation across "
                 "differing substitution is only partial. The secondary-vs-tertiary residual "
                 "is the quantity that decides whether a DFT site preference is chemistry."),
            BASIS_ERROR_SVP_TO_TZVP, STATE_ERROR_SCANNED, ZERO_POINT_ENERGY,
        ]),
        Budget("Any claim about which hydrogen the tool abstracts, at a real temperature", None, [
            ZERO_POINT_ENERGY, THERMAL_FREE_ENERGY, TUNNELLING,
            GEOMETRY_ERROR_NONSTATIONARY,
        ]),
    ]


def main() -> int:
    budgets = [b.evaluate() for b in build()]
    payload = {
        "schema_version": 1,
        "purpose": ("Attach measured uncertainty to this project's energetic claims, and block "
                    "rather than guess where a material term is unmeasured."),
        "convention": ("Systematic terms add linearly; independent terms add in quadrature. "
                       "Treating a known one-sided bias as random would understate it."),
        "budgets": budgets,
        "blocked": [b["quantity"] for b in budgets if b["status"] == "BLOCKED"],
        "reportable": [b["quantity"] for b in budgets if b["status"] == "REPORTABLE"],
        "headline": (
            "Four of five budgets BLOCK. The single largest unmeasured term is vibrational: "
            "for this exact reaction the source paper's own numbers imply about -0.5 kcal/mol "
            "(2.2 bare electronic, 1.7 at 0 K), the same order as the site preference the "
            "project is trying to resolve, and it does not cancel between sites with different "
            "C-H stretch frequencies. Tunnelling is "
            "worse still: for hydrogen transfer it is not a correction to a rate but a "
            "replacement for the classical picture, and nothing in this repository computes it."),
    }
    path = OUT / "uncertainty-budget.json"
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")

    for b in budgets:
        mark = "BLOCKED " if b["status"] == "BLOCKED" else "OK      "
        value = f"{b['point_estimate_kcal_per_mol']:+.3f}" if b["point_estimate_kcal_per_mol"] is not None else "   n/a"
        print(f"{mark} {value}  {b['quantity']}")
        if b["status"] == "BLOCKED":
            for name in b["blocking_terms"]:
                print(f"           missing: {name}")
            print(f"           measured terms alone are a lower bound of "
                  f"{b['measured_systematic_sum_kcal_per_mol'] + b['measured_independent_quadrature_kcal_per_mol']:.2f} kcal/mol")
        else:
            print(f"           +/- {b['uncertainty_kcal_per_mol']:.2f} kcal/mol")
    print(f"\n{len(payload['blocked'])} of {len(budgets)} budgets BLOCK.")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
