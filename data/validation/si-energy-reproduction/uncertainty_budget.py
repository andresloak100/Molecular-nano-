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
    "not computed; and the input it depends on, |nu| at a verified saddle, is unmeasured",
    "Hydrogen transfer over a thin barrier tunnels, and below the crossover "
    "temperature T_c = h c |nu| / (2 pi k_B) it is not a correction to a classical "
    "rate but a replacement for it. THE SIGN OF THE CONCLUSION DEPENDS ENTIRELY ON "
    "|nu|, WHICH THIS PROJECT HAS NOT MEASURED. Assuming a typical H-transfer value "
    "of 1500-2000i cm^-1 gives T_c = 344-458 K, so room temperature sits BELOW "
    "crossover and tunnelling dominates. But the only |nu| we actually have for this "
    "reaction is 259i cm^-1 (Temelso Table 1), giving T_c = 59 K, which would put "
    "room temperature far ABOVE crossover and make tunnelling a modest correction. "
    "A 2.2 kcal/mol barrier with a loose early transition structure (acceptor-H "
    "1.67 A) is physically consistent with shallow negative curvature, so 259i is "
    "not obviously wrong - but that structure carries three imaginary modes and is "
    "not a first-order saddle, so its 259i is not necessarily the reaction-coordinate "
    "frequency either. S1's verified saddles settle this directly. "
    "EXPERIMENT DOES NOT RESOLVE IT EITHER, AND MY EARLIER CLAIM THAT IT FAVOURED THE "
    "LOW BRANCH IS WITHDRAWN. Two arguments can be built from Opansky & Leone 1996 "
    "(C2H + CH4, 154-359 K, k = 1.2e-11 exp(-491/T)) and they point in OPPOSITE "
    "directions. (a) Form of fit: a Wigner factor varies 1.2-fold across that range at "
    "259i but 3.7-fold at 1500i, and a clean single exponential accommodates the former, "
    "not the latter - favours LOW |nu|. (b) Magnitude of the fitted value: starting from "
    "the 1.70 kcal/mol zero-point-corrected barrier, reproducing the measured apparent Ea "
    "of 0.976 REQUIRES substantial tunnelling, around 1000-1800i with no prefactor "
    "temperature dependence, and is unreachable within Wigner once a T^1 prefactor is "
    "included - favours HIGH |nu|. Neither is robust: (a) infers curvature from the "
    "authors' reporting choice rather than from residuals, and was transcribed from an "
    "abstract; (b) is acutely sensitive to the assumed prefactor exponent and to the true "
    "barrier, and uses Wigner outside its validity range at the cold end. The two "
    "arguments do NOT simply cancel, because (b) turns out to be self-undermining. "
    "The parabolic-barrier transmission diverges at u = 2 pi, and u = 2 pi is exactly "
    "T = T_c, so divergence and crossover are the same point. Argument (b) needs "
    "|nu| ~ 1648i, which puts T_c at 377 K and places the ENTIRE measured 154-359 K "
    "range below crossover, in the deep-tunnelling regime where apparent activation "
    "energy falls toward zero as T drops. That predicts a curved Arrhenius plot "
    "flattening at low temperature, which is not the clean single exponential the "
    "measurement reports. (b) therefore contradicts the data it was fitted to and is "
    "withdrawn as internally inconsistent, not merely imprecise. This does NOT repair "
    "(a), whose defect is independent and untouched. Net: weak support for the low "
    "branch, with a conditional bound of |nu| < 673 cm^-1 if the fit really is clean "
    "across the whole range, and a residual puzzle - the 1.70 to 0.976 gap still needs "
    "a non-tunnelling explanation, most plausibly thermal averaging over the reactant "
    "distribution, which generically depresses apparent Ea below a 0 K barrier for a "
    "low-barrier reaction. Status remains OPEN pending a verified |nu|. ""depend on its precision. The term stays unmeasured until a verified saddle supplies "
    "|nu| directly; it is no longer symmetric between the two branches.")


CURVATURE_PRECISION = Term(
    "precision of |nu| itself: single-step finite difference, unconverged", None, True, False,
    "stationary.py reports finite_difference_step_convergence_checked false; raised by S1",
    "Raised by S1 and it is the term that will outlive the other three. The crossover "
    "test needs |nu|, but |nu| is currently one finite-difference estimate with no step-size "
    "convergence check. S1 has queued a second Hessian at a different step. How much this "
    "matters depends entirely on where |nu| lands: at 259i the room-temperature boundary of "
    "1301 cm^-1 is a factor of five away, so the verdict survives being wrong by tens of "
    "percent; near 1301 a ten percent error flips which rate model applies.")

QUANTUM_AMPLITUDE_NOTE = (
    "Verified cross-lane result, A1 proposed and this lane confirmed. Positional spread uses "
    "sigma_quantum/sigma_classical = sqrt(x coth x) with x = hbar*omega/2kT, not Drexler's "
    "classical sqrt(kT/k). On the committed H2 Hessian at 4383.9 cm^-1 the ratio is 3.253, and "
    "the quantum spread is 87.51 mA at 298 K, 77 K and 4 K ALIKE, because it is pure zero-point "
    "motion; the classical formula gives 26.90, 13.67 and 3.12, wrong by a factor of 28 at 4 K. "
    "The ten-percent crossover frequency scales linearly with temperature: 336 cm^-1 at 298 K, "
    "87 at 77 K, 4.5 at 4 K, so cooling makes MORE modes quantum. For this design that is "
    "nonetheless favourable: a 30 N/m mount mode near 101 cm^-1 saturates low, and the margin "
    "widens from 21 spreads at 298 K to 43 at 4 K. Use the classical form only for the soft "
    "mount modes that set the tip spread, never for stiff bond modes.")


def crossover_temperature_kelvin(imaginary_frequency_cm: float) -> float:
    """T_c = h c |nu| / (2 pi k_B). Below T_c, tunnelling dominates the rate.

    Derived independently here and cross-checked against support session
    76190bf3's rate-physics lane; both give 343.5 K at 1500 cm^-1 and 458.0 K at
    2000 cm^-1. The formula is cheap and exact; the difficulty is entirely in
    obtaining a trustworthy |nu|.
    """
    planck, light_speed, boltzmann = 6.62607015e-34, 2.99792458e10, 1.380649e-23
    return planck * light_speed * float(imaginary_frequency_cm) / (2 * math.pi * boltzmann)


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
            GEOMETRY_ERROR_NONSTATIONARY, CURVATURE_PRECISION,
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
    payload["tunnelling_crossover"] = {
        "formula": "T_c = h c |nu| / (2 pi k_B); below T_c tunnelling dominates the rate",
        "kelvin_by_imaginary_frequency_cm": {
            str(nu): crossover_temperature_kelvin(nu) for nu in (259, 1000, 1500, 2000, 2500)
        },
        "room_temperature_kelvin": 298.0,
        "status": "UNRESOLVED",
        "why": (
            "Assuming a typical H-transfer |nu| of 1500-2000i puts room temperature BELOW "
            "crossover and tunnelling dominant. The only |nu| this project actually has for "
            "the reaction, 259i from Temelso Table 1, puts room temperature far ABOVE "
            "crossover and tunnelling minor. Those are opposite conclusions and the input "
            "distinguishing them is unmeasured. S1's verified saddles supply it."),
        "cross_checked_with": (
            "support session 76190bf3 rate-physics lane and S1; all three agree. S1 expresses the "
            "same boundary as a frequency, |nu| = 2 pi k_B T / (h c) = 1301 cm^-1 at 298 K, which "
            "reproduces here to 1301.4."),
        "boundary_frequency_cm_at_298K": 1301.4,
        "experimental_constraint": {
            "source": "Opansky & Leone 1996, C2H + CH4, 154-359 K, k = 1.2e-11 exp(-491/T)",
            "apparent_activation_energy_kcal_per_mol": 0.976,
            "wigner_kappa_variation_across_measured_range": {
                "259": 1.2, "1000": 2.8, "1301": 3.4, "1500": 3.7, "2000": 4.2},
            "RETRACTION": (
                "An earlier version of this file claimed the Arrhenius linearity favoured the low "
                "branch. That is withdrawn. A second argument from the same measurement, based on "
                "the magnitude of the fitted activation energy rather than the form of the fit, "
                "points the opposite way: reproducing +0.976 from a 1.70 zero-point-corrected "
                "barrier requires substantial tunnelling and therefore HIGH |nu|. Neither argument "
                "is robust and they do not compose. The question is open."),
            "magnitude_argument_required_nu_cm": {
                "Ea0=1.4, no prefactor T-dependence": 747,
                "Ea0=1.7, no prefactor T-dependence": 1648,
                "with a T^1 prefactor": "unreachable within Wigner"},
            "inference_form_of_fit_ONLY": (
                "A single exponential fitted cleanly over a 2.3-fold temperature span is "
                "compatible with the 1.2-fold Wigner variation at 259i and not with the "
                "3.7-fold variation at 1500i, which would curve the Arrhenius plot visibly. "
                "Evidence for the low branch, available without waiting for S1."),
            "consistency_chain": (
                "electronic +2.399 (in-house CCSD(T)) -> 0 K with ZPE +1.7 (Temelso Table 5) "
                "-> apparent Ea +0.976 (experiment). Monotone decreasing, the expected ordering."),
            "strength": "qualitative; Wigner is unreliable above kappa ~2, so the high-|nu| column is indicative only",
        },
    }
    payload["positional_spread_quantum_correction"] = QUANTUM_AMPLITUDE_NOTE
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
