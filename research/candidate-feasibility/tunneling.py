"""A barrier is not a rate, and for hydrogen at these temperatures it is not close.

Everything in this repository that touches kinetics computes a barrier height.
The unstated next step is that a barrier height gives a rate through
exp(-E/kT). For a transferred *hydrogen atom* that step is wrong, and not by a
little.

Hydrogen is light enough to tunnel through a barrier rather than over it. The
size of that effect is set by the barrier's curvature at the top - the
imaginary frequency - not by its height. A tall narrow barrier can be crossed
far faster than a short wide one, which is the reverse of the classical
intuition the rest of the pipeline is built on.

THE NUMBER THAT MATTERS is the crossover temperature

    T_c = hbar |omega*| / (2 pi k_B) = 0.2290 * omega*[cm^-1]  kelvin

Above it, tunneling is a correction to classical transition-state theory.
Below it, tunneling is the dominant mechanism and a classical rate is
qualitatively wrong, not merely inaccurate.

WHAT THIS REACTION'S FREQUENCY ACTUALLY IS, and a correction to an earlier
version of this file. It previously said hydrogen-transfer saddles "typically"
run 1000-2000i cm^-1 and concluded that room temperature sits below crossover.
That was a generic expectation for the reaction class asserted as a measurement
of this reaction. The measurement exists, in SOURCE_NOTES.md line 33, from the
source paper's own Table 1: the collinear methane transition structure has
imaginary frequencies of 259i, 50i, 50i cm^-1.

    259i  ->  T_c =  59.3 K.  At 298 K, kappa = 1.068: a 7% correction.
   1000i  ->  T_c = 229.0 K.  At 298 K, kappa = 3.63.
   1500i  ->  T_c = 343.5 K.  At 298 K the parabolic form diverges.

So on the only datum this repository has, room temperature is about five times
ABOVE crossover and tunneling is a small correction, not the mechanism. The
arithmetic in the earlier version was right; the input was not. Reproducing a
calculation verifies the calculation, not the claim it rests on.

BUT THE QUESTION IS GENUINELY OPEN, and this file should not be read as settling
it in the other direction either. Two arguments from the same source point
opposite ways. The measured 259i comes from a structure carrying THREE imaginary
modes, which is not a verified first-order saddle, at a geometry none of this
project's functionals owns. Against that, the forensics lane finds that
reproducing the experimental apparent activation energy from the zero-point
corrected barrier requires substantial tunneling, around 1648i. Neither argument
is robust. The honest state is that nobody knows which regime this reaction is
in, and S1's refined saddle will produce the first imaginary frequency on a
functional's own surface.

That is why this module surveys a range rather than committing to a value, and
why the useful output is the MAP from omega* to regime rather than a single
kappa. omega* is the number that decides the regime; report it beside any
barrier.

WHAT THIS MODULE COMPUTES, and what it deliberately refuses to.

Implemented, because both have closed forms this file can state and test:

  Wigner       kappa = 1 + u^2/24,  u = hbar omega* / (k_B T). The leading
               correction. Cheap, and only trustworthy while it is small.
  Bell         kappa = (u/2) / sin(u/2), the parabolic-barrier result, of
               which Wigner is the small-u expansion. It diverges exactly at
               u = 2 pi, which is T = T_c - the divergence is the parabolic
               approximation announcing its own failure, not a physical
               infinity.

Deliberately NOT implemented: the Eckart correction. It is the right tool
below T_c and its closed form is long enough that writing it from memory would
risk a plausible, wrong number - the exact failure mode this project has been
burned by three times today. It is named as the required refinement instead,
along with the instanton/small-curvature methods that supersede it.

A CONSEQUENCE FOR SELECTIVITY that is easy to miss. Tunneling does not just
scale the desired rate; it reweights the competition. Two channels with the
same barrier height but different widths tunnel at different rates, so below
T_c the narrower barrier wins regardless of which is lower. Any selectivity
argument in this repository that rests on comparing barrier *heights* is
therefore incomplete until the imaginary frequencies are compared too.

A TESTABLE PREDICTION falls out for free. Deuterium is twice the mass, so its
imaginary frequency is smaller by roughly sqrt(2) and it tunnels much less.
The resulting kinetic isotope effect is large - well above the classical
semiclassical ceiling of about 7 at room temperature - when tunneling
dominates. That is an experimental signature someone could actually measure,
and it is the kind of falsifiable claim this project is short of.

Run from the repository root:

    python research/candidate-feasibility/tunneling.py
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(LANE))

EVIDENCE = LANE / "evidence"

HBAR_J_S = 1.054571817e-34
BOLTZMANN_J_PER_K = 1.380649e-23
SPEED_OF_LIGHT_CM_PER_S = 2.99792458e10
KCAL_PER_MOL_TO_J = 4184.0 / 6.02214076e23

# The only imaginary frequency this repository has for this reaction, from the
# source paper's Table 1 via SOURCE_NOTES.md line 33. The structure carries
# three imaginary modes, so it is a transition-structure candidate rather than
# a verified first-order saddle, and this number inherits that status.
MEASURED_WAVENUMBER_CM = 259.0
MEASURED_WAVENUMBER_PROVENANCE = (
    "Temelso, Sherrill, Merkle & Freitas (2006) Table 1, collinear methane "
    "transition structure at UCCSD(T)/cc-pVDZ; modes 259i, 50i, 50i cm^-1. "
    "Three imaginary modes means this is not a verified first-order saddle, "
    "and no functional in this project owns that geometry."
)
# A counter-argument from the forensics lane: reproducing the experimental
# apparent activation energy from the zero-point-corrected barrier needs far
# more tunneling than 259i supplies. Carried so the survey is not read as
# settling the question.
EXPERIMENT_IMPLIED_WAVENUMBER_CM = 1648.0

SURVEY_WAVENUMBERS_CM = (259.0, 500.0, 800.0, 1000.0, 1200.0, 1500.0, 1648.0, 2000.0)
SURVEY_TEMPERATURES_K = (4.0, 77.0, 195.0, 298.15, 400.0, 500.0)


def angular_frequency(wavenumber_cm: float) -> float:
    """Convert an imaginary frequency in cm^-1 to rad/s (magnitude)."""
    return 2.0 * math.pi * wavenumber_cm * SPEED_OF_LIGHT_CM_PER_S


def crossover_temperature(wavenumber_cm: float) -> float:
    """T_c = hbar |omega*| / (2 pi k_B). Below this, tunneling dominates."""
    return HBAR_J_S * angular_frequency(wavenumber_cm) / (2.0 * math.pi * BOLTZMANN_J_PER_K)


def _u(wavenumber_cm: float, temperature_k: float) -> float:
    return HBAR_J_S * angular_frequency(wavenumber_cm) / (BOLTZMANN_J_PER_K * temperature_k)


def wigner_kappa(wavenumber_cm: float, temperature_k: float) -> float:
    """Leading-order correction. Only meaningful while it stays near 1."""
    return 1.0 + _u(wavenumber_cm, temperature_k) ** 2 / 24.0


def bell_kappa(wavenumber_cm: float, temperature_k: float) -> dict:
    """Parabolic-barrier correction, valid only above the crossover temperature.

    kappa = (u/2)/sin(u/2) diverges at u = 2 pi. That is exactly T = T_c, so
    the formula fails precisely where tunneling stops being a correction.

    This is the quantity other lanes call ``kappa_p``; same closed form, and
    the name is kept aligned so four lanes do not quote two symbols for it.
    """
    u = _u(wavenumber_cm, temperature_k)
    half = u / 2.0
    if half >= math.pi:
        return {
            "kappa": None,
            "valid": False,
            "reason": (
                "At or below the crossover temperature the parabolic-barrier "
                "result diverges. Tunneling is the dominant mechanism here and "
                "needs Eckart, instanton or small-curvature treatment; a "
                "classical rate is qualitatively wrong, not just inaccurate."
            ),
            "u": u,
        }
    return {"kappa": half / math.sin(half), "valid": True, "u": u}


def effective_barrier_reduction_kcal(kappa: float, temperature_k: float) -> float:
    """The barrier error you would have to make to fake this factor classically.

    Expresses a rate enhancement as the barrier lowering that would produce it,
    so it can be compared against this project's own accuracy claims.
    """
    return BOLTZMANN_J_PER_K * temperature_k * math.log(kappa) / KCAL_PER_MOL_TO_J


def isotope_effect(wavenumber_cm: float, temperature_k: float, mass_ratio: float = 2.0) -> dict:
    """Tunneling contribution to the H/D kinetic isotope effect.

    The imaginary frequency of a transfer mode scales roughly as 1/sqrt(mass),
    so deuterium's is smaller by about sqrt(2) and it tunnels much less. This
    is the tunneling factor only - zero-point differences also contribute to a
    real KIE and are not included, so it is a lower bound on the total and is
    not a predicted experimental KIE on its own.
    """
    deuterium_wavenumber = wavenumber_cm / math.sqrt(mass_ratio)
    light = bell_kappa(wavenumber_cm, temperature_k)
    heavy = bell_kappa(deuterium_wavenumber, temperature_k)
    ratio = None
    if light["valid"] and heavy["valid"] and heavy["kappa"]:
        ratio = light["kappa"] / heavy["kappa"]
    return {
        "hydrogen_wavenumber_cm": wavenumber_cm,
        "deuterium_wavenumber_cm": deuterium_wavenumber,
        "hydrogen_kappa": light.get("kappa"),
        "deuterium_kappa": heavy.get("kappa"),
        "tunneling_contribution_to_kie": ratio,
        "scope": (
            "Tunneling contribution only; zero-point differences are excluded, "
            "so this understates a measured KIE and is not one."
        ),
    }


def survey() -> dict:
    """Crossover temperatures and correction factors across plausible barriers."""
    rows = []
    for wavenumber in SURVEY_WAVENUMBERS_CM:
        t_c = crossover_temperature(wavenumber)
        temperatures = {}
        for temperature in SURVEY_TEMPERATURES_K:
            bell = bell_kappa(wavenumber, temperature)
            entry = {
                "temperature_k": temperature,
                "above_crossover": temperature > t_c,
                "wigner_kappa": wigner_kappa(wavenumber, temperature),
                "bell_kappa": bell.get("kappa"),
                "bell_valid": bell["valid"],
            }
            if bell["valid"] and bell["kappa"]:
                entry["equivalent_barrier_lowering_kcal_per_mol"] = (
                    effective_barrier_reduction_kcal(bell["kappa"], temperature)
                )
            else:
                entry["note"] = bell["reason"]
            temperatures[f"{temperature:g}K"] = entry
        rows.append({
            "imaginary_wavenumber_cm": wavenumber,
            "crossover_temperature_k": t_c,
            "room_temperature_is_below_crossover": 298.15 < t_c,
            "temperatures": temperatures,
            "isotope_effect_298K": isotope_effect(wavenumber, 298.15),
        })
    return {"rows": rows}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="tunneling-survey.json")
    args = parser.parse_args()

    EVIDENCE.mkdir(parents=True, exist_ok=True)
    out = EVIDENCE / args.out
    if out.exists():
        raise SystemExit(f"{out} exists; evidence is never overwritten. Choose --out.")

    data = survey()
    report = {
        "schema_version": 1,
        "question": "Is a classical barrier enough to get a rate for this reaction?",
        "answer_shape": (
            "Unresolved, and the resolution is one number. On the only measured "
            "imaginary frequency available (259i cm^-1) room temperature is five "
            "times above crossover and tunneling is a 7% correction. On the "
            "frequency implied by reproducing the experimental activation energy "
            "(~1648i) it is the dominant mechanism. Neither input is trustworthy: "
            "the first comes from a structure with three imaginary modes, the "
            "second is inferred rather than computed."
        ),
        "computes_chemistry": False,
        "inputs": "Imaginary frequencies are surveyed; the repository has no verified saddle of its own.",
        "measured_anchor": {
            "wavenumber_cm": MEASURED_WAVENUMBER_CM,
            "crossover_temperature_k": crossover_temperature(MEASURED_WAVENUMBER_CM),
            "kappa_at_298K": bell_kappa(MEASURED_WAVENUMBER_CM, 298.15).get("kappa"),
            "provenance": MEASURED_WAVENUMBER_PROVENANCE,
        },
        "experiment_implied_anchor": {
            "wavenumber_cm": EXPERIMENT_IMPLIED_WAVENUMBER_CM,
            "crossover_temperature_k": crossover_temperature(EXPERIMENT_IMPLIED_WAVENUMBER_CM),
            "basis": (
                "Forensics lane: the tunneling needed to reproduce the measured "
                "apparent activation energy from the zero-point-corrected "
                "barrier. Inferred from a rate, not computed from a surface."
            ),
        },
        "correction_history": (
            "An earlier version of this module asserted 1000-2000i as typical of "
            "the reaction class and concluded room temperature was below "
            "crossover. The arithmetic was right and the input was not. "
            "Corrected against SOURCE_NOTES.md line 33."
        ),
        "implemented": ["Wigner leading-order", "Bell parabolic-barrier", "crossover temperature", "tunneling contribution to H/D KIE"],
        "deliberately_not_implemented": {
            "eckart": (
                "Correct below the crossover temperature, but its closed form is "
                "long enough that reproducing it from memory risks a plausible "
                "wrong number. Named as the required refinement instead."
            ),
            "instanton_and_small_curvature": "Supersede Eckart below T_c; not implemented.",
            "zero_point_energy": "Changes the effective barrier independently of tunneling; not included here.",
        },
        "survey": data,
        "consequence_for_selectivity": (
            "Tunneling reweights competing channels by barrier WIDTH, not only "
            "height. Any selectivity argument resting on comparing barrier "
            "heights is incomplete until imaginary frequencies are compared."
        ),
    }
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    print("crossover temperature by imaginary frequency:\n")
    print(f"{'omega* (cm^-1)':>15} {'T_c (K)':>10} {'298 K below T_c?':>18}")
    for row in data["rows"]:
        print(f"{row['imaginary_wavenumber_cm']:15.0f} {row['crossover_temperature_k']:10.1f} "
              f"{str(row['room_temperature_is_below_crossover']):>18}")

    print("\ntunneling factor kappa (Bell), and the barrier error that would fake it:\n")
    print(f"{'omega*':>8} {'T (K)':>8} {'kappa':>12} {'= barrier lowering':>22}")
    for row in data["rows"]:
        if row["imaginary_wavenumber_cm"] not in (1000.0, 1500.0):
            continue
        for name, entry in row["temperatures"].items():
            if entry["bell_valid"]:
                print(f"{row['imaginary_wavenumber_cm']:8.0f} {entry['temperature_k']:8.1f} "
                      f"{entry['bell_kappa']:12.2f} "
                      f"{entry['equivalent_barrier_lowering_kcal_per_mol']:18.2f} kcal/mol")
            else:
                print(f"{row['imaginary_wavenumber_cm']:8.0f} {entry['temperature_k']:8.1f} "
                      f"{'DIVERGES':>12} {'below crossover':>22}")

    print("\ntunneling contribution to the H/D kinetic isotope effect at 298 K:\n")
    for row in data["rows"]:
        kie = row["isotope_effect_298K"]
        value = kie["tunneling_contribution_to_kie"]
        print(f"  omega* {row['imaginary_wavenumber_cm']:6.0f} cm^-1 -> "
              + (f"KIE_tunnel = {value:6.2f}" if value else "below crossover, not computable here"))

    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
