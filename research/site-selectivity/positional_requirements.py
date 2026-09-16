"""Turn a measured positional margin into a stiffness requirement on the mount.

`positional_uncertainty.py` answers "given this structure, how far does the tip
wander?".  This module answers the design question, which is the inverse: "given
that the margin is 2.495 Angstrom, how stiff must the mount be?".  A design tool
should state requirements, not only grade finished candidates.

Two results here are worth stating up front because they bear on the canonical
feasibility argument.

**1. Drexler's sigma = sqrt(k_B T / k) is a classical formula, and it
understates the spread wherever hbar*omega is not small against k_B T.**  The
quantum harmonic amplitude carries a factor sqrt(x coth x) with
x = hbar*omega / 2 k_B T.  Measured against a real DFT Hessian already in this
repository (H2 at PBE0-D3/def2-SVP, stretch 4383.9 cm^-1, in
data/validation/h2-integration/modes/), the bond-length spread is 87.3 mA
quantum against 26.9 mA classical at 298 K -- the classical value is 3.25 times
too small, and the quantum value is unchanged from 4 K to 298 K because it is
entirely zero-point motion.  So the classical formula does not merely miss a
correction there, it misses the whole effect.

**2. But that does not sink the argument, because a tip's positional spread is
dominated by its soft modes, and there the classical formula is fine.**  The
correction exceeds ten percent in sigma only above roughly 500 cm^-1 at room
temperature (`quantum_correction`, `crossover_wavenumber`).  A 10-100 N/m mount
mode moving tens of atomic masses sits near 100 cm^-1, where the correction is
under one percent.  The honest statement is therefore conditional: use the
classical formula for the soft mount modes that set the spread, and do not use
it for stiff bond modes, which contribute a zero-point floor it sets to nearly
zero.  This module computes both and reports which regime a mode is in, instead
of picking one and hoping.

**What the tail probability is and is not.**  The criterion used is crossing the
perpendicular bisector plane between the intended hydrogen and its nearest
competitor, which is a one-dimensional projection of the displacement, so the
probability is (1/2) erfc(d / (sigma sqrt(2))).  For a demanding target that
lands five or six sigma out, where a harmonic Gaussian is least trustworthy: the
real potential is not harmonic at such displacements, and the tail is the part of
a harmonic model to believe least.  **So the probability is not an error rate and
is not reported as one.**  What survives is the requirement it generates, and in
particular the *comparison* between a required stiffness and an achievable one.
That comparison is robust precisely when it is not close, which is the case worth
reporting.
"""

from __future__ import annotations

import math
from numbers import Real

from scipy.constants import Boltzmann, atomic_mass, electron_volt, hbar, speed_of_light
from scipy.optimize import brentq
from scipy.special import erfc, erfcinv

EV_PER_ANGSTROM2_TO_N_PER_M = electron_volt / 1e-20
_WAVENUMBER_TO_ANGULAR = 2.0 * math.pi * speed_of_light * 100.0


class PositionalRequirementError(ValueError):
    """The requested requirement is not well posed."""


def _positive(name: str, value) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value <= 0:
        raise PositionalRequirementError(f"{name} must be finite and positive")
    return float(value)


def _probability(name: str, value) -> float:
    value = float(value)
    if not math.isfinite(value) or not 0.0 < value < 0.5:
        raise PositionalRequirementError(f"{name} must lie strictly between 0 and 0.5")
    return value


def crossing_probability(sigma_angstrom: float, margin_angstrom: float) -> float:
    """One-sided probability of crossing the bisector plane at ``margin``.

    Harmonic equilibrium, one-dimensional projection. See the module docstring
    on why this is not an error rate.
    """
    sigma = _positive("sigma_angstrom", sigma_angstrom)
    margin = _positive("margin_angstrom", margin_angstrom)
    return 0.5 * float(erfc(margin / (sigma * math.sqrt(2.0))))


def required_sigma(margin_angstrom: float, target_probability: float) -> float:
    """Largest spread consistent with a target crossing probability."""
    margin = _positive("margin_angstrom", margin_angstrom)
    target = _probability("target_probability", target_probability)
    return margin / (math.sqrt(2.0) * float(erfcinv(2.0 * target)))


def quantum_correction(wavenumber_cm1: float, temperature_kelvin: float = 298.15) -> dict:
    """Ratio of quantum to classical amplitude for one harmonic mode.

    sigma_quantum / sigma_classical = sqrt(x coth x), x = hbar*omega / 2 k_B T.
    Tends to 1 for soft modes and grows as sqrt(x) for stiff ones.
    """
    wavenumber = _positive("wavenumber_cm1", wavenumber_cm1)
    temperature = _positive("temperature_kelvin", temperature_kelvin)
    omega = wavenumber * _WAVENUMBER_TO_ANGULAR
    x = hbar * omega / (2.0 * Boltzmann * temperature)
    ratio = math.sqrt(x / math.tanh(x))
    return {
        "wavenumber_cm1": wavenumber,
        "temperature_kelvin": temperature,
        "hbar_omega_over_2kT": x,
        "sigma_quantum_over_sigma_classical": ratio,
        "regime": (
            "zero-point dominated; the classical formula is not usable" if x > 3.0
            else "mixed; use the quantum amplitude" if x > 0.5
            else "near-classical; sqrt(kT/k) is adequate"
        ),
        "classical_formula_understates_sigma_by_percent": 100.0 * (ratio - 1.0),
    }


def crossover_wavenumber(temperature_kelvin: float = 298.15, tolerance: float = 0.10) -> dict:
    """Wavenumber at which the classical amplitude is wrong by ``tolerance``.

    Above this, sqrt(k_B T / k) understates sigma by more than the tolerance.
    """
    temperature = _positive("temperature_kelvin", temperature_kelvin)
    tolerance = _positive("tolerance", tolerance)
    target = 1.0 + tolerance

    def excess(wavenumber: float) -> float:
        return quantum_correction(wavenumber, temperature)["sigma_quantum_over_sigma_classical"] - target

    low, high = 1e-3, 1e6
    if excess(low) > 0 or excess(high) < 0:
        raise PositionalRequirementError("Crossover is not bracketed for this temperature")
    wavenumber = float(brentq(excess, low, high, xtol=1e-6, rtol=1e-12))
    return {
        "temperature_kelvin": temperature,
        "tolerance_fraction": tolerance,
        "crossover_wavenumber_cm1": wavenumber,
        "meaning": (
            "Above this wavenumber the classical sqrt(k_B T / k) amplitude understates the "
            f"root-mean-square displacement by more than {100.0 * tolerance:.0f} percent at this "
            "temperature. Below it, the classical formula is adequate."
        ),
    }


def required_stiffness(margin_angstrom: float, target_probability: float,
                       temperature_kelvin: float = 298.15, *,
                       effective_mass_amu: float | None = None) -> dict:
    """Stiffness the mount must supply to keep the tip inside ``margin``.

    The classical requirement k = k_B T / sigma**2 needs no mass.  The quantum
    requirement does, because the zero-point contribution depends on the mode's
    frequency and therefore on its effective mass; supply
    ``effective_mass_amu`` to get it.  For a soft mount mode the two coincide,
    and that agreement is itself the check that the classical answer is usable.
    """
    margin = _positive("margin_angstrom", margin_angstrom)
    target = _probability("target_probability", target_probability)
    temperature = _positive("temperature_kelvin", temperature_kelvin)

    sigma = required_sigma(margin, target)
    kt_ev = Boltzmann * temperature / electron_volt
    classical_ev = kt_ev / sigma**2
    result = {
        "margin_angstrom": margin,
        "target_crossing_probability": target,
        "temperature_kelvin": temperature,
        "required_sigma_angstrom": sigma,
        "margin_in_required_sigma": margin / sigma,
        "required_stiffness_classical_n_per_m": classical_ev * EV_PER_ANGSTROM2_TO_N_PER_M,
        "statistics_note": (
            "Classical equipartition. Valid when the mode setting the spread is soft; "
            "see quantum_correction and crossover_wavenumber."
        ),
        "not_an_error_rate": (
            "The target is a harmonic-equilibrium crossing probability for a one-dimensional "
            "projection, evaluated far into a Gaussian tail where a harmonic model is least "
            "reliable, and a driven assembly step does not sample an equilibrium distribution. "
            "Use the required-versus-achievable stiffness comparison, not this probability."
        ),
    }

    if effective_mass_amu is not None:
        mass = _positive("effective_mass_amu", effective_mass_amu) * atomic_mass
        target_variance = (sigma * 1e-10) ** 2

        def variance_gap(omega: float) -> float:
            x = hbar * omega / (2.0 * Boltzmann * temperature)
            variance = (hbar / (2.0 * mass * omega)) / math.tanh(x)
            return variance - target_variance

        # Variance falls monotonically with omega, so bracket and solve.
        low, high = 1e8, 1e18
        if variance_gap(low) < 0 or variance_gap(high) > 0:
            raise PositionalRequirementError(
                "No angular frequency reproduces the required spread for this mass; "
                "the zero-point amplitude alone may already exceed it"
            )
        omega = float(brentq(variance_gap, low, high, xtol=1e-6, rtol=1e-14))
        stiffness_si = mass * omega**2
        wavenumber = omega / _WAVENUMBER_TO_ANGULAR
        result.update(
            effective_mass_amu=effective_mass_amu,
            required_stiffness_quantum_n_per_m=stiffness_si,
            required_mode_wavenumber_cm1=wavenumber,
            quantum_over_classical_stiffness=(
                stiffness_si / result["required_stiffness_classical_n_per_m"]
            ),
            mode_regime=quantum_correction(wavenumber, temperature)["regime"],
        )
    return result


def assess(margin_angstrom: float, achievable_stiffness_n_per_m: float,
           target_probability: float, temperature_kelvin: float = 298.15) -> dict:
    """Compare a required stiffness against one a mount is believed to supply."""
    achievable = _positive("achievable_stiffness_n_per_m", achievable_stiffness_n_per_m)
    requirement = required_stiffness(margin_angstrom, target_probability, temperature_kelvin)
    required = requirement["required_stiffness_classical_n_per_m"]
    achieved_sigma = math.sqrt(
        (Boltzmann * temperature_kelvin / electron_volt)
        / (achievable / EV_PER_ANGSTROM2_TO_N_PER_M)
    )
    return {
        "requirement": requirement,
        "achievable_stiffness_n_per_m": achievable,
        "achievable_sigma_angstrom_classical": achieved_sigma,
        "stiffness_headroom_factor": achievable / required,
        "requirement_met": bool(achievable >= required),
        "margin_in_achievable_sigma": margin_angstrom / achieved_sigma,
        "verdict_basis": (
            "A necessary geometric condition on positioning only. Meeting it does not make the "
            "operation work: it says nothing about reaction barriers, competing chemistry, "
            "approach and withdrawal, or whether the mount is real. Rigid anchors make every "
            "computed stiffness an upper bound, so a real mount is softer than modelled."
        ),
    }
