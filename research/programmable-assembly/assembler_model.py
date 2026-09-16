"""Programmable (template-directed) molecular assembly: the fidelity model.

Motivation
----------
A programmable assembler — a machine that reads a program and builds the
corresponding structure, as the ribosome reads mRNA and builds a protein — may
be a more tractable route to custom molecular structures than rigid positional
mechanosynthesis. This module does NOT claim to design such a machine. It
supplies the one quantitative fact that decides whether any such scheme can
work: **error propagation**.

An assembler that adds N building blocks, each with per-step fidelity f
(probability the correct block is added correctly), yields correct full-length
product with probability f**N. This compounding is unforgiving, and it is the
reason the ribosome invests so heavily in accuracy. The whole feasibility
question therefore reduces to a chemistry question this repository already
answers: the per-step fidelity is set by the free-energy discrimination
DeltaDeltaG-doub-dagger between the correct and the competing reaction, exactly
the quantity the site-selectivity lane (A1) and the quantum backend compute.

So this layer sits ON TOP of the existing engine: `nanodesign` validates a
single bond-forming/breaking step and returns a barrier discrimination; this
model turns that discrimination into a per-step fidelity and a whole-structure
yield, and tells you how much better the chemistry (or how much proofreading)
a target structure demands.

Everything here is standard transition-state-theory and kinetic-proofreading
arithmetic (Eyring; Hopfield 1974; Ninio 1975). Units: energies in kcal/mol,
temperature in kelvin. No claim of chemical accuracy is made for any particular
assembler; this is the arithmetic that a validated per-step number feeds into.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

R_KCAL_PER_MOL_K = 1.987204258640832e-3  # gas constant, kcal/(mol*K)
ROOM_TEMPERATURE_K = 298.15


def rt(temperature_k: float = ROOM_TEMPERATURE_K) -> float:
    """RT in kcal/mol. ~0.5925 kcal/mol at room temperature."""
    if temperature_k <= 0:
        raise ValueError("temperature must be positive")
    return R_KCAL_PER_MOL_K * temperature_k


def fidelity_from_discrimination(ddg_doubledagger_kcal: float,
                                 temperature_k: float = ROOM_TEMPERATURE_K,
                                 proofreading_stages: int = 0) -> float:
    """Per-step fidelity from the correct-vs-competing barrier gap.

    Kinetic branching under transition-state theory with a shared prefactor:
    the correct and incorrect pathways branch with rate ratio
    exp(DeltaDeltaG / RT), so the intrinsic single-check fidelity is

        f0 = 1 / (1 + exp(-DeltaDeltaG / RT)).

    Kinetic proofreading (Hopfield/Ninio) applies the same discrimination over
    ``proofreading_stages`` additional independent checks, each driven
    irreversibly, so the error is raised to the (stages + 1) power:

        error = (1 - f0) ** (stages + 1),   f = 1 - error.

    ``proofreading_stages=0`` is a single check (no proofreading). The ribosome
    is well described by ~1 extra stage. Proofreading is not free: each stage
    discards intermediates and must be paid for with free energy (see
    ``driving_note``); this model reports the accuracy, not that budget.
    """
    if proofreading_stages < 0:
        raise ValueError("proofreading_stages must be >= 0")
    x = ddg_doubledagger_kcal / rt(temperature_k)
    # numerically stable logistic
    f0 = 1.0 / (1.0 + math.exp(-x)) if x > -700 else 0.0
    error0 = 1.0 - f0
    error = error0 ** (proofreading_stages + 1)
    return 1.0 - error


def product_yield(n_steps: int, fidelity: float) -> float:
    """Fraction of assemblies that are correct over their whole length: f**N."""
    if n_steps < 0:
        raise ValueError("n_steps must be >= 0")
    if not 0.0 <= fidelity <= 1.0:
        raise ValueError("fidelity must be in [0, 1]")
    return fidelity ** n_steps


def required_fidelity(n_steps: int, target_yield: float) -> float:
    """Per-step fidelity needed to reach ``target_yield`` over ``n_steps``."""
    if n_steps <= 0:
        raise ValueError("n_steps must be >= 1")
    if not 0.0 < target_yield <= 1.0:
        raise ValueError("target_yield must be in (0, 1]")
    return target_yield ** (1.0 / n_steps)


def required_discrimination(n_steps: int, target_yield: float,
                            temperature_k: float = ROOM_TEMPERATURE_K,
                            proofreading_stages: int = 0) -> float:
    """Barrier discrimination (kcal/mol) needed for a target whole-structure yield.

    Inverts the fidelity chain: from the required per-step fidelity, back out
    the intrinsic single-check error, then the DeltaDeltaG that produces it.
    This is the number to hand back to the chemistry: "your competing-pathway
    barrier must be at least this far above the intended one."
    """
    f = required_fidelity(n_steps, target_yield)
    error = 1.0 - f
    error0 = error ** (1.0 / (proofreading_stages + 1))
    if not 0.0 < error0 < 1.0:
        raise ValueError("degenerate error; check inputs")
    # f0 = 1 - error0 ; DeltaDeltaG = RT * ln(f0 / (1 - f0))
    f0 = 1.0 - error0
    return rt(temperature_k) * math.log(f0 / error0)


@dataclass(frozen=True)
class AssemblyProgram:
    """An ordered list of building-block/operation identifiers — the 'mRNA'.

    ``blocks`` are opaque labels; this model is chemistry-agnostic and models
    only step count and per-step fidelity. It deliberately holds no biological
    sequence semantics: it is a materials-assembly abstraction, not a design
    tool for any biological agent (see README scope).
    """

    blocks: tuple

    @property
    def length(self) -> int:
        return len(self.blocks)


@dataclass(frozen=True)
class Assembler:
    """A programmable assembler characterized by its accuracy, not its structure.

    ``per_step_discrimination_kcal`` is the correct-vs-competing barrier gap for
    a single addition — the number a validated quantum calculation supplies.
    ``proofreading_stages`` is how many extra kinetic-proofreading checks the
    mechanism performs (0 = none; the ribosome ~1).
    """

    per_step_discrimination_kcal: float
    proofreading_stages: int = 0
    temperature_k: float = ROOM_TEMPERATURE_K

    def step_fidelity(self) -> float:
        return fidelity_from_discrimination(
            self.per_step_discrimination_kcal, self.temperature_k, self.proofreading_stages)

    def run(self, program: AssemblyProgram) -> dict:
        """Report the honest outcome of assembling ``program`` with this assembler."""
        f = self.step_fidelity()
        n = program.length
        y = product_yield(n, f)
        return {
            "n_steps": n,
            "per_step_discrimination_kcal_per_mol": self.per_step_discrimination_kcal,
            "proofreading_stages": self.proofreading_stages,
            "temperature_k": self.temperature_k,
            "per_step_fidelity": f,
            "per_step_error": 1.0 - f,
            "correct_full_length_yield": y,
            "expected_defects_per_assembly": n * (1.0 - f),
            "interpretation": (
                "correct_full_length_yield is the fraction of assemblies that are "
                "perfect over all steps (f**N). It is a necessary condition for a "
                "usable programmable assembler, not a validated device metric. "
                "per_step_discrimination must come from a validated calculation; "
                "an assumed value here proves nothing about real chemistry."
            ),
        }


def driving_note() -> str:
    return (
        "Fidelity is necessary but not sufficient. Each step must also be made "
        "effectively irreversible by coupling to a free-energy source (the "
        "ribosome uses aminoacyl-tRNA activation and GTP hydrolysis). "
        "Proofreading stages consume additional free energy per discarded "
        "intermediate. A complete feasibility case must budget that driving "
        "energy alongside the accuracy this module computes."
    )
