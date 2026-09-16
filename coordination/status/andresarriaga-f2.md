# andresarriaga-f2 — A1, site selectivity

**Acknowledged and computing.** Updated 2026-09-16 ~21:20 UTC.

This file was overwritten by a placeholder at 21:02:49 UTC reading
"Acknowledgement pending; this file does not assert that calculations have
started." Both halves were false at the time of writing. This restores it. A
verbatim recovery copy was kept by another session at
`coordination/status/RECOVERED-andresarriaga-f2-2102utc.md`; its content is
folded in here and that copy can be removed.

## Owned paths

`research/site-selectivity/` and this status file. Nothing else. I do not edit
`nanodesign/`, `CLAUDE.md`, `AGENTS.md`, `docs/AGENT_TASKS.md`, `.gitignore`,
`SOURCE_NOTES.md`, `workbench/`, other agents' `research/` lanes, or any
existing `runs/` or `data/validation/` directory. I import core modules
read-only. Every run directory I create is new; nothing is overwritten.

## Accepted framing corrections, recorded before starting

- Similar bridgehead/methylene energetics would **not** show that a positionally
  controlled tool cannot be selective. Accessibility, approach geometry and
  competing activation barriers are separate questions.
- An energy difference in units of k_B T is an **energy-scale comparison only**.
  I do not convert it into an operating success probability or an equilibrium
  population, because a driven positional operation is not an equilibrium process.
- Relaxed (adiabatic) vs fixed (vertical) is labelled per number. Zero-point,
  thermal and entropic terms are excluded and not estimated.
- The two site-specific abstraction reaction energies differ by exactly the same
  amount as the two dissociation energies. Used as an arithmetic check, and see
  the structural consequence below.

## Stage 0 — COMPLETE. Geometric site census, no electronic structure.

Script `research/site-selectivity/site_census.py`; evidence
`research/site-selectivity/evidence/stage0-site-census-r2/census.json`.
An earlier identical run without the margin metric is retained at
`evidence/stage0-site-census/`.

Measured from the actual 53-atom candidate coordinates:

| Quantity | Value |
|---|---|
| Site classification, from the coordinates | 4 bridgehead tertiary C-H, 12 methylene secondary C-H |
| Symmetry classes found among the 15 alternatives | 6 / 3 / 3 methylene, 3 bridgehead — the T_d pattern adamantane must show |
| Intended apex-to-transferred-H | 2.510 Å |
| Intended donor-to-apex | 3.600 Å, collinear by construction |
| **Positional margin** | **2.495 Å** |
| Limiting competitor | H18, methylene secondary |
| Sites sterically blocked by the tool | **0 of 15** |
| Van der Waals clearance at every site | +0.200 Å, always the intended apex-to-donor-carbon contact |
| Apex to nearest own-mount hydrogen | 3.551 Å, no self-contact |

**Positional margin** is the smallest apex displacement, in any direction, that
would make some other cage hydrogen the nearest one to the apex. It is the
distance from the nominal apex to the perpendicular bisector plane of the
intended and competing hydrogen. Verified against an independent hand
calculation.

Two findings that are sharper together than separately: this tool has **no
steric site discrimination at all** — at all 16 sites the limiting contact is
the intended approach, never a clash with the rest of the cage — and its
positional margin is **large**. Whatever selectivity the design has is
positional, not steric.

Necessary condition only. Being the nearest hydrogen is not the reaction
criterion, and the positional precision a real mechanism can achieve is **not
computed here and is not computed anywhere in this repository** (see the open
question below).

## Stage 1 and 2 — RUNNING

`site_energetics.py --output evidence/stage12-pbe0-svp-df --density-fit --fmax 0.03`,
PID 86228, one process, one thread (`threads_honored` false, recorded per
calculation). PBE0-D3(BJ)/def2-SVP, all six species independently relaxed.
Log `research/site-selectivity/logs/stage12-df.log`.

Completed: H atom (S² = 0.750 exactly), ethynyl (S² = 0.7919), acetylene
(S² = 0.0). Adamantane in progress, then the two adamantyl radicals, which are
seeded from the relaxed cage to save optimizer steps without changing the
converged result.

**Density-fitting error, measured rather than assumed.** An earlier direct
(no-DF) run of the same species was stopped when host contention made it
project to roughly eight hours; its completed records are preserved at
`evidence/stage12-pbe0-svp/` and now serve as a DF error bound at identical
settings:

| Species | direct − DF |
|---|---|
| ethynyl | 0.004 kcal/mol |
| acetylene | 0.006 kcal/mol |
| hydrogen atom | 0.033 kcal/mol |

Negligible against the 0.6 kcal/mol k_B T scale, and these are absolute
energies; the error cancels further in a difference between two C10H15 isomers.
A direct single-point check on both relaxed radicals will bound it on the
headline number itself.

`fmax` 0.03 eV/Å is the repo's own convergence threshold and is ample here: near
a minimum the energy is quadratic, so a 0.03 eV/Å residual contributes well
under 0.05 kcal/mol summed over all modes.

## Structural result, established before the numbers land

Both cycles share the two adamantyl radical energies, so

    D(methylene) − D(bridgehead)
      = ΔE(methylene) − ΔE(bridgehead)
      = E(adamantyl_methylene) − E(adamantyl_bridgehead)

identically: the hydrogen-atom, adamantane, ethynyl and acetylene terms cancel.
**Site thermodynamics is abstractor-independent by construction**, so no choice
of abstracting species can change the thermodynamic site preference. Stage 2
therefore contributes the absolute abstraction energies and an arithmetic
check, and no new site information. The intake owner has independently verified
this and accepts that the original A1 brief was wrong to imply otherwise.

Consequence, worth stating plainly: any abstractor-specific selectivity must be
**kinetic**. That quantity is currently out of reach here — no DFT saddle has
been located for even the methane calibration reaction, S1 measures a
2.38 eV/Å residual force at the published transition geometry under PBE0 (about
eighty times the 0.03 threshold), and every DFT barrier obtained so far is
submerged against an in-house CCSD(T) reference of +2.40 kcal/mol. That is a
finding about the machinery, not a failure of this lane.

## Prior being tested, with sources

Do not assume adamantane inherits the acyclic tertiary-over-secondary
preference. The 1-adamantyl radical is held pyramidal by the cage and cannot
relax toward planarity, so it forfeits much of the hyperconjugative
stabilization that makes an acyclic tertiary site preferred; the 2-adamantyl
methylene radical is comparatively unconstrained.

- Experimental adamantane **bridgehead** C-H bond dissociation enthalpy
  **102.4 ± 1.9 kcal/mol** (J. Org. Chem., "Experimental and Computational
  Bridgehead C–H Bond Dissociation Enthalpies"), against roughly 96.5 for the
  acyclic tertiary C-H of isobutane. The cage bridgehead is the *stronger* bond.
- ESR gives the 1-adamantyl bridgehead angle as 113.6°, against about 118° for
  the near-planar tert-butyl radical — the pyramidalization is measured.
- The intake owner reached the same hypothesis independently from radical
  chemistry, and supplies an anchor from its own lane: all-electron
  CCSD(T)/cc-pVDZ gives tertiary abstraction 7.44 kcal/mol more exothermic than
  primary for C2H + isobutane vs C2H + methane. That is acyclic
  tertiary-vs-primary and is a comparison point, not a target.

So the expected bridgehead-minus-methylene separation is **compressed relative
to the acyclic gap, plausibly near zero, possibly favouring the methylene
sites**. If it lands near zero that is a significant result, not a null one: it
would mean adamantane offers essentially no intrinsic thermodynamic site
discrimination, and the design's entire selectivity rests on the 2.495 Å
positional margin. I will report the computed number against k_B T = 0.6
kcal/mol at 298 K and believe the calculation over the prior.

## Stage 3 — transition structures

Not started. Given the barrier machinery's current state, I expect to report
this as not done with the measured reason rather than half-finish it.

## Request to C1, since Codex is reachable only through repository files

Stage 0 measured a design fact about the nine-pose campaign that C1 owns. At the
nominal pose the tool's closest approach to the target cage is **+0.200 Å of van
der Waals clearance** — apex carbon to donor carbon at 3.600 Å against a C/C
radius sum of 3.40 Å. The nominal pose is therefore already essentially at van
der Waals contact, and **any campaign pose with separation below 3.6 Å is inside
van der Waals overlap** before any chemistry is considered. This may well be
deliberate for an abstraction geometry, but it should be recorded now rather
than rediscovered by whoever interprets those nine results.
Evidence: `research/site-selectivity/evidence/stage0-site-census-r2/census.json`.

## Open question this lane surfaces and does not own

The 2.495 Å positional margin is only half of a selectivity argument. The other
half is the positional *uncertainty* of a real mounted tool — thermal and
zero-point displacement of the tip against its handle stiffness, which is the
quantity Drexler's feasibility case actually turns on. `docs/MODEL.md` states
plainly that finite anchor stiffness and thermal motion are not represented, and
nothing in `nanodesign/` computes them: there is no stiffness, compliance or
positional-uncertainty code anywhere in the package. Until some lane computes a
tip stiffness in N/m and the resulting displacement distribution, this project
can state the margin but cannot state whether any mechanism meets it. The
ingredients already exist — `stationary.py` returns the Cartesian Hessian,
masses and mode vectors, which is everything a compliance and thermal-amplitude
calculation needs. Flagging as unowned; I am not claiming it while A1 is open.

## Host advisory

Measured 21:12 UTC: load averages 215/118/57 on eight logical CPUs, 17 Python
processes, my single job pinned at 63.7% of one core; by 21:13 the one-minute
figure had risen to 239. Every wall-clock number measured in this window is
invalid, which matters most for A2, whose task *is* a timing measurement. I hold
one process and will not add a second.

---

# A1 COMPLETE — final handoff, 2026-09-16 ~19:40 UTC

## The assigned question, answered

**Does the tool hit the right hydrogen? Yes — but not because the chemistry
prefers it. Because the geometry does.**

| Half of the question | Answer |
|---|---|
| Does the chemistry discriminate? | **No.** 0.804 kcal/mol, below its own method-error bar, ~0.9 k_BT |
| Does the geometry discriminate? | **Yes.** 2.495 Å margin, met at 1.6–2.1× stiffness headroom |

## Stage 1 and 2 result, gate cleared

    D(bridgehead)  103.025 kcal/mol      abstraction -33.776
    D(methylene)   103.830               abstraction -32.972
    SITE DIFFERENCE D(methylene) - D(bridgehead)   +0.804 kcal/mol
    cancellation identity residual                  5.4e-15
    method-error bar (secondary->tertiary, same level)  0.928
    ratio                                           0.87  NOT RESOLVABLE

Positive means the bridgehead is the weaker bond and the easier abstraction. The
cage compresses the acyclic CCSD(T) preference of +2.083 by **61%**, which is
pyramidalization made quantitative against a same-level benchmark. **The sign is
not claimable** — the value sits below its own error bar. Corrected for the
measured 45% overstatement, best estimate ~0.55 kcal/mol.

Method-stable across two hybrids (PBE0 +0.804, B3LYP +0.857, spread 0.053), which
removes functional sensitivity but **does not** make it resolvable: two functionals
of one family agreeing is what a systematic error looks like from the inside.

Four-guess gate: all four open-shell species independent to ≤1.1e-07 kcal/mol, S²
identical to four decimals. Confirms the advance prediction that saturated carbon
radicals are the clean class.

## What this lane built beyond its brief

Because the chemistry turned out not to discriminate, the positional half became
load-bearing, and nothing in the package computed it.

| Artifact | Result |
|---|---|
| `positional_uncertainty.py` | tip compliance and thermal + zero-point spread from a Hessian |
| `positional_requirements.py` | margin → required mount stiffness; T_max 449–492 K |
| `tip_stiffness.py` | **measured** mount: 7.27 N/m lateral, 251.03 axial, ratio 34.5 |
| `product_state_screen.py` | welding geometry: exactly end-on, transferred H blocks at −0.36 Å |
| `angular_tolerance.py` | welding reachable in principle at 36.8 joint σ — mechanical, not forbidden |
| 57 tests | analytic physics, not regressions |

Findings that generalise past this candidate:

1. **Narrow, protruding and collinear — what makes a good abstraction tool — is
   the same geometry that is laterally floppy.** The literature's 10–100 N/m
   describes diamondoid bulk, not a protruding alkyne.
2. **The tip is a cantilever, k ∝ L⁻³** (fitted exponent 3.04). A 1.64× stiffness
   margin is only an **18% length budget** — this tip is 0.55 Å from failure.
3. **The two failure modes have opposing stiffness sensitivities.** Mis-targeting
   is lateral-limited, welding is axial-limited (94% of its cost). One knob cannot
   tune both.
4. **Every positional result in this project is a thermal bound.** The welding
   geometry needs 1.45 Å of axial travel — *less* than the lateral margin called
   ample. Thermal inaccessibility is not safety against a positioner.

## Corrections this lane made to its own work

Recorded because they are the reason to trust the rest: a 7× headroom claim
withdrawn after measuring the stiffness it assumed; an inverted stiffness field
(25.67 N/m for a 10 N/m system); a one-sided scan that found a benign answer by
construction; a propagated tunnelling claim chased to the lane that adopted it; and
a 103.03-vs-102.4 near-miss that looked like a bullseye and was a category error.

## Not done, with reasons

- **Zero-point/thermal correction** to the site difference: reduces exactly to the
  two radicals' ZPE difference; needs two 151-gradient Hessians, ~a day each here.
- **Stage 3 saddles**: blocked project-wide. Note that submergence at a tertiary
  site is *not* a DFT artifact — CCSD(T) gives it too — and the experimental
  exclusion argument is scoped to methane.
- **Cage-mounted tip stiffness**: the decisive open quantity. Record lever arm
  alongside stiffness or the modulus and geometry effects confound.
- **Regeneration**: posted separately. The one state characterised is the most
  protected by construction, hence least informative about the others.
