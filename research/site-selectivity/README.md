# A1 — Does the tool hit the right hydrogen?

Lane owner: session `andresarriaga-f2`. Owned paths: this directory and
`coordination/status/andresarriaga-f2.md`. Core `nanodesign/` modules are
imported read-only and never edited.

Atomically precise assembly means abstracting one *specific* hydrogen. The
53-atom candidate targets a bridgehead hydrogen of adamantane — carbon 0,
hydrogen 10. Adamantane has 4 bridgehead tertiary C-H bonds and 12 methylene
secondary C-H bonds. This lane asks whether the other fifteen are competitive.

The question splits cleanly in two, and the answers point in opposite
directions:

| | Question | Status |
|---|---|---|
| Geometric | Are the wrong sites distinguishable by position? | **Answered.** Yes, with a 2.495 Å margin, and none is sterically blocked. |
| Chemical, thermodynamic | Does the chemistry prefer the right site? | **Running.** Pre-registered expectation: no — it favours the *wrong* sites by roughly 2 kcal/mol. |
| Chemical, kinetic | Do the barriers prefer the right site? | **Out of reach.** No DFT saddle exists for even the calibration reaction. |

## Stage 0 — geometric site census. Complete.

`site_census.py` → `evidence/stage0-site-census-r5/census.json`. No electronic
structure; this is measurement of the actual candidate coordinates. Earlier runs
r1 through r4 are retained unmodified; **r5 is the live artifact**. The margin
value has never changed across runs; later runs add the margin metric, unambiguous
key names, and the competitor shell grouping.

| Quantity | Value |
|---|---|
| Site classification, from coordinates alone | 4 bridgehead tertiary, 12 methylene secondary |
| Symmetry classes among the 15 alternatives | 6 / 3 / 3 methylene, 3 bridgehead |
| Intended apex-to-transferred-H | 2.510 Å |
| Intended donor-to-apex | 3.600 Å, collinear by construction |
| **nearest_hydrogen_margin** | **2.495 Å**, limiting competitor H18 |
| site_relocation_distances, nearest alternative | 5.200 Å, H14 |
| Sites sterically blocked by the tool | **0 of 15** |
| Van der Waals clearance, every site | +0.200 Å, always the intended apex-to-donor-carbon contact |
| Apex to nearest own-mount hydrogen | 3.551 Å |

The 6/3/3/3 split is the pattern adamantane's symmetry requires, which is an
internal check that the classifier and the geometry agree.

### Two distances, and which one is the margin

These are different quantities that land close together in this geometry, and
conflating them caused a factor-of-two disagreement downstream.

- **`nearest_hydrogen_margin` = 2.495 Å.** The smallest apex displacement, in
  *any* direction, that makes some other cage hydrogen the nearest one to the
  apex. Computed as the distance from the nominal apex to the perpendicular
  bisector plane of the intended and competing hydrogen, so it is a true minimum
  over directions rather than one chosen path. Closed form, verified by hand:
  (4.327² − 2.510²) / (2 × 2.489) = 2.4956.
- **`site_relocation_distances` = 5.200 Å.** The apex translation needed to
  re-aim the tool into the *identical idealized collinear geometry* over another
  hydrogen. A deliberate re-aiming distance.
- **2.489 Å** is the plain H18-to-H10 internuclear distance. A structural fact
  about adamantane, not a margin, and not quoted here.

**This lane quotes 2.495 Å as the selectivity margin**, because it is the
smaller and therefore conservative one: it bounds how far the apex may wander
before it is addressing a different atom at all. Quoting 5.200 Å as "the margin"
overstates the tolerance by a factor of two. An earlier census put the 5.200 Å
block under a parent key named `selectivity_margin`, which is what misled
readers; r3 renames it to `site_relocation_distances`.

### The competitors form four degenerate shells

| Shell | Sites | Apex displacement | Type | Hydrogens |
|---|---|---|---|---|
| 1 | 6 | 2.4950 Å | methylene secondary | 14–19 |
| 2 | 3 | 4.0664 Å | methylene secondary | 21, 23, 25 |
| 3 | 3 | 4.1968 Å | bridgehead tertiary | 11, 12, 13 |
| 4 | 3 | 4.8821 Å | methylene secondary | 20, 22, 24 |

Members agree to machine precision. Three of the four shells are methylene, so
**site type alone does not identify a shell** — "the methylene competitor" is
ambiguous. Reporting a single limiting hydrogen index previously caused three
separate readers to mistake an arbitrary tie-break for a distinguished atom, and
each named a different three-element subset of the six-fold shell 1; the census
now reports the tied set and the shell grouping.

**These shells are not a symmetry reduction.** They are a property of this tool
placement on an unrelaxed candidate, and exact degeneracy is precisely what an
idealized cage must produce. Relaxation, a different mount orientation or any
azimuthal preference would lift the ties by an amount this file cannot bound.
That is a different thing from the species-level symmetry used in stage 1, which
*is* rigorous — see below.

One consequence worth stating: the three remaining bridgehead hydrogens sit in
shell 3, farther than nine methylene hydrogens. **Every nearest competitor is the
opposite site type**, which per the pre-registration below is also the
thermodynamically favoured type. Geometry and chemistry point the same unhelpful
way — the closest wrong targets are also the easiest ones to abstract.

### What stage 0 establishes and what it does not

It establishes that the fifteen wrong hydrogens are *positionally*
distinguishable, and — the sharper half — that this tool has **no steric site
discrimination whatsoever**. At all sixteen sites the limiting contact is the
intended apex-to-donor-carbon approach, never a clash with the rest of the cage.
So the margin is a positioning **requirement**, not a steric barrier that
prevents mis-abstraction.

It does not establish selectivity. Being the nearest hydrogen is not the
reaction criterion, both molecules are held rigid at unrelaxed geometries, no
approach path or barrier is computed, and a necessary geometric condition is not
a demonstration.

### Incidental finding, routed to C1

At the nominal 3.6 Å pose the tool sits at +0.200 Å of van der Waals clearance —
essentially at contact. Any campaign pose with separation below 3.6 Å is inside
van der Waals overlap before any chemistry is considered. Possibly deliberate
for an abstraction geometry, but it should be recorded rather than rediscovered
when the nine-pose results are interpreted.

## Stages 1 and 2 — site energetics. Running.

`site_energetics.py` → `evidence/stage12-pbe0-svp-df/`.
PBE0-D3(BJ)/def2-SVP, density fitting, six species each relaxed independently to
`fmax` 0.03 eV/Å, one process, one thread. `species.py` builds the species; the
adamantane cage is sliced out of the 53-atom candidate itself so the isolated
molecule is constructed identically to the target cage.

### Exactly two radicals, for rigorous reasons

Isolated adamantane has T_d symmetry: its 4 bridgehead hydrogens form a single
symmetry orbit and its 12 methylene hydrogens form another. So there are exactly
**two** distinct adamantyl radicals. Deleting H14 rather than H16 does not give
two similar molecules to be sampled; it gives the *same* molecule in two
orientations, with exactly equal energy under any method.

`tests/test_species_symmetry.py` verifies this from the actual coordinates
instead of asserting the point group, by comparing sorted interatomic distance
spectra — invariant to translation, rotation, reflection and atom ordering. All
12 methylene radicals are congruent to within 1e-10 Å, all 4 bridgehead radicals
likewise, and a negative control confirms the two classes are *not* congruent, so
the reduction cannot have silently merged the two sites.

This is why stage 1 computes one representative per environment and is complete
rather than a screening shortcut, and it is why a proposal to cut stage 1 from
six relaxations to one had no saving available: its six relaxations are six
distinct chemical species, not six site samples.

Stage 1 is the adiabatic C-H dissociation energy at each site,
D = E(radical) + E(H) − E(adamantane). Stage 2 is the ethynyl abstraction
reaction energy, ΔE = E(C2H2) + E(radical) − E(C2H) − E(adamantane).

### The two cycles give the same site difference, identically

Both share the two adamantyl radical energies, so

    D(methylene) − D(bridgehead)
      = ΔE(methylene) − ΔE(bridgehead)
      = E(adamantyl_methylene) − E(adamantyl_bridgehead)

The hydrogen-atom, adamantane, ethynyl and acetylene terms cancel exactly. Two
consequences, one useful and one structural:

1. It is a free arithmetic self-check, and the code asserts it numerically.
2. **Site thermodynamics is abstractor-independent by construction.** No choice
   of abstracting species can change the thermodynamic site preference, so
   stage 2 contributes absolute abstraction energies and nothing new about
   preference. The original A1 brief expected otherwise; the intake owner has
   independently verified the cancellation and accepts the correction.

Therefore any abstractor-specific selectivity must be **kinetic** — and that is
currently unavailable in this project, not merely unmeasured by this lane. No
DFT saddle has been located for even the methane calibration reaction, S1
measures a 2.38 eV/Å residual force at the published transition geometry under
PBE0 (about eighty times the 0.03 convergence threshold), and every DFT barrier
so far is submerged against an in-house CCSD(T) reference of +2.40 kcal/mol.

### Pre-registered expectation, recorded before the number exists

Do not assume adamantane inherits the acyclic tertiary-over-secondary
preference. The 1-adamantyl radical is held pyramidal by the cage and cannot
relax toward planarity, so it forfeits the hyperconjugative stabilization that
normally makes a tertiary site preferred; EPR puts the bridgehead angle at about
113.6° against about 118° for the near-planar tert-butyl radical.

Gas-phase experiment plus G3 (Fattahi & Kass 2012, reported in a review and
independently verified by a peer session against Crossref and Europe PMC) gives
the adamantane **bridgehead** C-H BDE as **102.4 ± 1.9 kcal/mol**, "notably
higher" than the secondary C-H bonds at about **100.3**. An earlier 96.3 value
from the same review is superseded within it and is not used.

So the expectation is not that the gap collapses to zero but that it **inverts**:

> **`dissociation_difference` = D(methylene) − D(bridgehead) ≈ −2.1 kcal/mol,
> i.e. the methylene sites are the thermodynamically easier abstraction.**

Qualifications kept with the prediction: the experimental ±1.9 is comparable to
the 2.1 gap; BDEs are thermodynamic and not kinetic; and the transition-state
term that would favour a tertiary site is polar matching, which is conditioned
on an *electrophilic* HAT reagent. Ethynyl is a hot, essentially non-polar σ
radical and collects little of it — a mechanistic inference, not a sourced claim
about ethynyl. Strip that term and both remaining terms, thermodynamics and
sterics, favour the methylene sites.

If that holds, the honest statement is stronger than "no intrinsic
discrimination": **the intrinsic chemistry favours the wrong sites, and
positional control must overcome an adverse preference rather than a neutral
one.** The measured number governs over all of this.

For comparison, not as a target: all-electron CCSD(T)/cc-pVDZ gives tertiary
abstraction 7.44 kcal/mol more exothermic than primary for C2H + isobutane
versus C2H + methane. That is acyclic, and the cage is exactly what breaks the
analogy.

### Open-shell hazard, load-bearing rather than procedural

The site difference *is* the two radical energies, so a wrong SCF solution in
either one replaces the answer rather than degrading it — and it would look
entirely plausible.

This project's most expensive error came from exactly that: the ethynyl radical
has two converged, both-stable UHF solutions 8.7 kcal/mol apart, with PySCF's
default `minao` on the higher one. DFT was initially reported guess-independent;
that claim has since been **retracted as over-scoped** — it held for five
species at fixed geometries. A peer lane then found propynyl at
PBE0-D3(BJ)/def2-SVP with `minao` **11.2 kcal/mol above** the solution `atom`,
`huckel` and `1e` agree on. Open-shell DFT is therefore not safe by default.

`scf_guess_scan.py` scans `minao`, `atom`, `huckel` and `1e` at each relaxed
open-shell geometry and reports the spread. It runs once the radicals relax.

**Spin contamination must not be used to choose between solutions.** In all
three cases examined here the *cleaner* S² was the wrong, higher solution:

| Species | Method | Clean-looking guess | S² | Error |
|---|---|---|---|---|
| ethynyl | UHF | minao | 0.7591 | +8.76 kcal/mol |
| methane TS | UHF | 1e | 0.7595 | +22.92 kcal/mol |
| propynyl | PBE0-D3 | minao | 0.7521 | +11.20 kcal/mol |

A known limitation of the current setup: `PySCFCalculator` exposes no
initial-guess setting, so the relaxations run on PySCF's default. If the scan
finds a spread, the geometry was optimized on the wrong surface too and a
post-hoc scan cannot repair it; that would require an `init_guess` option in the
core calculator, which is Codex's to add, or a lane-local relaxation driver. The
scan is a gate on the result, not a formality.

### Density-fitting error, measured rather than assumed

An earlier direct run without density fitting was stopped when host contention
made it project to about eight hours. Its completed records are preserved at
`evidence/stage12-pbe0-svp/` and now bound the DF approximation at identical
settings: ethynyl 0.004, acetylene 0.006, hydrogen atom 0.033 kcal/mol.
Negligible against the 0.6 kcal/mol k_B T scale, and these are absolute
energies — the error cancels further between two C10H15 isomers.

`fmax` 0.03 eV/Å is ample: near a minimum the energy is quadratic, so a
0.03 eV/Å residual contributes well under 0.05 kcal/mol summed over all modes.

### Timing caveat

Every wall-clock number from this period is contaminated. The host has been at
load 187–274 on eight logical CPUs with roughly twenty sessions active, and it
is also paging hard after 39 days of uptime, so jobs stall on memory as well as
queue for cores. A controlled measurement elsewhere in the project found the
same benchmark call taking 6.1 s quiet and 63.4 s loaded, a 10.4× slowdown. No
correctness output is affected — energies, geometries, convergence and spin
diagnostics are all load-independent.

## Stage 3 — transition structures

Not attempted. Given the state of the barrier machinery above, the honest
deliverable is to say so with the reason rather than half-finish it.

## The other half of the selectivity argument

`positional_uncertainty.py`, with `tests/test_positional_uncertainty.py`
(18 passing).

A margin is only half an argument. The other half is the positional
*uncertainty* of a real mounted tip — thermal and zero-point displacement
against the handle's finite stiffness — which is the quantity Drexler's
feasibility case actually turns on. `docs/MODEL.md` states that finite anchor
stiffness and thermal motion are not represented, and nothing in `nanodesign/`
computes them: there is no stiffness, compliance or positional-spread code in
the package. Until something computes it, this project can state a 2.495 Å
margin but cannot say whether any mechanism meets it.

This module supplies it as analysis over a Hessian that has already been
computed, so it adds no electronic-structure cost. It consumes
`nanodesign.stationary.characterize_stationary_point` read-only.

What it computes:

- **Tip stiffness in N/m** from the *compliance*, the inverse Hessian restricted
  to one atom, which lets the rest of the structure relax as a real handle does.
  The clamped diagonal Hessian block is reported alongside, because it is always
  stiffer and the gap between them is the error that clamping hides.
- **Thermal-plus-zero-point positional spread.** Classically Cov = k_B T H⁻¹;
  quantum mechanically the mode sum with ⟨Q²⟩ = (ħ/2ω)coth(ħω/2k_BT). The
  quantum form matters: for a stiff diamondoid tip at 300 K, ħω ≫ k_B T, so
  zero-point motion dominates and equipartition *underestimates* the spread.
- **Pair-distance fluctuation**, which is translation- and rotation-free and so
  is defined without anchors. It is also the natural measure for a transfer
  coordinate such as donor-to-hydrogen.
- **A length ratio** against a measured margin, and deliberately nothing else.

Validation, all against closed-form results rather than recorded outputs:

- Drexler's canonical case: 10 N/m at 300 K gives σ = 0.2035 Å.
- 1 eV/Å² = 16.02176634 N/m, derived from SI constants.
- The classical mode sum equals k_B T H⁻¹ to within 1e-12 — two independent routes.
- Quantum reduces to classical in the high-temperature limit, and at 1 K the
  ratio is the analytic √(ħω/2k_BT) ≈ 40.
- Compliance-based stiffness matches the closed form t(t+2c)/(c+t) for a coupled
  pair, with the c→∞ and c→0 limits checked, and is always softer than the
  clamped block.
- A diatomic's bond-length fluctuation matches √((ħ/2μω)coth(ħω/2k_BT)) exactly.
- An end-to-end run through the *real* core Hessian contract, using a plain ASE
  Lennard-Jones calculator so it costs microseconds and does not compete for the
  host. This one caught a genuine error in an earlier version of the test, where
  a 298 K amplitude was compared against the T=0 zero-point formula and was
  wrong by a factor of 4 — the soft dimer mode is classical, not zero-point
  dominated.

Fail-closed behaviour, because these are the ways the calculation becomes
meaningless: a saddle or any nonpositive mode is refused, a soft near-zero mode
is refused with an explanation that difference coordinates are the alternative,
a non-stationary geometry is refused unless explicitly overridden, and a frozen
atom has no distribution. Soft modes excluded from a pair-distance sum have the
variance they would have contributed reported, so the exclusion is auditable
rather than invisible — a near-zero mode carries a near-divergent amplitude, and
relying on its overlap vanishing by symmetry is luck rather than method.

What it does **not** do. It returns no error rate, success probability or
reliability figure, and it never will from these inputs. Equilibrium harmonic
statistics do not describe a driven assembly step: the tool is pushed along a
reaction coordinate, the model contains no barrier, anharmonicity grows exactly
where displacements are large, and the anchors are rigid — so every stiffness
here is an upper bound and every displacement a lower bound. A real mount is
softer and the real spread is larger.

Next step for it, not yet run: a tip Hessian for the actual candidate. A full
Hessian on the 53-atom system is out of reach on this hardware, so the route is
a partial Hessian over the apex and its neighbours with the cage frozen, which
is a bounded upper bound on stiffness and is stated as such.

## Reproducing

```
python research/site-selectivity/site_census.py
python research/site-selectivity/site_energetics.py --output research/site-selectivity/evidence/<new-dir> --density-fit --fmax 0.03
python research/site-selectivity/scf_guess_scan.py --geometries <stage12 dir> --output <new-dir>
python -m pytest research/site-selectivity/tests/ -q
```

Every output directory must be new; the drivers refuse a directory named `runs`
and `site_energetics.py` only reuses an existing one under `--resume`, and then
only for species already recorded complete.
