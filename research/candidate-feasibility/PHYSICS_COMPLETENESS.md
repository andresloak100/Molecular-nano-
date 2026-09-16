# What stands between this repository and a scientifically complete design model

Written from the A2 lane, 2026-09-16, drawing on results from every active
lane. It is a map of what is missing, not a claim that anything here is done.

The repository can currently compute an electronic energy and its gradient for
a 53-atom candidate, reproduce a published high-level barrier, and refuse to
overwrite its evidence. That is a real foundation and it is roughly level 2 of
the twelve below. The gap between "runs quantum chemistry" and "predicts
whether a molecular machine works" is most of the ladder.

## The ladder

| # | Capability | State | Owner |
|---|---|---|---|
| 1 | Explicit structures, atom conservation, anchors | **done** | C1 |
| 2 | Electronic energy and gradient, provenance recorded | **done** | C1 |
| 3 | Correct electronic state (SCF solution scans) | **partial** | all lanes, no API primitive |
| 4 | Relaxed endpoint geometries for the real candidate | **not done** | blocked on cost |
| 5 | A located saddle with exactly one imaginary mode | **in progress** | S1 |
| 6 | Zero-point energy on the barrier | **not done** | unowned |
| 7 | Free energy at the operating temperature | **not done** | unowned |
| 8 | Quantum rate theory (tunneling) | **scaffolded** | A2 |
| 9 | Enumerated competing channels, not just site choice | **not done** | unowned |
| 10 | Positional / thermal error rate | **done**, single-mode | A2 |
| 11 | Mount stiffness and mechanical force coupling | **partial**, k measured, k_force open | rung-5, A1 |
| 12 | Tool regeneration: a closed cycle | **not even posed** | unowned |

Rung 11 moved today: A1 measured a mounted tip's lateral stiffness at 7.27 N/m
instead of quoting the literature's 10–100, and a mechanical-coupling lane
opened to compute the force needed to drive the abstraction. Feasibility there
reduces to one comparison, `k_ceiling >= max(k_positioning(T), k_force)`, of
which two of the three terms are now numbers.

Rung 10's remaining weakness is worth naming because every other correction
found today ran the safe way and this one does not: it treats the tip as a
single harmonic coordinate, while the true positional variance sums over all
modes as `Σ kT/λ_k`. That understates σ by an unquantified amount — small
against the present margins, not zero, and closable only with the Hessian that
rung 11 also needs.

Below the ladder sit two cross-cutting requirements — multireference validity
of the reference method, and prospective experimental validation — neither of
which is satisfied.

## The four gaps that would change conclusions, in order

### 1. Whether classical rate theory applies here is unresolved, and one number decides it

The transferred particle is a hydrogen atom, and hydrogen tunnels. How much is
set by the barrier's curvature, not its height, through the crossover
temperature `T_c = ħ|ω*| / (2π k_B) = 0.2290 × ω*[cm⁻¹]` kelvin. Above `T_c`
tunneling is a correction; below it, tunneling is the mechanism and a classical
rate is qualitatively wrong rather than merely inaccurate.

| ω* (cm⁻¹) | T_c (K) | Is 298 K below it? | κ at 298 K |
|---|---|---|---|
| **259 (measured)** | **59** | no, 5× above | **1.07** |
| 1000 | 229 | no | 3.63 |
| 1500 | 343 | yes | diverges |
| 1648 (experiment-implied) | 377 | yes | diverges |

**A correction to an earlier version of this document.** It asserted that
hydrogen-transfer saddles "typically" run 1000–2000i and concluded room
temperature sits below crossover. That was a generic expectation for the
reaction class stated as a measurement of this reaction. The measurement was in
`SOURCE_NOTES.md` line 33 the whole time: 259i, 50i, 50i cm⁻¹ from the source
paper's Table 1. On that datum tunneling is a 7% correction, not the mechanism.

**But the question is genuinely open in both directions.** The 259i comes from
a structure carrying three imaginary modes, which is not a verified saddle, at
a geometry no functional here owns. Against it, the forensics lane finds that
reproducing the experimental apparent activation energy from the zero-point
corrected barrier requires roughly 1648i. Neither anchor is trustworthy. S1's
refined saddle will produce the first ω* on a functional's own surface.

Two conditional consequences, which fire only if that saddle comes back stiff.
Tunneling would reweight competing channels by barrier *width*, making any
selectivity argument that compares only heights incomplete. And it would
predict a large H/D kinetic isotope effect — a falsifiable experimental
signature this project is short of.

What is needed: ω* from a verified saddle, reported as a headline beside the
barrier rather than buried in a mode table. `tunneling.py` has the crossover
map and the parabolic corrections and deliberately stops short of Eckart rather
than reproduce a long formula unverified.

### 2. Nobody has asked what positional precision the design needs — now answered

A peer lane established that this tool has no steric discrimination across
adamantane's sixteen C–H sites, and the external literature predicts a near-zero
thermodynamic preference between site types. If neither sterics nor
thermodynamics chooses the site, position does, and position at finite
temperature is a distribution.

Computed (`positional_control.py`, geometry plus equipartition, no quantum
time):

- The apex may move **2.495 Å in any direction** before a different hydrogen is
  nearer than the target. The easiest escape is not lateral but tilted about
  120° from the tool axis, down toward the equatorial methylenes.
- A 10⁻¹⁵ error rate therefore needs per-axis σ < 0.292 Å, i.e. a stiffness
  above **0.30 eV/Å² (4.8 N/m)** at 298 K, or 1.2 N/m at 77 K.
- 4.8 N/m is very soft — a couple of orders of magnitude below ordinary
  covalent stiffness. At 0.5 eV/Å² the mis-targeting probability is 10⁻²⁶.

**And then a correction, because that requirement was compared against nothing.**
Saying 4.8 N/m "is very soft" measured it against an imagined mount. The rung-5
lane has since bracketed what is actually buildable, and load geometry — not
material — spans two orders of magnitude, because a cantilever softens as the
cube of its length:

| Mount | Stiffness | P(wrong site), 298 K | P, 77 K |
|---|---|---|---|
| Bending cantilever, soft end | 2 N/m | **1.2 × 10⁻⁶** | 3 × 10⁻²⁵ |
| Bending cantilever, stiff end | 20 N/m | 3 × 10⁻⁶⁵ | ~0 |
| Axial strut | 130–400 N/m | ~0 | ~0 |
| Single C–C bond, axial (cap) | 450 N/m | ~0 | ~0 |

The buildable range *starts below the requirement*. A long handle worked in
bending misses the target by nine orders of magnitude at room temperature and
clears it comfortably at 77 K. So the constraint is real and has a shape:

> **Mount stiffly and short, or operate cold.** Either suffices; neither is
> optional if you lack the other.

Neither lane could have found this alone — one had the requirement with no
buildable range, the other the range with no requirement.

**The condition that carries this result.** It is a criterion about which
hydrogen is **nearest**, not which one **reacts**. Those coincide only if the
competing barriers are comparable, which is A1's open question — and below the
tunneling crossover they would decouple further, since a more distant site with
a narrower barrier could win on width. So the finding is "thermal wander is not
the binding risk, *given* that nearest implies reacting", never "positional
control is solved."

Two methodological warnings attached, both from errors made here and caught.
Ranking rivals by lateral distance from the tool axis — the obvious approach —
picks three hydrogens 1.452 Å off-axis and gives a frightening 0.73 Å margin;
those three sit at z = −3.66, on the far side of the cage, unreachable at any
offset. And using the purely lateral crossover, 2.878 Å, overstates the
allowance by 15%, because thermal displacement is three-dimensional and finds
the easiest direction rather than the sideways one.

### 3. The tool is thermodynamically hungry, which explains the selectivity problem

The forensics lane measured the abstraction this design depends on, for a
tertiary C–H with an ethynyl abstractor, at CCSD(T)/cc-pVDZ: **−32.23 kcal/mol**.
Independently, this lane's acetylene C–H bond strength of 136.7 kcal/mol
(electronic, no zero-point) against a literature adamantane bridgehead C–H
bond of roughly 96–102 kcal/mol puts the same quantity in the same region.

A reaction that far downhill has an early transition state and a low barrier,
which is good for the intended operation. The same driving force is why the
tool discriminates poorly: a reagent with 32 kcal/mol of appetite is not fussy
about which C–H it takes. The reactivity-selectivity relationship is being
observed here, not violated, and it ties together three otherwise separate
findings — no steric discrimination, near-zero thermodynamic site preference,
and selectivity resting entirely on position.

The design implication is a genuine trade-off nobody has written down: a
*less* reactive tip would discriminate better but might not react at all. That
is a design axis, and the repository currently explores exactly one point on it.

### 4. The tool cycle is not closed, and has not been posed

Mechanosynthesis requires *reusable* tools. After abstracting a hydrogen, this
ethynyl tip is acetylene-terminated and cannot abstract again. Recharging it —
removing that hydrogen and restoring the radical — is a second reaction, with
its own barrier, its own selectivity problem, and its own error rate. It may
well be the harder half.

Nothing in the repository models it. Every result so far concerns a single
half-operation on a fresh tool. A design that performs one operation and then
stops is not a machine, and the cost of the cycle is not the cost of the
abstraction times two.

## Smaller gaps that still matter

**Zero-point energy is missing from every barrier, and the reason it matters is
not the obvious one.** An earlier version of this document said the 4–5
kcal/mol zero-point energy of a C–H stretch is larger than the +2.4 kcal/mol
benchmark barrier it would correct. That compares the wrong things: the barrier
correction is the *difference* in total zero-point energy between the
transition structure and the reactants, and the C–H stretch is weakened rather
than destroyed at the saddle, so most of it cancels. Temelso measures what
survives for this exact reaction — 2.2 kcal/mol electronic against 1.7 at 0 K,
so −0.5. The correction is a fifth of the barrier, not double it.

It still matters, for a sharper reason: −0.5 kcal/mol is the same order as the
1.0–1.3 kcal/mol kinetic site selectivity this design depends on, and it need
not cancel between two sites whose C–H frequencies differ. The threat is that
it is comparable to the signal and site-dependent, not that it is large.

**Free energy is not electronic energy.** Bringing a tool and a substrate
together costs translational and rotational entropy; every number here is a
bare electronic energy at 0 K.

**Competing channels are unenumerated.** Site selectivity is one failure mode.
Tip dimerization, radical recombination, hydrogen migration on the product
radical and abstraction from the tool's own handle are not modelled at all.

**Mount stiffness is unmeasured.** It adds to the positional stiffness above
(making that conclusion conservative), but it is also what a real mechanical
design would be engineered around. It needs a Hessian.

**Multireference validity is unchecked.** The S² of 1.21 at the reference
transition structure suggests genuine multireference character, which
single-reference CCSD(T) and every DFT functional describe poorly. A peer lane
has offered a T₁-diagnostic that reads amplitudes already computed and
discarded.

**No prospective experimental validation exists.** Every number is either
in-house or a reproduction of one published supporting-information table. The
first external anchors appeared today: the adamantane bond-strength and
radical-pyramidalization literature, and this lane's acetylene bond strength
landing where experiment says it should.

## The honest summary

The repository is a competent and unusually careful electronic-structure
harness that has not yet computed a rate, a barrier on its own surface, or a
machine cycle. Its discipline about provenance and refusing to overclaim is
genuinely better than typical, and it is the reason the gaps above are visible
rather than hidden.

The single most valuable next result is S1's verified saddle, because it
unblocks items 6, 8 and most of 9 at once. The single most under-examined
assumption is that a barrier height is a rate, which for this reaction at these
temperatures is not true.
