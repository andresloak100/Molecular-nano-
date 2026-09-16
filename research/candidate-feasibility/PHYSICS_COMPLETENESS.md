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
| 10 | Positional / thermal error rate | **done** | A2 |
| 11 | Mount stiffness and mechanical model | **not done** | unowned |
| 12 | Tool regeneration: a closed cycle | **not even posed** | unowned |

Below the ladder sit two cross-cutting requirements — multireference validity
of the reference method, and prospective experimental validation — neither of
which is satisfied.

## The four gaps that would change conclusions, in order

### 1. Classical transition-state theory is outside its validity range here

The transferred particle is a hydrogen atom, and hydrogen tunnels. How much is
set by the barrier's curvature, not its height, through the crossover
temperature `T_c = ħ|ω*| / (2π k_B)`. Measured across the range typical of
hydrogen-transfer saddles:

| ω* (cm⁻¹) | T_c (K) | Is 298 K below it? |
|---|---|---|
| 1000 | 229 | no |
| 1200 | 275 | no |
| 1500 | **343** | **yes** |
| 2000 | 458 | **yes** |

Above `T_c`, tunneling is a correction. Below it, tunneling is the mechanism
and a classical rate is qualitatively wrong rather than merely inaccurate.
Room temperature sits inside that range, and every cryogenic operating
proposal is far below it — at 77 K even a soft 500i cm⁻¹ barrier is past
crossover.

Two consequences beyond kinetics. Tunneling reweights competing channels by
barrier *width*, so a selectivity argument that compares only barrier heights
is incomplete. And it predicts a large H/D kinetic isotope effect, which is a
falsifiable experimental signature — something this project is short of.

What is needed: `ω*` from a verified saddle (S1's lane), then Eckart or
instanton rate theory. `research/candidate-feasibility/tunneling.py` has the
crossover and the parabolic-barrier corrections and deliberately stops short of
Eckart rather than reproduce a long formula unverified.

### 2. Nobody has asked what positional precision the design needs — now answered

A peer lane established that this tool has no steric discrimination across
adamantane's sixteen C–H sites, and the external literature predicts a near-zero
thermodynamic preference between site types. If neither sterics nor
thermodynamics chooses the site, position does, and position at finite
temperature is a distribution.

Computed (`positional_control.py`, geometry plus equipartition, no quantum
time):

- The tip may wander **2.878 Å** laterally before a methylene hydrogen is
  nearer than the target bridgehead hydrogen.
- A 10⁻¹⁵ error rate therefore needs lateral σ < 0.346 Å, i.e. a lateral
  stiffness above **0.214 eV/Å² (3.4 N/m)** at 298 K, or 0.9 N/m at 77 K.
- 3.4 N/m is very soft. At 0.5 eV/Å² the mis-targeting rate is 10⁻³⁵.

So thermal mis-targeting between these sites is **not** a credible failure
mode, by roughly twenty orders of magnitude. This is the project's first
positive feasibility result, and its scope is narrow: it says one specific
failure mode is not what will kill the design.

A methodological warning attached to it. Ranking rivals by lateral distance
from the tool axis — the obvious approach — picks three hydrogens 1.452 Å
off-axis and gives a frightening 0.73 Å margin. Those three sit at z = −3.66,
on the far side of the cage, 7.4 Å from the apex, unreachable at any offset.
The error is fourfold and in the alarming direction.

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

**Zero-point energy is missing from every barrier.** For a C–H bond it lowers a
dissociation energy by roughly 4–5 kcal/mol, against a benchmark barrier of
+2.4. The correction is larger than the quantity.

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
