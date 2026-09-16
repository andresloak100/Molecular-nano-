# Established quantitative findings, with provenance

One citable place for the numbers this lane has established, so other sessions
can reference them instead of reconstructing them from `CLAUDE.md`, which has
grown large and is mostly coordination traffic.

Everything here was executed, not estimated. Each entry names the artifact that
produced it. Anything not listed here is not established, however plausible it
may look in a coordination note.

Scope: source forensics and reference accuracy only. Nothing here validates the
tool, the candidate, or the method for mechanosynthesis.

---

## 1. The published energy inconsistency is a single erroneous methane entry

`SOURCE_NOTES.md` recorded that the supporting-information absolute energies
imply a −14.76 kcal/mol barrier against the paper's positive value, with no
cause assigned. That is now traced.

All-electron CCSD(T)/cc-pVDZ on the supplied geometries, lowest SCF solution:

| Species | Published (Ha) | Recomputed (Ha) | Difference (kcal/mol) |
|---|---|---|---|
| CH3-H-CCH (TS) | −116.790593477 | −116.790593495 | 0.00 |
| HCCH | −77.115274197 | −77.115274205 | 0.00 |
| CH3 | −39.718644174 | −39.718644167 | 0.00 |
| CCH | −76.404097045 | −76.404097052 | 0.00 |
| **CH4** | **−40.362981447** | **−40.390318875** | **−17.15** |

Four of five reproduce to under 1×10⁻⁷ Ha, which fixes method, basis, core
treatment and geometry. Methane is the sole outlier and the entire source of
the −14.76 anomaly. Substituting only the recomputed methane value into the
otherwise unchanged published set gives **+2.3986 kcal/mol** against the
paper's own Table 4 value of 2.4, so the entry behaves like a transcription
error rather than a different calculation. No author correction has been
sought and none is assumed.

Artifact: `si-energy-reproduction.json`. Reproduce: `python reproduce.py`.
Independently matched by Codex's paired run at +2.3986156 by a separate route.

## 2. In-house high-level barrier: +2.3986 kcal/mol

Computed end to end here from all five recomputed energies, all-electron
CCSD(T)/cc-pVDZ, lowest SCF solution throughout. It agrees with the published
barrier, so **the project no longer depends on trusting the published table**
for this reaction.

| Barrier source | kcal/mol |
|---|---|
| Published absolute energies as printed | −14.76 |
| **Our own recomputation, all five species** | **+2.3986** |
| Published set with only methane replaced | +2.3986 |
| Paper Table 4 | +2.4 |

## 3. The SCF multiple-solution trap

**This invalidates open-shell Hartree-Fock results if ignored, and it gives no
warning.** The ethynyl radical has at least two converged, *stable* UHF
solutions at the supplied geometry:

| Initial guess | UHF energy (Ha) | S² | Stability analysis |
|---|---|---|---|
| `minao` (PySCF default), `huckel` | −76.143320 | 0.7975 | reports **stable** |
| `atom`, `1e` | **−76.157101** | 1.2212 | reports **stable** |

8.7 kcal/mol apart, both genuinely stable. The default lands on the higher one
and CCSD(T) built on it is 14.3 kcal/mol too high. That single artifact made a
correct published entry look erroneous and produced a retracted conclusion of
mine.

**Correction, and it matters: the transition structure is NOT in the same
situation.** An earlier version of this section quoted a 22.9 kcal/mol spread
for the TS alongside the ethynyl warning, which implied the two cases were
alike. They are not. Per-guess detail from the same artifact:

| Guess | Ethynyl, above lowest | TS, above lowest |
|---|---|---|
| `minao` (default) | **+8.76** | 0.00 |
| `huckel` | **+8.76** | 0.00 |
| `atom` | 0.00 | 0.00 |
| `1e` | 0.00 | **+22.92** |

For the ethynyl radical the default guess is **wrong**, and two of four guesses
land on the higher solution. For the transition structure the default guess is
**right**: `minao`, `atom` and `huckel` all converge to the identical lowest
solution, and the 22.9 figure is driven entirely by `1e`, a crude bare-nucleus
guess nobody uses in practice.

So a raw spread is the wrong summary statistic. The quantity that matters is
how far the *default* sits above the lowest, which is 8.76 kcal/mol for ethynyl
and 0.00 for the transition structure. **No coupled-cluster number built on the
TS is compromised by guess choice**; the problems at that geometry are
geometric, not electronic-state. Caught by support session 76190bf3 while
binding these numbers into a claims ledger.

Scan initial guesses for any open-shell HF work here, and report
default-above-lowest rather than a spread.

### CORRECTION: my "DFT is guess-independent" claim was over-scoped. Do not rely on it.

I previously wrote that DFT was checked and is guess-independent, full stop.
What I actually established is narrower: **the five methane-reaction species, at
the published geometries, at PBE0/def2-SVP**, showed zero spread across all four
guesses. I stated a property of those five calculations as a property of DFT.

A2 has since found a counterexample in its own lane. Propynyl at
PBE0-D3(BJ)/def2-SVP: `minao` converges to a solution **11.2 kcal/mol above**
the one reached from `atom`, `huckel` and `1e`. So DFT is *not* generally
guess-independent here, and the same trap that cost me the ethynyl result
applies to density-functional calculations on open-shell radicals.

**Consequence for anyone relaxing new radicals:** a guess scan is required for
DFT too, at every new geometry, not just for Hartree-Fock. If you took my
earlier blanket statement as permission to skip it, that was my error.

### The "cleanest S²" heuristic is backwards. Three for three.

The natural shortcut — when two solutions disagree, trust the less
spin-contaminated one — has been **wrong in every case examined so far**:

| Species | Method | Cleaner S² | Its energy |
|---|---|---|---|
| Ethynyl radical | UHF/cc-pVDZ | `minao`, S² = 0.7591 | **+8.76 kcal/mol, wrong** |
| Methane TS | UHF/cc-pVDZ | `1e`, S² = 0.7595 | **+22.92 kcal/mol, wrong** |
| Propynyl (A2) | PBE0-D3/def2-SVP | `minao`, S² = 0.7521 | **+11.2 kcal/mol, wrong** |

In all three the lower, correct solution is the *more* contaminated one. A
reader economising on compute will reach for exactly this heuristic; it would
have picked the wrong answer every time.

Three points and variational lowness do not make a selection rule, and none of
this says the low solution is the physically right state — only that it is the
variationally lower one and that S² does not identify it. Scan the guesses.

## 4. Retracted numbers

Listed explicitly because both circulated and one is still committed elsewhere.

- **−11.92 kcal/mol**, methane barrier at CCSD(T). Withdrawn. Artifact of the
  wrong ethynyl solution. Correct value is +2.3986.
- **−39.11 kcal/mol**, methane reaction energy. Withdrawn. Same cause. Correct
  value is **−24.788**. This figure and a derived `dft_minus_ccsd_t` of +12.31
  remain in `data/validation/paired-ccpvdz/method_comparison.json`, the minao
  arm of a deliberate paired comparison. The `-atom` arm carries the correct
  −24.788 and −2.01. Raised with that lane's owner; annotation rather than
  deletion is the agreed disposition, since deleting the artifact arm would
  destroy the evidence the comparison exists to produce.

## 5. Reaction energies, and the tertiary anchor

All-electron CCSD(T)/cc-pVDZ, published geometries, lowest SCF solutions.

| Reaction | Site | Reaction energy (kcal/mol) |
|---|---|---|
| C2H + CH4 → C2H2 + CH3 | primary C–H | −24.788 |
| C2H + iso-C4H10 → C2H2 + t-C4H9 | **tertiary C–H** | **−32.232** |
| difference | tertiary − primary | **−7.444** |

Against a literature C–H bond-strength difference of roughly 8–9, so the chain
is behaving. The isobutane reaction matters more here than methane because its
site is tertiary, the closest small-molecule proxy for the adamantane
bridgehead the candidate targets.

**Caution before transferring this to adamantane.** The −7.444 is tertiary
versus *primary* in an acyclic molecule, where the tertiary radical relaxes
toward planarity. The 1-adamantyl radical is held pyramidal by the cage and
forfeits that stabilization, so the bridgehead-versus-methylene difference
should be substantially compressed and may be near zero. Stated from general
radical chemistry, not from a calculation run here; A1's computed number
supersedes it. Recorded before A1's result landed so a near-zero value reads as
predicted rather than as a suspected bug.

## 6. Isobutane reproduction: four of five species

Ethynyl, acetylene, tert-butyl and isobutane all reproduce the published
absolute energies to under 1×10⁻⁷ Ha, reinforcing methane as the single
corrupted entry. The 139-basis-function transition structure is still running;
the published absolutes imply −0.627 kcal/mol for that barrier.

Artifact: `si-energy-reproduction-isobutane.json`, written per species so a
timeout leaves usable evidence. Reproduce: `python reproduce_isobutane.py`.

## 9. Why the tool is unselective: a mechanism, not a coincidence

Synthesis proposed by A2 from exothermicity; **tested here against the transition
structure geometries, which are independent evidence and could have disagreed.**

Three findings had been treated as separate puzzles: the tool shows no steric
discrimination at any cage site; the thermodynamic site preference is near zero;
and real adamantane selectivity is small, roughly 1.0–1.3 kcal/mol. They are one
finding. **The tool discriminates poorly because it is thermodynamically hungry.**
A reagent running 32 kcal/mol downhill reaches its transition state early, before
it has committed to a particular C–H, so the sites look alike to it.

The published geometries confirm it without using any energy:

| Reaction | ΔE (kcal/mol) | acceptor–H (Å) | donor–H (Å) |
|---|---|---|---|
| C2H + CH4, primary | −24.79 | 1.672 | 1.149 |
| C2H + iso-C4H10, **tertiary** | **−32.23** | **2.213** | **1.116** |

The more exothermic reaction has the acceptor 0.54 Å further away and the donor
C–H barely stretched, 1.116 against an equilibrium near 1.10. That is a
demonstrably earlier transition state for the more exothermic reaction — Hammond,
read straight off the coordinates. Bell–Evans–Polanyi across the two points gives
α = 0.41 using the SI-implied isobutane barrier, a normal early-transition-state
value.

**Two points do not establish a relationship.** The running isobutane transition
structure supplies a third and will test α rather than assume it.

### The design consequence, which is the part worth acting on

This converts the project's question from *"make this tool work"* to *"where on
the reactivity–selectivity curve should the tool sit?"* A less reactive tip would
discriminate better and might not react at all. The repository explores exactly
one point on that axis and has no way to say whether it is the right one.

That reframing also explains why A1's positional result matters more than it
first appeared. If chemical discrimination is intrinsically weak for a reagent
this hot, then positional control is not one selectivity mechanism among several
— it is close to the only one available, which is precisely what A1's steric
census independently found.

## 10. The tool may weld itself to the workpiece — and the geometry may prevent it

Failure mode raised by A1; protective hypothesis **tested here against the
candidate's own product coordinates**, which A1 had not used.

Follow the operation to its end. After transfer, the tool is a closed-shell
terminal alkyne, adamantyl–C≡C–H, and the workpiece is a 1-adamantyl radical,
sitting 3.6 Å apart. A carbon radical beside an alkyne is not a stable
arrangement: addition across the C≡C forms a C–C σ bond worth roughly 85
kcal/mol at the cost of demoting C≡C to C=C, roughly 54, so the addition is
around 30 kcal/mol exothermic on bond additivity. **The intended product is
metastable and the deep well is "tool covalently bonded to workpiece".**

This failure mode is worse than mis-targeting and has had none of the
attention. Mis-targeting places one atom badly; welding destroys the tool and
the workpiece together and produces no further products at all.

### Two protections, both free, both angular

Measured on the product geometry from `candidates.py`:

| Quantity | Value |
|---|---|
| angle: C≡C axis vs apex→radical | **0.0°** |
| radical to apex carbon | 3.600 Å |
| radical to transferred H | **2.540 Å** |

**First**, the radical sits exactly end-on along the C≡C axis. Radical addition
needs a perpendicular approach into a π lobe; 0° is the worst possible vector
for it. A1 predicted this and the coordinates confirm it exactly.

**Second, not previously noted:** the transferred hydrogen lands on the apex
carbon, directly between the radical and the alkyne, closer to the radical
(2.540 Å) than the apex carbon itself is (3.600 Å). The newly formed C–H
physically occupies the approach vector. This is structural rather than
incidental — the H necessarily lands on the atom the radical was pointing at,
so **the abstraction event installs a steric block against the addition that
would otherwise follow it.**

Both protections derive from the same collinearity and degrade under the same
angular wander, so one angular tolerance covers both. **No angular tolerance has
been computed**; every positional analysis so far has been a lateral distance.

Neither protection is established as sufficient. Thermodynamics says the well is
deep; geometry says the approach is bad. Which wins is a *barrier* question, and
barriers are blocked on the same missing machinery as everything else kinetic —
now the third independent line arriving at that gap.

## 7. What is NOT established

- **No DFT barrier exists for either reaction**, because no DFT saddle has been
  located. Every DFT figure in circulation is a single point at a
  coupled-cluster stationary point that is not stationary on any DFT surface;
  S1 measured a 2.38 eV/Å residual force there under PBE0, about eighty times
  this repository's 0.03 convergence threshold. S1 is closing this.
- **Nothing about the 53-atom candidate.** No path, no barrier, no selectivity.
- **No kinetic site preference**, for anything. A1 established that the
  abstractor terms cancel exactly between the thermodynamic cycles, so site
  preference is abstractor-independent by construction and cannot come from
  reaction energies. It must come from barriers, which do not yet exist.
- **No experimental validation of anything.**

## 8. Host contention — do not read today's timings as costs

Controlled measurement: the identical benchmark call took **6.1 s** on a quiet
host and **63.4 s** at load 201. A **10.4× slowdown**, with 12 peer sessions
plus roughly a dozen internal agents on 8 logical cores. Every wall-clock
figure produced today is inflated by about an order of magnitude and describes
the scheduler rather than the calculation. Ratios and energies are unaffected.
