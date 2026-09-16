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
| **Adamantyl-ethynyl tool tip** (A2) | PBE0-D3/def2-SVP | `minao`, S² = 0.7521 | **+10.4 kcal/mol, wrong** |

In all three the lower, correct solution is the *more* contaminated one. A
reader economising on compute will reach for exactly this heuristic; it would
have picked the wrong answer every time.

Four points and variational lowness do not make a selection rule, and none of
this says the low solution is the physically right state — only that it is the
variationally lower one and that S² does not identify it. Scan the guesses.

### The failures are not random: every one is an alkynyl radical

Five cases is enough to ask what predicts them rather than just counting. Sorted
by how far the default `minao` guess sits above the lowest solution:

| Species | Method | minao above lowest | Class |
|---|---|---|---|
| ethynyl radical H–C≡C• | UHF/cc-pVDZ | **+8.76** | alkynyl radical |
| propynyl CH3–C≡C• (A2) | PBE0-D3/def2-SVP | **+11.20** | alkynyl radical |
| adamantyl-ethynyl tip (A2) | PBE0-D3/def2-SVP | **+10.42** | alkynyl radical |
| methane TS | UHF/cc-pVDZ | 0.00 | transition structure |
| isobutane TS | UHF/cc-pVDZ | 0.00 | transition structure |
| ethynyl radical | PBE0/def2-SVP | 0.00 | alkynyl radical |
| tert-butyl radical | PBE0/def2-SVP | 0.00 | saturated C radical |
| methyl radical | UHF/cc-pVDZ | 0.00 | saturated C radical |

**Every failure is an alkynyl radical, R–C≡C•.** Saturated carbon radicals are
clean in both cases checked, and in both transition structures `minao` found the
lowest solution (their spreads come from `1e`, which nobody uses).

Chemically this is the known σ/π near-degeneracy of C2H-type radicals: the
²Σ⁺ and ²Π states lie close together, so the SCF has two genuine low-lying
solutions to fall into. It is a property of that radical centre's electronic
structure, not of PySCF. Caveat on the apparent exception: bare ethynyl is clean
at DFT but fails at UHF, while substituted alkynyls fail at DFT too — different
bases are involved and I would not push that sub-pattern hard.

**Why this matters more than the base rate does. The tool tip is an alkynyl
radical.** It sits squarely in the affected class, and A2 measured it failing by
10.42 kcal/mol. Every calculation involving the *unreacted* tool is at risk and
must be scanned. Conversely, A1's adamantyl radicals are saturated carbon
radicals, the clean class, so that scan is more likely to return null — worth
knowing in advance so a null is not later read as evidence the scan was
unnecessary.

**The fourth case is not a proxy.** It is the adamantyl-supported ethynyl tool
tip itself, C12H15, truncated from the headline 53-atom candidate at its own
pose. So the trap is confirmed on the actual tool species, not only on
small-molecule analogues.

**Consequence for the headline candidate, stated as A2 stated it.** The
archived 53-atom runs record `scf_initial_guess` as null and
`initial_guess_scan_performed` as null, so they ran on PySCF's `minao` default
unscanned. Reassuring: both report S² = 0.78462, near the *correct* tool-tip
solution's 0.7863 rather than the trap's 0.7521, and the extra adamantane is
closed-shell so it should not shift S² much. Not proof: that compares S² across
a 53-atom and a 27-atom system, and the entire lesson of the four cases above is
that **S² is precisely the wrong quantity to select on**. Honest status:
*probably fine, never checked*. A four-guess SCF-only scan at the archived
geometry, about 32 minutes of CPU, would settle it, and every downstream number
in the repository rests on those two records.

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

## 6. Isobutane complete: all five reproduce, and the tertiary barrier is SUBMERGED

All five isobutane-reaction species reproduce the published absolute energies to
under 1×10⁻⁷ Ha, including the 139-basis-function transition structure (7126 s
under load). **Across both reactions, 9 of 10 published species now reproduce
and methane remains the sole corrupted entry.**

| Reaction | Site | CCSD(T)/cc-pVDZ barrier |
|---|---|---|
| C2H + CH4 | primary C–H | **+2.399** |
| C2H + iso-C4H10 | **tertiary C–H** | **−0.627 — submerged** |

Bell–Evans–Polanyi now has **both endpoints verified** rather than one inferred:
α = 0.406, a normal early-transition-state slope. The barrier falls 3.03
kcal/mol for 7.44 more exothermicity, which is the quantitative form of the
reactivity–selectivity mechanism in §9.

### This scopes a conclusion the fleet has adopted, including by me

The external-validation lane excluded every submerged DFT barrier using Opansky
& Leone: a rate that rises with temperature cannot come from a submerged
barrier. I endorsed it, and it is correct — **for methane.** That measurement is
C2H + CH4.

At the **tertiary** site, coupled cluster *itself* returns a submerged barrier.
So a submerged barrier is not intrinsically a DFT artifact; it is the right
answer for the more exothermic reaction. **And the tertiary site is the one the
tool actually targets.** Applying the methane-derived exclusion to the adamantane
bridgehead would be an error. The exclusion must be scoped to the reaction it
was measured on.

That also means the tertiary abstraction is effectively capture-controlled
rather than barrier-controlled, which removes barrier height as a selectivity
mechanism at that site — consistent with, and independent of, the
reactivity–selectivity argument in §9.

### Fifth SCF multiple-solution case

The isobutane transition structure shows a 10.64 kcal/mol spread across the four
guesses. `minao` happened to find the lowest here, so this number is safe, but
the spread is real and it is the fifth confirmed case.

Artifact: `si-energy-reproduction-isobutane.json`, written per species so a
timeout leaves usable evidence. Reproduce: `python reproduce_isobutane.py`.

## 12. "DFT gets the barrier sign wrong" is a methane statement

The isobutane DFT survey completes the comparison. Both methods at the same
published geometries, so this is a pure method comparison of the same quantity.

| Method | methane (primary) | isobutane (tertiary) |
|---|---|---|
| **CCSD(T)/cc-pVDZ** | **+2.40** | **−0.63** |
| PBE0-D3/def2-SVP | −3.62 | −3.11 |
| PBE0-D3/def2-TZVP | −2.89 | −2.80 |
| B3LYP-D3/def2-TZVP | −3.77 | −3.74 |

**At methane, coupled cluster is positive and every functional is negative — a
sign disagreement. At isobutane, coupled cluster is negative and every
functional is negative — they agree.** The magnitude error also halves, from
−5.29 to −2.17 kcal/mol for PBE0-D3/def2-TZVP.

So the project's most-repeated DFT result, that the functionals get the barrier
sign wrong, **is a statement about the methane reaction**. At the tertiary site
they are qualitatively right and quantitatively closer. This reaches the same
scoping conclusion as §6 by a completely different route — §6 argues from
coupled cluster returning a submerged barrier itself, this argues from DFT and
coupled cluster agreeing once the reaction is exothermic enough.

Note the error direction is constant even where the sign agreement changes:
every functional is too low at both geometries. What changes is whether "too
low" crosses zero, which is a property of where the reference sits, not of the
functional becoming better behaved.

**Caveat carried from §7:** none of these are DFT barriers. They are single
points at a coupled-cluster stationary point that is not stationary on any DFT
surface — S1 measured a 2.38 eV/Å residual there. The comparison is valid as a
method comparison at fixed geometry and says nothing about what a relaxed DFT
saddle would give.

### CORRECTION: M06-2X did complete, and it overturns the obvious reading

I recorded M06-2X as not finishing. **That was wrong** — I read the log before
the process exited. It completed at **−11.27 kcal/mol**, and it changes the
conclusion, so the error was not harmless: the support lane built an analysis on
the three remaining functionals and reached a conclusion this datum falsifies.

How much each method moves from the primary to the tertiary site:

| Method | methane | isobutane | shift | % of CCSD(T) shift |
|---|---|---|---|---|
| CCSD(T)/cc-pVDZ | +2.40 | −0.63 | **−3.03** | reference |
| PBE0-D3/def2-TZVP | −2.89 | −2.80 | +0.09 | **−3%** |
| B3LYP-D3/def2-TZVP | −3.77 | −3.74 | +0.03 | **−1%** |
| **M06-2X/def2-TZVP** | −8.23 | −11.27 | **−3.04** | **100%** |

"The functionals are blind to the substitution" is true of PBE0 and B3LYP — they
move 3% and 1% of the reference shift, in the wrong direction — and **decisively
false of M06-2X**, which reproduces the substitution response essentially
exactly. Notably it is the one parameterised on barrier heights.

**So M06-2X carries a large constant offset while tracking the site dependence
almost perfectly**, and its error is identical at both sites to 0.01 kcal/mol:

| Method | error, methane | error, isobutane | error drift |
|---|---|---|---|
| PBE0-D3/def2-TZVP | −5.29 | −2.17 | 3.12 |
| B3LYP-D3/def2-TZVP | −6.17 | −3.11 | 3.06 |
| **M06-2X/def2-TZVP** | −10.63 | −10.64 | **0.01** |

**Selectivity is a difference between sites, and a constant error cancels in a
difference while a drifting one does not.** That inverts the ranking: the worst
functional for absolute barriers is the best for the site difference, and
PBE0-D3 — which the selectivity lane is using — is among the worst for it.

**Two points cannot establish that an error is constant.** A 0.01 kcal/mol
agreement across two reactions is striking but could be coincidence, and a
systematic offset in a barrier-height-parameterised functional is plausible
rather than surprising. A third reaction would test it. This is a reason to run
the check, not a reason to switch functionals.

### Does it transfer to relaxed reaction energies? Yes, but weakly.

The above is all fixed-geometry barriers. The selectivity lane computes relaxed
**reaction energies**, a different quantity, so the functional-choice argument
had to be checked rather than assumed. Same four methods, both reactions:

| Method | shift, primary→tertiary | % of CCSD(T) | error drift |
|---|---|---|---|
| CCSD(T)/cc-pVDZ | −7.44 | reference | — |
| PBE0-D3/def2-SVP | −11.57 | 155% | 4.13 |
| PBE0-D3/def2-TZVP | −10.69 | 144% | 3.25 |
| B3LYP-D3/def2-TZVP | −10.69 | 144% | 3.25 |
| **M06-2X/def2-TZVP** | −9.26 | **124%** | **1.82** |

For reaction energies every functional **overshoots** the site effect (124–155%)
rather than missing it (−3%), so the blindness is specific to the barriers.
M06-2X is still best for the site difference, with the smallest drift — but the
margin collapses from a factor of ~300 on barriers to a factor of ~1.8 here.

**Practical reading:** switching functional would cut the site-difference method
error from about 3.25 to about 1.82 kcal/mol. Real, but 1.82 still swamps a
signal of order 1 kcal/mol. **Not a reason to redo completed work**, and not a
rescue for the underlying resolution problem.

### CORRECTION: these rows are the wrong error bar for a one-step site comparison

I recommended the 3.25 figure to the selectivity lane as its method-error bar.
**That was wrong twice over**, caught by A1:

| Value | What it actually is | Why it was wrong |
|---|---|---|
| 3.25 | PBE0-D3/**def2-TZVP**, primary→tertiary | wrong basis *and* wrong comparison |
| 4.13 | PBE0-D3/def2-SVP, primary→tertiary | right basis, wrong comparison |
| **0.928** | PBE0-D3/def2-SVP, **secondary→tertiary** | **the correct one** |

That lane works at def2-SVP, so 3.25 was the wrong row — and the flattering one.
More importantly, primary→tertiary is a **two-step** substitution change, while
bridgehead-versus-methylene is tertiary-versus-secondary, **one step**. The right
bar is the §11 calibration, which was computed for exactly that comparison at
exactly that level, at their request. **I had the right number and pointed at a
different one.**

**The error is strongly non-uniform across substitution steps**, which sharpens
the point rather than merely rescaling it:

    primary -> tertiary   (two steps)   4.13
    secondary-> tertiary  (one step)    0.928
    implied primary -> secondary        3.20

Most of the method error lives in the step involving **methane**, which is the
only species here with no alkyl substituent at all — an outlier rather than one
point on a trend. A comparison that avoids methane avoids the bulk of the error.
(This decomposition assumes the two errors add across steps, which two
measurements cannot establish.)

**Consequence: the resolution verdict improves materially.** With a 0.928 bar a
site difference of order 2 kcal/mol is resolvable and one of order 1 is not.
With 3.25 or 4.13, nothing that lane could measure would survive. The direction
is unchanged — PBE0-D3 still overstates the difference by 45% and still flatters
selectivity.

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
