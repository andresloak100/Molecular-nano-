# Is the in-house CCSD(T) reference itself trustworthy?

Owner: `molecular-nano-6c` (cloud session; compute ran off-host, no Mac load).
Lane opened 2026-09-16 ~23:15 UTC. New files only; no archived record or core
module touched. This is a scientific-reliability check, so its reading is
offered to the forensics **deputy** and **lead** for sign-off before anything
here is promoted into shared results.

## Why this lane exists

The project's sharpest quantitative claim is that every DFT functional gets
this reaction's barrier wrong — submerged by several kcal/mol against an
in-house CCSD(T)/cc-pVDZ reference of **+2.40 kcal/mol**, evaluated on the
supplied collinear transition structure. That claim is only as good as the
reference. CCSD(T) is a **single-reference** method: it is reliable only when
one Hartree-Fock determinant dominates. Two warning signs were already on the
record but never quantified as a reliability statement — the transition
structure's UHF reference has S² = 1.21 (vs 0.75 for a clean doublet), and
multireference character was "suspected."

The standard coupled-cluster reliability diagnostics settle it, and they are
free byproducts of the CCSD amplitudes the reference already computes but
never saved (`highlevel.py` computes `t1, t2` and discards them). This lane
recomputes the five cc-pVDZ species on the archived geometries at the lowest
SCF solution (guess `atom`, the one the forensics lane validated) and reads
the diagnostics off.

## Result

| Species | role | T1 diag | D1 diag | max \|t2\| | HF S² |
|---|---|---:|---:|---:|---:|
| methane | closed-shell control | 0.0049 | 0.009 | 0.033 | 0.00 |
| acetylene | closed-shell control | 0.0110 | 0.028 | 0.090 | 0.00 |
| methyl radical | doublet, well-behaved | 0.0088 | 0.013 | 0.042 | 0.76 |
| **ethynyl radical** | reactant radical | **0.0836** | **0.153** | 0.104 | 1.22 |
| **transition structure** | the +2.40 anchor | **0.0638** | **0.153** | 0.106 | 1.21 |

Thresholds (literature, not invented): T1 > 0.02 (Lee & Taylor 1989) and
D1 > 0.05 (Janssen & Nielsen 1998) mark single-reference CCSD(T) as
questionable for closed-shell systems.

**The transition structure's T1 is 0.064 — three times the closed-shell
threshold and 5.8x the closed-shell controls. Its D1 is 0.153, three times
its threshold. The ethynyl radical is worse still (T1 = 0.084).** The two
closed-shell species and the ordinary methyl radical are all clean. So the
single-reference diagnostic is bad for exactly the two open-shell species the
+2.40 reference is built from, and good for everything else.

## What this does and does not mean

**Read honestly, because the diagnostic conflates two effects.** A UHF-based
T1 is inflated both by genuine near-degeneracy (multireference character) and
by spin contamination, and here S² = 1.21 shows heavy contamination. This
lane cannot cleanly separate the two from a single UHF reference — and saying
otherwise would repeat the kind of overclaim the results log has already had
to retract once. What is solid:

1. The in-house **UHF/UCCSD(T)** reference is evaluated on determinants that
   are strongly spin-contaminated and flagged by both single-reference
   diagnostics. It is a **lower-confidence anchor** than "gold standard"
   language implies — not a demonstrated wrong number.
2. The published paper used **ROHF-based RCCSD(T)**, which does not carry this
   spin contamination. That our UHF number nonetheless matched their +2.4 to
   two figures is reassuring but not a proof the UHF reference is sound; it
   could be partly cancellation.
3. This changes the confidence on the **yardstick**, not any DFT result. DFT
   being ~5 kcal/mol below a *low-confidence* +2.40 is a weaker statement than
   DFT being ~5 kcal/mol below a rock-solid one.

## What "scientifically perfect" needs here (the concrete gap)

To turn +2.40 from a lower-confidence anchor into a defensible reference:

1. **ROHF-based RCCSD(T)** on the same geometries, to remove spin
   contamination and compare against the UHF/UCCSD(T) number. Agreement would
   largely vindicate the reference; divergence would move it. Cheap; runnable
   off-host here.
2. **A multireference benchmark** (CASSCF then NEVPT2 or CASPT2) on the TS and
   ethynyl radical, with an active space chosen from the natural-orbital
   occupations. This is the only thing that actually resolves near-degeneracy
   vs contamination. More expensive but still small-molecule; a candidate for
   this container or a GPU trial (G1).
3. Until one of those exists, shared docs should call +2.40 an in-house
   single-reference estimate with a flagged reliability caveat, not a gold
   standard.

## Files

- `diagnose_reference.py` — recomputes the five species and extracts T1/D1,
  largest t2, and S². Reproducible; stdlib + PySCF; reads archived geometries
  read-only.
- `reference-reliability.json` — machine-readable result with thresholds,
  per-species records, the TS-vs-control reading, and explicit limits.

## Boundaries

Read-only against `data/reference/` geometries; writes only in this lane; no
archived CCSD record or core module modified. Interpretation is submitted for
deputy/lead sign-off, not self-promoted into the results section.
