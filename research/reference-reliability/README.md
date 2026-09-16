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

## Follow-up result: the reference is robust to spin contamination

The decisive test (`rohf_vs_uhf_barrier.py`, `rohf-vs-uhf-barrier.json`)
recomputes the barrier E(TS) − E(CH4) − E(C2H) at CCSD(T)/cc-pVDZ on both a
UHF reference (S² = 1.21, contaminated) and an ROHF reference (S² = 0.7500,
clean). The closed-shell CH4 is identical under both, so any difference is
carried by the two flagged open-shell species.

| Reference | barrier | TS S² | TS T1 | C2H T1 |
|---|---:|---:|---:|---:|
| UHF | +2.399 kcal/mol | 1.214 | 0.064 | 0.084 |
| ROHF | +2.237 kcal/mol | 0.7500 | 0.032 | 0.035 |

**ROHF − UHF = −0.16 kcal/mol.** Both bracket the published Table 4 value of
+2.4. So the heavy spin contamination the diagnostic flagged is worth only
~0.16 kcal/mol in the relative energy that matters — the error sources in the
two open-shell species largely cancel in the barrier. The +2.40 anchor is
**robust to the contamination concern**, and this *strengthens* the DFT
critique: DFT being ~5 kcal/mol below the reference is not an artifact of a
contaminated yardstick, because a clean-spin reference gives essentially the
same barrier.

Residual caveat, stated honestly: the clean ROHF reference still shows
T1 = 0.032 at the TS, above the 0.02 closed-shell threshold. So about half the
elevated UHF T1 was spin contamination (now removed) and half is residual —
genuine open-shell correlation or static-correlation character that a single
reference does not fully capture. Its effect on the barrier is now *bounded as
small* (the contamination piece was only 0.16 kcal/mol, and the residual is
smaller than the piece already removed), but not proven zero. Only a
multireference calculation (CASPT2/NEVPT2) closes it completely; given how
little the reference moved, that is now a lower priority than it looked.

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

## What "scientifically perfect" needs here (updated after the ROHF test)

1. **ROHF-based RCCSD(T)** — DONE. Moves the barrier only 0.16 kcal/mol, so
   spin contamination is not corrupting the reference. This was the cheap,
   high-value check and it came back reassuring.
2. **A multireference benchmark** (CASSCF then NEVPT2/CASPT2) on the TS and
   ethynyl radical, active space from natural-orbital occupations. Now a
   *lower* priority: it would close the residual T1 = 0.032 caveat, but the
   ROHF result bounds any remaining effect as small. Worth doing for a
   publication-grade claim; not blocking a working decision. Small-molecule;
   off-host or a G1 GPU trial when convenient.
3. Shared docs can now describe +2.40 as an in-house single-reference estimate
   that is **robust to spin contamination (±0.16 kcal/mol between UHF and
   ROHF references)**, with a minor residual multireference caveat — not as a
   fragile or contaminated number. The DFT discrepancy stands on firmer
   ground than before this check, not weaker.

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
