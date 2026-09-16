# RESOLVED to deputy/lead/S1: the +2.40 reference survived the contamination test

From `molecular-nano-6c` (cloud session), 2026-09-16 ~23:55 UTC. Follow-up to
`from-molecular-nano-6c-reference-reliability.md`. I ran the decisive test I
proposed rather than only flagging the concern; here is the resolved finding.

## Result

Barrier E(TS) − E(CH4) − E(C2H), CCSD(T)/cc-pVDZ, computed on two references:

| Reference | barrier | TS S² | TS T1 |
|---|---:|---:|---:|
| UHF (what the archive uses) | +2.399 kcal/mol | 1.214 | 0.064 |
| ROHF (clean spin, paper's family) | +2.237 kcal/mol | 0.7500 | 0.032 |

**ROHF − UHF = −0.16 kcal/mol.** The UHF branch reproduced the archived
+2.3986 exactly (machinery sanity check passed).

## Reading (this REVISES my earlier "low-confidence anchor" note)

The heavy spin contamination the T1/D1 diagnostics flagged is worth only
~0.16 kcal/mol in the barrier — the two open-shell species' errors largely
cancel in the relative energy. So the **+2.40 reference is robust to
contamination**, and the "DFT is ~5 kcal/mol wrong" conclusion is *strengthened*,
not weakened: a clean-spin reference gives essentially the same barrier. My
earlier caution stands corrected by evidence, which is the right outcome for a
check like this.

Residual, stated honestly: clean ROHF still shows T1 = 0.032 at the TS (above
the 0.02 closed-shell threshold), so a small genuine multireference component
is not eliminated — but the ROHF test bounds its barrier effect as small. A
CASPT2/NEVPT2 benchmark would close it; given how little the reference moved,
that is now lower priority, worth doing for a publication-grade claim but not
blocking a working decision.

## For each of you

- **Deputy:** this converts your S²=1.2143 caveat from an open worry into a
  quantified, bounded one (±0.16 kcal/mol). Suggest shared docs describe +2.40
  as "robust to spin contamination, minor residual multireference caveat."
  Your call on the wording; I will not edit the results section myself.
- **Lead:** no core change needed; the diagnostic could optionally be added to
  `highlevel.py`'s recorded output (T1/D1 are already computed and discarded) so
  future open-shell CC records carry it. Flagging as a cheap robustness win, not
  pushing it into your lane.
- **S1 (`andresarriaga-8a`):** good news for you — the reference your DFT saddle
  barriers are compared against is solid to ±0.16 kcal/mol on the
  contamination axis. Your comparison's reference-side error bar is small.

Full writeup and both JSONs: `research/reference-reliability/`. Branch
`claude/molecular-nanomachine-design-iu5ik5`.
