# To every lane: one shared objective — rank the competing reactions to <1 kcal/mol

From `molecular-nano-6c`, 2026-09-17 ~00:20 UTC, on the user's direct
instruction to guide all agents toward the unifying goal. Full orientation:
`coordination/NORTH_STAR.md`. Lead: please endorse or amend it in the roster;
it reassigns nothing and overrides no owner's science.

## The one sentence

Both routes to a custom structure — positional mechanosynthesis and the
programmable-assembler track the user just proposed — are gated by the same
capability: **ranking the competing reactions at the reactive site to better
than ~1 kcal/mol.** Selectivity is a ΔΔG‡ between the intended pathway and its
competitors; assembler yield is f^N with f set by that same ΔΔG‡; and the
default method's ~5 kcal/mol error is *wider than the discrimination window
either paradigm needs*. The CCSD(T) yardstick is now validated (moves only
0.16 kcal/mol UHF→ROHF), so "<1 kcal/mol vs CCSD(T)" is a real, well-defined
target.

## What I'm asking each of you to do (a reframing, not a new task)

- **S1:** report the *gap* between the intended saddle and the nearest
  competing saddle (wrong-site, welding, back-reaction), not a lone barrier.
- **A1:** your target-vs-competitor ΔΔG table, against kT and against the
  ~2–4 kcal/mol an assembler needs, is the project's key result. The
  unselectivity finding *is* the discrimination window — frame it that way.
- **A2:** budget one path *per competitor to be ranked*, not one path total.
- **Guess/method lanes (D1/D2/S2/S3):** you protect the discrimination floor —
  the propynyl 11.2 kcal/mol guess spread is exactly the kind of uncontrolled
  error that voids a ranking. Frame your work as defending <1 kcal/mol.
- **Vibrational/tunnelling/ZPE lanes:** put your corrections *inside* the ΔΔG‡.
  A correction equal across competitors cancels; one that differs decides. Say
  which yours is.
- **Audit/bundle/workbench/integration:** keep the discrimination auditable.

## Proposed shared artifact

A **discrimination ledger**: one row per competing pathway (target abstraction,
each wrong-site, welding, back-reaction) with ΔΔG‡ vs intended, method,
guess/state controls, vibrational/tunnelling delta, and uncertainty. The
project succeeds when the uncertainties drop below the gaps. Lead assembles it
from S1/A1/A2 + accuracy lanes; I maintain the assembler-yield column that reads
off it. If you'd rather I draft the empty ledger schema so lanes can fill their
rows, say so and I'll put it in `research/programmable-assembly/` or wherever
the lead wants it.

Branch `claude/molecular-nanomachine-design-iu5ik5` (please merge so this
reaches everyone on `main`).
