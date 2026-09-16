# North star: the one capability that gates the whole project

Written 2026-09-17 ~00:20 UTC by `molecular-nano-6c` on the user's direct
instruction to *guide all the agents toward this*. It is an orientation for
every lane, offered to the lead (Codex root) to endorse or amend in
`ROSTER.md`/`AGENT_TASKS.md`. It reassigns no lane and overrides no owner's
science; it says what all the lanes are collectively *for*, so the twenty of us
pull in one direction.

## The thesis, in one sentence

**Every path to a custom molecular structure — the positional-mechanosynthesis
track this repo is built on, and the programmable-assembler track the user just
proposed — is gated by one capability: ranking the competing reactions at the
reactive site to better than ~1 kcal/mol.** Everything else is in service of
that, or is downstream of it.

## Why this is the gate (the numbers, from our own evidence)

1. **Selectivity is a discrimination, not a barrier.** The tool must prefer the
   target hydrogen over ~15 competitors (and over welding to the workpiece, and
   over the back-reaction). That preference *is* the free-energy gap ΔΔG‡
   between the intended pathway and its nearest competitor. A1 has already found
   the tool essentially unselective on the intrinsic energetics — i.e. the gap
   is small. A small gap is exactly what a <1 kcal/mol method is needed to
   resolve.

2. **Programmable assembly obeys the same law.** `research/programmable-assembly`
   shows an N-block assembler yields correct product with probability f^N, and
   the per-step fidelity f is set by that same ΔΔG‡. A 100-unit structure needs
   ~4 kcal/mol of clean discrimination per step (or ~2 with kinetic
   proofreading). The ribosome hits ~10⁻⁴ error with ~2.7 kcal/mol + one
   proofreading stage. Same currency: barrier discrimination.

3. **The default method cannot yet pay in that currency.** PBE0-D3(BJ) is
   ~5 kcal/mol off on the calibration barrier and gets the sign wrong. That
   error is *larger than the entire discrimination window* either paradigm
   needs. A method that cannot rank competing pathways to ~1 kcal/mol cannot
   design either an assembler or a mechanosynthetic sequence.

4. **The yardstick is now trustworthy, so the target is real.** The
   reference-reliability lane showed the CCSD(T) +2.40 reference moves only
   0.16 kcal/mol between contaminated-UHF and clean-ROHF references. So "get
   within 1 kcal/mol of CCSD(T)" is a well-defined, defensible target, not a
   chase after a wobbly number.

**Conclusion: the project's single most valuable output is a method (or
correction scheme) that ranks the competing pathways at the reactive site to
<1 kcal/mol, and every barrier we report should be a ΔΔG‡ against its real
competitors, not an absolute number in isolation.**

## What this means for each lane (reframing, not reassigning)

- **S1 (saddle search):** a located saddle is necessary but not the deliverable.
  The deliverable is the *gap* between the intended abstraction saddle and the
  nearest competing saddle (wrong-site abstraction, tool–workpiece welding, back
  reaction). Report ΔΔG‡ between pathways, each on its own verified saddle.
- **A1 (selectivity):** you own the heart of the north star. Frame the
  unselectivity result as a discrimination window: target-vs-each-competitor
  ΔΔG, against kT (0.59 kcal/mol) and against the ~2–4 kcal/mol an assembler
  needs. That table is the project's key result.
- **A2 (feasibility/cost):** the real budget is not one path — it is one path
  *per competitor that must be ranked*. Cost the discrimination, not a single
  barrier. Your reduced-handle fidelity work feeds directly into whether the
  discrimination survives a tractable model.
- **Method/guess lanes (D1, D2, S2, S3, state-scan):** guess dependence is a
  >1 kcal/mol error source (propynyl: 11.2 kcal/mol between guesses at PBE0-D3).
  That directly threatens the <1 kcal/mol target. Frame your work as *protecting
  the discrimination floor*: any uncontrolled spread larger than ~1 kcal/mol
  invalidates a ranking.
- **Vibrational/tunnelling/ZPE lanes:** ZPE and tunnelling corrections are
  several kcal/mol and can *differ between competing pathways* (different
  transition-state frequencies, different H-transfer tunnelling). They belong
  inside the ΔΔG‡, not only in the absolute barrier. A correction that is equal
  for all competitors cancels; one that differs is decision-relevant. Say which.
- **Reference-reliability + programmable-assembly (this session):** the yardstick
  is validated; next is whether *DFT itself* can hit <1 kcal/mol on ΔΔG (not
  absolute energies), functional by functional, against CCSD(T). And the
  assembler model turns any validated ΔΔG‡ into a concrete yield. Both consume
  the discrimination.
- **Audit / bundle / workbench / integration lanes:** keep doing what makes the
  discrimination *auditable* — provenance, evidence bundles, honest labels. A
  ranking nobody can reproduce is not a result.

## The shared artifact this suggests

A single **discrimination ledger**: one row per competing pathway at the
reactive site (target abstraction, each wrong-site abstraction, welding, back
reaction), each with its ΔΔG‡ vs the intended pathway, the method, the
guess/state controls, the vibrational/tunnelling delta, and an uncertainty. The
project succeeds when that ledger has entries whose uncertainties are below the
gaps. That is the concrete form of "scientifically perfect" for this problem.
Proposed owner: the lead, assembling from S1/A1/A2 + the accuracy lanes; I will
maintain the assembler-yield column that reads off it.

## What this does not change

Honesty rules stand: converged ≠ validated; no absolute-barrier claim without a
verified saddle; guess scans for open-shell work at every geometry (DFT
included — propynyl proved it); deputy/lead sign-off before a number enters
shared results. The north star sharpens what we compute; it does not lower the
bar on how we report it.
