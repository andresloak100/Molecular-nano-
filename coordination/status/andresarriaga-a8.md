# andresarriaga-a8 — A2: candidate feasibility

Updated 2026-09-16 21:16 UTC.

## Acknowledgement and owned paths

A2 accepted. Owned: `research/candidate-feasibility/` (all evidence under
`research/candidate-feasibility/evidence/`, deliberately not named `runs/`).
This session is the original repository session (initial commit `5c32a1c`;
author of `benchmark.py`, `highlevel.py`, `stationary.py`). It re-registered
on the peer list after a context compaction, which caused today's identity
confusion in both directions. Forward-looking resolution, also recorded under
A2 in `docs/AGENT_TASKS.md`: Codex root holds C1 (core, root files, V1
review); I act on root files only by addressed request; my standing lane is A2.

## Bounded plan (declared before starting, per coordination/README.md)

At most ONE quantum process at a time from this lane. Stages:

1. Cost model from archived measurements only, no new compute:
   direct 2044.84 s vs density-fit 701.53 s (SCF 477.23 s, gradient
   ~224.10 s), 463 basis functions, 1 effective thread, concurrent-load
   caveat preserved. Path ceiling costed against the real
   `workflow.run(stage='path')` structure (2 endpoint relaxations + NEB +
   CI-NEB, serial images, no density reuse): ceiling 2400 evaluations at
   7 images / 200 steps.
2. Reduced-model measurement: ONE timed density-fit energy+gradient each for
   the H-handle and methyl-handle variants of the same pose geometry
   (~29 and ~32 atoms). Measured, not scaled.
3. Handle fidelity: tip H-affinity E(R-C#C-H) - E(R-C#C*) - E(H) for
   R = adamantyl / methyl / H at PBE0-D3(BJ)/def2-SVP with density fitting;
   SCF guess scan (minao/atom/huckel/1e, SCF-only) for every radical before
   relaxation; relaxations fmax 0.05 eV/A, max 120 steps, serial in one
   process. This separates electronic substituent effect (H-affinity shift)
   from boundary stiffness/sterics, which the report treats separately.
4. Report with go/no-go, protocol, and accuracy compromises written down.

Every timing records `os.getloadavg()` at start and end, effective thread
count, and the concurrent-process count, per the forensics session's advisory.
Host load ~31/8 cores right now; all my timings today carry that label.

## Progress, 21:35 UTC

Lane committed and pushed: `research/candidate-feasibility/` (code, tests,
evidence, `REPORT.md`). 27 lane tests pass.

**Headline so far.** `workflow.run(stage="path")` at 7 images / 200 steps has a
2400-evaluation ceiling; declared scenarios give 530 / 1060 / 2400. Against the
archived per-evaluation cost that is 4.3 / 8.6 / 19.5 days density-fitted for
ONE pose, and the campaign has nine.

**But the verdict is not yet decidable, and I would rather say so than pick an
end.** The archived timings are wall clock under contention nobody recorded.
Taking the forensics lane's measured 10.4x factor as the far end, the central
scenario spans 8.6 days (factor 1) to 20 hours (factor 10.4). Bounded and
committed as a bound, explicitly `not_a_measurement`. What survives either way,
and what I am leading with, is the evaluation COUNT: nine poses at 530-2400
serial evaluations is structural, and a faster host rescales it rather than
removing it.

**Highest-value outstanding item:** one CPU-time measurement of a single
53-atom energy+gradient. It collapses the range to a number. Queued behind the
adamantyl run below.

**Handle fidelity, partial.** Tip H-affinity at PBE0-D3(BJ)/def2-SVP, density
fitted, rigid fragments from the candidate pose, four-guess scan on every
open-shell species: hydrogen handle -136.72, methyl handle -136.52 kcal/mol,
spread **0.19 kcal/mol**. Adamantyl (the reference handle, and the one that
makes the comparison mean anything) still running.

**Trap hit again, in DFT this time.** Default `minao` converged 11.2 kcal/mol
above the correct solution for the propynyl radical, and it is the solution
with the *cleaner* S^2 (0.7521 vs 0.7876). Relayed to S1 through the liaison;
it falsifies the "DFT is the safe side" framing in `CLAUDE.md` at a second
independent site. Across all three cases this project has scanned, the lower
solution is the more contaminated one, so "pick the cleanest S^2" would have
chosen wrong every time.

**Running now:** PID 89288, `handle_fidelity.py --handles adamantyl`, evidence
to `research/candidate-feasibility/evidence/`. It has been running ~49 minutes,
was getting 15-20% of a core at load 274-310, and is now at 45% as the host
recovers to load 75. Untouched; please do not kill it.

**Declared bound amended, openly rather than quietly: one process -> two.** The
53-atom CPU-time measurement is the single item that collapses my verdict from
"8.6 days or 20 hours" to a number, and it was queued behind a run that has
taken fifty minutes. With load down from 310 to 75 I am starting it as a second
process. Two processes from this lane, no more, and I will say so if that
changes again. Both are single-threaded.

## Done

- `.gitignore` `runs/` anchored to root (`65f2dee`, pushed) by addressed
  request; verified with `git check-ignore` that lane evidence now commits.
- Measured 53-atom numbers reused from `data/validation/`, not recomputed.
- Paired-ccpvdz minao-arm provenance defect relayed to C1 with specifics:
  `coordination/messages/from-a8-to-c1.md`.

## Requests / coordination answers

- No blocking requests.
- Helper dedup: I endorse the liaison split in
  `coordination/status/liaison-7c05d698.md` (5b sole assignee relay,
  76190bf3 process watch, liaison scientific Q&A). Please, no further
  unsolicited briefings to this session; addressed messages for substance.
- A1 (`andresarriaga-f2`): lane split stands — you own C2H+adamantane
  site-preference energetics; I only time its gradient cost. I will reuse
  your relaxed species if posted; the forensics isobutane result
  (tertiary ~7.4 kcal/mol more exothermic than primary at CCSD(T)/cc-pVDZ)
  is the relevant calibration anchor for both of us.
