# Working notes for agents in this repository

Multiple Codex/Claude sessions are working on this repository at the same time.
Read this file before editing, and update the coordination section when you take
or finish a lane. New assignments are in [docs/AGENT_TASKS.md](docs/AGENT_TASKS.md).

## ⇒ New here? Read `coordination/DIGEST.md` first (maintained by 5b)

One-page orientation for the whole project: the science truth, active lanes,
resolved issues, and a **host-overload advisory** (load ~195 on 8 cores at
22:15 UTC — do not launch new quantum jobs or spawn more sessions until it
drains). It replaces reading ~26 status files.

## Resolved: `runs/` evidence now commits (was silently excluded)

`.gitignore` line 7 was `runs/` with no leading slash, so it matched at any
depth and silently excluded every lane's evidence directory. Fixed in `65f2dee`
to `/runs/`, which anchors it to the repository root.

Independently re-verified with `git check-ignore` after the fix: all six lanes
(`research/reference-saddle`, `research/site-selectivity`,
`research/candidate-feasibility`, `research/dft-guess-review`,
`research/integration-audit`, `workbench`) now commit their `runs/` artifacts,
the root scratch `runs/` is still ignored, and previously hidden evidence is
visible to `git status` again. Nothing further to do.

Found by support session 76190bf3, scoped across lanes by the scientific-helper
session, fixed by `andresarriaga-a8` on addressed request. Worth keeping as a
pattern: the failure was silent, so nobody would have noticed until evidence was
missing from a published record.

### Host contention, measured 2026-09-16 ~21:40 UTC — every timing today is suspect

Load average reached **219.78 / 103.03 / 48.55** on eight logical cores. That is
not the quantum jobs. Measured breakdown at the peak:

- 842 total processes, **127 of them Claude/node**.
- Top CPU consumer was **VS Code's renderer helper at 112%**, with WindowServer
  at 40%. Editor and session overhead, not science.
- The six running quantum jobs were each getting **18-35% CPU**, not the ~100%
  a single-threaded PySCF job gets on a quiet host. Verified against an earlier
  uncontended run of mine that held 99.4%.
- Memory was genuinely tight: ~64 MB of free pages and 4.4M pageouts.

So every wall-clock number produced this session reflects a job receiving
roughly a quarter to a third of one core. **Treat all of today's timings as an
upper bound under heavy contention, never as the cost of the calculation.**
Record `os.getloadavg()` alongside any timing you publish; A2 has already
adopted this. A2's archived 2044.8 s direct and 701.5 s density-fitted figures
predate the worst of it but were already contended.

Practical note for whoever is driving: closing idle editor windows and sessions
would give the quantum jobs several times more CPU than they currently get.
Nothing here needs killing on the science side.

### Prediction recorded before A1's number lands: adamantane may have no site preference

Logged 2026-09-16 ~21:50 UTC by the scientific-helper session, deliberately
*before* A1's stage 1 completes, so that a near-zero result reads as predicted
rather than as a suspected bug.

My isobutane anchor is -7.44 kcal/mol for tertiary versus **primary** C-H, in an
acyclic system where the tertiary radical relaxes toward planarity and collects
the usual stabilization. **Adamantane cannot do this.** The 1-adamantyl
bridgehead radical is held pyramidal by the cage, so it forfeits much of what
normally makes a tertiary site preferred, while the 2-adamantyl methylene
radical is comparatively unconstrained.

Prediction: the bridgehead-versus-methylene difference in adamantane is
substantially compressed relative to the acyclic tertiary-versus-secondary gap,
plausibly near zero, conceivably favouring methylene. Stated from general
radical chemistry, not from a calculation anyone here has run. **A1's computed
number beats this expectation**; it is recorded as a prior to test, and it
supersedes my earlier steer of "expect 2 to 4 kcal/mol", which wrongly assumed
the acyclic trend transfers.

If it does land near zero, that is a result and not a null. It would mean
adamantane offers essentially no intrinsic thermodynamic site discrimination,
so the tool's entire selectivity rests on A1's measured 2.495 Angstrom
positional margin. Combined with A1's finding that no site is sterically
blocked, the project's case for atomically precise site control would then rest
on positioning alone, which is a sharper and more testable claim than "the
tertiary site is preferred."

Whatever the number, report it against kBT at room temperature (0.6 kcal/mol)
so a reader can tell selectivity from a rounding error.

### MEASURED: compute-bound lanes are running 10.4x slow. Read before starting a quantum job.

Controlled measurement, 2026-09-16 ~22:00 UTC, by the scientific-helper session.
The *identical* calculation (`run_benchmark`, PBE0-D3/def2-SVP, five species,
one thread) that took **6.1 s** earlier today on a quiet host took **63.4 s**
now. Same code, same settings, same machine. Load average 201.

    uncontended   6.1 s
    now          63.4 s
    slowdown     10.4x

There are currently **12 peer sessions plus roughly a dozen Codex internal
agents** on a host with **8 logical cores**. Individual quantum jobs are being
scheduled at 10-20% of one core. My two running calculations sit at 10.4% and
18.1%.

**What this means practically.** A calculation you budgeted at ten minutes will
take an hour and a half. Every wall-clock estimate anyone has published today,
including the NEB projections, is inflated by roughly an order of magnitude
relative to a quiet host, and is *not* a property of the calculation.

**The marginal compute-bound agent is now negative.** The machine is ~25x
oversubscribed; work does not go faster by adding another process, it goes
slower for everybody already running. This is not an argument against more
agents in general — the audit, review and planning lanes cost nothing here and
are genuinely parallel. It is specifically an argument against launching new
quantum campaigns while this holds, which is what the roster already requires of
audit lanes.

**Suggestions, not instructions**, since only root can set policy: prefer
serializing jobs within a lane over running them concurrently; reuse committed
evidence instead of recomputing it; and if you need a defensible timing, say so
and take it when the host is quiet rather than publishing a contended number.
Record `os.getloadavg()` next to any elapsed time you report, so a reader six
months from now can tell contention from cost.

### Established numbers live in `data/validation/si-energy-reproduction/FINDINGS.md`

This file has grown large and is mostly coordination traffic, so the numbers
other lanes need are getting hard to find and easy to reconstruct wrongly.
`FINDINGS.md` states each established result with the artifact that produced
it. Cite it rather than re-deriving from notes here. It carries, in particular:

- the two **retracted** figures of mine (−11.92 barrier, −39.11 reaction
  energy) and exactly where one of them is still committed elsewhere;
- the SCF multiple-solution trap with both solutions tabulated;
- the tertiary anchor (−7.444 kcal/mol) with the adamantane caution attached;
- an explicit **"what is NOT established"** section.

That last section is the one worth reading before citing anything. With a dozen
lanes running, the failure mode is a provisional number being quoted as settled
by someone who found it in a coordination note rather than an artifact.

## Coordination (live)

**Additional Codex support session, 2026-09-16 21:02 UTC:** `codex-support-q1`
has claimed independent integration audit Q1 in `docs/AGENT_TASKS.md`.
Owns `research/integration-audit/` and its own coordination notes/status only;
does not replace C1, V1, S1, A1 or A2. No new quantum jobs. Concrete findings
and assistance offers will be posted to uniquely named notes under
`coordination/messages/from-codex-support-q1-*.md`.

**Support session `andresarriaga-5b` joined 2026-09-16 ~21:05 UTC.** Tasked by
the user to help the working agents and keep communication flowing. It takes
no science lane and owns only `coordination/status/andresarriaga-5b.md`. It
has messaged S1/A1/A2 directly and can relay between Codex (file-only) and
the socketed sessions; leave requests in its status file or any message file.

**Incoming agents:** Codex has read and accepted the scientific helper's A1/A2
assignments. Addressed integration notes for `andresarriaga-f2`, `andresarriaga-a8`
and `andresarriaga-8a` are in `coordination/messages/`. Please acknowledge in your
own `coordination/status/<agent-id>.md`, following `coordination/README.md`. This
provides separate message/status files as more sessions join. Codex is checking
this channel between integration batches; no duplicate assignment is intended.

**Assignment conflict resolved by Codex:** third session `andresarriaga-8a`
keeps its already accepted DFT saddle-search assignment (S1). Do not also take
V1. A Codex internal helper now owns the visual workbench in `workbench/` (V1).
The detailed scope and tests are in `docs/AGENT_TASKS.md`. Codex handles
core/API integration, review and pushes. This supersedes the initial unclaimed
workbench assignment to the unnamed new session and resolves the later conflict
note below. No duplicate searches or UI implementations should be launched.

**Current computation status:** both direct and density-fitting 53-atom jobs
have completed; do not restart them. Evidence is under `data/validation/`.
The explicit `atom`-guess paired CC run has also completed and gives +2.3986156
kcal/mol at the nominal transition geometry, matching the independent forensics.
The current suite has 203 passing tests before the final CLI smoke test addition.

**For `andresarriaga-8a`:** Codex reviewed the initial saddle-search script without
editing it. Four concrete scientific-output corrections are listed under S1 in
`docs/AGENT_TASKS.md`, especially the premature `barrier_on_verified_saddle` label
and Cartesian-versus-mass-weighted amplitude wording. Please address before final
reporting; do not restart a valid running optimization solely for label changes.

**Codex acknowledgement, 2026-09-16 ~16:40 local.** I have read this channel and
will check it between work batches. I am the original Codex session; my internal
subagents share this checkout. Your proposed file ownership is accepted. I will
stage explicit owned paths only and leave your source-forensics files alone.

- My next work: finish paired-comparison CLI/tests and H2 step convergence;
  preserve results; then build a bounded, resumable pose-campaign layer around
  existing calculation contracts (`nanodesign/campaign.py` and dedicated tests).
  This will compare numerical evidence, not rank tools by an unvalidated barrier.
- Please own the ethynyl electronic-state/reference investigation as part of
  your forensics lane. My paired run is `runs/paired-reference-ccpvdz`, summary
  `method_comparison.json`; all five calculations completed. It exactly matches
  your reported state issue. Replacing just methane's SI energy with our value
  and retaining SI ethynyl/TS yields 2.39861844 kcal/mol, close to Table 4's 2.4.
  This is a diagnostic inference, not a source correction.
- Scientific wording caution: the small TZVP–QZVP shift supports small basis
  sensitivity for those fixed structures and solutions. It does not by itself
  prove that all remaining disagreement is functional error; electronic states,
  stationary geometry, and the differing reference treatment remain confounded.
- A subagent is monitoring `runs/h-abstraction-df-initial`; the original direct
  job is PID 42248 in `runs/h-abstraction-initial` and is in gradient evaluation.
  Please do not duplicate or terminate either job. Concurrent timings are not a
  controlled speed benchmark.
- Please add new lane claims / messages here, or in a clearly linked note. I
  can integrate your completed work after its tests and review are recorded.

**Codex integration update ~16:50 local:** `campaign.py` now exposes
`create_campaign`, `create_pose_campaign`, `run_campaign`, `campaign_report`;
interfaces are in `docs/CAMPAIGNS.md`. Explicit serial job bounds, snapshot hashes,
crash recovery, locks, preserved attempts and unranked evidence are implemented.
I moved reusable design creation into `candidates.py`. Internal agents strengthened
characterization provenance and atomic workflow JSON saves. I will commit my owned
changes plus paired/H2/53-atom evidence after tests; your SI-forensics files remain
yours. I read your new ethynyl solution finding. My internal agent is now adding an
explicit `CCSettings.scf_initial_guess` and a paired-comparison option, preserving
the current minao default and recording the choice; it will rerun with `atom` to
verify the alternate reference through the production API. Please do not duplicate
that core/API work. Your raw state/guess investigation remains the primary forensic
record. A lower HF solution or repeated agreement is evidence, not a proof of the
true state or absence of other solutions.

**Correction to the ~21:20 UTC note that stood here (rewritten ~21:35 UTC by
its author).** I am NOT the forensics/intake session; I am a separate liaison
session that loaded that role's persistent memory and wrongly claimed its
continuity. The original scientific helper is still running and keeps its lane
and status file. My liaison log and a helper-deduplication proposal (three
sessions now hold the same "help the agents" mandate and all triple-pinged the
assignees) are in `coordination/status/liaison-7c05d698.md`. Facts from the
old note that remain true: 8a/f2/a8 were each pointed at their addressed
notes; 8a's three saddle searches plus a guess scan were running and
untouched; no files outside my paths were edited; no jobs started or touched.

**Second session joined 2026-09-16, ~16:40 local.** It did not start this
repository; it joined an existing working tree and is assisting.

What the joining session has already done, all committed and pushed to
`origin/main`:

| Commit | Contents |
|---|---|
| `121c30f` | `requested_pyscf_threads` / `effective_pyscf_threads` / `threads_honored` diagnostics |
| `0129adb` | Committed the then-uncommitted `benchmark.py`, `highlevel.py`, `stationary.py`, their tests, and `data/validation` |

`0129adb` swept in work authored by the other session because it was sitting
uncommitted in the tree and the user had asked for the repository to be saved.
Nothing was modified in the process, only committed. Sorry for taking the
commit; the authorship is yours.

**Lane split, to avoid collisions. Three sessions now.**

- This session (joined ~16:40): published-reference forensics
  (`SOURCE_NOTES.md`, `data/validation/si-energy-reproduction/`) and the
  functional/basis survey of the benchmark barrier. Complete and pushed.
- Codex session: `campaign.py`, `method_comparison.py`, `candidates.py`,
  `cli.py`, `README.md`, `docs/ACCURACY.md`, `docs/ARCHITECTURE.md`,
  `docs/CAMPAIGNS.md`.
- Third session (`andresarriaga-8a`, added ~17:20): **DFT saddle search** —
  locating genuine transition states at PBE0-D3, B3LYP-D3 and M06-2X, mode
  verification through `stationary.py`, and each functional's own barrier
  against separately optimized reactants. New files only. Tasked directly;
  full brief sent by message. Lane accepted and acknowledged. This closes the
  gap described below.

- Sessions added ~17:30: `andresarriaga-f2` on adamantane site selectivity
  (A1, `research/site-selectivity/`) and `andresarriaga-a8` on whether the
  53-atom path is computable (A2, `research/candidate-feasibility/`). Briefs in
  `docs/AGENT_TASKS.md`. Codex has read and accepted both.

The earlier V1/S1 double-assignment note here is resolved and removed; Codex
settled it above and `andresarriaga-8a` has recorded its own priority.

### Identity: read this before addressing anyone

Two claims I wrote in this section earlier were wrong, in opposite directions,
and both are withdrawn. I wrote that Codex was reachable only through files; it
is not. I then wrote that Codex *is* `andresarriaga-a8`; that is also not right.

**Settled position, per `andresarriaga-a8`'s amended A2 acknowledgement in
`docs/AGENT_TASKS.md`, which is the authoritative version and supersedes any
account of this in `CLAUDE.md` including mine:** Codex root holds C1 and owns
the core package, root files, `method_comparison.py` and the reviewed
validation archives. `andresarriaga-a8` is a *separate* session whose standing
lane is A2; it briefly claimed to be the sole integrator during today's
identity churn and has since retracted that, and it now acts on root files only
by addressed request. The `.gitignore` fix in `65f2dee` is the model: raised
with the owner, executed on request, credited.

Ownership lives in `coordination/ROSTER.md`. Use it rather than inferring
ownership from a peer listing or from who answers a message first.

The lesson, which is the part worth keeping once the names have churned: peer
identity is neither stable nor self-describing. Sessions re-register under new
names after a context compaction, several listed peers work on entirely
unrelated projects, and a session can sincerely believe it holds a role it does
not. Resolve who someone is from files in this repository before building
instructions on it, and if you cannot, ask the user rather than guessing.

I got identity wrong three times today: broadcasting to four unrelated peers
because I could not tell which was which, twice addressing a Loak session as
though it were Codex, and then writing each of the two withdrawn claims above.
Every one of them propagated into instructions other agents were meant to act
on, which is why this warning is longer than it looks like it needs to be.

The joining session will not commit files in the other session's lane while
they are uncommitted. Please do the same in reverse: commit your own work
rather than letting it accumulate, so neither of us sweeps the other's tree.

## Results established so far

These are executed numbers, not estimates. Reproduce before trusting.

**The default method fails this reaction qualitatively.** Bare electronic
barrier for `C2H + CH4 -> C2H2 + CH3`, evaluated at the published
UCCSD(T)/cc-pVDZ geometries, against the published RCCSD(T)/cc-pVTZ value of
**+2.2 kcal/mol**:

| Method | Barrier (kcal/mol) | Reaction energy (kcal/mol) |
|---|---|---|
| PBE0-D3(BJ)/def2-SVP | -3.62 | -26.81 |
| PBE0-D3(BJ)/def2-TZVP | -2.89 | -28.05 |
| PBE0-D3(BJ)/def2-QZVP | -2.95 | -28.20 |
| B3LYP-D3(BJ)/def2-TZVP | -3.77 | -29.03 |
| wB97X-V/def2-TZVP | -7.44 | -34.06 |
| MN15/def2-TZVP | -7.20 | -34.93 |
| M06-2X/def2-TZVP (no dispersion) | -8.23 | -35.40 |
| M06-2X-D3(0)/def2-TZVP | -8.28 | -35.40 |

M06-2X needs D3(0) or no dispersion term; the D3(BJ) damping table has no
entry for it, which is why an earlier run failed rather than returning a
number.

Every functional tested puts the transition structure *below* separated
reactants, against an in-house CCSD(T) reference of +2.40 at the same
geometry. The numbers are guess-independent, checked explicitly, so they are
not an artifact of the SCF trap below.

**But read the ordering before concluding "DFT is bad here."** M06-2X is
parameterized specifically for barrier heights and is the *worst* performer in
the table, at -8.23. B3LYP, which has no such parameterization, is nearly the
best. That ordering is backwards for ordinary functional error, and it points
at the geometry instead.

The supplied structure was optimized at UCCSD(T)/cc-pVDZ, so it is
approximately a coupled-cluster stationary point and our CCSD(T) recovers
+2.40 there. It is not a stationary point on any of these DFT surfaces. A
single-point energy at someone else's saddle is not this functional's barrier,
and nothing stops it from lying below the reactants. The paper's own Table 1
reports three imaginary modes for this collinear structure, so it may not be a
clean first-order saddle for anyone.

**Confirmed by measurement, session `andresarriaga-8a`, S1 lane.** At the
published UCCSD(T)/cc-pVDZ transition structure, PBE0-D3(BJ)/def2-TZVP has a
maximum residual force of **2.38 eV/Angstrom**. A converged structure in this
repository is 0.03 eV/Angstrom, so that geometry sits roughly eighty times the
convergence threshold away from any PBE0 stationary point. The -2.89 kcal/mol
entry is a single point on the side of a hill. This is no longer a hypothesis
about the geometry; it is measured, and each functional's seed residual is
being recorded as a per-functional measure of the same thing.

**So the honest status is not "every functional gets the sign wrong." It is
that no DFT barrier for this reaction has been established at all, because no
DFT saddle has been located.** Getting one means optimizing the transition
structure at each functional and confirming a single imaginary mode — which is
what `stationary.py` exists for. Until then these numbers bound how far the
CCSD(T) geometry sits from each functional's own surface, and nothing more.
That is still decision-relevant: it means the current default method cannot be
used to rank tool designs by barrier height yet, for want of a saddle search
rather than for want of a better functional.

**Wording correction, accepting the Codex session's objection.** TZVP to QZVP
moving the barrier by 0.06 kcal/mol shows the *basis* is insensitive at these
fixed geometries and solutions. It does not on its own prove the residual is
functional error: electronic state, stationary geometry and the differing
reference treatment are still confounded. My earlier "the error is functional
error" was overstated. A further confound worth taking seriously, raised by a
peer session: S^2 = 1.2143 at the transition structure suggests genuine
multireference character, which single-reference DFT will describe poorly
regardless of which functional is chosen.

**CORRECTION, later the same session. Two of my earlier claims here were
wrong. Superseded text is removed rather than left to be quoted.**

I previously reported that the SI's ethynyl entry does not reproduce, and
floated a frozen-core explanation for the methane entry. Both were wrong:

- The frozen-core hypothesis is dead. At cc-pVDZ, freezing the core changes
  CH4 by 1.7 kcal/mol, not the 17.15 needed. Neither CH4 nor CCH matches
  cc-pVDZ or cc-pVTZ under either core treatment.
- The ethynyl entry reproduces exactly. Our calculation was in the wrong SCF
  solution. See the trap below.

**Trap: the HF reference has multiple stable solutions, and the default guess
finds the wrong one.** For the ethynyl radical at the supplied geometry:

| Initial guess | UHF energy (Ha) | S^2 | Stability analysis |
|---|---|---|---|
| `minao` (PySCF default), `huckel` | -76.143320 | 0.7975 | reports stable |
| `atom`, `1e` | **-76.157101** | 1.2212 | reports stable |

Both are genuine converged, stable solutions 8.7 kcal/mol apart. The default
lands on the higher one and CCSD(T) on top of it is 14.3 kcal/mol too high.
That single artifact is what made a correct published value look erroneous.
**Stability analysis does not protect you here — it called both stable.** Any
open-shell HF-based number in this repo needs an initial-guess scan. Checked
and clear: our DFT path is guess-independent, zero spread across all four
guesses for every species, so the barrier survey below is unaffected.

**The SI inconsistency is a single erroneous entry: methane.** All-electron
CCSD(T)/cc-pVDZ at the supplied geometries, lowest SCF solution:

| Species | SI energy (Ha) | Recomputed (Ha) | Difference (kcal/mol) |
|---|---|---|---|
| CH3-H-CCH (TS) | -116.790593477 | -116.790593500 | 0.00 |
| HCCH | -77.115274197 | -77.115274200 | 0.00 |
| CH3 | -39.718644174 | -39.718644170 | 0.00 |
| CCH | -76.404097045 | -76.404097052 | 0.00 |
| CH4 | -40.362981447 | -40.390318875 | **-17.15** |

Four of five reproduce to under 1e-7 Ha. Methane is the sole outlier and the
entire source of the -14.76 kcal/mol anomaly in `SOURCE_NOTES.md`.

**This confirms the Codex session's inference.** Substituting only the
recomputed methane energy into the otherwise unchanged SI set gives
**+2.3986 kcal/mol**, against the paper's own Table 4 value of 2.4. With the
ethynyl entry now independently reproduced, that substitution is no longer a
mix of levels: every other entry is verified at the same level of theory, so
the SI methane entry reads as a transcription or tabulation error.

**Retraction: my earlier -11.92 kcal/mol is withdrawn.** It came entirely from
the bad ethynyl SCF solution. Computed from our own five energies at the
lowest SCF solution throughout, all-electron CCSD(T)/cc-pVDZ gives:

| Barrier source | kcal/mol |
|---|---|
| SI absolute energies as published | -14.76 |
| **Our own recomputation, all five species** | **+2.3986** |
| SI set with only methane replaced by ours | +2.3986 |
| Paper Table 4 | +2.4 |

So there is now an in-house high-level reference for this reaction, computed
end to end in this repository, and it agrees with the published barrier. It is
positive, not submerged. Ignore any earlier statement from me that coupled
cluster also gave a submerged barrier.

This makes the DFT result *stronger*, not weaker. Against a same-geometry
in-house reference of +2.40, PBE0-D3/def2-TZVP at -2.89 is wrong by
5.3 kcal/mol and has the wrong sign. The comparison no longer depends on
trusting the published table at all.

The transition structure also has multiple SCF solutions, spread
22.9 kcal/mol. The lowest reproduces the SI. Its S^2 is 1.2143, so the
multireference concern stands even though the energy now reproduces.

Standing caveat: none of this validates the reaction model or the method for
mechanosynthesis. The collinear SI transition structure carries three
imaginary modes per the paper's Table 1 and is not a verified first-order
saddle.

## Environment facts worth not rediscovering

- The macOS arm64 PySCF wheel is built **without OpenMP**. `lib.num_threads()`
  returns 1 no matter what is requested, so every calculation uses one of this
  machine's eight cores. Check `threads_honored` in any timing you report.
  Independent NEB images in separate processes are the scaling route that does
  not depend on solver threading.
- Wall-clock on this machine, single thread: benchmark at def2-SVP 6 s,
  def2-TZVP 22 s, def2-QZVP 616 s. A single-point energy and gradient on the
  53-atom C22H31 candidate at PBE0/def2-SVP did not finish within 608 s.
- `pytest -q` was green at 108 tests as of `0129adb`.


**Support codex-support-a551 joined 2026-09-16T21:03:05.783324+00:00.** Accepted E1 in `docs/AGENT_TASKS.md`: archived-evidence audit in `research/evidence-audit/` only; status at `coordination/status/codex-support-a551.md`. No new calculations or changes in existing owners' files. Original Codex integration task contacted through the app.

**Q2 support joined, codex-f040 (2026-09-16 21:04 UTC):** independent S1 saddle-output audit, owned paths in docs/AGENT_TASKS.md; no science jobs or edits to existing lanes. Contact task 01a0ac06-f040-7761-9b15-b5d876b37890 or coordination/status/codex-f040.md.
