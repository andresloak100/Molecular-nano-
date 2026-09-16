# Scientific helper (forensics / intake) status

Updated: 2026-09-16 21:20 UTC. Same role as the "original scientific helper" /
"joining session ~16:40" in CLAUDE.md; now operating from a fresh session at the
user's request ("help these agents with whatever they need"). Peers can reply to
my inter-session messages directly; Codex should reply here or in CLAUDE.md.

## Lane

Unchanged: source forensics (`SOURCE_NOTES.md`,
`data/validation/si-energy-reproduction/`) and scientific intake. Complete and
pushed; no new claims. I am not starting quantum jobs.

## Actions this session (2026-09-16 ~21:15 UTC)

- Read CLAUDE.md, AGENTS.md, docs/AGENT_TASKS.md, coordination/* in full before
  touching anything. No files in other lanes edited.
- Sent addressed inter-session messages to all three assignees, since none had
  posted a status file yet:
  - `andresarriaga-8a` (S1): pointed at coordination/messages/andresarriaga-8a.md
    and the four S1 review corrections in docs/AGENT_TASKS.md; confirmed the
    53-atom jobs are complete and must not be restarted; confirmed the +2.3986
    kcal/mol in-house CC reference now reproduced through the production API.
  - `andresarriaga-f2` (A1): pointed at its brief and Codex's integration note;
    restated the three non-negotiable constraints (1-thread PySCF, open-shell
    SCF guess scans, converged != validated) and the adiabatic/vertical and
    kBT-wording requirements; asked for acknowledgement in its status file.
  - `andresarriaga-a8` (A2): pointed at its brief, the measured direct/DF
    timings (2044.84 s / 701.53 s, 463 basis functions, overlapping-workload
    caveat), and the 2400-evaluation path cost ceiling; asked for
    acknowledgement and a bounded first experiment in its status file.
- Verified running processes before messaging: three `saddle_search.py`
  processes (pbe0 / b3lyp / m06-2x) plus `guess_scan_fixed_geometry.py`, all
  session 8a's. Untouched. Load advisory (4 of 8 cores busy, single-thread
  each) included in every message.

## Requests to Codex

None blocking. I will relay any assignee request that needs a `nanodesign/`
change into docs/AGENT_TASKS.md under the relevant task, per the documented
protocol. I am deliberately not committing `coordination/` files; the directory
is yours and uncommitted — fold this file into your next coordination batch
commit, or tell me here to commit only my own status file.

---

## Note: two sessions are writing to this file

Appended 2026-09-16 ~21:25 UTC by the original scientific-helper session, which
is **still running** (it assigned A1/A2, ran the forensics, and holds PID 57526).
The section above was written by a different session that understood itself to
hold this role. I am merging rather than overwriting, because clobbering another
agent's status is exactly what the protocol forbids, and the relay work logged
above is real and useful.

Proposed split so we stop colliding on one path: the other session keeps the
relay/liaison log above, and should move it to its own status file. Scientific
results for this lane go below. If the user intended a handover, I am happy to
stand down from the lane, but the handover did not reach me and I have live
results and a running job, so I am continuing until told otherwise.

Correction to the section above, since it will otherwise mislead: **Codex is
not file-only.** It is `andresarriaga-a8` [e3c1ee] and is reachable by direct
message. It dropped off the peer listing and re-registered after a context
compaction. I originally wrote the file-only claim into `CLAUDE.md` and it was
wrong; `CLAUDE.md` now carries the correction.

## Scientific results in this lane

Pushed: `121c30f`, `0129adb`, `e12a476`, `203dc84`.

**The published SI inconsistency is one bad entry: methane.** All-electron
CCSD(T)/cc-pVDZ on the supplied geometries reproduces four of five labeled
energies to under 1e-7 Hartree; methane alone is off by -17.15 kcal/mol, and
substituting the recomputed value gives +2.3986 kcal/mol against the paper's
Table 4 value of 2.4. Independently matched by Codex's paired run at +2.3986156.

**SCF trap.** The ethynyl radical has at least two converged, *stable* UHF
solutions 8.7 kcal/mol apart; PySCF's default `minao` guess finds the higher
one and `stability()` calls both stable. CCSD(T) on it is 14.3 kcal/mol high.
The supplied transition structure admits solutions spread 22.9 kcal/mol. Scan
guesses for any open-shell HF work. DFT is guess-independent here (checked).

**Retraction, and a number still live in a committed file.** My earlier
methane reaction energy of -39.11 kcal/mol was computed on the wrong ethynyl
solution. Corrected value is **-24.79**. The stale -39.11, and a
`dft_minus_ccsd_t` of +12.31, remain in
`data/validation/paired-ccpvdz/method_comparison.json`; the `-atom` arm has the
correct -24.79 and -2.01. That is the intended output of a paired comparison,
but the minao arm records `scf_initial_guess: null` and omits the
electronic-state caveat that the atom arm carries, so the wrong-solution number
is the one a reader can quote without knowing what it is. Raised with Codex;
relevant to V1, which is briefed to render both arms.

**Isobutane extension, in progress (PID 57526).** Same reproduction applied to
`C2H + iso-C4H10 -> C2H2 + t-C4H9`, which matters more than methane here
because its abstraction site is *tertiary* — the closest small-molecule proxy
for the adamantane bridgehead C-H the 53-atom candidate targets.

All four cheap species reproduce to under 1e-7 Hartree: ethynyl (guess `atom`,
spread 8.76), acetylene, tert-butyl (101 bf, 153 s), isobutane (106 bf, 54 s).
Only the 139-basis-function transition structure is outstanding. This
reinforces methane as the single corrupted entry.

Recomputed isobutane reaction energy **-32.23 kcal/mol** against methane's
corrected -24.79, so tertiary abstraction is 7.44 kcal/mol more exothermic than
primary. Literature C-H bond strength difference is roughly 8-9, so the chain
is behaving. Directly relevant to A1: the adamantane bridgehead is tertiary.

Concurrent-load caveat: several sessions compute on the same eight cores, so
none of these timings are controlled benchmarks.
