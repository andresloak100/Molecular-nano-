# Project organization — lead agent, deputy, and operating rules

Recorded 2026-09-16 ~21:45 UTC, revised ~22:20 UTC, by session
`molecular-nano-6c` (cloud, remote execution environment) on the user's direct
instruction: *"assign the lead agent and give structure to the other agents to
perform better."* This file carries that user directive. It changes reporting
and process only; **no accepted lane is reassigned.**

## Lead agent: Codex root (C1)

Codex root is the lead. The reasons are on the record: it is the only agent
with cross-cutting ownership (core package, CLI, root docs, tests), it reviews
and integrates every other lane's output, it settled the assignment conflicts
and the identity churn, it runs its own internal helpers, and — since this was
first recorded — it independently created and maintains the ownership roster.
Formalizing what is already true beats inventing a new hierarchy mid-flight.

**Lane ownership lives in [`ROSTER.md`](ROSTER.md), which Codex root
maintains. That file is authoritative for who owns what; this file does not
duplicate it.** The first revision of this file carried its own roster table;
it was superseded within the hour by root's more complete one, which is
exactly how this should work.

**Lead responsibilities:** owns the roster and `docs/AGENT_TASKS.md` as the
assignment boards; reviews and integrates lane deliverables into core/shared
files; first stop for conflicts (the user is final arbiter); may query any
lane whose status file goes silent for a full work cycle. Lane owners commit
and push their own lane paths per the roster's write boundaries and root's
commit-handoff protocol; anything touching core or another agent's files goes
through the lead.

## Scientific-validity deputy: original scientific helper

The original scientific/forensics session (owner of `SOURCE_NOTES.md` and
`data/validation/si-energy-reproduction/`) is deputy for scientific validity.
Every major correction so far — the SCF multiple-solution trap, the SI methane
entry, the barrier-wording fixes, the contention measurements — came from or
through this lane. **Sign-off rule:** any number or claim promoted into shared
docs (`README.md`, `docs/ACCURACY.md`, the results section of `CLAUDE.md`)
needs the deputy's or the lead's explicit review noted in a status file.

## Operating rules

Each rule exists because its absence already cost this project something.

1. **Status discipline.** One status file per agent under
   `coordination/status/`, UTC-timestamped at the start and end of every work
   batch. Silent lanes are queryable by the lead.
2. **Identity.** Peer listings are not identity; sessions re-register after
   compaction and can sincerely misstate their own role (three
   misidentifications in one day, plus a helper that wrongly claimed lane
   continuity). The roster and status files are the identity record; resolve
   who someone is from files before addressing or assigning, and ask the user
   when you cannot.
3. **Commit cadence.** The `runs/` gitignore fix landed (`65f2dee`); re-add
   previously ignored evidence and verify with `git status --ignored`. Commit
   your own lane promptly, respecting root's explicit commit-handoff protocol
   in the roster, so nobody sweeps anyone's tree again.
4. **Host contention is measured, not hypothetical.** The Mac reached load
   219 on 8 cores; an identical benchmark ran **10.4x slower** contended, and
   quantum jobs get 10–35% of one core. Per the standing CLAUDE.md note: do
   not launch new quantum campaigns on the shared host while this holds;
   serialize within a lane; record `os.getloadavg()` next to every published
   timing. Audit/review/planning lanes remain cheap and genuinely parallel.
5. **Off-host compute exists — use it.** See below. The marginal
   compute-bound agent is negative *on the Mac*; it is positive in an
   isolated cloud container.
6. **Claims discipline.** Converged is not validated; no "barrier" language
   without a verified first-order saddle; failures are results; per-guess SCF
   scans for open-shell HF. Deputy sign-off per above.
7. **Escalation.** Blocked on another lane → write the request under the
   relevant task board entry and your status file. No response within one
   work batch → escalate to the user. Never edit another lane's files to
   unblock yourself.

## Off-host compute channel (`molecular-nano-6c`)

This session runs in an isolated cloud Linux container (4 cores, 15 GiB RAM)
that shares nothing with the Mac — no load contribution, no memory pressure —
and has the project's first **working OpenMP PySCF** (2.14.0;
`lib.num_threads()` honors 4, verified). Portability is established: the full
suite passes on Linux and the PBE0-D3(BJ)/def2-TZVP benchmark reproduces the
Mac values to <0.01 kcal/mol on different OS/arch/BLAS.

Live lane (user-directed): `research/offhost-compute/` — uncontended,
multi-threaded timings of the exact archived 53-atom single points (direct
and density-fitting), to give A2 and C2 a clean per-gradient cost the
contended Mac cannot produce, and to test how far off-host offloading moves
C2's 19.5/57-day path projections. Root: please add this lane to the roster
or object in your status file.

To delegate a bounded job: leave a spec in `coordination/messages/` addressed
to `molecular-nano-6c` (script or exact CLI invocation, settings, output
directory, budget) and have the user ping the cloud session. Results return
as pushed evidence with full provenance. The container is ephemeral: nothing
is real until pushed, so results are pushed on completion.

## Branch note

This session must push to `claude/molecular-nanomachine-design-iu5ik5`, not
`main`. It merges `origin/main` into that branch to stay current; the lead
should merge the branch back (or the user can authorize a direct push to
`main`) so Mac-side agents see these files without the user relaying.
