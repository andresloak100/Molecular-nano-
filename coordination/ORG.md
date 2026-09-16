# Project organization — lead agent and reporting structure

Recorded 2026-09-16 ~21:45 UTC by session `molecular-nano-6c` (cloud, remote
execution environment) on the user's direct instruction: *"assign the lead
agent and give structure to the other agents to perform better."* This file
carries that user directive. It changes reporting and process only; **no
accepted lane is reassigned**, per the standing rule in
`coordination/README.md`.

## Lead agent: Codex root (C1)

Codex root is the lead. The reasons are already on the record: it is the only
agent with cross-cutting ownership (core package, CLI, root docs, tests), it
already reviews and integrates every other lane's output, it has settled the
two assignment conflicts to date (V1/S1, and the A2 self-assignment), and it
runs its own internal helpers. Formalizing what is already true beats
inventing a new hierarchy mid-flight.

**Lead responsibilities:**

1. Owns `docs/AGENT_TASKS.md` as the single assignment board. New lanes and
   task changes go through the lead (scientific intake may draft assignments,
   as it did for A1/A2, with the lead's acknowledgement).
2. Reviews and integrates lane deliverables into core/shared files. Lane
   owners still commit and push their **own lane paths** directly; anything
   touching the core package, root docs, `.gitignore`, or another agent's
   files goes through the lead.
3. First stop for conflicts between agents. The user is the final arbiter.
4. May query any lane whose status file has gone silent for a full work
   cycle, and reassign a lane only if its owner confirms it has stopped or
   the user directs it.

## Scientific-validity deputy: original scientific helper

The original scientific/forensics session (owner of `SOURCE_NOTES.md` and
`data/validation/si-energy-reproduction/`) is deputy for scientific validity.
Every major correction so far — the SCF multiple-solution trap, the SI methane
entry, the barrier-wording corrections — came from this lane or was verified
by it. **Sign-off rule:** any number or claim promoted into shared docs
(`README.md`, `docs/ACCURACY.md`, `CLAUDE.md` results section) needs either
the deputy's or the lead's explicit review noted in a status file.

This also settles the current collision on
`coordination/status/scientific-helper.md`: the **still-running original
session keeps the lane and that file**; the later relay session takes the
liaison role below under its own status file, as already proposed inside that
file. Endorsed here so both can stop negotiating and work.

## Roster and reporting lines

| Agent | Role / lane | Reports to | Owned paths |
|---|---|---|---|
| Codex root | **Lead**; C1 integration | User | `nanodesign/`, `cli`, root docs/tests, `.gitignore` |
| Codex `stationary_check` (+ viewer helper) | V1 visual workbench | Codex root | `workbench/` |
| `andresarriaga-8a` | S1 DFT saddle search | Lead (science: deputy) | `research/reference-saddle/` |
| `andresarriaga-f2` | A1 site selectivity | Lead (science: deputy) | `research/site-selectivity/` |
| Codex (as `andresarriaga-a8`) | A2 candidate feasibility | — (lead's own lane) | `research/candidate-feasibility/` |
| Original scientific helper | **Deputy**; forensics | User / lead | `SOURCE_NOTES.md`, `data/validation/si-energy-reproduction/` |
| `andresarriaga-5b` | Liaison / relay | Lead | `coordination/status/andresarriaga-5b.md` |
| `codex-support-q1` | Q1 integration audit | Lead | `research/integration-audit/` |
| `codex-support-a551` | E1 evidence audit | Lead | `research/evidence-audit/` |
| `codex-f040` | Q2 S1-output audit | Lead | per `docs/AGENT_TASKS.md` |
| `molecular-nano-6c` (this session, cloud) | Remote compute, portability, org relay | User / lead | `coordination/ORG.md`, own status/messages |

## Operating rules

Each rule exists because its absence already cost this project something.

1. **Status discipline.** One status file per agent under
   `coordination/status/`, updated with a UTC timestamp at the start and end
   of every work batch. Silent lanes are queryable by the lead.
2. **Identity.** Peer listings are not identity — three misidentifications
   happened in one day. The status file is the identity record; resolve who
   someone is from files before addressing or assigning to them.
3. **Commit cadence and the `runs/` fix.** The `.gitignore` fix **has
   landed** (`65f2dee`, `runs/` → `/runs/`), so lane evidence under
   `research/*/runs/` now commits. Every lane owner must `git add` their
   previously ignored run artifacts and verify capture with
   `git status --ignored` before citing any run as saved. Commit your own
   lane promptly so no one sweeps your uncommitted tree again.
4. **Job registry.** Before launching any quantum job, record its directory,
   PID, and thread count in your status file, and check the other status
   files first. Leave at least 2 of the Mac's 8 cores free. Declare a bounded
   job count before any campaign (standing rule, now lead-enforced).
5. **Claims discipline.** Converged is not validated. No "barrier" language
   without a verified first-order saddle; failures are results; per-guess SCF
   scans for any open-shell HF number. Deputy sign-off per above.
6. **Escalation.** Blocked on another lane → write the request under the
   relevant task in `docs/AGENT_TASKS.md` and your status file. No response
   within one work batch → escalate to the user (directly or via liaison).
   Never edit another lane's files to unblock yourself.
7. **Handoffs.** Exact method/settings, artifact paths, limitations, and a
   ready-for-review flag, per `coordination/README.md`.

## Remote compute offer (new capacity)

`molecular-nano-6c` runs in an isolated cloud Linux container that — unlike
the shared Mac — has a **working OpenMP PySCF build**: PySCF 2.14.0,
`lib.num_threads()` honors requests (4 threads, 4 cores, 15 GiB RAM).
Verified on 2026-09-16: all 197 non-quantum + 7 quantum tests pass on Linux,
and the PBE0-D3(BJ)/def2-TZVP fixed-geometry benchmark reproduces the Mac
values to <0.01 kcal/mol (barrier −2.89, reaction energy −28.05 kcal/mol;
33 s wall on 4 threads). Different OS, arch, and BLAS — so this doubles as an
independent cross-platform reproduction of the benchmark table.

To delegate a heavy or queue-blocking job (QZVP points, 53-atom evaluations,
guess scans): leave a spec in `coordination/messages/` addressed to
`molecular-nano-6c` — script or exact CLI invocation, settings JSON, output
directory, and a bounded budget — and have the user ping the cloud session.
Results come back as pushed evidence with full provenance. Caveat: the
container is ephemeral; nothing is real until pushed, so delegated results
are pushed immediately on completion.

## Where this session's work lives

This session is required to push to branch
`claude/molecular-nanomachine-design-iu5ik5`, not `main`. The lead should
fetch and merge that branch (or the user can tell this session to push to
`main` explicitly) so the Mac-side agents see this file. Until merged, the
user is the relay.

## Suggested lead queue (advisory, lead may reorder)

1. S1 saddle results — the first genuine DFT barrier, gating everything.
2. A1 stages 1–2 — the selectivity question nobody was on.
3. A2 verdict with numbers.
4. V1 review and integration.
5. The `data/validation/paired-ccpvdz/method_comparison.json` minao-arm
   labeling fix raised by the deputy (stale −39.11 / +12.31 quotable without
   caveat) — small, and it removes a live footgun for V1.
