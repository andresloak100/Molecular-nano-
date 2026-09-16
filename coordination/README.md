# Agent messages and status

This directory is the shared communication channel for sessions that cannot address each other directly. Read `AGENTS.md`, `CLAUDE.md` and `docs/AGENT_TASKS.md` first. The original scientific helper handles incoming scientific-agent assignments; Codex handles core integration and review. Preserve an already accepted task rather than assigning the same person twice.

For the expanded team, consult the authoritative `coordination/ROSTER.md`. Codex tasks can also send direct messages to root task `01a0abd2-e250-7a62-a6d6-91f406c345fd`. Ownership changes need one root decision; post progress in your own status file, not by simultaneously rewriting shared task boards.

Current lanes:

| Agent | Lane | Owned implementation |
|---|---|---|
| Codex root | C1, integration and review | Core package, CLI, root docs/tests |
| Codex `stationary_check` + its viewer helper | V1, local visual evidence workbench | `workbench/` |
| `andresarriaga-8a` | S1, DFT saddle search | `research/reference-saddle/` |
| `andresarriaga-f2` | A1, intrinsic site energetics and selectivity limits | `research/site-selectivity/` |
| `andresarriaga-a8` | A2, candidate cost and reduced-model feasibility | `research/candidate-feasibility/` |
| Original scientific helper | Source forensics and scientific intake | `SOURCE_NOTES.md`, `data/validation/si-energy-reproduction/` |

Assignments A1/A2 were posted by the scientific intake owner. Their presence in this table does not imply they have acknowledged or started. Each agent should create and maintain `coordination/status/<agent-id>.md` with its acknowledgement, owned paths, current work, UTC update time, running process/output directories, latest tests, and any requests. Use one file per agent to avoid concurrent edits to a single board. Brief result handoffs should include exact method/settings, artifact paths, limitations and whether the work is ready for review.

Read addressed notes under `coordination/messages/` before reporting. Reply in your own status file. Do not copy another agent's calculations or terminate its processes. Reuse completed evidence. Declare a bounded job count and optimization limits before a new expensive campaign; memory settings are not hard process limits.

Measured host facts on 2026-09-16: eight logical CPUs, 24 GiB RAM. Tested PySCF OpenMP calls expose one effective thread; this is not a universal statement about all linked libraries. Multiple processes share these finite resources, so report concurrent workload caveats with timings.
