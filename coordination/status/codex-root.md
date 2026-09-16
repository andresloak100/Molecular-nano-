# Codex integration status

Updated: 2026-09-16 21:47 UTC.

C1/root task: `01a0abd2-e250-7a62-a6d6-91f406c345fd`.

Current user priority is scientific accuracy for the supported H-abstraction
operation. See `coordination/ROSTER.md` and the root scientific-priorities note.
S1/A1/A2/source-forensics retain existing chemical calculations; root has launched
no new quantum jobs. None of the code/test handoffs validates a molecular tool.

Reviewed integration ready for publication:

- S2 starting-guess survey and B1 evidence transport are in the CLI. Explicit
  charge/spin and preserved subset/failure reporting prevent implicit state
  assumptions. Bundle integrity is separate from scientific validity.
- S3/S4 CLI integration found ASE can clamp invalid negative frame indices and
  leak an exception for missing positive frames. Root fixed explicit frame
  selection from captured bytes; tests now cover binary/text/compressed/database
  frames. 24 survey CLI checks pass after repair; 19 bundle CLI checks pass.
- Focused production integration initially passed 204 tests. After frame and
  snapshot changes, 127 reader/characterization/survey/campaign/audit checks pass.
  These are software checks with synthetic force models, not chemical evidence.
- Root added a separate exact input_snapshot_sha256 for the rendered
  characterization frame, preserving the original source structure_sha256.
  N1 and H1 now pass actual mocked producer/consumer integration; N1 has 54
  checks, H1 has 23 total cases and H2 retained its historical 36-check receipt
  plus a fresh two-case metadata/reuse receipt. Prototype remains outside the CLI.
- E2 force/Hessian reconstruction and G2 strict offline GPU protocol are frozen
  handoffs. No GPU adapter, actual GPU equivalence or acceleration is established.
- P2 completed the operation-specific validation protocol, with P3 independent
  source/code review and a hash-bound receipt. Scientific gates remain open.
- Final new-command integration: 186 tests passed, including survey/bundle APIs,
  both new command families and the older paired-comparison CLI. No quantum jobs.
- Prior published foundation: 358 core tests and remote CI passed at 9e7ac9f;
  workbench/planner b550e54 and offline GPU-readiness b2657c9 passed remote CI.
  Workbench remains available at http://127.0.0.1:8765 (server session 45775).

Root will publish only explicit reviewed paths using an isolated Git index and
compare-and-swap commit update. Other agents' live evidence, staging and source
notes remain theirs. Completed supporting lanes should freeze until a concrete
next scientific need is assigned; do not add broad audits or duplicate jobs.
