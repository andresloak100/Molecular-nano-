# Codex integration status

Updated: 2026-09-16 21:20 UTC.

- C1 owner: core package/CLI, integration docs, reviewed validation data and tests.
- Core checkpoint `9e7ac9f` pushed; GitHub CI passed 358 core tests. Includes explicit DFT guesses, paired provenance checks, exact-byte geometry loading, campaign history integrity, setup failure recording and reconstructible Hessian force evidence.
- Both 53-atom direct and density-fitting jobs are complete and archived. No Codex-owned quantum jobs are currently running.
- V1 is complete and handed off by `stationary_check`/`viewer_ui`: 25 backend tests, 6 independent import-boundary regressions and JS syntax pass. Root browser review confirmed exact structure measurements, endpoint changes, nine pending poses, actual calculation values, historical annotations and uncomputed path labeling; no browser errors. Server: http://127.0.0.1:8765, session 45775.
- Combined published research audits pass 113 tests plus 10 subtests; Q2 saddle audit is independently published. Core/archived evidence remains scientifically unvalidated.
- S1/A1/A2 scientific work is separate and live. A2's accepted status resolves the identity confusion: `andresarriaga-a8` owns candidate feasibility, not C1 integration. Original scientific helper retains the isobutane run and forensics paths.
- Next-milestone assignments have been sent and accepted: S2 fixed-geometry guess surveys, S3 independent review, H1 resumable characterization prototype, E2 stationary evidence verifier, B1 portable evidence bundles and B2 independent review. Exact paths are in the roster. No new quantum jobs were authorized for these software tasks.
- Current integration: publish the finished workbench, C2 cost planner, docs and CI. G1 GPU readiness is a separate in-progress handoff. New source files from next-milestone agents are not part of this finished checkpoint until reviewed and tested.
