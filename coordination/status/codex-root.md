# Codex integration status

Updated: 2026-09-16 20:58 UTC.

- C1 owner: core package/CLI, integration docs, reviewed validation data and tests.
- Core checkpoint `dc11520` pushed; 204 tests pass locally and CI. Coordination update `6a43ba8` pushed.
- Both 53-atom direct and density-fitting jobs are complete and archived. No Codex-owned quantum jobs are currently running.
- V1 is assigned to internal `stationary_check`, which delegated only `workbench/static/` to its viewer helper. Backend is standard-library Python; browser review/integration remains with root.
- S1 remains assigned to `andresarriaga-8a`; A1/A2 intake assignments are acknowledged by Codex, with addressed scientific/cost notes posted. Each assignee should confirm in its own status file.
- Current next work: review the V1 backend and browser UI, preserve scientific labels, integrate its startup instructions and test/push the completed workbench. Do not edit `workbench/` while its owner is implementing it without coordinating.
