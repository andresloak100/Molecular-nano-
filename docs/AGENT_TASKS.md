# Shared agent assignments

The user explicitly asked Codex to assign work to the additional agent on 2026-09-16. This is the shared handoff channel; acknowledge a task here before starting and name your owned files. All sessions share one checkout. Do not reset, overwrite or commit another agent's uncommitted files. Stage explicit owned paths. Keep result provenance and scientific status visible.

## V1 — New joining agent: local visual design workbench

**Assignment posted; acknowledgement pending. Own `workbench/` only.** Codex is integrating core Python modules, CLI, root documentation and tests; the other scientific session owns source/reference forensics. Do not modify their files without coordinating here.

Build a useful, working local visual workspace around the actual model and recorded calculations. This is a structure/evidence inspection tool, not a fictional operating nanomachine animation.

Deliver these connected views:

1. **Atomic structure:** an interactive 3D viewer of real XYZ/extended-XYZ coordinates. Show atom indices and element colors; distinguish fixed anchors, substrate and tool; highlight donor carbon, transferred H and acceptor carbon from the design metadata. Rotation, zoom, atom selection and pair distances must work. Preserve Å units and the actual supplied geometry.
2. **Pose study:** show the nine entries in `examples/pose-campaign`, their separation/offset settings and pending status. Allow selection of a pose and its initial/final structures. Read the schema in `docs/CAMPAIGNS.md` and `nanodesign/campaign.py`. Support a user-selected local campaign directory through a documented, bounded mechanism.
3. **Calculation evidence:** inspect the completed 53-atom direct and density-fitting records, including energy, maximum free force, spin diagnostics, anchor loads, method/settings, elapsed time and remaining validation. A single-point energy is not a barrier or a relaxed design. Show the documented timing caveat.
4. **Reference comparison:** show the actual paired `minao` and `atom` runs from `data/validation/paired-ccpvdz*`, making the different electronic solutions and unverified state explicit. Use recorded numbers and links to evidence. Path and vibration views may load existing artifacts when available; otherwise display “not computed.”

Technical boundaries:

- Keep the app, its dependencies, startup instructions, and focused tests inside `workbench/` so it can be developed independently. Prefer a small maintainable implementation. Do not add framework/dependency changes to root project files.
- Bind any local server to loopback only. Do not expose arbitrary filesystem paths; constrain reads to the configured project/campaign root, including symlink/path traversal checks. This assignment is read-only: do not launch quantum jobs, mutate saved evidence or implement hardware control.
- Never populate missing scientific values with fabricated data, success percentages, “best design” labels or implied physical validation. Numerical completion and scientific validation must remain separate.
- Explain bond rendering: if inferred from distance or taken from nominal metadata, label it as such. It does not establish electronic bond order.
- No deployment, cloud spending or external publication is needed. No need to touch the long-running calculation machinery.

Acceptance checks:

- Start locally with a documented command and test the interface in a browser.
- Confirm the 53-atom view really contains 22 C and 31 H, with fixed indices `[7,8,9,35,36,37]` and transfer indices donor 0, H 10, acceptor 26.
- Verify the default pose's donor-to-apex distance is 3.6 Å and changing selected pose changes the displayed coordinates correctly.
- Confirm all nine sample campaign rows show pending, and archived completed calculations display their actual status without promoting the candidate to validated.
- Check empty/missing/failed results, invalid import paths, and no network dependency for scientific data. Report tests and any remaining limitations.
- Post changed file paths and a short handoff here. Commit only `workbench/` once tested, or leave an explicit ready-for-Codex-review note. Coordinate pushes if another session is pushing.

**Acknowledgement / progress:** awaiting the new agent.

## S1 — Existing scientific helper: reaction-specific validation

Retain ownership of `SOURCE_NOTES.md` and your source-forensics scripts/artifacts. Your initial-guess investigation has been integrated as explicit production controls; do not duplicate `highlevel.py` or paired-comparison API changes.

Next useful independent task: characterize whether the supplied eight-atom nominal transition structure can be refined to a DFT stationary point at the stated method, including force residual, displacement-step sensitivity and mode directions. Work in a new `research/reference-saddle/` directory and preserve run artifacts separately. Set a bounded calculation budget and report actual outcomes, including failure to find a first-order saddle. Do not call a fixed-geometry energy difference a DFT barrier. Coordinate any need to change core `stationary.py`/`workflow.py` with Codex.

This task is suggested for claim when your current source writeup is complete; do not launch duplicate searches until you acknowledge it here.

## C1 — Codex integration owner

Own core package modules, root CLI/docs/tests, campaign orchestration and the reviewed validation archives. Finish/push the current tested checkpoint, review V1's read-only interfaces and results, then integrate its launcher without conflicting edits. Maintain the distinction between an executed research model and a validated molecular-machine design.
