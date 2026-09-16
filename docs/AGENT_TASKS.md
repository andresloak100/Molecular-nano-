# Shared agent assignments

The user explicitly asked Codex to assign work to the additional agent on 2026-09-16. This is the shared handoff channel; acknowledge a task here before starting and name your owned files. All sessions share one checkout. Do not reset, overwrite or commit another agent's uncommitted files. Stage explicit owned paths. Keep result provenance and scientific status visible.

## V1 — Codex internal helper: local visual design workbench

**Assigned to Codex internal helper `stationary_check`. Own `workbench/` only.** The newly added external session already accepted the saddle-search lane, so it does not own this task. Codex is integrating core Python modules, CLI, root documentation and tests; the other scientific session owns source/reference forensics. Do not modify their files without coordinating here.

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

**Acknowledgement / progress:** dispatched by Codex to its internal helper.

## S1 — Third session `andresarriaga-8a`: reaction-specific validation

**Already accepted through the scientific helper's direct handoff.** Keep that task and its original brief. Do not take V1 or duplicate `highlevel.py`/paired-comparison API changes. The original scientific helper retains `SOURCE_NOTES.md` and its source-forensics artifacts.

Next useful independent task: characterize whether the supplied eight-atom nominal transition structure can be refined to a DFT stationary point at the stated method, including force residual, displacement-step sensitivity and mode directions. Work in a new `research/reference-saddle/` directory and preserve run artifacts separately. Set a bounded calculation budget and report actual outcomes, including failure to find a first-order saddle. Do not call a fixed-geometry energy difference a DFT barrier. Coordinate any need to change core `stationary.py`/`workflow.py` with Codex.

This clarifies and records the already accepted task; it is not a second search request. Preserve the existing run if you have started. Post file ownership, actual calculations and remaining limitations here or in the linked live coordination note.

**Codex read-only integration review of the initial `saddle_search.py`:**

- Please rename `barrier_on_verified_saddle` and `saddle_verified_single_imaginary_mode` to candidate/evidence language. One negative mode plus an approximate transfer direction does not verify connectivity or electronic state. Include `stationary_within_force_tolerance` in any candidate acceptance check; preserve `transition_state_verified: false` while IRC/equivalent connectivity is absent.
- The mode vectors returned by `stationary.py` are Cartesian displacements normalized by the mass metric. Your sum of squared Cartesian norms is a Cartesian displacement share, not a mass-weighted amplitude share. Label it accurately or multiply each squared displacement by its mass. The 25% cutoff is a declared heuristic, not a proof of reaction identity.
- If only two of four SCF guesses converge, agreement between those two does not establish a completed four-guess check. Expose missing guesses and label agreement as applying only to the converged subset.
- High HF S² can flag a poor or broken-symmetry reference, but does not by itself diagnose multireference character or establish how every DFT functional must behave. Preserve spin diagnostics without promoting that inference to a finding.

These are review requests for the file owner; Codex has not edited your research implementation or interrupted your calculations.

## C1 — Codex integration owner

Own core package modules, root CLI/docs/tests, campaign orchestration and the reviewed validation archives. Finish/push the current tested checkpoint, review V1's read-only interfaces and results, then integrate its launcher without conflicting edits. Maintain the distinction between an executed research model and a validated molecular-machine design.

---

## Intake for sessions added 2026-09-16 ~17:30

The user asked the scientific/forensics session to take agent intake for the two sessions added at this time. Tasks A1 and A2 below are written against work that is currently unclaimed and that does not touch C1, V1 or S1 files. Acknowledge under your task before starting, and name the files you own.

Three constraints apply to both tasks and are not negotiable, because each has already cost this project a wrong answer or a wasted run.

1. **The PySCF build on this machine has no OpenMP.** `lib.num_threads()` returns 1 regardless of the `threads` setting, so every calculation uses one of eight cores. Check `threads_honored` in the diagnostics before quoting any timing. Parallelism means separate OS processes, nothing else.
2. **Open-shell Hartree-Fock has multiple converged, *stable* SCF solutions here.** For the ethynyl radical the default `minao` guess sits 8.7 kcal/mol above the solution found from `atom` or `1e`, and PySCF's stability analysis calls both stable. This produced a wrong published conclusion earlier today. Scan initial guesses for any open-shell HF work. DFT was checked and is guess-independent at the geometries tested, but verify rather than assume at new geometries.
3. **A converged number is not a validated one.** Do not report a fixed-geometry energy difference as a barrier, do not rank anything by an unvalidated quantity, and record failures as results. Every run directory must be new so nothing silently overwrites earlier evidence.

### A1 — session `andresarriaga-f2`: does the tool hit the right hydrogen?

Own `research/site-selectivity/`. Do not edit core package modules.

**How to reach the core-module owner.** Codex owns `nanodesign/` and is *not* reachable by inter-session messaging; it does not appear in a peer listing and has no socket. It communicates only by writing in this repository: this file, `CLAUDE.md`, and `coordination/messages/`. So if you need a change inside `nanodesign/`, write the request in this file under your task and acknowledge in `coordination/status/<your-agent-id>.md` per `coordination/README.md`. If you get no response and are genuinely blocked, escalate to the user rather than waiting or editing core modules yourself. Do not assume a session in a peer listing owns anything here; most of the listed peers are working on unrelated projects, and guessing has already produced one misrouted round of coordination.

This is the question the whole repository is built to answer and nobody is on it. Atomically precise assembly means abstracting one *specific* hydrogen. Adamantane has 4 bridgehead (tertiary) C-H bonds and 12 methylene (secondary) C-H bonds. The 53-atom candidate targets a bridgehead H, carbon index 0 and hydrogen index 10 in `nanodesign/candidates.py`. If the methylene sites are energetically competitive, then positional control alone does not deliver selectivity, and that is a first-order problem for the whole design premise.

Work in increasing cost, reporting each stage before moving on:

1. **Intrinsic site preference.** Compute the C-H bond dissociation energy at a bridgehead site and at a methylene site: E(adamantyl radical) + E(H atom) - E(adamantane), for both radicals, at a consistent DFT level and basis. Relax each species. Report the difference between sites. Adamantane is 26 atoms and tractable; the isolated H atom is a doublet with `spin=1`.
2. **Reaction preference with the actual abstracting species.** Repeat as ethynyl abstraction reaction energies: C2H + adamantane -> C2H2 + adamantyl, for both sites. This removes the free H atom, which is a poor proxy.
3. **Only if 1 and 2 are done and time remains,** attempt the two transition structures. Expect this to be expensive and say so rather than half-finishing it.

Report the site preference in kcal/mol with the method stated, and state plainly whether it is large enough to matter at room temperature. A thermodynamic preference is not a kinetic one; label which you computed. Note that our DFT has not reproduced a correct barrier for the calibration reaction, so treat DFT site *energetics* as screening, and say so.

### A2 — session `andresarriaga-a8`: is the 53-atom candidate reaction path actually computable?

Own `research/candidate-feasibility/`. Read-mostly: do not restart or kill the existing 53-atom jobs, and do not edit core modules.

The repository's headline candidate is the 53-atom adamantane-plus-ethynyl-tool system, and no reaction path has ever been run for it. Before anyone tries, somebody needs to establish whether it is computable on this hardware and, if not, what the smallest faithful substitute is. Deliver a go/no-go with numbers, not an opinion.

1. **Measure, do not estimate.** Time a single energy-plus-gradient evaluation for the 53-atom system at PBE0/def2-SVP, with and without density fitting. Completed records exist under `data/validation/`; reuse them rather than recomputing where they answer the question. Report basis-function count and seconds per gradient with the effective thread count attached.
2. **Cost out the path.** A CI-NEB needs roughly (images - 2) gradients per optimizer step, plus two endpoint relaxations, each of which needs its own optimizer steps. Use the measured per-gradient cost and a defensible step count to give a wall-clock range for the default 7-image, 200-step settings. State your assumptions.
3. **Give the verdict.** If it is infeasible, say so with the number, then find the smallest model that preserves the chemistry being tested and cost that instead. The obvious reduction is the tool handle: the second adamantane cage exists to make the tip a realistic mounted species, and replacing it with a much smaller group tests whether it changes the local chemistry. Quantify what the reduction costs in fidelity rather than asserting it is fine.
4. **Check the cheaper knobs honestly.** Density fitting, a smaller basis for screening, and a looser force convergence all trade accuracy for speed. Report what each buys and what it costs on a system where we can afford to check both ways.

The useful deliverable is a protocol somebody can actually run this week, with its accuracy compromises written down. "Needs a GPU cluster" is an acceptable conclusion if the numbers support it, but it must come with the numbers.
