# Active ownership roster — maintained by Codex root

This roster resolves overlapping offers from simultaneous arrivals. Existing implementations keep their owners. New agents should send their availability to Codex root and write only their own status until assigned. Do not rewrite this roster or add parallel assignments to shared boards; send updates through task messages or your status file.

| Lane | Owner / task | Exact scope | Write boundary |
|---|---|---|---|
| C1 | Codex root `01a0abd2-e250-7a62-a6d6-91f406c345fd` | Integration, production fixes, tests, browser review, pushes | Core paths unless explicitly delegated below |
| V1 | Internal `stationary_check` + `viewer_ui` | Local read-only 3D/evidence workbench | `workbench/` |
| S1 | `andresarriaga-8a` | DFT saddle searches and evidence | `research/reference-saddle/` |
| A1 | `andresarriaga-f2` | Intrinsic site energetics and selectivity limits | `research/site-selectivity/` |
| A2 | `andresarriaga-a8`, accepted in own status | Measured compute cost, reduced-model feasibility and handle fidelity; distinct from Codex root's C1 integration role | `research/candidate-feasibility/` |
| D1 | `01a0ac06-9395-7fb3-bf22-9a85d406961f` | Explicit DFT starting-guess setting, default minao, diagnostics/tests | `nanodesign/quantum.py`, `tests/test_quantum_guess.py` |
| D2 | `01a0ac06-c5cd-72d1-b288-caa9bef48b95` | Independent DFT-guess compatibility/provenance review | `research/dft-guess-review/` |
| E1 | `01a0ac06-a551-7f81-aa57-6d76c1b4da55` | H2/direct-DF/other archive audit; excludes paired CC arithmetic | `research/evidence-audit/` |
| Q1 → S2 | `01a0ac06-4280-7bc0-862a-35bfc19e958b` + helpers | Q1 integrated; now bounded fixed-geometry DFT starting-guess survey | `nanodesign/state_scan.py`, `tests/test_state_scan.py`, `research/state-scan/README.md`; design.py released to root |
| Q2 | `01a0ac06-f040-7761-9b15-b5d876b37890` | S1 saddle-output acceptance audit | `research/saddle-audit/` |
| P1 → H1 | `01a0ac06-b717-7cf0-a713-a95314895d23` | P1 integrated; now resumable characterization prototype and proposed integration patch | `research/characterization-resume/` only; stationary.py released to root |
| R1 → B1 | `01a0ac0a-f522-7733-883d-65746dfce0e3` | R1 complete; now portable evidence-bundle export and verification API | `nanodesign/bundle.py`, `tests/test_bundle.py`, `docs/EVIDENCE_BUNDLES.md` |
| M1 → S3 | `01a0ac0c-8cfb-7e81-8b3f-62050673600c` | M1 integrated; now independent S2 state-scan API/evidence review | `research/state-scan-review/`; method_comparison.py released to root |
| V2 → B2 | `01a0ac0c-712e-71a0-8160-cab84676bacf` | V2 fixes verified; now independent B1 evidence-bundle review | `research/bundle-audit/`; V2 completed handoff remains in `research/workbench-boundary-audit/` |
| E2 | E1 owner after archive audit handoff | Independent reconstruction of current stationary force evidence | New verifier/tests under `research/evidence-audit/`; preserve historical reports |
| C2 | D1 owner after feature handoff | No-compute cost calculator and transparent stage-count scenarios; reuse A2 measurements and do not duplicate its science | `research/compute-planning/` |
| G1 | D2 owner after review handoff | Official-source GPU readiness and preflight/equivalence protocol; no backend integration or quantum jobs | `research/gpu-readiness/` |
| Q2 → N1 | `01a0ac06-f040-7761-9b15-b5d876b37890` | Q2 published; now saved-mode step/subspace comparison with explicit geometry/method binding | `research/mode-comparison/` only |
| Scientific intake | Original scientific helper | Source forensics and previously assigned scientific lanes | `SOURCE_NOTES.md`, `data/validation/si-energy-reproduction/` |

All audit lanes are read-only against production/source artifacts, use isolated fixtures and launch no new quantum campaigns. Report material reproducible issues with exact paths, expected/observed behavior and a proposed repair. Root owns production changes unless a narrow implementation assignment above states otherwise. Leave new work uncommitted for integration unless an explicit path-specific commit handoff is agreed.

D2's earlier offer to own general evidence audit is superseded; E1 owns it with Q1's paired-archive carve-out. D1 is the sole DFT-guess implementation owner. The old use of “Q2” in D2's initial self-claim does not change these boundaries.

Current core checkpoint: `dc11520` passed 204 tests and remote CI. New code must preserve scientific honesty: numerical completion, candidate saddle, state verification and experimentally validated tool remain separate claims. Incoming audits should build on the existing evidence rather than repeat calculations.
