# Active ownership roster — maintained by Codex root

Updated 2026-09-16. Root task: `01a0abd2-e250-7a62-a6d6-91f406c345fd`.
Read [the current scientific priorities](messages/from-root-scientific-priorities.md).
Secondary scientific summaries also need the concrete
[A1/A2/digest corrections](messages/from-root-to-science-summary-owners.md).
The user prioritizes accurate predictions for one supported-tip H-abstraction
operation. Software verification, chemical calibration and physical operation
remain separate claims. No lane has validated a molecular assembler.

Existing implementations keep their owners until an explicit handoff. New
arrivals should send availability to root and write only their own status until
assigned. Do not rewrite this roster or add parallel assignments to historical
boards. Do not duplicate quantum jobs or stop another owner's processes.

## Scientific and integration responsibilities

| Lane | Owner / task | Exact scope and write boundary |
|---|---|---|
| C1 | Codex root | Production integration, CLI, root docs/tests, review and publication; delegated paths below remain exclusive until handoff |
| S1 | `andresarriaga-8a` | Existing reference DFT saddle searches, solution follow-up and corrected interpretation; `research/reference-saddle/` |
| A1 | `andresarriaga-f2` | Existing intrinsic site energetics and geometric selectivity limits; `research/site-selectivity/` |
| A2 | `andresarriaga-a8` | Existing reduced-handle electronic fidelity and measured compute cost; `research/candidate-feasibility/`; distinct from root/C1 |
| Scientific intake | Original scientific helper | Source forensics and existing isobutane work; `SOURCE_NOTES.md`, `data/validation/si-energy-reproduction/` |
| P2 | `01a0ac0a-f522-7733-883d-65746dfce0e3` | Completed and released: source-backed operation validation protocol; `docs/VALIDATION_PROTOCOL.md` |
| P3 | `01a0ac0c-712e-71a0-8160-cab84676bacf` | Completed and released: independent P2 scientific-contract review; `research/validation-protocol-review/` |
| S3 CLI | `01a0ac0c-8cfb-7e81-8b3f-62050673600c` | Completed and released: 24 integration tests for root's three survey commands; `tests/test_cli_state_scan.py` |

## Electronic and derivative evidence milestone — published

These assignments supersede renewed availability offers. No duplicate lanes.
Root owns optional backend capture wiring in `nanodesign/quantum.py` and CLI/API
integration. Existing S1/A1/A2 calculations continue; no new SCF/gradient
campaigns or additional agent sessions are authorized for these software lanes.
E3/E4/X1/F1/F2/L1/L2/G3 and capture-lifecycle handoffs are complete, frozen and
released to root. S5 survey integration review is also complete: 21 cases pass after the
failed-gradient SCF-context binding repair; its test file is released to root.
Completed owners should wait for a concrete finding, not open another lane.
Published commit `c66be23` passed remote CI `35156590649`, including 720 core
tests and every new research/review suite. Four new tiny actual solver capture
cases ran on the remote runner; the chemical validity gates remain open.

| Lane | Owner / task | Exact write scope |
|---|---|---|
| E3 electronic snapshots/comparison | f522 | `nanodesign/electronic_state.py`, `tests/test_electronic_state.py`, `docs/ELECTRONIC_EVIDENCE.md` |
| E4 independent snapshot/math review | 712e | `research/electronic-state-review/` |
| X1 explicit cross-AO-overlap bridge | f040 | `research/ao-overlap-bridge/`; tiny <=4 atom/<=10 AO integral-only verification allowed, no SCF/gradient jobs |
| Capture lifecycle review | 8cfb | `research/electronic-capture-integration-review/`; read-only recommendations then synthetic hook tests |
| F1 saved energy/force consistency | 4280 | `research/force-energy-consistency/`; prefer actual H1 checkpoint adapter |
| F2 independent derivative review | 9395 | `research/force-energy-review/` |
| L1 local harmonic residual correction | b717 | `research/local-relaxation-diagnostic/`; no nonlinear energy-error bound or actual optimization |
| L2 independent residual/curvature review | a551 | `research/local-relaxation-review/` |
| G3 actual CPU numerical-field mapping | c5cd | `research/gpu-readiness/CAPTURE_FIELDS.md` only; no GPU adapter |

All snapshots and comparisons expose missing state/branch evidence and retain
scientific validation as unestablished. No automatic lowest-energy state choice,
perfect-model claim, numerical sensitivity-as-error-bound, or operating-rate
inference is permitted. Use exact source bindings and real mocked-producer
integration checks. Coordinate pair contracts directly, then send root a bounded
handoff. Read-only reviews do not grant permission to edit producer modules.

## Completed handoffs under root integration

These sources are frozen. Owners remain available for concrete findings; new
feature work requires an exact root assignment. No broad repeated audits, new
sessions or speculative calculation campaigns are needed.

| Lane / owner | Handoff | Current boundary |
|---|---|---|
| B1 / f522 | `nanodesign/bundle.py`, `tests/test_bundle.py`, `docs/EVIDENCE_BUNDLES.md`; 108 synthetic checks | Released to root; independent B2 has 11 checks and a synthetic roundtrip |
| B2 / 712e | `research/bundle-audit/` | Released to root; owner now P3 |
| Bundle CLI / internal `benchmark_validation` | `tests/test_cli_bundle.py`; 19 checks | Released to root |
| S2 / `01a0ac06-4280-7bc0-862a-35bfc19e958b` | `nanodesign/state_scan.py`, `tests/test_state_scan.py`, `research/state-scan/README.md`; 34 synthetic checks | Released to root |
| S3 / 8cfb | `research/state-scan-review/`; 24 independent checks, two findings repaired | Released to root; owner now CLI tests |
| H1 / `01a0ac06-b717-7cf0-a713-a95314895d23` | `research/characterization-resume/`; 23 checks including snapshot metadata, unapplied production proposal | Research prototype only; production workflow unchanged |
| H2 / `01a0ac06-9395-7fb3-bf22-9a85d406961f` | `research/characterization-resume-review/`; 36 independent checks plus two-case updated metadata/reuse receipt | Released to root |
| E2 / `01a0ac06-a551-7f81-aa57-6d76c1b4da55` | New stationary-force verifier/tests and E2 README section under `research/evidence-audit/`; 59 checks | Released to root; preserve historical E1 files and other owner's claims-ledger relocation |
| N1 / `01a0ac06-f040-7761-9b15-b5d876b37890` | `research/mode-comparison/`; 54 checks including actual mocked producer integration | Released to root; explicit snapshot precision and exact force-reference geometry binding |
| G2 / `01a0ac06-c5cd-72d1-b288-caa9bef48b95` | `research/gpu-readiness/`; 229 checks including strict v2 protocol | Released to root; no GPU adapter, execution or measured parity |

## Published foundation

`3eff5d5` published state surveys, portable evidence bundles, force/Hessian
reconstruction, mode comparison and the validation protocol. Remote CI
`35154502109` passed, including 547 core tests. These are software checks and
do not certify a molecule, reaction or assembled machine.

`9e7ac9f` passed 358 core tests and remote CI. `b550e54` published the local
workbench, computation planner and existing audits, also with passing CI.
`b2657c9` published offline GPU-readiness tooling with 98 checks; no GPU execution
or speedup was measured. The V1 workbench and earlier D1/D2/E1/Q1/Q2/P1/R1/M1
implementation/review lanes are complete. Their historical briefs do not grant
overlapping ownership of released core files.

All research audits use isolated fixtures and leave source evidence unchanged.
Report material reproducible issues with paths, expected/observed behavior and a
proposed repair. Leave new work uncommitted for integration unless root agrees
to a path-specific publication. Root uses an isolated Git index to avoid
capturing another agent's staging. Each owner maintains only their own status.
