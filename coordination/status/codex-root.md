# Codex integration status

Updated 2026-09-16. Root task: `01a0abd2-e250-7a62-a6d6-91f406c345fd`.

The electronic and derivative evidence milestone is published at `c66be23` and
passed remote CI `35156590649` on 2026-09-16 at 22:15 UTC. Existing S1/A1/A2
chemistry stays with its owners; no new
local SCF or gradient campaign was launched by this milestone. None of the
software checks establishes chemical accuracy or an operating molecular tool.

## Reviewed changes

- Optional current-call converged-SCF snapshots in the CPU backend, separate
  from gradient/whole-call acceptance. Failure and interruption clear results;
  saved SCF evidence remains correctly identified. Independent lifecycle review
  and 23 core synthetic cases cover cache restoration and secondary failures.
- E3 snapshot export/read/validation and same-geometry spin-resolved occupied
  subspace/density comparison: 92 owner tests and 35 independent E4 tests.
  Review repaired incomplete settings, an actual import defect and a metric
  compatibility error that could yield direction-dependent comparisons.
- S2 opt-in capture freezes the option and binds local artifact bytes, settings,
  coordinates and call identity. S5 adds 21 actual producer/serializer checks
  with synthetic numerics and real installed imports. Failed-gradient snapshots
  now also bind available SCF reference and AO-count facts.
- The read-only state-compare CLI has 31 independent checks. Survey preparation
  exposes --capture-electronic-state without launching work; two additional CLI
  option checks supplement the previous 24 cases.
- X1 cross-AO research bridge: 37 tests including tiny real one-electron integrals;
  explicit expanded basis, pair/frame binding and metric-relative reconstruction
  checks. It is not integrated automatic branch tracking.
- F1 saved total-energy/force consistency: 42 owner plus 39 independent tests;
  actual mocked H1 integration, unequal displacement formula, finite arithmetic
  and energy-rounding warnings. No physical error bound or branch certificate.
- L1 local quadratic residual correction: 55 owner plus 25 independent tests;
  Cartesian Hessian/force reconstruction, soft-mode and conditioning refusal,
  separate signed saddle contributions. No actual relaxation or nonlinear bound.
- G3 documents actual CPU fields and the distinction between retained SCF grids
  and unavailable transient force grids. No GPU adapter or parity claim.

## Verification and publication

The full local non-quantum suite passed 706 cases before the two final CLI flag
cases were added; the CLI suite verifies those separately. Root also reran all
243 new independent/research cases (E4, lifecycle, X1, F1/F2 and L1/L2), all passing.
Remote CI passed all **720 core tests**, including four new tiny actual
RKS/UKS x direct/DF capture cases and the eight earlier quantum cases. No local
SCF/gradient work was run for the new cases. All research, workbench and prior
integration steps passed as well. Exact run:
https://github.com/andresloak100/Molecular-nano-/actions/runs/35156590649

Prior published checkpoint `3eff5d5` passed CI `35154502109`, including 547 core
checks. The workbench remains at http://127.0.0.1:8765.

Root publishes explicit reviewed paths with an isolated Git index, preserving
other agents' live evidence and staging. Completed lanes are frozen. Additional
scientific interpretations in A1/A2/the secondary digest require the concrete
corrections in messages/from-root-to-science-summary-owners.md; they are not
validated merely because their files appear elsewhere in this shared repository.
