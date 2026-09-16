# P1 — reaction-path and force-contract audit

Owner: support task `01a0ac06-b717-7cf0-a713-a95314895d23`.
Review performed on 2026-09-16. This directory contains **synthetic software
tests**, not molecular calculations or evidence of chemical feasibility.
The core integration owner C1 applied the job-status repair and authorized this
support task to repair `nanodesign/stationary.py` with new regression tests in
`tests/test_stationary_resolution.py`. Existing scientific artifacts were left
unchanged; no quantum jobs were launched.

## Confirmed defect: setup failures leave a running record

`nanodesign/workflow.py`, lines 121–126 in the original reviewed version, wrote both
input snapshots and constructs the initial/final calculators before entering
its `try`/`except`/`finally` block. An exception or user interruption there
escaped with `result.json` still reporting `status: running`, without the
terminal error, validation, or elapsed-time fields. The process has already
terminated. This does not create a false barrier, but makes the job record
misrepresent its lifecycle and loses failure provenance.

`test_path_contracts.py` reproduces an `OSError` during the final snapshot,
a `KeyboardInterrupt` at the same point, and calculator initialization failure.
The reproduction deliberately leaves the JSON destination writable; no code
can guarantee a final persisted status if that destination becomes unwritable.

`setup-failure.patch` records the fix C1 has now applied: move those operations
into the existing error boundary. Before that fix, the audit had **2 passing
and 3 failing cases**. All five now pass on the shared implementation. The patch
is retained as the handoff record and should not be reapplied.

## Repaired defect: unresolved steps can erase an unstable mode

A synthetic H atom with mass 1 and Hessian `diag(-1, 2, 3)` has one negative
mode. At the default 0.005 Å step the original code recovered that result using
seven force evaluations. At accepted steps of 1e-15 Å and below, ASE reused the
baseline forces for every displacement: only one evaluation occurred, yielding
a zero Hessian and no resolved negative modes. Large coordinate magnitudes can
also round a requested displacement to zero or strongly distort its size.

The repair preflights every signed displacement before force work. It rejects
steps invisible to ASE's position cache and offsets that do not represent the
requested step within the declared one-part-per-million relative tolerance.
That tolerance guards coordinate representation; it is not a force-error or
chemical-accuracy claim. A runtime check also rejects an attached calculator
with a coarser cache instead of accepting its stale forces. Existing calculator
geometry, results and diagnostics are restored on success or failure.

## New force evidence for independent reconstruction

C1 also requested a prospective repair for E1's finding that archived Hessians
could not be reconstructed from saved displaced forces. New successful results
include `finite_difference_evidence` with reference coordinates, full baseline
forces, and the full force array for each signed displacement. Each record
identifies the atom, axis, requested and actual offset, displaced coordinate,
and diagnostic calculation call ID where available. Length and force units,
ordering, and force scope are explicit; fixed-anchor forces remain included.

For column j, select the free-atom force components and compute
`H[:, j] = -(F_plus - F_minus) / (2 * step_angstrom)`. The raw matrix yields the
asymmetry diagnostics; `(H + H.T) / 2` yields the reported symmetric Hessian.
The independent regression uses a deliberately nonconservative synthetic field
so a nonzero asymmetry must survive reconstruction and JSON serialization.

This adds evidence to future complete characterizations. It does not recover
missing forces from earlier runs or establish chemical accuracy, reaction
connectivity, or a validated transition state.

## Passing behavior

Using a known analytic double-well and real ASE/FIRE/NEB execution, the added
checks confirm:

- A seven-image path preserves both anchors exactly and recovers the known
  0.2 eV synthetic barrier.
- Opposing 2 eV/Å anchor loads remain visible individually despite a zero
  global sum; external holding-force signs are the negatives of internal loads.
- Raw anchor forces remain in the result, while the free-force norm excludes
  them. Numerical path completion keeps `design_validated: false`.
- An unconverged endpoint stops before NEB, saves a failed status and free-force
  evidence, and does not report a candidate barrier.

The existing workflow and stationary suites also passed **32 tests** before
the proposed change. These checks exercise software contracts only and do not
calibrate the electronic method or validate a molecular design.

## Reproduce

From the repository root:

```sh
.venv/bin/python -m pytest -q --disable-warnings research/path-contract-audit/test_path_contracts.py
.venv/bin/python -m pytest -q --disable-warnings tests/test_workflow.py tests/test_stationary.py tests/test_stationary_resolution.py research/path-contract-audit/
```

The combined focused suite passed **52 tests** on the repaired shared tree.
No expected-failure markers hide the original defects. ASE 3.29.0 emits NumPy
deprecation warnings in this environment.

## Handoff

- Setup failure finding and validated patch sent directly to C1 and Q1/4280.
- Q1/4280 identified the initial setup-boundary concern; this audit supplied
  the independent reproduction, focused patch, and validation.
- C1 applied the workflow patch. The authorized stationary implementation and
  new regression file are ready for C1's final integration review.
- The independent mode helper checked unequal masses, coupled Hessians,
  frozen anchors, cache restoration and force-evidence reconstruction.
