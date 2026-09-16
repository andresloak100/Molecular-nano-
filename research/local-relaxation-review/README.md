# L2 independent local quadratic diagnostic review

Review completed 2026-09-16 for L1's research diagnostic. **25 independent
checks passed**, with no material defect reproduced in the reviewed scope.
All numerical examples are synthetic. This receipt does not validate a
molecular electronic state, nonlinear energy accuracy, a transition state,
or a physical mechanosynthesis operation.

## Evidence and scope

`test_independent.py` independently compares the reported displacement with
`numpy.linalg.solve(H, -g)` and evaluates the quadratic polynomial directly.
It checks coupled Cartesian coordinates, force/gradient sign reversal,
rotation invariance, positive and negative saddle energy changes, exact
cancellation with nonzero displacements, soft and ill-conditioned directions,
index mismatch, symmetry and explicit control validation. An unforced soft
mode still prevents full-space inversion; no silent subspace claim is accepted.

`producer_fixture.py` drives the actual H1 checkpoint producer with an analytic
force backend. Its three-atom structure has a fixed middle atom, unequal masses,
a nonzero free gradient, coupled free-coordinate curvature and large support
forces. Thirteen synthetic acquisitions produce the unmodified saved evidence
consumed by L1 through N1 and E2. Checks cover full-precision geometry behind
rounded snapshots, correct free/all-atom mapping with zero anchor shifts,
mass-independent Cartesian estimates, byte-preserved evidence, and refusal of
snapshot, geometry, force and settings tampering. Missing legacy force evidence
cannot produce an estimate. A separate nonconservative synthetic force matrix
reconstructs successfully but is refused by the raw-Hessian asymmetry gate.

The explicit force-tolerance counterexample has a gradient of 0.029 eV/angstrom
along a curvature of 0.0001 eV/angstrom². Its quadratic stationary displacement
is 290 angstrom and signed change is -4.205 eV. L1 marks the oversized step
inconclusive and establishes no energy-error bound. These are analytic fixture
values, not a prediction for a molecular structure.

## Reproduction

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q research/local-relaxation-review/test_independent.py
```

Final receipt: **25 passed in 2.85 s**. ASE emitted 1,260 existing NumPy shape
assignment deprecation warnings; no test failure. No electronic-structure jobs,
production edits, actual relaxation or source-artifact rewrites were performed.

## Exact reviewed sources

| Source | SHA-256 |
|---|---|
| `research/local-relaxation-diagnostic/diagnose.py` | `2dc7a9ebe0c367ecf1a97780efd1021d400d711f770ebdee4a957a08b4c7fc05` |
| `research/characterization-resume/prototype.py` | `7429161b9be1bfd2c512a1e2580de1bb126cfe199f5bb405935e61dd61c32187` |
| `research/evidence-audit/verify_stationary.py` | `589c94b7510a6bc5df17bcc93d79277c808b3acfea8aecac2cb5dcb7eafea972` |
| `research/mode-comparison/compare_modes.py` | `2ad76f953a41f324f457afe3c0c386d20add0d8bb8de5a2e19f66fdc4449936f` |

The source owner declared diagnose.py frozen before the final pass; its hash
matched the earlier 24-check pass. The independent review owns only this
folder and its own coordination status. L1 owns its implementation and additional
producer tests. Changes are left uncommitted for C1's integration.

## Interpretation limits

The formula is a signed stationary Newton estimate in the recorded free
Cartesian space. For an indefinite Hessian it need not lower energy. Positive
and negative mode contributions can cancel despite large coordinate shifts.
Cartesian curvature is not a mass-weighted frequency. Matching saved inputs,
reconstructing forces and successful mock flags do not prove electronic branch
continuity or force accuracy. Caller-selected curvature, conditioning,
asymmetry and displacement limits are diagnostic gates, not calibrated
physical tolerances. No real quantum result was evaluated in this review.
