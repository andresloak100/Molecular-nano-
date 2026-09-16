# Local quadratic stationarity diagnostic — L1

This research tool uses saved residual forces and a saved Cartesian Hessian to
estimate the displacement to the **stationary point of their local quadratic
model**. It performs no actual relaxation, electronic calculation, or change to
source evidence. A result is not a rigorous energy-error bound, a corrected
barrier, or evidence of an operating molecular tool.

A small maximum force alone can conceal a large displacement along a soft
direction. At a saddle, uphill and downhill energy contributions can also
cancel. This tool reports both effects explicitly. It complements N1's
step/mode comparison and E2's saved-force reconstruction; it does not implement
another optimizer, vibration solver, or energy-versus-force check.

## Use with an existing workflow or H1 record

From the repository root:

```sh
.venv/bin/python research/local-relaxation-diagnostic/diagnose.py \
  PATH_TO_COMPLETED_CHARACTERIZATION/result.json \
  --curvature-floor 0.001 --max-condition 100000 \
  --max-atom-shift 0.1 --max-relative-asymmetry 0.001 \
  --expected-index 0 --output NEW_DIAGNOSTIC.json
```

These example thresholds illustrate syntax only. They have **not** been
calibrated for any molecular operation. Every threshold is required; there
are no implicit accuracy targets. Curvature is in eV/Å², maximum atom shift is
in Å, and condition/asymmetry limits are dimensionless. Index 0 requests a
positive-curvature model; index 1 requests a first-order saddle model.

`--structure` supplies a nonstandard snapshot path, still subject to its
recorded hash. Omit `--output` for JSON on standard output. An existing output
is never overwritten. Exit 0 means an estimate was produced within the
declared diagnostic limits; exit 2 covers inconclusive, refused, unavailable,
or invalid input/output. Neither exit code certifies scientific accuracy.

## What the estimate means

In the existing free Cartesian coordinate order, use `g = -F_free` and

```text
E_quadratic(q + delta) = E(q) + g^T delta + 1/2 delta^T H delta
delta = -inverse(H) g
Delta E_quadratic = E_quadratic(q + delta) - E(q)
                  = -1/2 g^T inverse(H) g
```

These are the full Newton stationary step and its local quadratic energy
prediction. Geometry optimizers assess and constrain such local models using
additional evaluations and trust-radius logic; see the primary
[geomeTRIC optimization documentation](https://geometric.readthedocs.io/en/latest/how-it-works.html)
and [transition-state documentation](https://geometric.readthedocs.io/en/latest/transition.html).
This diagnostic does not evaluate the predicted new position and does not
establish a region where the quadratic model is accurate.

The code diagonalizes the saved **Cartesian** Hessian. It does not mistake
mass-weighted frequencies or mass-normalized mode vectors for a Euclidean
Cartesian eigensystem. Masses are checked for input consistency in the adapter;
they do not enter the Cartesian Newton formula.

For each Cartesian eigenpair `(lambda, v)`, with `a = v^T g`, the signed energy
contribution is `-a²/(2 lambda)`. Positive-curvature contributions lower the
quadratic energy; negative-curvature contributions raise it toward a saddle.
The report includes both sums and the sum of their absolute mode contributions.
That nonnegative sum is only a diagnostic scale: it is not a physical energy
error bound. Individual components within degenerate eigenspaces depend on
the eigensolver's basis; total displacement and grouped sums are invariant.

The displacement includes all three components of every recorded free atom.
The reported maximum is the largest atom-vector norm; RMS uses the number of
free atoms, not the number of Cartesian components. Adapter results additionally
give the complete atom order with exact zero displacement at frozen anchors.
Anchor forces do not enter the free-coordinate Newton step. No coordinates are
exported as an optimized geometry. Nonlinear constraints are outside this API.

## Refusal and interpretation

- **`refused`:** any absolute curvature at or below the caller floor, excessive
  absolute spectral condition number, excessive raw Hessian asymmetry, invalid
  evidence, or nonfinite numerical operation. No inverse is taken after a
  soft/conditioning refusal, and `estimate` remains null. Even a soft mode with
  zero projected gradient is refused; the tool does not assume it is external.
- **`inconclusive`:** an invertible quadratic model has a different index than
  requested or predicts a maximum atom shift above the caller limit. The full
  unmodified estimate remains visible. It is never clipped into acceptance.
- **`estimated`:** a full-rank quadratic estimate meets those declared checks.
  Nonlinear validity, actual stationarity and electronic-state continuity
  remain unverified. In particular, a small linearized gradient residual only
  shows that the quadratic linear system was solved.
- **`unavailable`:** the source lacks the full saved-force evidence needed for
  reconstruction. Historical records are preserved without invented forces.

The pure numerical API requires an explicitly symmetric Hessian and reports
the unsigned spectral condition `max(abs(lambda))/min(abs(lambda))`.
Soft modes are not deleted, sign-flipped, regularized, or inverted with a
pseudoinverse. No translation/rotation removal is applied, including in an
anchored system. Raw antisymmetry is checked separately in the file adapter;
symmetrization does not repair noisy or electronically inconsistent forces.

## Evidence binding and API

`diagnose_quadratic(H, g, *, curvature_floor_ev_per_angstrom2,
max_condition_number, max_atom_shift_angstrom, expected_index)` is the pure
numerical API. `H` and `g` are finite JSON-style lists, with shapes
`(3*n_free_atoms, 3*n_free_atoms)` and `(n_free_atoms, 3)`. It labels its inputs
as unbound arrays. Invalid shape/units-as-values/controls raise `DiagnosticError`;
valid but unsuitable numerical models return the statuses above.

`diagnose_saved_result(path, *, structure_path=None,
max_relative_hessian_asymmetry, **numerical_options)` performs the actual
workflow/H1 integration. It reuses, without editing:

1. N1 `load_saved_result` and `checked_modes` for source/snapshot hashes,
   element/mass/order/free/frozen consistency, coordinate convention, explicit
   units and declared quantum settings.
2. E2 `verify_stationary` for the complete baseline/displacement force stencil,
   same-reference coordinate checks, gradient sign, Hessian reconstruction and
   saved numerical metadata.

The exact saved force-reference coordinates must bind to the snapshot exactly
or through its declared eight-decimal serialization. The original source hash
may refer to a selected frame of a binary trajectory and is distinct from the
exported snapshot hash. The report preserves both, exact reference coordinates,
metadata limitations, verifier results and implementation/dependency hashes.
Older N1 default-guess interpretations, if encountered, remain explicit notes;
they do not establish the electronic state. Hashes bind bytes, not physical
authenticity. E2 does not verify the electronic event log or state continuity.

This adapter intentionally refuses a standalone summary that lacks workflow
provenance; the pure API is available for explicitly unbound numerical work.
Files are subject to N1's 64 MiB limit and calculations to 600 free coordinates.
The adapter loads N1/E2 and NumPy/ASE support, with no production quantum
module imports. Its report is an assessment of the bytes loaded at that time;
it does not lock live producer directories.

## Verification and limits

**55 L1 checks pass:** 45 closed-form cases and 10 real producer/consumer
integration checks. The latter call the actual production workflow and H1
checkpoint producer with only the calculator replaced by analytical forces.
They cover nonzero gradient, coupled curvature, unequal masses, large anchor
forces, selected binary input frames, exact coordinates behind text rounding,
source immutability, no added calculator calls, altered snapshot/force rejection,
missing evidence, and refusal to overwrite source files. These are software
fixtures, never molecular predictions. L2 independently reviews the formulas
and actual H1 integration in `research/local-relaxation-review/`.

```sh
.venv/bin/python -m pytest -q research/local-relaxation-diagnostic
```

The next scientific evidence remains branch continuity, converged force and
curvature sensitivity, actual optimized structures and forward/backward
connectivity. This tool supplies a more informative local diagnostic while
those questions remain open. It cannot turn a harmonic estimate into a
calibrated chemical or operating accuracy claim.
