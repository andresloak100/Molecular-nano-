# Saved vibration-step comparison (N1)

This read-only tool compares two completed characterization records at the
same geometry and declared electronic inputs. It starts no quantum jobs,
changes no input records, and gives no convergence or transition-state
certificate. It complements E2's independent force-to-Hessian reconstruction
and H1's resumable characterization rather than replacing either.

## Use

From the repository root, using its existing environment:

```sh
.venv/bin/python research/mode-comparison/compare_modes.py \
  data/validation/h2-integration/modes/result.json \
  data/validation/h2-integration/modes-half-step/result.json

.venv/bin/python -m pytest -q research/mode-comparison
```

Use `--output NEW_FILE.json` to create a new report; existing files are never
overwritten. `--left-structure` and `--right-structure` supply explicit
structure snapshots when they are not named `input.extxyz` beside each result.
Their bytes must match `input_hashes.input_snapshot_sha256` when that key is
present. `structure_sha256` identifies the original input source, which may be
a binary or multi-frame trajectory and is not interchangeable with the saved
snapshot. For legacy records without the snapshot key, comparison is supported
only when the original source hash also exactly matches the snapshot bytes.
An explicit invalid snapshot key is rejected without falling back to the source
hash. Historical source-only records that do not meet that narrow condition
remain unsupported; their provenance is not invented or rewritten.
`--degeneracy-gap-cm1` sets the declared spectral grouping threshold (default
10 cm⁻¹). The threshold is a grouping choice, not a physical accuracy target.

The CLI exits **0** when a comparison was produced, and **2** when the inputs
are not comparable. Neither exit status means that a molecular design is
validated. Each file is limited to 64 MiB and the comparison to 600 free
Cartesian coordinates. No production quantum modules are imported.

## Binding before comparison

The file adapter accepts completed workflow `schema_version: 1`,
`stage: characterize` results with a `stationary` object. It verifies the
one-frame snapshot hash and explicit Å/eV/cm⁻¹ units, accepting the current
optional `force: eV/angstrom` unit; binds element order,
coordinates, isotope masses, free/frozen atoms and calculation settings; and
reconciles design/result/diagnostic quantum settings. An explicit diagnostic
starting guess cannot contradict those settings. Recorded force-reference
coordinates, when present, must match the snapshot exactly or after the known
eight-decimal ASE serialization. Exact full-precision force-reference coordinates
then become the primary cross-record geometry binding. Snapshot rounding and its
maximum difference are reported explicitly; differences between two full-precision
references are never hidden by the rounded snapshots. Historical records without
that reference remain bound only at the saved snapshot's precision.

The comparison requires identical coordinate frames and atom/free-index order.
It does not automatically align, rotate or permute molecules. That restriction
prevents a silent change of coordinate basis. Only the characterization's
displacement step may differ. Method, charge, spin, starting guess, grid,
electronic thresholds or known software-version differences reject a
step-only comparison. Missing software metadata is exposed as a limitation.

The one historical normalization is the schema-1 omitted `scf_initial_guess`:
the documented old default was `minao`. Every such inference is recorded in
`metadata_notes_by_source` and `legacy_or_missing_metadata`. Explicit null,
unsupported guesses, contradictory diagnostics and other missing state inputs
are rejected. This is input provenance, not a claim that two calculations
reached the same electronic solution.

Bare matrices or S1 summaries without the workflow's full geometry/method
binding are intentionally not accepted by the file adapter. H1 has confirmed
that its successful workflow-shaped output fits this adapter when the caller
supplies `input_context.design`. Its emitted result now includes the distinct
snapshot hash and workflow units; absence of design metadata is not inferred
away. `compare_bound(left, right)` is available for consumers that already
hold explicit geometry/masses/settings, but callers are responsible for the
provenance of that binding. It independently checks dimensions, index
partitions, valid inputs and mass-orthonormality before comparing.

## What is compared

For Cartesian mode rows V and diagonal atomic mass matrix M, the mass-weighted
rows are W = V√M. The full matrix `|W_left W_rightᵀ|²` reports squared overlaps,
independent of each eigenvector's arbitrary sign. Both complete bases must be
orthonormal in that mass metric first.

Frequencies are sorted for reporting. Adjacent **spectral ranks** are joined
when the gap in either spectrum is at most the declared threshold. This is
connected grouping: a chain can span more than the threshold, so every
group's frequency range is exposed. Inside a multi-mode group, the singular
values of the cross-basis overlap give principal-angle cosines. Their squares
measure subspace overlap without forcing individual matches between rotated
degenerate eigenvectors. A group covering the entire coordinate space is
explicitly flagged as directionally uninformative: complete orthonormal bases
always span that space.

Frequency shifts are differences by sorted rank, **not assertions of physical
mode identity**. The full overlap matrix reveals mode-order exchanges.
Raw negative-frequency counts and resolved imaginary-mode counts are separate;
a small negative-to-positive change inside the recorded resolution tolerance
remains unresolved. Identical displacement steps are explicitly identified and
do not provide a step-variation check.

## Executed evidence and verification

[The saved H₂ comparison](h2-saved-step-comparison-v3.json) reads the original
0.003 Å and 0.0015 Å records without recomputing chemistry. Its stretch-frequency
change is **−0.0895837007 cm⁻¹**. Both records have zero resolved imaginary modes
and five unresolved modes. The low-frequency five-dimensional group is compared
as a subspace; this does not certify those unresolved modes as stable.

The report includes source/structure hashes and the comparator implementation
hash. Legacy missing displaced forces remain visible; E2 reports their raw
Hessian reconstruction as unavailable. Original evidence has not been amended.

**54 focused tests pass.** Synthetic fixtures include a known 45° basis rotation
within a degenerate pair (individual squared overlaps 0.5, subspace overlaps
1 and 1), 30° leakage outside that group (principal squared overlaps 1 and
0.75), sign flips, unequal masses, spectral-order swaps, threshold crossings,
full-space grouping, invalid inputs and state/method provenance mismatches.
These fixtures are algebra tests, never electronic-structure results.

Interpretation remains bounded: two steps reveal observed sensitivity, not an
asymptotic error estimate. Agreement cannot establish state identity, endpoint
connectivity, chemical accuracy or a manufactured molecular machine.

The initial `h2-saved-step-comparison.json` is retained as an earlier comparator
snapshot. The v2 report retains the geometry-precision refinement. The v3 report above
adds distinct source/snapshot hash provenance and matches the current implementation hash. H1 confirmed that its exact checkpoint
manifest geometry and stationary force reference retain full precision; its XYZ
is a convenience snapshot and may use the same ASE text rounding.

## Real-producer integration regressions

`test_producer_workflow.py` calls the actual `workflow.run_characterization`
with only the production calculator replaced by analytical test forces. It
selects a nontrivial frame from a binary multi-frame trajectory, preserving
full-precision reference coordinates while verifying the distinct emitted
snapshot hash, current force units, source preservation and two-step comparison.
It also proves that an unsupported legacy source-only hash and an explicitly
wrong snapshot hash cannot be silently accepted.

`test_producer_h1.py` calls H1's actual `run_characterization_checkpoint` with
an analytical calculator and real caller-supplied design context. Its unmodified
successful producer records load and compare; source/context/evidence remain
unchanged and edited snapshots are rejected. These are software integration
fixtures, not quantum calculations or physical calibration. The tests depend
on the H1 prototype supplied in the same repository.
