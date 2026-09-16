# Independent S1 saddle-evidence review

Q2 support, `codex-f040`, 2026-09-16. This directory contains a read-only
audit, cheap synthetic regression fixtures, and snapshots of existing evidence.
It starts no quantum jobs and changes no S1 calculation files. It does not
certify a transition state or a molecular assembly tool.

## Findings delivered to S1 and C1

1. **New-geometry electronic ambiguity is already observed.** In the saved
   M06-2X-D3(0)/def2-TZVP relaxed ethynyl scan, all four guesses converge but
   span **2.208784977 kcal/mol**. `1e` gives -76.59667371933023 Hartree;
   `atom` gives -76.59315379638618 Hartree. The prior agreement at the
   published geometry does not cover this optimized geometry. Carry the
   electronic ambiguity alongside geometric saddle candidacy, and explicitly
   label the optimized transition structure's scan as not performed. Selecting
   the lowest energy does not establish state identity. If the saddle energy
   were held fixed, changing just the reactant reference would change the
   energy difference; a consistent saddle on the intended state is still needed.

2. **The transfer-direction heuristic accepts negligible axial noise.**
   At the real seed coordinates, a synthetic mass-normalized displacement of
   only H3 in direction `[0, 1, 0]` is rejected, while `[1e-14, 1, 0]` is
   accepted. The latter is effectively transverse: its two bond-length
   derivatives are about +/- 1e-14. Reversing the mode sign or rotating the
   structure and displacement together preserves this outcome. Opposite
   derivative signs do not measure meaningful transfer participation; for
   nonzero derivatives, the added antisymmetric-dominates test is redundant
   with the opposite-sign test. Use a declared normalized projection magnitude
   and an uncertain region, retaining the heuristic's limits. This is a
   **synthetic adversarial probe**, not an observed computed vibration.

3. **Completed electronic evaluations are undercounted.** Each of the six
   completed reactant optimization records has one more distinct
   `calculation_completed` call ID than its `electronic_evaluations` field:
   methane 6 versus 5, ethynyl 5 versus 4. The optimizer engine history excludes
   the final energy/force call. Rename the history count to optimizer requests
   or separately report the actual completed calculator calls. Existing logs
   suffice; no rerun is needed.

4. **A failed Hessian can leave an apparently unfinished record.**
   The calls to `characterize` and `analyze_transfer_mode` in S1's `main` are
   outside the optimizer failure handler. A displaced-geometry SCF failure
   exits before recording a characterization failure or study completion.
   Save an explicit stage failure and retain optimizer results. This is a
   source-review finding, not an observed failure of the current jobs.

The earlier fixes to mass weighting, force stationarity, candidate language
and converged-subset wording are present in the source reviewed here.
The running jobs loaded an earlier source revision: their current summaries
still use the older guess-scan fields. Keep the original records intact and
record postprocessing/version provenance separately. No optimizer-convergence
bug was found: installed geomeTRIC raises when its iteration limit is exhausted.

## Reproduce

From the repository root with its existing environment:

```sh
.venv/bin/python -m pytest -q research/saddle-audit
.venv/bin/python research/saddle-audit/audit.py
.venv/bin/python research/saddle-audit/reproduce_review.py
```

`audit.py` discovers S1's `saddle_study.json` files even under ignored run
directories. Explicit file paths can also be supplied. It checks saved energy
arithmetic, atom balance when recorded, scan completeness and spread, gradients,
frequency thresholds, dimensions, and mass-orthonormal mode vectors. It reports
pending/missing evidence separately. `--output NEW_FILE.json` refuses to replace
an existing file. Normal stdout output is also available. The script exits
successfully after generating a report even when findings exist; consumers
must inspect `findings` and missing evidence. No passing audit establishes
connectivity, electronic-state identity, model accuracy or manufacturing.

`reproduce_review.py` imports only S1's definitions and calls its algebraic mode
analysis using synthetic directions at the actual saved seed coordinates. It
also counts saved completed calculator events. It never calls the quantum
calculator. These fixture directions and assigned frequencies must never be
used as scientific results.

## Recorded review snapshot

- [Saved-record audit](observed-audit-final-20260916.json): three incomplete
  saddle studies, all with completed reactants and no saved transition result
  at inspection time; the M06-2X ethynyl ambiguity is retained explicitly.
- [Mode probes and evaluation counts](review-reproductions-20260916T2108Z.json):
  nine synthetic direction/transform probes, six saved optimization/event
  comparisons, and SHA-256 hashes of inspected inputs/source.
- Focused tests passed **44/44**. They use synthetic records, not
  substitute electronic-structure calculations.

Snapshots are point-in-time observations of active work. A source-file hash
records the code inspected, not proof of which source version an already
running process loaded. The audit does not infer whether a process is alive
from the absence of a terminal study record. Multi-file snapshots are not an
atomic capture of running jobs.

The final audit also binds the Hessian coordinates and free-atom set to the
recorded optimized transition geometry, and checks every recorded guess energy
even if a guess manifest omits it. The initial snapshot is retained alongside
the final snapshot; the final report includes the auditor implementation hash.

## Independent projection now implemented

`mode_evidence.py` supplies a read-only second opinion on the recorded mode.
For the bond-distance coordinate q = r(H,A) - r(H,D), gradient g, Cartesian mode
v and diagonal mass matrix M, the score is
`(g·v)^2 / [(vᵀ M v)(gᵀ M⁻¹ g)]`. It lies between zero and one and is invariant
to mode sign/amplitude and rigid coordinate rotation/translation. The default
0.25 overlap threshold is explicitly a screening choice, not calibrated chemical
accuracy or a transfer probability. Hydrogen displacement share and the two
bond-distance derivatives are reported separately: carbon motion can also
change this coordinate.

The audit now attaches this independent projection to valid geometry-bound S1
characterizations with reaction indices, and flags disagreements with S1's
existing transfer heuristic. Geometric saddle candidacy stays separate from
this directional check and from unverified endpoint/electronic-state identity.

[New synthetic comparison](review-projection-reproductions-20260916.json) records
S1's answers alongside the independent score. At the actual seed coordinates,
H-only axial motion scores 0.959728326; a transverse direction with 1e-14 axial
noise scores approximately 9.60e-29, below the declared threshold. These are
synthetic direction probes, not computed vibrational results.

[Current saved-record audit](observed-audit-with-projection-20260916.json) includes
hashes for both auditor and projection implementation. Earlier snapshots remain
unchanged. All 44 focused tests passed, including mass-metric invariance, invalid
inputs, geometry binding, misleading transfer flags, and unchanged input records.
