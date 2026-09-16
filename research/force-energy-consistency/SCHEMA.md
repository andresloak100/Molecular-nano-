# F1 normalized record schema, version 1

Public interfaces in `consistency.py`:

```python
analyze_stencils(document, *, force_tolerance_ev_per_angstrom=None) -> dict
read_stencil_document(path) -> dict
from_h1_checkpoints(paths) -> dict
```

The analyzer is pure and does not mutate its argument. The two readers do not
acquire calculations. `read_stencil_document` rejects duplicate JSON keys,
nonfinite numbers and inputs larger than 64 MiB. Native checkpoint reading uses
H1's existing validated loader. It imports producer definitions but never calls
an acquisition function. No files named by a normalized source ID are opened.

The document is an object with integer `schema_version: 1` and a nonempty list
`stencils` (at most 1800). Extra source metadata is preserved in the document's
canonical digest. Each stencil contains:

| Field | Type and meaning |
|---|---|
| `atom_index`, `axis` | Integer atom index and Cartesian axis 0, 1 or 2 |
| `requested_minus_offset_angstrom` | Finite negative JSON number |
| `requested_plus_offset_angstrom` | Finite positive JSON number |
| `reference` | Baseline evaluation; null explicitly means unavailable |
| `minus`, `plus` | Evaluation or null for an absent/unacquired slot |

Every evaluation contains these fields, including explicit nulls where allowed:

| Field | Type and meaning |
|---|---|
| `context` | Calculation identity object defined below |
| `positions_angstrom` | Finite N×3 full-precision coordinates in unchanged atom order/frame |
| `source_sha256` | Lowercase 64-digit SHA256 of captured source artifact, or null if unavailable |
| `source_record_id` | Nonempty record/slot identifier within that artifact, or null |
| `calculation_call_id` | Nonempty originating calculation ID, or null |
| `status` | `completed`, `failed`, `interrupted` or `not_run` |
| `scf_converged`, `gradient_completed` | Explicit booleans |
| `energy_ev` | Finite total potential energy or null; never SCF-only if forces include dispersion |
| `forces_ev_per_angstrom` | Finite N×3 raw forces or null; the baseline is required for a comparison |
| `state_evidence_id` | Nonempty opaque per-geometry identifier or null; only a declaration |

Required `context` fields:

- `atomic_numbers`: a nonempty integer list (1–118, at most 600 atoms).
- `fixed_indices`: unique, in-range integer indices. All others are free.
- `pbc`: exactly `[false, false, false]`.
- `units`: exactly `{"length":"angstrom","energy":"eV","force":"eV/angstrom"}`.
- `force_scope`: exactly `"raw_unconstrained"`.
- `energy_scope`: exactly `"total_potential_energy"`, meaning the total potential
  whose negative derivative is the supplied force, including configured D3 once.
- `quantum_settings`: all explicit current `QuantumSettings` fields: `charge`,
  `spin`, `xc`, `basis`, `dispersion`, `grid_level`, `conv_tol`, `max_cycle`,
  `threads`, `memory_mb`, `density_fit`, `scf_initial_guess`. No defaults are
  filled in. Complete extra recorded settings participate in exact matching.
- `software_versions`: mapping of version names to strings or explicit nulls.
- `resolved_numerics`: explicit quadrature/auxiliary-basis policy object or null.
  An empty/missing policy is not evidence of equivalent realized numerics.

Scientific context fields explicitly absent/null make the study unavailable;
malformed types, wrong units or contradictory nonmissing bindings make it invalid.
Missing required structural schema keys are malformed input. Numeric booleans,
numeric strings and nonfinite values are rejected. Repeated contexts match by
canonical JSON with sorted keys and compact separators, retaining numeric types.
All stencils must refer to the same exact baseline geometry and context.

The requested positions are calculated by one binary64 addition to the baseline
coordinate. Saved displaced positions must match that result exactly, with no
other coordinate changes. Actual offsets must retain the requested sign and
agree within relative error 1e-6. Duplicate represented stencils are invalid.
No coordinate registration, unit guessing or reconstruction from rounded XYZ is
performed.

An evaluation call may recur when the measured payload is exactly identical,
even under another captured checkpoint hash after a legitimate resume. A call ID
associated with different context, geometry, energy, forces or state metadata is
contradictory. Pure-API source hashes are declared bindings, not checked external
files. The adapter's source digest is obtained from the exact native snapshot
validated by H1 under a shared lock.

For H1, original total Hartree values and conversion factors accompany the
evidence under `source_energy_conversion`. These conversions are explicit;
missing values are not replaced, and different known factors are refused.
The adapter retains failed-attempt details and source input context under
`source_checkpoints`. Every planned stencil is emitted, including missing slots.
Do not treat the number of missing slots as a promise about future execution.

Output `comparisons` include the signed residual `F_FD − F_analytical`, its
absolute value, literal threshold comparison, actual steps, interpolation
weights, rounding sensitivity and bound source/call/state identifiers. An
unavailable slot is kept in `unavailable_stencils`; missing global context leaves
the study unavailable with an explanation in `findings`. Contradiction invalidates
the whole study so an earlier partial success is not returned as accepted evidence.

Only successful comparable stencils contribute to `step_sensitivity`, whose
`stencil_indices` identify the subset. No comparison ranks geometries, selects
an electronic state or certifies the declared scientific model.
