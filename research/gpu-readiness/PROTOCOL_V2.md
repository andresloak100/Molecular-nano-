# Version 2: actual input files and resolved numerical choices

G2 extends the recorded-number comparison into a checked evidence chain. It
does not implement a GPU backend, run a calculation, or authenticate hardware
execution. The original `protocol/`, `local-preflight.json` and `verification.json`
remain unchanged. `protocol-v2/` is a new, unexecuted protocol with the same five
small cases. All actual numerical metadata is explicitly unavailable at planning
time; no default is promoted into a measured setting.

## Commands

```sh
.venv/bin/python research/gpu-readiness/gpu_readiness.py plan --output /tmp/new-protocol-v2
.venv/bin/python research/gpu-readiness/gpu_readiness.py verify --plan /tmp/new-protocol-v2/plan.json
.venv/bin/python research/gpu-readiness/equivalence.py compare \
  --plan /tmp/new-protocol-v2/plan.json --cpu /path/cpu.json --gpu /path/gpu.json
```

Planning and verification execute no chemistry. `verify` can succeed while
`resolved_numerics_available` and `numerical_parity_passed` remain false. Only
future execution records can supply the missing measured configuration.

The default `compare` requires schema version 2. The deliberately explicit
`--records-only` option retains the historical v1 arithmetic check; its output
never accepts parity evidence or claims to have verified files or resolved
numerics. The low-level `compare_records` Python function has the same restricted
v1 scope. Use `protocol.compare_protocol_files` for v2 acceptance.

## Checks before any accepted comparison

1. Read the plan and all its declared members through a pinned directory,
   descriptor-relative opens and no-follow flags. Canonical relative paths only;
   symlinks, directory escapes and non-regular files are rejected. The current
   implementation requires the POSIX facilities available on macOS/Linux; it
   refuses a weaker fallback on platforms lacking them.
2. Recompute every coordinate, source-record and attribution digest from bytes
   actually read. Verify membership and agreement with each case's digests.
   Limits are 128 members, 16 MiB per file and 64 MiB total protocol bytes.
3. Parse one finite XYZ frame, check ordered elements and Å declaration, and
   compare complete requested settings and declared changes with the preserved
   source record. Original-machine paths remain provenance strings; they are
   never followed or required to exist.
4. Require each CPU/GPU record to hash the exact full plan bytes with
   `plan_sha256`. Semantically identical but rewritten plan bytes are a different
   execution input and cannot reuse those records silently.
5. Require complete resolved numerical metadata on both records and reject
   missing, malformed or conflicting values. Only then apply the existing
   convergence, force/energy/spin, D3, synchronization and tolerance checks.

Original files are never rewritten during comparison. Hashes establish internal
consistency, not authorship or a trusted external reference. A producer able to
fabricate every record can also fabricate matching digests; the checker does not
claim execution attestation.

## Execution record additions

Retain the required fields in the main README, set `schema_version: 2`, and add
the exact plan digest plus `resolved_numerics`:

| Section | Required content |
|---|---|
| `precision` | `energy`, `forces`, `density`, each explicitly `float64` |
| `orbital_basis` | `sha256` of the actual expanded basis definition used; matching a basis-name string alone is insufficient |
| `auxiliary_basis`, with DF | `status: resolved` and `sha256` of the actual expanded auxiliary basis |
| `auxiliary_basis`, without DF | `status: not_applicable`; no invented basis hash |
| `quadrature` | `level`, per-element `atom_grid`, `pruning`, `radial_method`, `becke_scheme`, `point_count`, `canonicalization`, `signature_sha256` |

`atom_grid` maps every element present in the case to positive integer
`radial` and `angular` counts. The quadrature level must match the requested
level. Recipe identities must be explicit names; `auto`, `default`, `unknown`,
`unresolved`, null and empty strings are unavailable evidence. Boolean values
are not accepted as counts.

The future producer must hash the actual expanded orbital and auxiliary basis
definitions under the same documented serialization convention on both backends.
That convention must bind atom order, shell definitions, exponents, contraction
coefficients and units, not just a library label. This protocol can check supplied
basis digests for validity and agreement, but does not reconstruct absent basis
data or declare it authenticated.

## Exact grid identity without an ordering artifact

Use `resolved_numerics.canonical_grid_signature(rows)` on the actual quadrature
grid. Each row contains `(x, y, z, weight)`, with coordinates in Bohr and weights
in Bohr³. This retains point–weight pairing. The canonicalization identifier is
`paired-xyz-weight-f64le-sorted-v1`:

- Convert finite numeric values to IEEE-754 binary64 and normalize signed zero.
- Sort the complete four-value rows lexicographically; preserve duplicates.
- Hash the little-endian unsigned 64-bit row count followed by little-endian
  binary64 rows, with SHA-256.

Reordering rows does not change identity. Reassigning weights, changing a value
or adding/removing a row does. There is no quantization or hidden tolerance.
Grids differing only by floating-point roundoff can yield close energies yet fail
this exact-grid-identity requirement. Report that mismatch; do not reinterpret it
as an established chemical error or silently loosen the requirement.

Optional raw `grid_points_sha256` and `grid_weights_sha256` are checked as hash
strings but are not compared, because independent raw ordering can differ.
The canonical paired-row signature and actual point count govern identity.
As with basis hashes, the checker cannot authenticate an absent raw grid from
producer-supplied digest claims alone.

## Output and evidence limits

The strict output separates `input_files_verified`, `resolved_numerics_verified`,
`record_validation_passed` and `numerical_parity_passed`. A fully passing supplied
case sets `parity_evidence_accepted`; this is not acceptance of a complete GPU
adapter or a chemical model. State/ground-state/scientific validation stay false.
Missing resolved values produce `unavailable:` messages, bad types or structures
produce `malformed:`, and conflicting resolved choices produce `mismatch:`.

An actual GPU adapter still needs review of derivative response support,
dispersion inclusion, precision, transfer/synchronization and raw anchor forces,
followed by real calculations. No GPU execution or new energy/force result was
produced by G2.
