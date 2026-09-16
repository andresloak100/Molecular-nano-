# H1 — resumable finite-difference characterization prototype

Owner: `codex-b717`, assigned by C1. All implementation in this directory is a
research prototype. It has not been added to the production CLI and has not
launched quantum calculations. Every executed acquisition test uses an explicitly
identified synthetic force calculator.

The prototype preserves completed force calculations when a characterization
fails or is interrupted. An explicit resume reads those records, checks that the
inputs still match, and acquires only missing forces in a **new output directory**.
The original attempt remains byte-for-byte unchanged. Final curvature analysis
uses the existing production `stationary.py` through a read-only force replay
adapter; it does not introduce a second implementation of Hessian mathematics.

## Executable interface

`prototype.py` exposes:

```python
run_characterization_checkpoint(
    atoms, output,
    settings=effective_quantum_settings,
    input_context=source_provenance,
    resume_from=None,
    step=0.005,
    fmax=0.03,
    frequency_tolerance=20.0,
    max_free_coordinates=120,
)
```

The supplied `atoms` must have the exact production `PySCFCalculator` attached
with the same `QuantumSettings`. That backend constructs a new mean field for
each geometry. Arbitrary stateful calculators are deliberately unsupported;
tests replace the type guard with an explicit controlled mock. Requested quantum
settings include the selected SCF starting guess. Input context must be supplied
explicitly; a production caller must obtain source hashes from the same byte
snapshots it parsed. The prototype also pins the actual in-memory geometry.

To inspect a saved checkpoint without running a solver:

```sh
.venv/bin/python research/characterization-resume/prototype.py /path/to/attempt
```

Inspection checks internal integrity and the request/geometry binding. Actual
resume additionally compares the entire manifest with the new invocation.

## Saved contract

Every attempt contains `input.extxyz`, `result.json`, `force_checkpoint.json`
and a persistent `.checkpoint.lock`. The checkpoint uses a versioned envelope:

```text
{payload: {...}, sha256: SHA256(canonical_JSON(payload))}
```

Canonical JSON sorts keys, uses compact separators, forbids NaN/Infinity and is
encoded as UTF-8. The payload contains:

- `manifest`: ordered atomic numbers, full coordinates and masses, fixed atoms,
  cell/PBC and initial atomic properties; complete quantum and characterization
  settings; caller source context; package versions and implementation hashes;
  explicit units and the exact baseline/atom/axis/signed-displacement plan.
- `manifest_sha256`, attempt identity, status and timestamps.
- `records`: a complete, ordered prefix of accepted force evaluations. Each has
  the exact request, full unconstrained force array, successful quantum
  diagnostics, calculation call ID and originating attempt ID.
- `resume_source`: path, original snapshot hash, attempt ID and original status.
- Failure details outside the reusable records, including failed quantum
  diagnostics when an acquisition fails.

The primary geometry is the full-precision JSON in
`result.checkpoint_manifest.geometry.positions_angstrom` (also in the checkpoint
manifest); the stationary force-reference JSON uses the same positions. The
`input.extxyz` file is a convenience rendering: this ASE installation writes
eight decimal places, so it can differ by up to its text-rounding envelope.
Reuse and exact cross-result identity must bind the JSON geometry, not that
rounded rendering. N1 identified and confirmed this consumer boundary.

New successful results separately bind the exact emitted snapshot bytes with
`result.input_hashes.input_snapshot_sha256`. The caller's original
`input_hashes.structure_sha256` keeps its source-file meaning, including a
multi-frame or binary source; it is never replaced with the snapshot hash.
The original caller context remains unchanged in `input_context` and the
checkpoint manifest. Top-level `result.units` explicitly declares length,
energy, force and frequency units to match the workflow consumer contract.

Each successful force record is validated and atomically persisted before the
next force request. Baseline forces are saved **before** checking the free-force
stationarity threshold, preserving anchor loads even for a rejected baseline.
Temporary-file replacement, file fsync and directory fsync keep a previous
complete snapshot available during a failed write. A failure to persist stops
acquisition. A hard termination between force completion and durable saving may
require repeating that last unsaved evaluation; the prototype does not promise
exactly-once solver execution.

On resume, the whole source snapshot is checked before creating output or
starting work. Changed geometry, masses, fixed atoms, settings, starting guess,
step, source context, versions or implementation fingerprints reject reuse.
Unknown/duplicate requests, invalid arrays, failed convergence diagnostics,
call-ID contradictions and digest mismatches also reject reuse. Invalid source
bytes are never repaired or overwritten automatically. The digest detects
integrity changes; it is not an authentication signature or proof that the
recorded forces are chemically accurate.

Exclusive output locks and shared source locks reject an active source worker.
After worker death, the operating system releases its lock; a remaining running
snapshot can be read without rewriting that source status. Locking currently
uses `fcntl` on macOS/Linux. Completed reuse does not inject cached forces into
the live quantum backend or invent electronic-log events. The backend's original
ASE cache and diagnostics are restored on success or failure.

## Counts and consumer compatibility

Successful results retain `result.stationary`, `quantum_settings`, `input_hashes`
and the existing `finite_difference_evidence` shape. A separate checkpoint block
reports planned/completed records, reused records, new acquisition calls and
lineage. New acquisition calls include a failing call; failed calls are excluded
from completed/reusable records. The planned record count minus the accepted
prefix length counts **missing record slots**, not work guaranteed to execute.
The saved baseline must pass the stationarity guard before further acquisitions;
for example, a nonstationary baseline can leave twelve missing slots while an
unchanged resume performs zero new calls and rejects the run again. Replaying
mathematics is separate from these solver counts.

For the executed three-free-coordinate fixture:

| Attempt | Durable records at exit | Source records reused | New acquisition calls |
|---|---:|---:|---:|
| Interrupted on third call | 2 | 0 | 3, including the interrupted call |
| Resume into new directory | 7 | 2 | 5 |
| Resume an already complete run | 7 | 7 | 0 |

E2's independent `research/evidence-audit/verify_stationary.py` successfully
verified a complete resumed synthetic result without a schema adapter. That
checks saved-array consistency; electronic state, transition-state connectivity,
chemical calibration and molecular-machine feasibility remain unverified.

## Validation and handoff

```sh
.venv/bin/python -m pytest -q --disable-warnings research/characterization-resume/test_checkpoint.py
```

The initial focused suite has **21 passing tests**: partial resume and original
source preservation, exact missing-call counts and origin IDs, incompatible
inputs, corrupted/failed records, nonstationary baseline retention, and cache
restoration. H2/9395 independently reviews interruption boundaries, source locks
and remaining-work counts in `research/characterization-resume-review/`; E2 owns
the algebra verifier and N1 owns comparisons between mode calculations. H2's
**36 independent checks passed** against the frozen prototype, including
baseline/plus/minus failures, fifteen identity mismatches, ten malformed-record
cases, active locks, complete replay and persistent atomic-write failure. H2
found no implementation defect and clarified the missing-slot wording above.

The later C1/N1 integration follow-up adds two snapshot-provenance regressions
and updates the partial-resume check. All three focused cases pass: actual
multi-frame source versus exported snapshot, full-precision JSON preservation,
untouched caller provenance, explicit units, and zero-call completed reuse.
The original 36-check H2 receipt describes its earlier frozen source revision;
the separate updated metadata smoke receipt records two passing existing checks
(complete zero-call replay and chained partial-resume lineage), with unchanged
source hash `7429161b9be1bfd2c512a1e2580de1bb126cfe199f5bb405935e61dd61c32187`.
N1 also ran two actual-H1 producer/consumer integration tests successfully,
without editing emitted JSON: distinct source/snapshot hashes, explicit units,
full-precision reference binding and tampered-snapshot rejection are covered.

`production-hook.patch` is an **unapplied proposal** for C1 review. It requires
first promoting the reviewed prototype to `nanodesign/characterization_checkpoint.py`
and reviewing its imports, backend restriction and implementation fingerprints
for that production location.
It adds opt-in Python workflow parameters only and preserves the existing
default path; CLI integration is a separate C1 decision. Prototype checkpoints
will intentionally fail the implementation-fingerprint check after promotion
unless a separately reviewed migration is added. No production file has been
modified by H1.
