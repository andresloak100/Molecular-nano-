# Fixed-geometry DFT starting-guess surveys

`nanodesign.state_scan` compares explicitly requested SCF starting guesses at
one frozen molecular geometry. It addresses the multiple-solution issue found
in this project's open-shell calculations without selecting a ground state or
calling equal energies the same electronic state.

The module uses the production `QuantumSettings` and `PySCFCalculator`. It does
not change the quantum engine, optimize coordinates, run orbital-stability
analysis, build reaction paths, or identify a preferred guess. Development and
regression checks use synthetic calculators; no new scientific calculation has
been executed to validate this survey workflow.

## API

The installed command line exposes the same bounded workflow:

```sh
# Preparation and inspection only; no electronic calculations.
nanodesign state-scan-create examples/h-abstraction/initial.xyz \
  --settings examples/h-abstraction/design.json --out runs/tip-guess-survey
nanodesign state-scan-report runs/tip-guess-survey

# Starts one real energy-plus-gradient calculation; run when capacity permits.
nanodesign state-scan-run runs/tip-guess-survey --max-jobs 1
```

`--settings` is required and must explicitly include `charge` and `spin`, either
in a settings object or a design's `quantum` object. Use `--image` for a trajectory
frame and `--guesses atom minao` to request a smaller ordered subset. The run
command returns 2 when any saved attempt failed, was interrupted, or was
abandoned; a successful bounded partial run returns 0 with pending guesses
explicit. The read-only report command returns 0 for internally consistent
evidence even when that evidence describes failed attempts. Interruptions return
130. Missing or contradictory evidence returns 2; nothing is repaired on report.

```python
from nanodesign.quantum import QuantumSettings
from nanodesign.state_scan import (
    create_state_scan, run_state_scan, state_scan_report,
)

# Input preparation only: this does not invoke the quantum solver.
plan = create_state_scan(
    "data/reference/methane_ethynyl_ts.xyz",
    "runs/methane-ts-guess-survey",
    QuantumSettings(charge=0, spin=1, xc="pbe0", basis="def2-svp"),
    image=-1,
    guesses=("minao", "atom", "1e", "huckel"),
)

# Inspection is read-only, including when another worker is running.
report = state_scan_report("runs/methane-ts-guess-survey")

# Calling this starts ONE real energy-plus-gradient calculation.
# Run only when compute capacity is available and this work is intended.
report = run_state_scan("runs/methane-ts-guess-survey", max_jobs=1)
```

The source path can be an ASE-readable single structure or trajectory. `image`
selects an integer frame, with -1 meaning the last frame. Single atoms are
supported. Supply `QuantumSettings` explicitly: charge and spin come from those
settings, not from coordinate-file metadata. Spin is N-alpha minus N-beta.
Nonperiodic, finite geometry, supported elements (H through Kr), sensible atom
separation, and electron/spin consistency are checked before starting a survey.

There must be one to four distinct guesses from `minao`, `atom`, `1e`, and
`huckel`, supplied in the desired order. That order is retained in reports.
Every new output directory must be absent. Creation captures the original file
bytes, selected geometry/frame and all settings. Editing or deleting the
original source later does not change the planned calculation. Changing any
saved input or the immutable plan invalidates the survey.

`run_state_scan` allows an integer `max_jobs` from 1 through 4 and evaluates
pending guesses serially. A failed guess consumes a job. Subsequent invocations
skip every attempted guess, including failures and interruptions. To repeat a
guess, create a new survey; this preserves previous evidence. There is no
automatic retry, method fallback, or resumption of SCF orbitals.

The solver computes both energy and analytical forces for each guess. Source
constraints do not mask these forces, and every atom stays at its frozen
position. A single point need not have small forces to be a valid single-point
evaluation. The job count is not a wall-time or memory limit; fixed `max_cycle`
and the solver's advisory memory setting still apply.
Each attempt records host load averages before/after its elapsed time; those
times are observations under concurrent load, not controlled speed benchmarks.

## Evidence and interpretation

| Artifact | Meaning |
|---|---|
| `plan.json` | Immutable guess order, canonical settings, exact source hash, frame and geometry hashes, units |
| `inputs/<source filename>` | Exact original bytes, including any compression and every source trajectory frame |
| `scan.json` | Plan digest and one attempt entry per requested guess |
| `attempts/<guess>/result.json` | Per-guess geometry/settings bindings, energy/raw forces, diagnostics or preserved failure |
| `attempts/<guess>/electronic.jsonl` | Production solver event log when the calculation reaches electronic work |
| `.scan.lock` | Retained advisory execution lock; do not remove during work |

The selected frame is identified by its index, exact source-file hash, and a
hash of its parsed coordinates and element order. Each attempt also hashes its
full settings, including that attempt's starting guess. Reporting verifies
every prior result with a stored checksum, not only the latest attempt.
These hashes establish internal consistency, not independent authentication of
the original computation.

Ordinary failures are saved and execution can continue to another pending
guess within the same explicit job limit. KeyboardInterrupt/SystemExit save an
interrupted attempt and propagate to the caller. After a process dies, the next
execution recovers an already-saved terminal result or marks an unfinished
attempt `abandoned`; it preserves partial bytes and never silently retries.
Read-only reporting does not perform that recovery. Concurrent executions of
the same survey are rejected by a macOS/Linux advisory lock.

Three separate flags prevent incomplete evidence from becoming a successful
state search:

- `all_requested_attempts_finished`: every requested guess has reached a
  terminal status, which may include failure, interruption or abandonment.
- `all_requested_guesses_converged`: all requested evaluations have finite
  energy/raw forces plus matching SCF and gradient success evidence.
- `complete_four_guess_check`: all four distinct supported guesses meet that
  convergence contract. Successful evaluation of a smaller subset is not a
  completed four-guess check.

`status="finished"` describes attempt bookkeeping only. Failed, interrupted,
abandoned and pending guesses remain explicitly listed. A numerically accepted
evaluation requires matching full settings and guess diagnostics, coherent
energy units, finite raw forces, and coherent determinant spin diagnostics.
All attempts use fresh calculators; no cached source-file energy or another
guess's calculation is reused.

`energy_spread` identifies its contributing guesses and describes only that
converged subset. It is null with fewer than two accepted evaluations, since
one result cannot establish agreement across guesses. eV and Hartree follow the
recorded solver conversion; kcal/mol uses exact SI constants and the
thermochemical calorie. The conversion is included in each report. Existing
historical benchmark values are not rewritten.

No energy grouping, lowest-guess selection, state identity, ground-state
verification, reaction barrier, or operating reliability is inferred. The
`ground_state_verified`, `electronic_state_identity_verified` and
`method_accuracy_validated` flags remain false. The per-calculator record still
describes one starting guess; survey completeness is a separate report property.

## Verification and handoff

```sh
.venv/bin/python -m pytest -q tests/test_state_scan.py
```

The production regression set passed 34 tests using synthetic solver
responses only. It covers bounded continuation and order, all-four versus
subset labels, preserved ordinary failures and interrupts, immutable input and
result checks, fresh calculations despite cached source energies, selected
trajectory frames and atomic inputs, invalid controls/geometries, and refusing
existing outputs. It also checks execution-lock exclusion, required completed
result checksums and strict JSON types in saved geometry/settings/frame bindings.
An independent S3 review passed 24 additional synthetic contract checks under
`research/state-scan-review/`. Both reported integrity findings were repaired;
the before/after verification receipts preserve that review history. These
checks validate software behavior only; no quantum calculations were launched.

S2 handed this module, its production tests and this README to C1 for integration
on 2026-09-16. C1 owns subsequent changes and any CLI/UI wiring.
Existing S1/A1/A2 scientific scans and
running jobs are unchanged; this module is a reusable option for subsequent
explicitly budgeted work, not a request to repeat them.
