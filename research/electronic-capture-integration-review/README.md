# Electronic capture: production lifecycle review

Owner: support `codex-8cfb`, assigned by C1 on 2026-09-16. This directory is a
bounded integration review. C1 owns `nanodesign/quantum.py`; E3/f522 owns the
snapshot writer/comparator and E4/712e independently reviews its mathematics.
No electronic calculations or production edits are performed here.

## Scope and decision

The original proposal captured only a completed energy-and-force evaluation.
C1 subsequently proposed retaining a converged SCF reference when its gradient
fails. **That is useful evidence if its narrower scope is explicit.** An orbital
snapshot is not an accepted gradient, a completed calculation, an identified
physical electronic state, or a ground-state certificate.

Recommended hook contract:

1. Start each actual `calculate` attempt with empty results and fresh diagnostics
   before fallible ASE setup. Preserve the immutable capture option separately
   from physical `QuantumSettings`.
2. Use the current, local final RKS/UKS mean-field object, including any density-
   fitting wrapper. Capture under the existing PySCF lock only after convergence,
   finite SCF energy and finite spin diagnostics. Do not make another SCF call.
3. Save a detached, validated snapshot to a new, immutable per-call destination.
   Bind the actual molecule, resolved AO representation, settings and `call_id`.
   E3's writer owns the byte/dimension bounds and safe persistence contract.
4. Store the reference, checksum, call ID and SCF-only phase in `diagnostics`,
   initially with whole-evaluation acceptance false. The file must not claim
   completed forces merely because the SCF calculation converged.
5. Continue gradient, dispersion, finite-value and force-shape checks. Preserve
   a captured SCF-only reference if any later stage fails; clear usable energy
   and force results. Keep the actual failure cause and stage identifiable.
6. Publish usable results and whole-evaluation acceptance only after completion
   logging and thread restoration succeed. Logging or restoration failures must
   not leave a successful current-call acceptance flag beside empty results.
7. Cache reads reuse the same call and reference. They must not generate another
   snapshot or implicitly recalculate to obtain previously disabled capture.

The root-owned hook now implements SCF-only retention with an immutable
`electronic_state_directory` constructor option outside `QuantumSettings`.
Export failure is immediate and fail-closed: no gradient is started after a
requested export fails. Its diagnostics distinguish `capture_error` from the
SCF/gradient flags. Full acceptance is set after thread restoration and the
completion event.

## Concrete lifecycle risks found

**Hessian cache restoration.** `stationary.py` saves and restores calculator
`atoms`, `results`, and a deep copy of `diagnostics`. An additional standalone
`calculator.electronic_snapshot` attribute would retain the last displaced
geometry's orbitals after the baseline energy/forces are restored. An immutable
reference within diagnostics participates in the existing restoration contract.
All displacement artifacts may remain on disk as separately identified evidence.

**Setup outside cleanup.** At the reviewed baseline, `super().calculate` runs
before fresh diagnostics and outside the failure boundary. ASE directory setup
can fail there. If capture is only reset later, a previous call's reference can
remain exposed as though it describes the new failed attempt.

**Interrupt and completion order.** The reviewed handler catches `Exception`,
not `KeyboardInterrupt` or `SystemExit`, and assigns results before the completion
event. Capture integration should clear usable results on all failed/interrupted
paths while re-raising the original interruption. A retained SCF-only reference
must still show that the whole evaluation was not accepted.

**Primary failure precedence.** Current ordinary scientific errors are preserved
when failure logging raises `OSError`. Extend this principle to the new capture
path: a secondary log/cleanup failure must not replace an already observed
scientific or interruption cause. A restoration error from a `finally` block can
otherwise replace the original exception. Record secondary failures separately.

**Export errors are different from physics errors.** If explicitly requested
capture fails before the gradient begins, immediate fail-closed behavior is a
small, coherent policy, provided diagnostics identify capture/export failure.
Alternatively root may continue the scientific evaluation and report the export
failure separately; then a later scientific failure takes precedence. Do not
invent gradient completion or reuse a prior snapshot under either policy.

## Independent synthetic verification plan

Once root supplies the final hook, replace PySCF solver factories and force
evaluation with small deterministic mocks; never run an SCF or gradient solver.
Exercise the real production calculator and, where appropriate, its real E3
serializer against those mocks.

| Case | Required observable |
|---|---|
| Capture disabled | No capture-specific work/import requirement or artifact |
| RKS/UKS, direct/density fitted | Current final solver, explicit settings, geometry and call ID reach capture |
| Energy then forces, unchanged atoms | One solve, one snapshot, one call ID |
| Geometry change or explicit recalculation | New call/reference; previous artifact preserved |
| Failed/nonfinite SCF or spin | No accepted SCF reference or force result |
| Failed gradient/dispersion or invalid forces | SCF-only artifact retained, whole-call acceptance false, no results |
| Capture writer failure | Declared export policy, no stale prior reference, exact failure retained |
| Completion/failure-log failures | No falsely accepted call; primary failure survives |
| Interruption and thread cleanup | Original interruption propagates; no usable results or false acceptance |
| ASE setup failure after prior success | New failure cannot expose previous capture as its own |
| Hessian characterization success/failure | Baseline cache and its reference restored together; displaced references stay distinct |

These are software lifecycle tests. They cannot validate orbital continuity,
chemical accuracy, or the physical assembly operation.

## Verified lifecycle subset

Run from the repository root:

```sh
.venv/bin/python -m pytest -q research/electronic-capture-integration-review/test_capture_lifecycle.py
```

**10 checks passed in 1.20 seconds** on 2026-09-16. These deliberately complement
the integration owner's broader `tests/test_quantum_capture.py` coverage:

- The real characterization workflow restores its baseline capture reference,
  diagnostics, atom cache, energy and forces after either a completed Hessian or
  a displaced-gradient failure. Every displaced snapshot remains identified by
  its own call ID and different geometry; the baseline artifact is unchanged.
- A gradient `RuntimeError`, `KeyboardInterrupt` or `SystemExit` retains its
  identity when thread restoration also fails, including the additional case
  where failure logging raises an I/O error. Secondary errors are recorded,
  energy/forces are cleared, and the SCF-only artifact remains unaccepted for the
  whole evaluation.

- Two added tests exercise the **actual E3 capture, writer and reader through
  the production hook**, using synthetic solver arrays and a synthetic AO metric.
  Completed RKS and gradient-failed UKS cases retain the exact call ID, solver
  coordinates, full settings, byte checksum and SCF-only scope. Saved evidence
  stays detached after the solver's coefficient buffer is changed.

The initial eight cases replace both PySCF and the E3 serialization boundary;
their saved files explicitly say `synthetic_capture_boundary_only`. The two
additional cases use the real E3 schema and persistence implementation, with
producer version marked `synthetic-solver-test-only`. Only temporary directories
receive these files. No SCF, gradient solver or AO integral is evaluated, and
these checks do not independently establish the orbital comparison mathematics
or physical state identity. Source hashes observed after the final focused run:

- `quantum.py`: `ae92c852ac9989542b3c0cdbe2861a4d2ffe18ef426f1aa1b3d9869ab2489ee4`
- `electronic_state.py`: `f7e277e66b334807a6ba5641559271afbeeac5e4c21c436c80e3c59dc6dbcb3d`
- `stationary.py`: `efa54d1ceb79f7ca2c6e166f89047e002ec8d1bb7904c0e400883a19a5d85a82`
- `test_capture_lifecycle.py`: `0b65a219ccc6278c07b1e7e9046bd86610f54062633b181d4da08d89afc84510`

The initial two roundtrip failures were caused by the test loader omitting its
package context; correcting that fixture fixed them without production changes.
Root and E3 received the final outcome and precise scope. The bounded review is
ready for integration, with no open finding. No production implementation was
edited by this reviewer. The machine-readable receipt is `verification.json`.
