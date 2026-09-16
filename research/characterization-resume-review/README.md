# H2 independent characterization resume review

Owner: `codex-9395`, assigned by C1. H1 owns the reviewed prototype in
`research/characterization-resume/`; this directory contains independent tests
and their recorded execution. No production code, source evidence or quantum
calculation was changed or run by this review.

**Result: 36 tests passed; no unresolved implementation defect found within
the tested interruption, reuse, identity and remaining-call scope.** Final run
started 2026-09-16 at 21:33:54 UTC and pytest reported 7.69 seconds. The wrapper
including interpreter startup took 10.09 seconds. See `test-output.txt` and
`verification.json` for the exact command, output and source fingerprints.
All five reviewed implementation/test files were unchanged during that run.
Warnings came from ASE's existing NumPy shape-assignment deprecation.

The reviewed H1 prototype SHA256 is
`84e2154c066c69277ff09d2287e7d7e7226bf8efb7191ccec2942fcad7144326`.
An internal read-only reviewer also checked the final durability, failure
lineage and accounting cases and reported no material unresolved finding.

## Reproduce

From the repository root with project dependencies installed:

```sh
.venv/bin/python -m pytest -q research/characterization-resume-review/test_resume_contract.py
```

`fixtures.py` defines an explicit analytic software fixture with three atoms,
one fixed atom and unequal masses. Six free Cartesian coordinates require
one baseline and twelve signed displacements. The fixture supplies controlled
interruptions, independent evaluation IDs and a nonzero raw force on the fixed
atom. It produces no electronic-state or molecular-chemistry evidence. Tests
replace H1's exact-backend type guard with this fixture; no real PySCF solver
is invoked. Saved checkpoints live only in pytest temporary directories.

## What was verified

- Keyboard interruption and ordinary calculation failure at baseline, first
  plus and first minus: only completed records are reused; failed diagnostics,
  failed call IDs and failure types remain visible outside reusable records.
- Exact new backend evaluation counts and independently constructed visited
  coordinate sequences, including chained failed resumes and zero-call replay
  of a completed checkpoint. Original record values, call IDs and attempt IDs
  survive reuse, including raw loads on the fixed atom.
- Source directory bytes and caller geometry/calculator cache are preserved.
  A locked source worker rejects reuse; a released-lock durable `running`
  checkpoint can be resumed without rewriting the original status.
- Fifteen changed-input cases reject before output creation or force work:
  positions, elements, masses, anchors, cell, charges, magnetic moments,
  quantum method settings, initial guess, source context, step, force tolerance,
  frequency tolerance, coordinate bound and recorded software versions.
- Ten corrupted-checkpoint cases reject before force work: truncated JSON,
  wrong digest, reversed displacement sign, duplicate slot, failed record,
  missing force row, Boolean force value, contradictory settings/call ID and
  false completion. Semantic corruptions have recomputed envelope checksums,
  so these checks exercise more than byte-integrity validation.
- A nonstationary baseline is saved before rejection. Reusing that same
  baseline performs zero new evaluations and stops at the same guard.
- A persistent failure of atomic checkpoint replacement after the first
  displacement preserves the previous durable baseline and stops acquisition.
  Resume correctly repeats the unsaved displacement and completes the remaining
  stencil. Temporary files are cleaned up. This is an injected write failure,
  not a physical power-loss or filesystem crash experiment.

## Remaining calls and failure costs

For a stationary baseline and six free coordinates, the complete stencil needs
13 successful force records. A failed call costs an attempt but is not reusable.

| Interrupted call | Durable successful records | Calls already attempted | New calls to complete after resume | Total attempts across both runs |
|---|---:|---:|---:|---:|
| Baseline | 0 | 1 | 13 | 14 |
| First plus | 1 | 2 | 12 | 14 |
| First minus | 2 | 3 | 11 | 14 |
| Already complete | 13 | 13 | 0 | 13 |

Two failed attempts followed by completion total 15 calls in the executed
chained fixture. Failure of durable saving can also require repeating an
already calculated but unsaved force. Exactly-once solver execution is not
claimed; accepted durable records determine reuse.

For `d` free coordinates and `k` accepted durable records, `1 + 2*d - k`
counts missing stencil slots. It describes new evaluations to finish only if
the baseline guard and all future calculations succeed. The rejected-baseline
test has 12 missing slots but unchanged resume permits **zero** new calls.
H1 accepted this distinction and updated its consumer documentation without
changing the reviewed prototype.

The existing C2 planner remains unchanged. A saved per-call timing proxy could
be multiplied by the conditionally executable missing count, but retries,
unsaved work and changed runtime conditions prevent that from being a fixed
cost or deadline. This review measured no quantum runtime or speedup.

## Scope limits and handoff

This review tests software recovery and provenance consistency. It does not
establish electronic-state continuity, force accuracy, transition-state
connectivity, experimentally validated chemistry or a molecular machine.
E2 owns independent saved-force algebra and N1 owns mode/subspace comparisons.
H1 remains a research prototype; production promotion changes implementation
fingerprints and needs the integration owner's corresponding review.

H2 files are handed to C1 uncommitted: this report, `fixtures.py`,
`test_resume_contract.py`, `test-output.txt` and `verification.json`.

## Focused follow-up after snapshot metadata repair

On 2026-09-16 at H1/C1's request, H2 ran only the existing complete-resume and
chained partial-resume cases against fresh checkpoints after H1 separated
the emitted snapshot hash from the caller's original source hash and added
result units. **Both tests passed in 1.32 seconds**, retaining exact call counts
and original-source bytes. The reviewed prototype's new SHA256 is
`7429161b9be1bfd2c512a1e2580de1bb126cfe199f5bb405935e61dd61c32187`.
All five implementation/test hashes were unchanged during this focused run.

`metadata-smoke-output.txt` and `metadata-smoke-verification.json` record this
new check separately. The earlier 36-test output and receipt remain unchanged
and apply to their recorded earlier prototype hash. H2 did not repeat the
full suite or expand scope; H1 owns the direct metadata-field regression tests.
The two new receipt files and this addendum are released to C1 for integration.
