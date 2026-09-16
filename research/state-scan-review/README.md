# S3 independent state-scan review

Owner: support task `01a0ac0c-8cfb-7e81-8b3f-62050673600c` (`codex-8cfb`).
C1 assigned this directory and the owner's status only. S2 task
`01a0ac06-4280-7bc0-862a-35bfc19e958b` owns the production module and tests.
This review starts no quantum calculations and changes no source evidence.

## Review contract

S2's proposed standalone-geometry interface is:

- `create_state_scan(structure_path, output, settings, image=-1, guesses=...)`:
  preserve exact source bytes, selected frame, and supplied quantum settings.
- `run_state_scan(directory, max_jobs=1)`: run an explicit serial bound of
  pending guesses, each with a fresh calculator and identical fixed geometry.
- `state_scan_report(directory)`: inspect saved input/attempt evidence without
  starting calculations or changing calculation status.

The proposed initial scope has up to four distinct supported guesses and no
retry within a scan. A new scan is an explicit repeat, preserving the old one.
The implementation now uses this interface. The review below uses an independent
synthetic ASE calculator through the public API; no production files are edited.

The independent review will target:

1. Source edits after planning, selected trajectory frames, caller settings,
   and per-guess isolation. Only the declared starting guess may vary.
2. Failure and interruption evidence, including failure before energy exists
   and an exception after earlier guesses have completed.
3. Pending versus attempted versus converged subsets; fewer than two accepted
   energies cannot establish a spread. A failed guess is not an equal energy.
4. Tampered snapshots/old attempts, repeated run calls, and read-only reporting.
5. Scientific claims: a lower observed energy is neither a selected production
   guess nor proof of ground-state or electronic-state identity.

Tests will use synthetic calculator outputs in temporary directories, clearly
separate from scientific results. No research scripts will be imported because
some current scripts execute calculations at import time. Concrete findings,
source hashes, executed commands and final outcomes are recorded below.

## Final outcome

**All 24 independent checks pass** after S2 repaired the two findings below.
The final [verification receipt](verification-after-repair.json) records stable
source hashes throughout the run, exit code 0, and 1.82 seconds of test time.
The checked production module SHA256 is
`fea9e2fbf09a12180638b75d07f824e5bd783f3f5bffd6126758523034567c24`.
Warnings are existing ASE/NumPy deprecations. No quantum calculations ran.

Completed attempts now require a valid recorded digest; present invalid digests
are rejected. Canonical JSON binding checks distinguish boolean/numeric aliases
in the frozen frame, settings and geometry. Recovery of legitimately unfinished
attempts still preserves the original bytes and does not repeat a guess.
No open S3 blocker remains within this bounded review.

## Reproduced findings, now repaired

The first complete independent run had **18 passing checks and six failing
regression cases**. All six failures belonged to two concrete findings reported
directly to S2 and C1. S2 implemented both repairs; S3 verified them independently.

1. **Completed attempts need a recorded result digest.** Removing the ledger's
   `result_sha256`, setting it to null, or setting it to an empty string still
   permits the reader to accept the completed result. A completed result should
   require a valid digest. Running attempts must retain their separate crash
   recovery behavior, since a process can stop before writing a terminal digest.
2. **Semantic input binding cannot rely on Python equality alone.** A saved
   frame `True` compares equal to planned frame `1`; a saved charge `False`
   compares equal to `0`; a coordinate `False` compares equal to `0.0`. The
   reader accepts these malformed record fields alongside unchanged declared
   settings/geometry hashes. These probes refresh only the outer result file
   digest to isolate semantic checks; they are synthetic writer-corruption
   cases, not claims that signed or authentic data were defeated. Validate the
   actual binding types/canonical values and their hashes.

[Before-repair receipt](verification-before-repair.json) captures all 24 cases,
the command, output, return code, and source hashes before and after execution.
The positive checks include binary trajectory capture, explicit settings,
fresh calculators, unmasked fixed-atom forces, backend geometry mutation,
failure budgets, interruption, older-attempt tampering, subset interpretation,
and preservation of completed/partial/missing crash evidence during recovery.

The initial test probe omitted a required Hartree-to-eV diagnostic; that fixture
was corrected before the recorded verification. It was not a production defect.

## Repeat the review

From the repository root, using a new output file:

```sh
.venv/bin/python research/state-scan-review/verify.py \
  --output research/state-scan-review/verification-next.json
```

The runner executes only this directory's mock tests, with a 60-second process
timeout. It records source hashes and flags concurrent changes. Its exit code
is the test return code, or 2 when inspected files changed during execution.
The destination must be new; earlier verification records remain unchanged.

This tests software provenance and bookkeeping. It does not establish quantum
accuracy, electronic-state identity, the ground state, or a functioning tool.
