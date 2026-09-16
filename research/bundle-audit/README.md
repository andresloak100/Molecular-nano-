# B2 independent evidence-bundle review

B2's final check **passes**: 11 independent tests plus a byte-complete synthetic
roundtrip of nine files / 630 bytes. Source hashes were unchanged during the run.
The exact command, output, reviewed source hashes and per-file roundtrip digests
are preserved in `verification-01.json`.

Reviewed B1 module SHA-256:
`0627420167d4528935a2f35cca570176792675392b68bbf9b0876a1e8ee7f9c1`.
No new material implementation finding remains from this bounded review.

## Ownership and scope

Assigned by C1 to task `01a0ac0c-712e-71a0-8160-cab84676bacf`.
B1 (`01a0ac0a-f522-7733-883d-65746dfce0e3`) owns `nanodesign/bundle.py`,
its production tests and documentation. B2 owns this directory and its own
coordination status only. B1 and B2 communicated directly about the manifest,
path boundaries, publication failure semantics and final frozen source hash.
No B2 production edits or commits were made.

All fixtures are synthetic files created in temporary directories. No quantum
job, real-archive copy, network request or source-evidence modification occurred.

## Reproduce

From the repository root:

```sh
.venv/bin/python research/bundle-audit/run_review.py
```

To preserve another receipt, supply a new path:

```sh
.venv/bin/python research/bundle-audit/run_review.py --output research/bundle-audit/verification-next.json
```

The runner refuses to overwrite an earlier receipt. It checks that reviewed
source files stayed unchanged during the run. Its imports and checks use only
Python's standard library; it does not call a solver.

## Executed checks

- Exact byte-preserving, read-only roundtrip, including failed and pending
  attempts, empty/hidden files, UTF-8 names/text, opaque binary data and explicitly
  synthetic XYZ coordinates. The moved bundle verifies independently of the
  original source location.
- Output targeting the source through a symbolic-link alias is rejected without
  adding output inside the source.
- Same-size source changes after capture and selected source replacement by
  either another directory or a symlink are rejected before publication.
- Interrupted payload writes, interrupted pending-manifest writes and failed
  final-manifest linking do not produce a verifiable partial bundle.
- A failure removing the pending manifest after publication leaves an extra
  entry that verification rejects. This is explicitly distinguished from
  failures before publication; no unconditional cleanup is expected.
- Payload changes after their initial verification read, manifest replacement
  between read and inventory, and duplicate JSON keys are detected.
- The additional receipt roundtrip saves every fixture's original/final length
  and digest, confirms source bytes remain unchanged, and retains
  `scientific_validation="not_assessed"`.

B1 separately owns the ordinary production test suite, including count/size
limits, missing/extra entries, invalid metadata, ordinary traversal, symlinks and
nonregular files. B2 does not claim those as additional independent tests here.

## Fixes and limitations

During coordination, B1 itself identified and repaired destination-parent
pinning and added a recheck of written payload bytes/inventory. The final B2
run used that patched module. These are not represented as new B2 findings.

The manifest is an integrity record, not source authentication, an atomic
cross-file snapshot, scientific validation or proof that all external referenced
files are included. The fixtures establish software behavior for the tested
cases only. Actual chemical accuracy and manufacturing feasibility remain
outside this review.

Review deliverables: `fixtures.py`, `test_verification.py`, `test_mutation.py`,
`run_review.py`, this README, and `verification-01.json`.
