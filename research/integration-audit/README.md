# Independent integration audit (Q1)

This support lane checks production input/campaign handling and the two archived
paired DFT/CC comparisons. It does not change quantum methods, launch a research
campaign, or certify the molecular design. C1 owns production fixes and integration.
Current ownership is in `coordination/ROSTER.md`.

C1 subsequently delegated the narrow input-byte consistency repair to Q1 in
`nanodesign/design.py` with `tests/test_design_snapshot.py`. The loader now
captures the design document and both coordinate files once, parses those bytes,
and hashes those same bytes. Text, binary, compressed and filename-only ASE
formats retain their format handling. The targeted design/campaign/workflow suite
passed **75 tests, with 1 deselected, in 5.05 s** after this repair. C1 repaired
the symlink-copy and historical-attempt integrity issues in `campaign.py`.
All **7 independent regression probes pass**, with zero expected failures, after
those repairs. The original observed failures remain documented in
`input_review.md`. Production path ownership returns to C1 for integration.

## Baseline

The existing suite selected by `.venv/bin/python -m pytest -q -m 'not quantum'`
passed **197 tests, with 7 deselected, in 10.35 s** on 2026-09-16. This selection
still includes some very small actual quantum tests whose functions have no
`quantum` marker; it is not a fully solver-free run. No new quantum campaign or
expensive calculation was started by Q1. This baseline does not cover every
failure mode, and is not scientific validation.

## Scope and handoffs

- Input/campaign correctness: `input_review.md` and
  `test_input_findings.py` contain the focused review and reproductions.
- Paired archive consistency: `audit_paired_evidence.py` independently checks
  recorded arithmetic and provenance with the Python standard library;
  `test_paired_evidence.py` exercises malformed or changed evidence.
- The preliminary workflow setup-failure concern was transferred to P1/b717
  when C1 clarified ownership. P1 confirmed that snapshot-write or calculator
  initialization failures can leave a finished attempt marked `running` and
  handed three regression cases to C1. The reproduction and fix belong to P1/C1,
  not this lane; see `research/path-contract-audit/test_path_contracts.py`.

## Paired archive result

Both archives passed. `paired_evidence_report.json` preserves the generated
machine-readable audit; no source evidence was rewritten. The older schema has
19 warnings for absent initial-guess/state fields; the newer `atom` archive has
none. These are distinct from consistency failures, of which there were zero.

| Archived run | Recomputed CC nominal TS-relative energy (kcal/mol) | Audit |
|---|---:|---|
| `paired-ccpvdz` (historical default-guess run) | -11.9238733722 | Consistent; some state/guess metadata unavailable |
| `paired-ccpvdz-atom` | +2.3986156144 | Consistent; recorded unverified-state flags preserved |

These quantities are not independently verified activation barriers. The older
archive itself does not record its initial guess; its label/default history does
not substitute for missing evidence. This audit imports no production chemistry
code and runs no quantum calculations. Five bounded test methods pass, including
temporary-copy corruption and missing-record cases.

Run from the repository root (Python 3.11+):

```sh
python3 research/integration-audit/audit_paired_evidence.py
python3 -m unittest discover -s research/integration-audit -p test_paired_evidence.py -v
```

The auditor writes JSON to stdout and exits nonzero on inconsistency. It can
also accept one or more archive directories as positional arguments. The pinned
historical unit-conversion convention is recorded in its output. Internal hash
agreement is not independent authentication of the original computations.

Audit results describe internal evidence consistency only. Nominal transition
structure relative energies remain fixed-geometry quantities. A missing field in
a historical schema is unavailable evidence, not proof that the corresponding
scientific check passed. Original calculation archives are read-only.

Files remain uncommitted for C1 integration, per the ownership roster. Shared
assignment-board edits and any other agent's implementation are outside this
handoff. See `coordination/status/codex-support-q1.md` for live status.
