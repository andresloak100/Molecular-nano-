# DFT initial-guess compatibility review

Owner: D2 / `codex-c5cd`. Implementation belongs to D1 / `codex-9395`;
the campaign compatibility change belongs to C1. This directory contains
independent regression checks and a read-only review, not molecular results.

## Finding and impact

Adding `QuantumSettings.scf_initial_guess = "minao"` preserves the previous
solver default but extends `asdict(settings)`. Before the corresponding reader
fix, `campaign._result` compared that expanded dictionary directly with saved
`quantum_settings`. Historical records lack the new key. Consequently completed,
failed and nonconverged attempts became unreadable, and even resuming a fully
completed campaign failed before it could recognize that no work remained.
Pending campaigns still loaded because they had no result settings to compare.

`reproduction-unfixed.json` records this failure against the actual added field,
without simulation. It includes exact source hashes and shows that every fixture
file remained unchanged. `reproduction-before.json` records the earlier
process-local serialization simulation; use the unsimulated record as the direct
reproduction. None of these checks launches a quantum calculation.

## Required compatibility contract

For comparison in memory only, an absent historical `scf_initial_guess` means
`"minao"`. All other saved fields must still match strictly. In particular:

- Missing guesses must not match an explicitly requested `"atom"`.
- An explicitly invalid guess (including null) must not become a default.
- Missing basis, spin, functional, threads or density-fitting controls remain
  errors, as do unknown saved fields.
- Original result, design and plan bytes and their existing hashes remain intact.
- Hash verification must still reject any later modification to a saved result.
- Reading or resuming completed evidence must not launch another calculation.
- Readable numerical evidence remains unvalidated scientifically.

## Reproduce

From the repository root, with the existing environment:

```sh
.venv/bin/python research/dft-guess-review/reproduce_legacy.py
.venv/bin/python -m pytest research/dft-guess-review/test_campaign_compatibility.py -q
```

The optional `--simulate-added-default` switch is only for reproducing the
serialization extension on code that predates the new field. It patches
serialization in the isolated process; it never edits the source files.

Fixtures are synthetic software documents constructed in temporary directories.
They are explicitly labeled, and their energies and forces carry no scientific
meaning. The tests exercise the real campaign reader/resumption and file-integrity
checks with the solver replaced by a guard that fails if called. They cover
backward compatibility, not DFT accuracy or state identity.

## Handoff

The feature-only code produced 5 failures and 15 passing checks in the first
20-case run. The suite then added five explicit invalid-type cases (25 total).
After C1's narrow reader fix, **all 25 checks pass**. The independent four-case
reader probe also passes and preserves every fixture file. `verification.json`
records the executed command, source hashes and output. The internal provenance
review additionally ran 70 existing mocked/analytical tests, all passing; see
`provenance-review.md` for scope and the optional paired-reader hardening finding.
General archived-evidence arithmetic, reaction-path physics and the visual
workbench are owned by other agents and are outside this review.
