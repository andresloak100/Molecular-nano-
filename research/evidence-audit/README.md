# Independent saved-evidence audit (E1)

The completed audit passes four groups of consistency checks on the existing
SVP/TZVP reference, H2 vibration, and 53-atom direct/density-fitting archives.
This verifies recorded arithmetic and provenance links, not chemical accuracy,
electronic-state identity, an actual reaction barrier, or an operating machine.

Run from the repository root:

```sh
.venv/bin/python research/evidence-audit/audit_evidence.py
.venv/bin/python -m pytest -q research/evidence-audit/test_audit_evidence.py
```

The checker uses NumPy only for small archived matrices. It imports neither
`nanodesign` nor a quantum solver and launches no calculations. The existing
`data/validation/verify_density_fitting_comparison.py` is reused for its
direct/DF arithmetic checks. Run without Python `-O`, which would disable that
existing checker's assertions.

`--validation-root PATH` audits a separate archive copy. Reads reject paths or
symlinks escaping that root. `--output NEW_FILE.json` optionally saves a new
report; it refuses to overwrite an existing file. Any failed group makes the
overall status `failed` and returns a nonzero exit code. The report includes
hashes of inspected files, the audit script, and the reused DF checker.

## Saved outcome

The deliverable is [audit-2026-09-16-final.json](audit-2026-09-16-final.json).
An earlier passing checkpoint is retained separately as
`audit-2026-09-16.json`; its implementation predates the additional diagnostic
and event consistency checks. Use the final report for this handoff.

| Group | What was checked | Result |
|---|---|---|
| PBE0/def2-SVP and def2-TZVP | 16 manifest checksums; ten species records against summaries; SI atom order and bohr-to-angstrom conversion; state/settings, finite energies/forces, energy units, spin arithmetic and relative energies | Pass |
| H2, both displacement steps | Original input hashes and relaxed coordinates; Hessian eigenvalues; mass-normalized mode orthogonality and eigenvector residuals; signed frequencies; classification, forces, completion events and step summary | Pass |
| 53-atom direct and DF | Original input hashes, 22 C + 31 H and six anchors; interpretation checksums; documented reproduction-input changes; existing independent energy/force comparison | Pass |

The H2 stretch change is **−0.0895837007 cm⁻¹** between 0.003 and 0.0015 Å
displacements. Both spectra retain five unresolved modes. The largest saved
mass-weighted mode-equation residual is **2.13 × 10⁻¹⁴** in the numerical
matrix units. These are algebraic and step-sensitivity checks on H2 only.

All **22 focused tests pass**. They modify temporary copies to establish that
the checker rejects altered coordinates, energies, Hessians, modes, frequencies,
step summaries, false convergence flags, missing/duplicated event evidence,
unsupported validation labels, nonfinite values, missing files and escaping
paths. One test also confirms failure exit status and report overwrite refusal.
No source calculation record is changed by the tests.

## Findings and limits communicated to the owners

1. **Missing raw force evidence in historical H2 archives.** Each mode log
   documents thirteen completed SCF-plus-gradient calls, but individual displaced
   coordinates/forces and their call-ID mapping were not saved. The Hessian's
   derivation and its reported pre-symmetrization asymmetry therefore cannot be
   reconstructed from these artifacts. Passing this audit checks the saved
   Hessian and downstream algebra. C1 assigned prospective force-evidence
   recording to the P1 stationary-workflow owner; historical archives remain
   unchanged. This omission is not evidence that the Hessian is wrong.
2. **Small historical unit-conversion difference.** Benchmark exports use
   0.0433641153087705 eV per kcal/mol; current exact SI definitions give
   0.043364104241800934. The difference is about 0.255 ppm, changing the two
   nominal reference relative energies by less than 0.000001 kcal/mol. Both
   historical reproduction and current-SI recalculation are shown in the report.
   This has no material bearing on the recorded scientific discrepancies.
3. **Historical omissions remain explicit.** The direct run predates event and
   effective-thread instrumentation. Its original design JSON bytes are not in
   the portable archive, so the recorded original design checksum was not
   reconstructed from its embedded object. Portable reference manifests omit
   the original machine-specific `source_directory`; their retained original
   manifest hashes were likewise not reconstructed. Available file hashes and
   documented reproduction inputs were checked instead.
4. **Integrity is not authenticity.** Agreement between a file and a saved
   checksum does not independently certify who produced it or the physical
   correctness of the computation. Raw finite-difference force evidence is
   absent in these archives, and no experimental evidence was introduced.

## Coordination and ownership

Owner: `codex-support-a551`, E1. Files in this directory and
`coordination/status/codex-support-a551.md` only. C1 owns integration and core
changes. Q1 (`research/integration-audit/`) independently audits the paired
CC archives; this directory deliberately excludes them. DFT initial-guess
implementation belongs to 9395, its compatibility review to c5cd, and path/
stationary changes to b717/P1. No S1/A1/A2/V1 files or ongoing jobs were modified.

Ready for C1 review and integration; no standalone quantum run or cloud work
is required to use the checker.

Related handoffs received: Q1's
[`paired_evidence_report.json`](../integration-audit/paired_evidence_report.json)
reports both paired archives consistent, with historical metadata omissions
explicitly marked. P1 reports that future displacement evidence now includes
units, atom/axis/offset order and optional diagnostic call IDs. P1 subsequently
reported a passing JSON-roundtrip reconstruction test covering both the raw
asymmetry and symmetric Hessian, with 52 combined tests passing; its files are
stable awaiting C1 integration. E1 did not independently rerun that separate
suite. Those changes are owned and validated in their respective lanes.
