# SCF-solution provenance audit of committed evidence

Run: 2026-09-16 ~21:45 UTC, support session 76190bf3, at commit `6a43ba8`
plus uncommitted working-tree docs. Tool: `audit_scf_provenance.py` (read-only,
stdlib, exits nonzero on findings). Scope: every git-tracked `*.json`.
Requested by the scientific-helper session after its ethynyl SCF-trap
discovery. This report records findings; no evidence files were edited, since
`data/validation/` is C1-owned.

## Verdicts

| File | Verdict |
|---|---|
| `data/validation/paired-ccpvdz-atom/` (all 6 files) | OK_LABELED, correct solutions, explicit `scf_initial_guess: atom` throughout |
| `data/validation/si-energy-reproduction/si-energy-reproduction.json` | OK_SCANNED, four-guess scan recorded per species, lowest solution selected |
| `data/validation/paired-ccpvdz/ccsd_t/ethynyl_radical.json` | **UNLABELED + ARTIFACT SOLUTION** (HF S^2 = 0.7591) |
| `data/validation/paired-ccpvdz/method_comparison.json` | **UNLABELED + MIXED SOLUTIONS** |
| `data/validation/paired-ccpvdz/ccsd_t/methane_ethynyl_ts.json` | UNLABELED (but on the lowest known solution, S^2 = 1.2143) |
| `data/validation/paired-ccpvdz/ccsd_t/methyl_radical.json` | UNLABELED (S^2 = 0.7619, no known competing solution) |
| All other tracked JSON | no open-shell CC content, out of scope |

## The corrupted quantities

`paired-ccpvdz/method_comparison.json` combines the ethynyl radical on the
**high/artifact UHF solution** (S^2 = 0.7591) with the transition structure on
the **lowest known solution** (S^2 = 1.2143). Any quantity involving ethynyl
in that file is therefore a cross-solution artifact:

| Quantity | minao arm (artifact) | atom arm (correct) |
|---|---|---|
| CCSD(T) reaction energy | -39.11 kcal/mol | -24.79 kcal/mol |
| DFT minus CCSD(T) reaction energy | +12.31 kcal/mol | -2.01 kcal/mol |
| ethynyl DFT-minus-CC gap | -3.58 eV | -2.96 eV |

The per-species DFT energies and the closed-shell CC energies in the minao arm
are unaffected (DFT is guess-independent here, checked previously). The TS,
methyl, methane and acetylene CC energies match the atom arm and are fine; the
defect is confined to ethynyl-derived quantities plus the missing provenance.

## Why the file cannot currently defend itself

The minao arm records **no** `scf_initial_guess` field anywhere (the atom arm
records it three times per species). A reader has no way to tell which
solution produced the numbers without recomputing S^2 fingerprints. That is
the auditable-provenance gap, distinct from the wrong-solution defect.

## Recommendations (for C1, which owns the archive)

1. Keep the minao arm — the pair is valuable precisely because it documents
   the trap — but annotate it: explicit `scf_initial_guess: "minao"` (or a
   correction note) and a warning that ethynyl-derived quantities are
   artifact-solution numbers, with a pointer to the atom arm.
2. Make `scf_initial_guess` a required recorded field for any open-shell
   CC archive going forward. The `CCSettings.scf_initial_guess` work already
   underway covers the API side; this covers the evidence side.
3. Optional: run `audit_scf_provenance.py` in the test suite so a regression
   fails loudly. Exit code is already CI-shaped. Fingerprint table must be
   extended when new open-shell species (adamantyl, t-butyl) enter evidence.

## Limitations

- Fingerprints cover only ethynyl and the methane-ethynyl transition
  structure at the published cc-pVDZ geometries, tolerance 0.005.
- The audit reads committed files only; uncommitted evidence (e.g. the
  running isobutane reproduction) is not judged here, though its incremental
  file already records full four-guess scans and correct selections.
- Absence of a fingerprint match is not proof the lowest solution was found;
  it only means no *known* wrong solution was matched.
