# R1: quickstart reproducibility handoff

Owner: support task `01a0ac0a-f522-7733-883d-65746dfce0e3` (`codex-f522`).
C1 assigned this lane directly. Only this directory and the owner's status file
were edited. Root owns production/documentation repairs and integration.

## Executed result

**34 checks passed; zero failed:** 26 actual CLI subprocess invocations and eight
artifact/output checks. Executed 2026-09-16 21:09:36–21:10:37 UTC in temporary
working directories outside the checkout, using the existing installed Python
3.14.7 environment. No quantum jobs were launched and no solver result/event
files appeared in the temporary work. The checked source hashes were unchanged
between the start and finish of this audit.

The raw commands, exit codes, stdout/stderr, environment package versions and
source hashes are preserved in [results.json](evidence/attempt-01/results.json).
This checks user-facing command behavior, not chemical accuracy.

| User flow | Observation |
|---|---|
| Installed command and `python -m nanodesign` | Help works outside the checkout. |
| All ten subcommands | Help exits successfully without launching calculations. |
| README candidate → check | Three expected files; 53 atoms, C22H31, 163 electrons, spin 2S=1, six anchors; `design_validated: false`. |
| Bundled candidate | Loads and checks successfully by absolute path. |
| README generated pose campaign | Nine pending poses, documented plan/state/snapshot files; no calculations. |
| Campaign from an existing design | One pending pose. |
| Copied bundled campaign | Reports nine pending poses after relocation. |
| Archived-result audit | Reports missing evidence and keeps `design_validated: false`; numerical archive arithmetic was not re-audited. |
| Missing command/files | Expected nonzero exit and useful error output. |
| Reusing candidate/campaign directory | Refused with a nonzero exit. |

To repeat from the repository root, use a **new** evidence destination:

```bash
.venv/bin/python research/quickstart-audit/audit_quickstart.py \
  --out research/quickstart-audit/evidence/attempt-02
```

The script accepts `--executable /absolute/path/to/nanodesign` when needed. It
allows help for calculation commands but permits execution only of candidate,
check, audit, campaign-create and campaign-report. Every command has a 30-second
timeout. It leaves the source checkout and scientific archives untouched.

## Documentation repairs completed by C1

**Closed after inspection at 2026-09-16 21:12 UTC:** C1 repaired both findings
in the root README. The characterization command now consistently uses the
generated design (current line 74), and the output/log description is scoped by
command (current line 81). The original observations are retained below for
traceability. The CLI probes were unchanged and did not need to be repeated for
these documentation repairs.

These two findings came from source/documentation inspection; no expensive
calculation was run to demonstrate them. References below refer to the source
hashes captured in attempt-01 and may move when C1 edits the README.

1. **Keep the generated design throughout the example.** README line 63 invokes
   `characterize examples/h-abstraction/design.json`, while the preceding
   candidate and relaxation commands use `designs/h-abstraction/design.json`.
   The default inputs agree today. A user changing the generated design's
   method or anchors would not carry that choice forward into characterization.
   `run_characterization` loads settings and constraints from its supplied
   design, so this is a provenance pitfall rather than a current default-path
   failure. Repair: use `characterize designs/h-abstraction/design.json`.

2. **Name command-specific output records.** README line 70 says post-start
   failures are in `result.json` and events in `electronic.jsonl`, immediately
   after examples for several command types. That is correct for `calculate`
   and `characterize`; `benchmark.py` writes `benchmark.json` and
   `electronic/<species>.jsonl`, and `method_comparison.py` writes
   `method_comparison.json` with DFT evidence under `dft/`. Suggested wording:
   “For `calculate` and `characterize`, inspect `result.json` and
   `electronic.jsonl`. A benchmark records its summary in `benchmark.json`
   and per-species events in `electronic/`. A paired comparison records
   `method_comparison.json`, DFT evidence under `dft/`, and coupled-cluster
   records under `ccsd_t/`.”

The documented relaxation output `structure.extxyz` was checked against the
writer and is correct; no repair is needed there.

## Limits and coordination

- This used an already installed environment. A fresh dependency install,
  another platform/Python version, performance, actual calculations, and browser
  behavior were not tested.
- Core input/campaign integrity belongs to Q1, path contracts to P1, saved
  arithmetic to E1/Q1, DFT guess compatibility to D2, and workbench implementation
  to V1. This audit intentionally stays at the user-entrypoint boundary.
- The full test suite was not run. `pytest -m 'not quantum'` is not an adequate
  no-solver guard here because some real-calculation tests lack that marker.
- Initial lane acknowledgement, preliminary findings, completed check counts and
  confirmed closure of both repairs were sent directly to C1. The internal `offline_inventory` helper reviewed
  documentation read-only and launched no computations.
- Leave these files uncommitted for C1's explicit path-specific integration.
