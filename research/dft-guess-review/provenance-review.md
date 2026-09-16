# D2: DFT initial-guess provenance integration review

Reviewed 2026-09-16 by internal reviewer `codex-c5cd-provenance`, for parent
`codex-c5cd`. Ownership: this report and the reviewer's own coordination status
only. No core/research-owner files edited and no quantum calculations run.
The checkout was live: `QuantumSettings.scf_initial_guess` was added by its owner
while this review was running. Final observations include that change atop
HEAD `67ab3bdbf06895e5f5c62e9eed785136879674aa`.

## Outcome

New workflow, characterization and benchmark records carry the effective DFT
guess automatically through the existing complete settings serialization. The
new calculator diagnostics also explicitly record one guess and leave scan,
ground-state and electronic-state verification false. No writer migration is
needed in the reviewed core paths.

One concrete validation gap is exposed by the extra setting: the paired-method
DFT evidence gate does not validate it, and accepts inconsistent or invalid
explicit values. This is a pre-existing partial-settings validation issue made
more consequential by the new electronic-solution control; it is not evidence
that production has used the wrong guess. The paired workflow always creates fresh DFT calculations with its default
guess, so the demonstrated altered-record cases are hardening suggestions, not
historical compatibility breaks or observed production errors. There is no
benchmark or paired-method resume loader whose old dictionary equality is newly
broken. Campaign
compatibility is reviewed separately by the parent.

## Related hardening suggestion, reproduced without quantum work

At `nanodesign/method_comparison.py:82`, `_paired_input` reads each DFT species'
`quantum_settings`. Lines 83-88 check basis, XC, dispersion, nominal charge/spin
and completion flags, but do not compare `scf_initial_guess` with the expected
DFT setting or with `quantum_diagnostics.settings`. The comparison writer records
`dft_settings` at line 128, which could therefore disagree with accepted species
metadata if a backend or future archive reader supplies conflicting evidence.
The strict CC settings equality at line 179 does not protect the DFT side.

A read-only probe loaded
`data/validation/paired-ccpvdz-atom/dft/benchmark.json`, deep-copied its dictionaries,
and called `_paired_input(copy, Path('data/validation/paired-ccpvdz-atom'),
'methane', 'cc-pvdz')`. Each of the following passed:

| Variant | Top-level/diagnostic guess | Species guess | Accepted |
|---|---|---|---|
| Historical record unchanged | absent | absent | yes |
| Explicit conflict | `minao` | `atom` | yes |
| Explicit null | `minao` | `null` | yes |
| Invalid explicit value | `minao` | `invented` | yes |

Minimal reproduction (only reads archived files):

```python
import copy
import json
from pathlib import Path
from nanodesign.method_comparison import _paired_input

root = Path('data/validation/paired-ccpvdz-atom')
dft = json.loads((root / 'dft/benchmark.json').read_text())
for value in ('atom', None, 'invented'):
    altered = copy.deepcopy(dft)
    altered['quantum_settings']['scf_initial_guess'] = 'minao'
    record = altered['species']['methane']
    record['quantum_settings']['scf_initial_guess'] = value
    record['quantum_diagnostics']['settings']['scf_initial_guess'] = 'minao'
    _paired_input(altered, root, 'methane', 'cc-pvdz')  # currently accepts all
```

Suggested owner follow-up: validate explicit DFT guess values and check
consistency against expected settings and diagnostics, while giving historical
missing fields a declared legacy policy. Do not normalize `None`, an invalid
value, or an explicit disagreement into `minao`. Geometry/basis pairing and
electronic-state identity should remain distinct regardless of this check.

## Writer and reader map

| Path | Existing behavior | New-field impact |
|---|---|---|
| `nanodesign/design.py:54`, `:90` | Constructs `QuantumSettings(**data['quantum'])`; returns raw design plus original file hashes | Omitted legacy input gets the new default at runtime; raw input remains unchanged. Invalid explicit null is rejected by the new settings validator. |
| `nanodesign/workflow.py:117`, `:126`, `:154` | Saves `asdict(settings)` and uses the same immutable settings for initial/final/path calculators | Effective guess is saved for new runs and shared by every image. Raw `design.quantum` can still omit the field, showing the difference between input and effective settings. |
| `nanodesign/workflow.py:233`, `:243`, `:245` | Characterization saves full quantum settings and the baseline diagnostics | New guess is present in the outer characterization result and baseline diagnostics, including failed/interrupted output paths. |
| `nanodesign/benchmark.py:183`, `:193`, `:220`, `:254` | Replaces only species charge/spin; saves full base/species settings and final diagnostics | Non-default DFT guess survives species replacement and is archived in new species and aggregate evidence. No fixed list drops the field. |
| `nanodesign/method_comparison.py:103`, `:111`, `:117`, `:128` | Separate `cc_initial_guess` option; DFT is constructed at its defaults and saved under `dft_settings` | `--cc-initial-guess atom` continues to affect CC only; DFT now records `minao`. This interface intentionally does not expose a DFT-guess option. Do not infer matched guesses from the CC label. |
| `nanodesign/cli.py:93` | Benchmark `--settings` reads either a settings dictionary or `design.quantum` | New explicit guess flows through ordinary JSON input; old omitted input remains loadable. |
| `nanodesign/cli.py:88` | `audit` reads JSON and calls numerical/scientific status gates | No dataclass reconstruction or strict settings-key equality; old artifacts remain readable. This is not a settings-provenance validator. |

## Stationary/vibration provenance boundaries

`nanodesign/stationary.py:170-177` returns `settings` containing only Hessian
step, force/frequency tolerances, and coordinate limit. It never serializes the
attached calculator's quantum settings. This is appropriate for its generic ASE
calculator interface, but the dictionary is not a self-contained DFT provenance
record. Preserve its outer `workflow.run_characterization` result when exporting
it. A bare `stationary` object cannot establish which SCF guess was used.

The S1 helper `research/reference-saddle/saddle_search.py:382-391` similarly
returns a characterization plus geometry/timing, and its outer summary records
quantum settings at line 548. Its raw-PySCF guess-scan function is a separate
calculation path; setting `QuantumSettings.scf_initial_guess` does not change its
explicit loop of four guesses. This review did not assess or modify S1's science
or acceptance criteria.

`nanodesign/quantum.py:336-350` writes electronic events without complete
settings; before the field addition this was lines 317-331. Events alone therefore
do not identify an initial guess. The current `run_characterization` provenance
scope (`workflow.py:246`) correctly says the saved diagnostics cover the initial
geometry and logs cover later displaced calculations. Keep the enclosing run
settings when extracting those logs. Adding a selected guess does not establish
that all displaced calculations followed one electronic solution.

## Historical absence and display

These completed records were sampled; each lacks a DFT `scf_initial_guess` in
its saved `quantum_settings` or `dft_settings`:

- `data/validation/h2-integration/modes/result.json`
- `data/validation/h-abstraction-direct-initial/result.json`
- `data/validation/pbe0-tzvp/benchmark.json`
- `data/validation/paired-ccpvdz-atom/method_comparison.json`

Do not rewrite them to look as though an explicit guess had been recorded.
Display “not recorded” for the archived field; an inference from the old backend
default should be separately labeled. The suffix `paired-ccpvdz-atom` denotes the
CC guess and is not proof of an `atom` DFT calculation.

Read-only UI observation for V1: `workbench/static/app.js:306-326` retains full
settings in its expanded diagnostics, so new workflow guesses are inspectable
there, but the main method list at line 312 does not show the DFT guess.
`referenceCard` at lines 339-352 displays only the CC guess and does not expose
`dft_settings`. This is an optional clarity improvement for V1, not a request to
modify its files from D2.

## Scientific claim review

The added diagnostics (`nanodesign/quantum.py:160-167`) explicitly keep
`initial_guess_scan_performed`, `ground_state_verified`, and
`electronic_state_identity_verified` false and explain the one-guess scope.
Existing benchmark flags at `benchmark.py:259`, paired-comparison flags at
`method_comparison.py:132-138`, and characterization validation at
`workflow.py:271-272` remain conservative. No added control should upgrade those
claims or treat one converged result as an electronic-state or reaction-path
verification.

## Validation and limits

- Runtime settings probe after the owner added the field: omitted input ->
  `minao`; explicit `minao` -> `minao`; explicit `atom` -> `atom`; explicit null
  rejected. The caller's legacy input dictionary was not mutated.
- The paired-evidence probe above performed only JSON reads, deep copies, XYZ
  parsing, and hash checks. It did not write archive files or calculate energies.
- Existing focused test modules were run with their analytical/mocked backends:
  `tests/test_characterization_workflow.py`, `tests/test_benchmark.py`, and
  `tests/test_method_comparison.py`. Result: 70 passed in 1.85 seconds; 165 existing ASE/NumPy deprecation warnings.
- These current test fixtures generally construct settings using the live
  dataclass; the characterization fixture at line 60 and paired fixture at line
  31 automatically adopt new defaults. They are not historical missing-field
  regression tests. The explicit guess assertion in method-comparison tests is
  currently CC-only. A future non-default DFT propagation test should use `atom`
  so a forgotten passthrough cannot hide behind the default.
- No quantum calculations, SCF state checks, or scientific energy reproduction
  were performed by this review. Findings concern code/data contracts only.
