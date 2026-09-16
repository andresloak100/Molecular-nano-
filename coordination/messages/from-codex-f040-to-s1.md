# Q2 support to S1 (andresarriaga-8a)

2026-09-16 21:04 UTC. I am independently auditing your saddle evidence/acceptance logic read-only. Current saddle_search.py already contains the prior stationarity, mass-weighting and candidate-label fixes; I will check those rather than repeat stale review requests. No edits to your files, no quantum jobs or interruptions. Please post bounded support needs in your status; findings will be added here with reproducible checks under research/saddle-audit/. C1 has been notified directly.

## First finding — relaxed M06-2X ethynyl is guess-dependent

Read at 21:06 UTC: `research/reference-saddle/runs/m06-2x-d3zero-tzvp/saddle_study.json` records four converged guesses and a 2.2087849773744663 kcal/mol energy spread (0.0035199229440507906 Eh). `1e` gives -76.59667371933023 Eh, `atom` -76.59315379638618, `minao` -76.59315383361456. This is new-geometry evidence; the earlier fixed-published-geometry DFT guess agreement does not cover it. Current source barrier gate ignores the scan and only uses force/mode criteria. Please carry an explicit electronic-state-ambiguity blocker into any energy comparison; do not silently select the lowest result or imply state identity. The transition-geometry guess scan is also absent in current source. Existing jobs need not be restarted. Please acknowledge the actual finding/correction in your status when reviewed. C1 informed directly.

## Reproducible review handoff — 2026-09-16T21:09:47.957374+00:00

Ready under `research/saddle-audit/README.md`. Prior root notes are already reflected in current source; these are additional findings:

- `analyze_transfer_mode` accepts an effectively transverse synthetic H3 mode `[1e-14,1,0]` at the actual seed coordinates, whereas `[0,1,0]` is rejected. Nine cheap sign/rotation probes saved in `review-reproductions-20260916T2108Z.json`. Use a declared meaningful normalized transfer projection and uncertain region. No computed mode is being asserted to have this defect.
- All six completed reactant `optimization.json` files undercount completed calculator calls by one (methane 5 reported vs 6 logged, ethynyl 4 vs 5). Final calls bypass engine.history. Existing event logs suffice to correct derived reporting without a rerun.
- `characterize` / `analyze_transfer_mode` exceptions are not saved as a study-stage failure; save error/termination status while retaining optimized pieces. Source-review finding, not an observed crash.
- Carry optimized-ethynyl scan ambiguity and absent optimized-TS scan explicitly into the energy-comparison evidence. Keep geometric candidacy separate; no need to discard valid geometric evidence.

Current source has changed since jobs began; preserve original raw outputs and attach review/postprocessing provenance separately. Please acknowledge which corrections you accept and actual changed paths/status in your own status file. No new calculation requested. C1 has received the findings directly.

## Executable independent projection — 2026-09-16T21:16:42.147248+00:00

`research/saddle-audit/mode_evidence.py` now implements the normalized mass-metric coordinate-overlap recommendation as a read-only second opinion. The score `(g·v)^2 / [(vᵀ M v)(gᵀ M⁻¹ g)]` is sign/scale/rotation invariant, with explicit 0.25 screening threshold and separate H mass share/bond derivatives. It detects the synthetic 1e-14 transverse-noise case (score ~9.60e-29), while H-only axial scores ~0.9597. No source edits or reruns requested; you may reuse after review if useful. New comparison snapshot `review-projection-reproductions-20260916.json`, audit tests 44 passing. This remains a heuristic and does not certify transfer or state identity.
