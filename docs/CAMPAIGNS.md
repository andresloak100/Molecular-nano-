# Pose campaign interfaces

`nanodesign.campaign` is the orchestration boundary for evaluating a finite family of poses. It uses the same real quantum calculation path as individual jobs. It currently enumerates explicitly supplied candidates; it does not learn a potential, optimize an unknown tool topology, or establish physical accuracy.

## Python entry points

| Function | Purpose |
|---|---|
| `create_campaign(designs, output, stage=..., state=..., fmax=..., steps=..., images=...)` | Validate and snapshot 1–100 existing designs into a new campaign. No solver is invoked. |
| `create_pose_campaign(output, separations, offsets, settings=None, **controls)` | Generate the finite Cartesian product of separations and lateral offsets for the supported H-abstraction candidate, then snapshot it. |
| `run_campaign(directory, max_jobs=1, retry_incomplete=False)` | Start at most the specified number of serial attempts. Completed jobs are skipped. |
| `campaign_report(directory)` | Return evidence rows in input order, with counts and explicit scientific interpretation. |

The entry points return JSON-compatible dictionaries. CLI equivalents are `campaign-create`, `campaign-run` and `campaign-report`. `--settings` for a generated grid accepts a quantum settings object or a design containing `quantum`; existing designs retain their own matching settings.

## Persistent contract

- `plan.json`: schema version, explicit calculation controls, quantum settings, stable pose identifiers, relative design paths and input hashes. It is immutable after creation.
- `designs/pose-NNNN/`: exact coordinate snapshots, runnable `design.json`, and the original `source-design.json` for provenance. Atom indices are unchanged.
- `campaign.json`: the plan hash and each pose's attempt ledger. Updates use atomic replacement.
- `runs/pose-NNNN/attempt-NNNN/`: ordinary workflow outputs, including `result.json`, input geometries and electronic events. Attempt directories are never reused.
- `.campaign.lock`: an operating-system advisory lock, automatically released on process exit. Do not remove this file during work.

Pending poses have no attempts. Attempts may be running, completed, not converged, failed, interrupted or abandoned. A worker that dies without writing a terminal result is marked abandoned on the next run. If the quantum job finished but the campaign worker died before updating its ledger, its saved matching result can be recovered without repeating the calculation. Explicit retries create new attempts from preserved inputs. All prior attempts remain inspectable.

Only one worker may execute a given campaign at a time. Different campaigns may run in separate processes, subject to the machine's resources. There is no distributed scheduler or GPU routing. The maximum-job setting bounds launched attempts, not runtime, memory, SCF iterations across the campaign, or the total calculations within a reaction-path job. SCF and optimization limits remain explicit in each job's settings.

## Scientific interpretation

Matching element order and settings makes basic comparisons more interpretable, but neither proves that electronic states match nor removes pose-dependent strain or interaction energy. Reports show reaction-basin and topology screens beside numerical quantities. A completed single-point calculation may still be an unsuitable starting geometry. Missing values remain absent/null instead of becoming favorable scores.

`design_validated` and `automatic_ranking` remain false. The interface can support a future graphical workspace without inventing reliability percentages or hiding numerical failures. A future automated optimizer will require a calibrated objective, competing pathways, state checks and uncertainty treatment in addition to this bookkeeping layer.
