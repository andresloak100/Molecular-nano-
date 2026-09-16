# Off-host compute lane

Owner: `molecular-nano-6c` (cloud Linux container, isolated from the shared
Mac). This lane exists to do two things the contended Mac cannot: (1) produce
uncontended, multi-threaded timings of the exact archived 53-atom single
points, and (2) verify that the production stack is numerically portable to a
different OS, architecture and BLAS. It launches no new science: every
calculation here is a byte-identical replay of an already-archived input.

## Environment

- Linux x86_64, 4 logical cores, 15 GiB RAM. Zero contribution to the Mac's
  load.
- PySCF 2.14.0 **with working OpenMP**: `lib.num_threads()` honors requests
  (the Mac's arm64 wheel caps at 1). This is the first environment on the
  project where `threads_honored` can be true; the G1 inspector confirms it
  in `gpu-readiness-offhost/offhost-preflight.json`.
- Full test suite green on Linux (197 non-quantum + 7 quantum).

## Result 1 — the production stack is numerically portable

Both 53-atom single points reproduce the archived Mac records essentially
exactly, across different OS/arch/BLAS:

| Run | Energy Δ (off-host − Mac), eV | Max free-force Δ, eV/Å |
|---|---|---|
| Density fitting | 9.8e-11 | 1.9e-11 |
| Direct | −8.7e-11 | 2.0e-13 |

These deltas are ~10 orders of magnitude below chemical significance and far
below the archive's own DF-vs-direct differences (−0.0065 eV, 4.7e-4 eV/Å). A
result computed on this container can be trusted to match one computed on the
Mac. This is what makes off-host offloading safe to build on.

## Result 2 — uncontended timings, and what they say about A2's verdict

The exact archived single points, replayed serially on an idle host, 4 OpenMP
threads (load logged in `timing-log.txt`):

| Run | Off-host wall (4 threads) | Archived Mac wall | Mac/off-host ratio |
|---|---:|---:|---:|
| Density fitting | 568.0 s | 701.5 s (1 contended thread) | 1.23x |
| Direct | 789.2 s | 2044.8 s (threads not recorded) | 2.59x |

Read these carefully — they are different hardware, so the ratio mixes core
speed, thread count and the Mac's contention, and is not a controlled
benchmark of either machine.

The load-bearing observation for A2: **the density-fitting run barely sped up
(1.23x) despite going from 1 to 4 threads on a fresh host.** If the archived
DF run had been heavily contended (say the 10.4x factor measured elsewhere at
load 201), four uncontended threads would have crushed it by far more than
1.23x. They did not. So the archived 701.5 s DF number was close to that
hardware's genuine cost, and the *days* end of A2's feasibility range — not
the 20-hour end — is the live one. The direct run's larger 2.59x gain is
consistent with it having been genuinely single-threaded (no OpenMP) on the
Mac, so four working threads help it more; contention need not be invoked to
explain it.

This is diagnostic, not final: wall-clock on 4 threads still confounds core
speed with threading. The single-thread CPU-time run below is what separates
them cleanly, and is the measurement A2 named as its top outstanding item.

## Result 2b — single-thread CPU time resolves A2's range toward "days"

A2 named a single-thread CPU-time measurement as the item its verdict hangs
on. Done (`cputime-summary.json`), DF 53-atom energy+gradient, one process:

- **Total CPU: 1807 s** (1748 user + 59 system), wall 1585 s, peak RSS 2.6 GB,
  1 major page fault (not memory-starved). cpu/wall ≈ 1.14 (BLAS spun slightly
  above one core in the gradient), so this is a clean, contention-independent
  cost.
- The archived Mac DF run was **701 s** wall at one thread. This container's
  single core took **1585 s** wall / 1807 s CPU for the same work — i.e. the
  Mac's single core is ~2.3x faster than this box's.

**Interpretation for A2, decisive:** if the archived Mac 701 s had been heavily
contended, an uncontended fast M-series core would have finished far under
701 s. Instead this (slower) container needed 1585 s single-threaded and only
matched the Mac by using four cores. So the archived 701 s was **near its
true uncontended cost**, not a 10x-inflated number. The DF path projection
therefore sits at the *days-per-pose* end of A2's range, not the 20-hour end.
The contention factor that dominates the shared host's *other* timings did not
materially inflate this particular DF record.

## Result 3 — projection arithmetic re-run on off-host cost

Re-running C2's evaluation-count arithmetic (2,407 evaluations for the default
7-image / 200-step / four-stage scenario, unchanged) with the off-host
per-evaluation cost, alongside the archived Mac cost:

| Path scenario source | Days (one pose) |
|---|---:|
| Density fitting, Mac archived | 19.54 |
| Density fitting, off-host 4-thread | 15.83 |
| Direct, Mac archived | 56.97 |
| Direct, off-host 4-thread | 21.99 |

Off-host 4-thread execution roughly halves the direct-DFT ceiling (57 → 22
days) and trims the DF ceiling (19.5 → 15.8 days). Both remain measured in
**days per pose**, and the campaign has nine poses. Offloading to this class
of container helps, but does not by itself move a full reaction path into
"run it overnight." The structural cost — 530 to 2,400 serial evaluations per
pose — is what a faster machine rescales but does not remove. Machine json:
`offhost-comparison.json`.

## Files

- `design-53atom-{df,direct}-4t.json` — replay designs (threads 4).
- `design-53atom-df-1t.json` — single-thread CPU-time design.
- `run_53atom_timings.sh`, `run_1thread_cputime.sh` — drivers (loads logged).
- `runs/53atom-*/` — completed evidence, byte-identical inputs to the archive.
- `analyze_offhost.py`, `offhost-comparison.json` — read-only comparison and
  projection re-run.
- `gpu-readiness-offhost/offhost-preflight.json` — G1's inspector on this box.
- `timing-log.txt`, `cputime-log.txt` — wall/CPU logs with load averages.

## Boundaries

Read-only against the archive; no core edits; no writes to any archived
record. The archived Mac evidence remains the record of what ran there; this
lane records what the same inputs do on a different, uncontended machine.
