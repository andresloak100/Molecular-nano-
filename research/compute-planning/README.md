# Computation planning from saved measurements

This read-only calculator estimates serial computation scenarios for the existing
53-atom candidate. It reads saved evidence and performs arithmetic. It launches
no quantum jobs and does not change the evidence, schedule work, or choose a
design. A2 owns the separate reduced-handle measurements and fidelity research
in `research/candidate-feasibility/`.

With seven images and **200 assumed position updates in each of four stages**,
the present workflow needs about **2,407 new energy-and-force evaluations** under
normal cache reuse. Holding the recorded starting-geometry cost constant gives
**19.54 days with density fitting or 56.97 days with direct DFT**. These are
conditional projections from uncontrolled, overlapping measurements. They are
not completion deadlines, confidence intervals, or guarantees of convergence.

## Use

From the repository root, with Python 3.11 or later:

```sh
python3 research/compute-planning/estimate.py
python3 research/compute-planning/estimate.py --images 7 --steps 50
python3 research/compute-planning/estimate.py --max-free-coordinates 141
```

The script needs only Python's standard library. JSON goes to standard output;
it includes source paths/hashes, measurement settings, stage counts, hours/days,
and interpretation limits. `--project-root` can point to another checkout holding
the same archive layout. Missing, failed, changed, or nonfinite timing evidence
causes an explicit error; no substitute timing is invented.

`--steps` is an **assumed actual update count in each stage**, not a new optimizer
control. The production workflow's one `--steps` argument sets a separate maximum
for each stage, which can converge earlier or fail. Library function
`path_counts(images, endpoint_steps, neb_steps, climbing_steps)` supports different
assumed stage counts; `endpoint_steps` applies to each of the two endpoints.
Zero assumed updates still require initial geometry evaluations.

## Measured inputs

Both archived runs use PBE0-D3(BJ)/def2-SVP on identical initial coordinates,
neutral doublet C22H31: **53 atoms, 463 basis functions**, six fixed atoms
`[7, 8, 9, 35, 36, 37]`. The single-point structures are unrelaxed and do not
establish a reaction barrier or a functional tool.

| Record | Full single-point elapsed | Calculator elapsed | Requested / recorded effective PySCF threads |
|---|---:|---:|---|
| [Direct](../../data/validation/h-abstraction-direct-initial/result.json) | 2,044.842914 s | 2,044.835078 s | 2 / not recorded |
| [Density fitting](../../data/validation/h-abstraction-df-initial/result.json) | 701.531096 s | 701.440628 s | 2 / 1 |

Projections use **full single-point-run elapsed per new geometry**, including SCF,
the analytical gradient, and small workflow overhead. Gradient-stage time alone
would omit the SCF work required at a new geometry. The difference between the
two timing fields is retained rather than silently changing the timing basis.

The direct run predates thread instrumentation: its missing effective-thread
field stays `null`. Later checks on the same tested PySCF build observed one
effective thread. The repository records an eight-logical-CPU, 24-GiB host in
[coordination notes](../../coordination/README.md). Other processes shared it;
neither the observed timing ratio nor these projections establish a controlled
DF speedup. See the original [comparison caveat](../../data/validation/h-abstraction-df-initial/comparison.json).

## Stage accounting

Let `m = images - 2`, `R` be updates in **each** endpoint relaxation, `P` ordinary
NEB updates, and `C` climbing-image updates. With unchanged caches between stages:

```text
initial endpoint:      R + 1
final endpoint:        R + 1
ordinary NEB:          m × (P + 1)
climbing NEB:          m × C
total:                2R + 2 + m × (P + C + 1)
```

Each optimizer has a step-zero force request. Ordinary NEB requires a new initial
evaluation for every interior image. Climbing NEB starts on the same image
objects and unchanged coordinates, so its initial quantum results are cached;
only the NEB force projection changes. For seven images and 200 updates per stage:

- **2,400** counts updates only and omits initialization.
- **2,407** includes initialization and the expected five-image cache reuse.
- **2,412** counts both band starts without crediting that cross-stage reuse.

The last two are accounting scenarios, not statistical bounds. Restarts,
additional state checks, repeated calculations, cache loss, different geometry
costs, or errors can put real runtime outside these values. Early convergence,
failed endpoint/NEB guards, or unchanged image coordinates can reduce work.

FIRE logging/convergence checks, NEB force projection and final summaries request
energies/forces more often than these counts. The quantum calculator supplies
**both energy and force in one evaluation** and ASE normally caches repeated
reads. IDPP initialization temporarily uses its own distance objective and
restores the quantum calculators; it does not normally launch DFT. IDPP and
optimizer/I/O overhead are not separately timed here. These details were checked
against [workflow.py](../../nanodesign/workflow.py),
[quantum.py](../../nanodesign/quantum.py) and installed ASE 3.29.0.

| Assumed updates per stage | New evaluations with cache reuse | Direct projection, days | DF projection, days |
|---:|---:|---:|---:|
| 20 | 247 | 5.85 | 2.01 |
| 50 | 607 | 14.37 | 4.93 |
| 100 | 1,207 | 28.57 | 9.80 |
| 200 | 2,407 | 56.97 | 19.54 |

## A complete vibration check is additional work

Six whole-atom anchors leave `3 × (53 - 6) = 141` free Cartesian coordinates.
Central differences require **283 force requests in the Hessian routine**:
one baseline plus two displacements per coordinate. No external-mode subtraction
reduces this Cartesian evaluation count.

The production characterization wrapper creates a fresh calculator and first
checks the initial force. Its baseline is reused by the Hessian routine, giving
**284 force-array requests but 283 new quantum evaluations** for a complete run.
At a constant saved single-point cost, that projects to **55.15 hours (DF)** or
**160.75 hours (direct)**. A second independent displacement-step check requires
another full stencil under the same assumptions. See
[stationary.py](../../nanodesign/stationary.py) and
[run_characterization](../../nanodesign/workflow.py).

The default `max_free_coordinates=120` **rejects this 141-coordinate candidate
before any quantum calculation**. If the user explicitly raises the production
limit, a failing initial-force guard costs one calculation and stops the stencil.
The calculator's `--max-free-coordinates` option only reports the guard's effect;
it changes no production limit. Stencil precision and other input guards can
also reject work. A completed Hessian still needs mode, connectivity, step-size,
and electronic-state interpretation.

## CPU and GPU boundary

The implemented adapter runs CPU PySCF; band images run serially and changed
geometries do not reuse a saved SCF density. The tested build ignores additional
PySCF thread requests because OpenMP is unavailable; this does not describe every
linked numerical library. Independent processes still share CPU and memory.
This calculator does not divide time by the host's core count.

[Architecture documentation](../../docs/ARCHITECTURE.md) identifies GPU4PySCF on
NVIDIA/CUDA as a possible future route. The repository has no integrated or tested
GPU adapter or parallel NEB scheduler. No GPU speedup, GPU memory requirement,
hardware purchase recommendation, or CPU/GPU numerical equivalence is inferred
here. New backend measurements must cover the required energies **and forces**.
Changing the tool handle or reactive region changes the physical model; reuse
of these 53-atom timings does not certify a smaller model's fidelity.

## Checks

```sh
.venv/bin/python -m unittest discover -s research/compute-planning/tests -v
```

Verified on 2026-09-16: **11 tests passed in 1.928 seconds**, including the
workflow integration case. The estimator also ran with `python -S`, confirming
that it requires no installed scientific packages. The only recorded warning
was an existing ASE/NumPy deprecation during the integration check.

Tests cover initialization/cache accounting, characterization guards, rejected
inputs, archived measurements, missing thread evidence, tamper/failure handling,
and unchanged input bytes. One optional integration test uses production workflow
transitions and real ASE caching/IDPP with a synthetic analytic counting fixture
and prescribed optimizer updates. It is a software check, not a chemistry
simulation; it launches no quantum jobs. Without the project's ASE/NumPy
environment, that one integration check is skipped.
