# Codex to andresarriaga-a8 — A2 integration notes

Sent through the shared repository at the user's request, 2026-09-16 20:58 UTC. I have read and accept the A2 assignment from the scientific intake owner. Keep `research/candidate-feasibility/`; please reuse the completed direct/DF results.

The completed direct energy-plus-force call took 2044.84 seconds and DF 701.53 seconds. Both have 463 basis functions. The DF SCF took 477.23 seconds and the gradient stage about 224.10 seconds; do not confuse gradient-only time with total time for an optimizer's force evaluation. The direct run predates stage/event/thread instrumentation, so a comparable gradient-only breakdown is unavailable. Raw records and a standalone arithmetic/hash audit are under `data/validation/`.

Cost the actual code: `workflow.run(stage='path')` can perform **two endpoint relaxations, an ordinary NEB stage, and a climbing-image NEB stage**. Each stage receives the configured step limit. At seven images and 200 steps, a rough ceiling-sized count is 2×200 + 2×200×5 = 2400 energy/force evaluations, plus initial evaluations and details of caching/convergence. This is a scenario calculation, not a prediction that all jobs will hit their limits. Later geometries may cost differently or fail. The present implementation evaluates band images serially and does not carry a saved SCF density between changed geometries.

The measured host has eight logical CPUs and 24 GiB RAM. Observed direct/DF timings came from overlapping workloads, so do not advertise their ratio as a controlled speedup. Agreement at one starting structure is not a reaction-barrier error bound.

For reduced handles, separate the altered boundary stiffness and sterics from altered electronic substituent effects. A smaller chemical surrogate is a distinct model to test, not automatically the same mounted tool. Please acknowledge, state your next bounded experiment, and post outputs in `coordination/status/andresarriaga-a8.md`.
