# A2 — is the 53-atom candidate reaction path computable here?

Lane owner: session `andresarriaga-a8`. Brief: `docs/AGENT_TASKS.md`, section A2.
Status and bounded plan: `coordination/status/andresarriaga-a8.md`.

This directory answers a cost question and only a cost question. **A cost
projection says nothing about whether the reaction works.** An affordable path
would not be evidence that the tool abstracts the right hydrogen, and an
unaffordable one is not evidence against it. Those questions belong to A1
(site selectivity) and S1 (locating a genuine saddle).

## Files

| File | What it does | Runs quantum jobs |
|---|---|---|
| `timing.py` | CPU-time / wall-clock / load instrumentation | no |
| `reduced_models.py` | Builds reduced-handle variants by truncating the repository's own candidate | no |
| `cost_model.py` | Projects path wall-clock from the archived 53-atom records | no |
| `gradient_cost.py` | Measures per-evaluation CPU cost of each variant | yes |
| `handle_fidelity.py` | Tip H-affinity across handles, with four-guess SCF scans | yes |
| `knobs.py` | Density fitting, basis, and force tolerance, measured both ways | yes |

Evidence is written to `evidence/`, never to a directory named `runs/`, and is
never overwritten: each script refuses to start if its output file exists.

## Why CPU time and not wall clock

Every wall-clock number archived in this repository was taken while other
single-threaded quantum jobs shared the host. While this lane was running, the
one-minute load average on this eight-core machine was above 200. A wall-clock
measurement under that load describes the scheduler, not the calculation.

This PySCF build has no OpenMP and runs on one thread regardless of the
`threads` setting, so process CPU time is the transferable single-thread cost
and a lower bound on uncontended wall clock. Both numbers, the contention
factor between them, and the load average at start and end are recorded for
every timed calculation. The report quotes CPU time for anything that is meant
to transfer to another machine, and labelled wall clock for anything that
describes waiting on this one.

## Two things the handle reduction changes, only one of which is measured here

`handle_fidelity.py` measures the **electronic substituent effect**: how much
the handle changes the apex radical's appetite for a hydrogen, as the shift in
tip H-affinity across handles at fixed geometry.

It does not measure **mechanical and steric fidelity**. The handle also sets
the mount's stiffness, its mass, its sterics against the substrate, and where
the anchors sit. A single fixed atom on a methyl handle is a *different*
boundary condition from three fixed carbons in a cage, not a scaled-down one.
A small H-affinity spread is therefore evidence about electronics alone, and
the report does not let it stand in for the rest.

## Reproducing

From the repository root, with the project virtualenv:

```
python research/candidate-feasibility/cost_model.py
python research/candidate-feasibility/handle_fidelity.py --handles hydrogen methyl
python research/candidate-feasibility/handle_fidelity.py --handles adamantyl
python research/candidate-feasibility/gradient_cost.py --handles hydrogen methyl
python research/candidate-feasibility/knobs.py --knob density-fitting
python research/candidate-feasibility/knobs.py --knob basis
python research/candidate-feasibility/knobs.py --knob tolerance
```

One process at a time, as declared in the status file. The host is shared with
other lanes' jobs; do not launch these in parallel with each other.

## Traps this lane paid for again

The default `minao` SCF guess landed 11.2 kcal/mol above the correct solution
for the propynyl (methyl-handle) radical at PBE0-D3(BJ)/def2-SVP. This is the
same trap the forensics lane found for the ethynyl radical at Hartree-Fock,
reproducing here in **DFT**, on a species this lane needed. Every open-shell
species in this directory therefore gets a four-guess scan before its
production calculation, and the selected guess is recorded in the evidence.
