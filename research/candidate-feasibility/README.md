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
| `timing.py` | CPU-time / wall-clock / load / page-fault instrumentation | no |
| `reduced_models.py` | Builds reduced-handle variants by truncating the repository's own candidate | no |
| `cost_model.py` | Projects path wall-clock from the archived 53-atom records | no |
| `positional_control.py` | Thermal mis-targeting: the tolerance the design needs, and whether it has it | scan only |
| `tunneling.py` | Whether a classical barrier can become a rate; the ω* → regime map | no |
| `gradient_cost.py` | Measures per-evaluation CPU cost of each variant | yes |
| `handle_fidelity.py` | Tip H-affinity across handles, with four-guess SCF scans | yes |
| `knobs.py` | Density fitting, basis, and force tolerance, measured both ways | yes |

Two documents: `REPORT.md` is the A2 cost deliverable;
`PHYSICS_COMPLETENESS.md` maps what stands between this repository and a model
that could actually predict whether a molecular machine works.

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

The default `minao` SCF guess lands well above the correct solution for the
open-shell species this lane needs, at PBE0-D3(BJ)/def2-SVP:

| Species | minao error | minao S² | correct S² |
|---|---|---|---|
| Ethynyl, HC≡C· | none — all four guesses agree | 0.7913 | 0.7913 |
| Propynyl, CH₃-C≡C· | **+11.20 kcal/mol** | 0.7521 | 0.7876 |
| Adamantyl-ethynyl, C₁₂H₁₅· | **+10.42 kcal/mol** | 0.7521 | 0.7863 |

The last row is the tool tip of the repository's headline candidate. Combined
with the forensics lane's Hartree-Fock cases, the pattern is four for four:
**the wrong, higher solution is the one with the cleaner S²**. The ideal
doublet value is 0.75, so the trap looks like the better answer by exactly the
diagnostic a reader would reach for. Do not select on S².

Every open-shell species here therefore gets a four-guess scan before its
production calculation, and the selected guess is recorded in the evidence.

Note for whoever uses the archived 53-atom records: they were run on the
default guess, with `scf_initial_guess` and `initial_guess_scan_performed` both
null. Their S² of 0.78462 sits near the correct region rather than the trap's,
so they are probably fine — but that is an inference across two different
systems, not a check, and S² is the wrong thing to infer from.

## What this lane found beyond cost

Two results that came out of asking what the cost was *for*:

- **Positional control has enormous margin.** The apex may move 2.495 Å in any
  direction before a different hydrogen is nearest, which needs only 4.8 N/m of
  stiffness at 298 K for a 10⁻¹⁵ error rate. Thermal wander is not the binding
  risk — *given* that the nearest hydrogen is the one that reacts, which is
  A1's open question.
- **Whether a classical barrier can become a rate is unresolved**, and one
  number settles it: the imaginary frequency at the saddle. On the repository's
  only measured value (259i cm⁻¹) tunneling is a 7% correction; on the value
  implied by the experimental activation energy (~1648i) it is the mechanism.
