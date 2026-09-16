# Validation contract: supported-tip hydrogen abstraction

**Status: research protocol; no chemically validated tool or operating reliability is established.**
Evidence snapshot: 2026-09-16. This document specifies the observations needed to
support a design decision, rather than treating successful software execution as
that decision. Live lane results may advance after this snapshot; link a new
dated assessment to new evidence instead of silently upgrading these conclusions.

## 1. Define the operation before evaluating it

The first calibration is the isolated reaction
`C2H· + CH4 → C2H2 + CH3·`. It tests a small reaction and the calculation pipeline.
It is not the supported tool, a diamond surface, or a complete assembly cycle.

The current tool candidate is the constructed, nonperiodic **C22H31, 53-atom,
neutral doublet** in [the example design](../examples/h-abstraction/design.json).
Zero-based indices identify donor carbon **0**, transferring hydrogen **10**,
and accepting ethynyl apex carbon **26**. Atoms **7, 8, 9, 35, 36, 37** are fixed;
the first three belong to the target and the last three to the tool handle.
The initial and final files contain the same atoms. The six fixed atoms give
141 free Cartesian coordinates. They represent an idealized boundary, not a
measured mount stiffness. Generated bond lists and product coordinates are
hypotheses to test, not observed bond orders or relaxed products.

For an operating claim, define the controlled mount pose `λ(t)`, temperature,
environment, loading/withdrawal schedule and uncertainty in pose. Success means
H10 is removed from C0 and retained at C26, the intended target radical and tool
product survive withdrawal, and specified non-target bonds remain intact.
Regeneration is a separate operation with an explicit hydrogen destination and
material/energy accounting. Hydrogen removal alone does not demonstrate carbon
incorporation or manufacture of a molecular machine.

### Required decision inputs

All entries below are **unset** for the present candidate. The scientific owner
must specify values and a rationale before a calculation can pass an operating
accuracy or reliability gate. Software defaults do not supply these values.

| Decision input | Required declaration |
|---|---|
| Intended use | Screening, ranking, quantitative prediction, or experimental operating specification; exact chemistry/model domain |
| Operating conditions | Temperature range, environment, charge/spin preparation, applied fields if any, pose trajectory, dwell time, loading rate |
| Mechanical uncertainty | Physical mount model and distribution/range of position, orientation, compliance and thermal motion |
| Success criterion | Allowed final products and damage; maximum error/no-operation probability; cycle count and required confidence |
| Numerical accuracy | Budgets for the particular relative energies, forces, geometries and mode classification that drive this decision |
| Model accuracy | Allowed changes under basis/reference/cluster/boundary/thermal-model refinement; defensible reference uncertainty |
| Resource bound | Assigned operator, bounded jobs/steps, effective threads, memory assumptions, stop/review conditions |

A chosen force cutoff such as `0.03 eV/Å`, SCF tolerance `1e-9 Ha`, grid level 3,
or unresolved-frequency cutoff `20 cm^-1` is an algorithm setting. None is an
energy error bound or an operating tolerance. Force-to-energy conversion would
require additional information about displacement and curvature; a flat region
can satisfy a force cutoff without fixing the barrier to useful precision.

## 2. Keep the observables separate

| Observable | Definition and permissible interpretation |
|---|---|
| Nominal reference energy | `ΔE_nom = E(published nominal TS coordinates) − E(CH4) − E(C2H)` with each coordinate file, state and method specified. A fixed-geometry diagnostic, even when a legacy key says `barrier`. |
| Reaction energy | `ΔE_rxn = ΣE(products) − ΣE(reactants)` on specified states and geometries. State whether structures were optimized. This alone does not give a rate. |
| Conditional electronic barrier | `ΔE‡(λ) = E(first-order saddle; λ) − E(reactant minimum; λ)` after stationary-point and connectivity checks on the same electronic surface and boundary model. A separated-reactant reference must be reported separately. |
| Site H-removal energy | `D_i = E(target radical_i) + E(H) − E(target-H)` with consistent method and relaxation policy. Differences between sites measure thermodynamics, not competing transition rates. |
| Mechanical response | Raw force vectors, residual forces on free coordinates, anchor group forces and torques about a declared origin; changes with pose and boundary model. Constraint-masked anchor forces cannot measure mount load. |
| Thermal activation/rate | A specified free-energy barrier and kinetic model, including ensemble, nuclear effects and transmission/recrossing assumptions. It is not the fixed-geometry electronic energy above. |
| Operational reliability | Probability of the specified outcome over the complete pose/time protocol and uncertainty distribution, including no reaction and competing channels; eventually compared with measured outcomes. |

Never combine a saddle from one method with a reactant from another under an
unqualified barrier label. A composite energy protocol is possible only if its
geometry, state and correction definitions are explicit and its error is assessed.
Use the same units throughout a difference, retaining original values and the
conversion used. An energy below separated reactants can coexist with a barrier
above a pre-reactive complex; a failed saddle search proves neither barrierlessness
nor impossibility.

## 3. What the saved evidence currently establishes

These are observations from existing files, not new calculations performed for
this protocol. Preserve original archives; interpretation belongs in dated,
hash-bound sidecars or a new assessment.

| Evidence | Measured observation | Limit on the claim |
|---|---|---|
| [Original paired cc-pVDZ archive](../data/validation/paired-ccpvdz/method_comparison.json) and [interpretation](../data/validation/paired-ccpvdz/interpretation.json) | Nominal relative energies: DFT −3.557388, CCSD(T) −11.923873 kcal/mol. | The schema-1 archive omits initial-guess/state-verification metadata; MINAO attribution is retrospective. The CC value is not an accepted barrier. |
| [Explicit atom-guess paired archive](../data/validation/paired-ccpvdz-atom/method_comparison.json) | At the same archived coordinates and cc-pVDZ basis: DFT −3.557388, CCSD(T) +2.398616; difference −5.956004 kcal/mol. | CC explicitly used `atom`; archived DFT settings still lack an explicit guess. State identity is unverified. The difference is not a calibrated DFT error bar. |
| [Ethynyl CC record](../data/validation/paired-ccpvdz-atom/ccsd_t/ethynyl_radical.json) | UHF reference `S² = 1.221177`, versus 0.75 for a pure doublet. | This describes the reference determinant, not correlated CC spin purity or proof of multireference character. |
| [SVP](../data/validation/pbe0-svp/benchmark.json) and [TZVP](../data/validation/pbe0-tzvp/benchmark.json) | PBE0-D3(BJ) nominal relative energies −3.624839 and −2.890575 kcal/mol; nominal-geometry maximum forces 1.657059 and 2.383279 eV/Å. | These supplied geometries are not stationary on those DFT surfaces. The +0.734264 kcal/mol shift is fixed-geometry basis sensitivity, not a total error estimate. |
| [SI reproduction](../data/validation/si-energy-reproduction/si-energy-reproduction.json), [forensic conclusions](../data/validation/si-energy-reproduction/FINDINGS.md) | Four source energies reproduce within 1e-7 Ha using the investigated references; methane remains discrepant. | A proposed source transcription explanation is an inference. Keep the source unchanged. |
| [53-atom direct result](../data/validation/h-abstraction-direct-initial/result.json), [DF result](../data/validation/h-abstraction-df-initial/result.json), [comparison](../data/validation/h-abstraction-df-initial/comparison.json) | PBE0-D3(BJ)/def2-SVP, grid 3, 463 basis functions. DF−direct energy −0.006521680 eV; maximum atomic force-vector difference 0.000473513 eV/Å. DF maximum free force 0.958058 eV/Å. | Completed unrelaxed single points only. One geometry does not bound barrier/force errors elsewhere. DF recorded requested/effective threads 2/1; the older direct record lacks those fields. Concurrent elapsed times are not a speed benchmark. |
| [A1 geometric census](../research/site-selectivity/evidence/stage0-site-census-r3/census.json) | Four bridgehead and twelve methylene C–H sites; 2.494963 Å is the minimum apex displacement to an intended/competing-H perpendicular-bisector boundary with target atoms fixed. | A static nearest-site boundary distance, not an operating positioning tolerance, reaction probability or demonstrated chemical selectivity. |
| [A2 reduced-handle evidence](../research/candidate-feasibility/evidence/handle-fidelity-hydrogen-methyl.json) | Rigid-geometry H-affinity spread 0.194414 kcal/mol for the recorded H/methyl handles at DF PBE0-D3(BJ)/def2-SVP. | The adamantyl shift is null here. This does not establish full-handle electronic equivalence, sterics, stiffness or mounting fidelity. |

The paired comparison uses RHF/RCCSD(T) for closed-shell species and
UHF/UCCSD(T) for open-shell species, all-electron cc-pVDZ with no frozen orbitals.
It must not be relabeled as the paper's different restricted-open-shell/basis
comparison. The source paper itself reports an index-3 nominal methane structure
at UCCSD(T)/cc-pVDZ and discusses additional bending modes and symmetry constraints;
its tabulated nominal energy is not a certificate for our desired first-order
saddle. [Temelso et al., 2006, Table 1 and §3.1](https://vergil.chemistry.gatech.edu/static/pdfs/temelso_2006_11160.pdf)

Source coordinates and the recorded bohr-to-Å conversion are in
[provenance.json](../data/reference/provenance.json). The SI hash there binds the
original source, including its limitations. Reference reaction indices in S1 are
**donor 4, H 3, acceptor 2**, not the 53-atom indices above.

S1 owns ongoing DFT saddle research. Its mutable, unreleased working records are
not incorporated into this published evidence snapshot; no completed S1 saddle
and connectivity assessment is asserted here. Obtain an owner-released evidence
snapshot before updating those conclusions. Agreement between electronic guesses
at supplied coordinates must be reassessed at relaxed geometries. Missing terminal
results say nothing by themselves about process liveness.

## 4. Evidence and decision format

For each gate below, record an independent decision:
`supported_for_declared_scope`, `inconclusive`, `failed_declared_criterion`, or
`not_run`. These are proposed scientific-review labels, not implemented CLI fields.
An unset acceptance target means **inconclusive for that decision**, even when a
useful numerical observation is available. Do not collapse the gates into a
global `validated: true` flag.

Each record must identify the claim/domain; source and coordinate hashes; atom
order; charge and `spin = Nα−Nβ`; method, dispersion, basis, auxiliary basis and
resolved grid policy; constraints; initial guess and electronic branch evidence;
software versions; raw energies/forces; convergence and failures; and the exact
acceptance criterion with its rationale. Record requested and effective resources
and host contention when available. Missing historical fields remain missing.
Include the reviewer, date, dependencies, unresolved alternatives and next action.
Hashes and independent arithmetic checks establish consistency, not authenticity
of a physical claim or completeness of the model.

## 5. Gates for a defensible prediction

### G1 — Electronic surface and reference suitability

At each distinct reactant, product, candidate saddle and representative path
geometry, compare explicit starting guesses under identical numerical controls.
Retain every requested attempt, including nonconvergence. Report the converged
subset's spread; agreement of that subset cannot stand in for missing guesses.

Then inspect orbital stability, spin diagnostics, occupations and real-space
spin/density character. Follow electronic branches continuously across geometry,
including **every positive and negative force displacement used for a Hessian**,
as well as path/IRC images. Save the wavefunction/density information needed to
evaluate continuity. Switching independently to the lowest returned SCF solution
at each point can splice surfaces and corrupt gradients or curvature.
Internal/external SCF stability tests concern allowed orbital variations around a
solution; they do not prove the physical ground state. [PySCF SCF documentation](https://pyscf.org/user/scf.html)

**Acceptance:** a documented, reproducible branch assignment and reference model
suitable for the claimed domain, with unresolved branches either excluded by
evidence or retained as distinct alternatives. Persistent ambiguity leaves energy
ranking, Hessian interpretation and downstream kinetics inconclusive. Where
single-reference treatment is questionable, a specialist must define and converge
an appropriate correlated/active-space comparison; an arbitrary diagnostic cutoff
or active space does not close this gate. [PySCF multiconfigurational methods](https://pyscf.org/user/mcscf.html)

**Current gap:** production diagnostics explicitly report stability unchecked and
do not save orbital/density-overlap evidence. S2 supplies fixed-geometry guess
surveys, not state verification. Reconstructible forces and stable numerical modes
therefore do not by themselves establish a single electronic surface.

### G2 — Numerical convergence of the decision observables

Use matched coordinates and the G1 branch to isolate individual effects. Tighten
SCF controls; compare analytic forces with independent central energy differences
at several resolvable displacement sizes; examine force differences and residuals.
Repeat selected comparisons with finer DFT grids and explicitly recorded pruning,
radii and quadrature choices. A grid level alone is not the complete quadrature
specification. [PySCF DFT documentation](https://pyscf.org/user/dft.html)

Compare direct and density-fitted calculations, basis enlargement and relevant
auxiliary-basis choices at endpoints, the saddle region and difficult intermediate
geometries. Include combined tighter settings to detect coupled sensitivities.
First compare fixed-geometry observables; then reoptimize relevant stationary
structures under the proposed final settings and recompute the actual differences
used for the decision. Check near-contact basis effects and the definition of
separated versus interacting fragments rather than assuming cancellation.

**Acceptance:** changes in each declared observable, and any unresolved numerical
residual, satisfy its predeclared budget over the sampled domain with no unexplained
state change. Report the tested settings/domain and observed changes, not an
unproved global bound. Do not add correlated error components in quadrature or
treat a two-basis plateau as proof of the complete-basis limit.

### G3 — Stationary points and reaction connectivity

Optimize reactants/products with the same Hamiltonian and constraint definition
used for the path. Record raw forces and free-coordinate residuals. A converged
NEB/climbing image provides a saddle candidate; its maximum energy is not yet a
verified barrier. [ASE NEB documentation](https://docs.ase-lib.org/ase/neb.html)

For the candidate, converge force residuals, displacement size and the relevant
mass-weighted Hessian spectrum. Require one resolved negative mode in the physical
free-coordinate space, with a direction consistent with the intended reaction;
require minima at the relevant endpoints. Investigate soft/extra negative modes
under tighter numerics and state checks rather than deleting them with a convenient
frequency cutoff. Distinguish Cartesian displacement from mass-weighted amplitude.
The current code does not project translation/rotation; identify external modes
explicitly for a free molecule. A constrained 141-coordinate tool Hessian describes
only that fixed-boundary model, not the full support's vibration spectrum.

Follow both signs of the reaction direction using an IRC or demonstrably converged
downhill connection under the same method, basis and constraints, finishing with
endpoint optimization. Verify complete atom mapping, resulting chemical identities
and electronic branch continuity, not only the donor-H/apex distances. Save every
image, energy, raw force, path step and convergence criterion. The ORCA manual
describes the same-method forward/backward IRC requirement; citing it does not mean
an ORCA backend exists here. [ORCA IRC documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/irc.html)

**Acceptance:** first-order stationary character survives numerical refinement and
both connections reach the intended basins. Otherwise label a candidate or an
inconclusive search. For a genuinely barrierless proposal, establish a continuous,
resolved downhill route on the relevant surface and separately assess capture and
dynamics; absence of a located maximum is insufficient.

### G4 — Chemical calibration and transfer to the supported candidate

First resolve G1–G3 for methane/ethynyl. Compare the same observable against an
independently assessed reference, keeping electronic energies, ZPE-corrected
energies, enthalpies and experimental activation energies separate. Compare both
fixed-coordinate methods and each method's own verified stationary structures;
these answer different questions. Document reference/basis/correlation dependence
and source uncertainty rather than declaring any converged CCSD(T) result exact.

Next test chemistry nearer to the proposed operation: a tertiary donor, the actual
target cage and progressive ethynyl-handle substitution. Preserve radical geometry
and spin diagnostics. Methane agreement alone cannot calibrate cage rigidity,
nearby contacts or the supported handle; a constant methane correction must not be
applied to the 53-atom path without independent transfer evidence.

**Acceptance:** discrepancies for a defined calibration set and its reference
uncertainties meet the declared use-specific budget, and transfer tests cover the
chemical environments used in the prediction. A set used to tune a method must be
distinguished from held-out validation. The supported candidate stays outside the
validated domain until its relevant transfer tests pass.

### G5 — Competing chemistry and positional selectivity

Build an explicit channel inventory: other accessible C–H sites; competing C–C
addition or exchange; H back-transfer; premature tip passivation; cage/handle bond
damage; reconstruction; and relevant electronic-state changes. Identify how each
channel was sought and what could remain undiscovered. Use A1's site energies to
prioritize searches, not to exclude kinetically accessible competitors.

For relevant poses, compute competing pathways using the same G1–G4 standards.
Evaluate the neighborhood of positions/orientations from the declared mechanical
uncertainty, including approach and withdrawal. The nearest-site boundary margin
is only a geometric screening descriptor. Neither an exothermic target reaction nor a
larger bond dissociation energy elsewhere establishes kinetic exclusivity.

**Acceptance:** the target outcome meets the declared selectivity criterion after
including competing-channel uncertainty in a justified kinetic/operating model.
Until then report channel-specific evidence and unsearched alternatives; do not
convert a geometric margin into a nanometer/Å positioning specification.

### G6 — Handle, boundary and mechanical convergence

Enlarge target and handle models and move anchors farther from the reactive center
while preserving the intended physical pose. Compare relaxed geometry, electronic
character, interaction/barrier energies, raw support loads and reaction paths.
Rigid-geometry H-affinity convergence tests electronic substitution only.

Replace perfectly fixed anchors with a physically specified finite-compliance
model when predicting a mounted device. Determine pose-dependent force and torque
response, cross-coupling and stability; converge the model size/boundary treatment
for these observables. Distinguish displacement-controlled and force-controlled
conditions, and account for external work under the selected control protocol.
Fixed anchors can suppress rearrangements or instabilities available to a real
mount. Loading history matters; dynamic-force theory illustrates why one static
force is not a time-independent failure threshold, but supplies no carbon-tool
parameters here. [Evans and Ritchie, 1997](https://pubmed.ncbi.nlm.nih.gov/9083660/)

**Acceptance:** changes under boundary/model refinement and physical pose
uncertainty meet the declared budgets, using independently justified support
properties. No finite-stiffness or mount-validation evidence currently closes this
gate. A2's electronic handle survey must not be promoted to that role.

### G7 — Temperature, dynamics and a complete cycle

Define the statistical ensemble consistent with the mount and environment.
Evaluate zero-point contributions, soft/hindered/anharmonic motion and hydrogen
nuclear quantum effects at the conditions of interest. Do not apply free ideal-gas
translational/rotational entropies to attached fragments. ASE offers distinct
thermochemistry models; selecting one does not establish that it describes this
device. [ASE thermochemistry documentation](https://docs.ase-lib.org/ase/thermochemistry/thermochemistry.html)

A rate model needs a justified dividing surface, sampling and treatment of
recrossing/tunnelling where relevant. If constant independent irreversible channel
hazards `k_i` are defensible during a fixed dwell, then the elementary model gives
`P_target(t) = (k_target/K) × [1 − exp(−Kt)]`, `K = Σ_i k_i`, and no-event
probability `exp(−Kt)` (with the zero-rate limit taken continuously). This formula
is a conditional model definition, not a prediction from our saved energies.
Reversibility, moving poses, memory or changing states require a different model.

Specify and assess approach, dwell, withdrawal, product survival and regeneration,
including the destination of removed H and tool wear. Current fixed-anchor NEB
represents a path at one imposed pose; moving anchors through a cycle requires
additional controlled-pose calculations or an appropriate dynamics implementation.

**Acceptance:** a converged, calibrated kinetic/dynamical treatment covers the
declared conditions and full cycle, includes competing outcomes/no event, and
meets the error and confidence target. Repeated-cycle reliability cannot be
obtained by assuming independent cycles without checking regeneration and damage.
No operating temperature, rate or success probability is currently established.

### G8 — Experimental validation of the stated operation

Specify an actual substrate, apex composition and mount that can be prepared and
characterized. Predeclare observables that distinguish site-specific H removal,
tip hydrogen capture, collateral damage and tip changes. Use appropriate spatial,
chemical and, where useful, isotope evidence; include no-tip/nonreactive-tip,
position and approach/retraction controls. Measure uncertainty in position,
temperature, timing and detection/classification errors.

Compare out-of-sample outcomes and their distributions to predictions, including
failed and ambiguous trials. Predeclare sample size, uncertainty method and
replication across preparations; zero observed failures does not establish zero
failure probability. Verify regeneration/repeated operation independently of a
single transfer event. Atomic H manipulation demonstrated on silicon is useful
experimental precedent, but it does not validate this ethynyl/diamondoid design.
[Huff et al., Atomic White-Out, 2017](https://arxiv.org/abs/1706.05287)

**Acceptance:** measured chemical identities, operating behavior and uncertainty
satisfy the declared claim under the tested conditions. Extrapolation to other
tips, surfaces or assembly operations remains a new validation task.

## 6. Owners and executable next steps

Assignments below identify existing evidence owners, not permission to expand
their scopes. [The roster](../coordination/ROSTER.md) and addressed C1 assignments
govern writes. Unassigned science needs an explicit new scope from C1.

| Need | Existing lane / boundary |
|---|---|
| Source/reference forensics | Scientific intake; `SOURCE_NOTES.md` and `data/validation/si-energy-reproduction/` |
| Reference stationary search | S1; `research/reference-saddle/`; Q2 acceptance audit in `research/saddle-audit/` |
| Fixed-geometry guess surveys | S2 production API, S3 independent review; state identity remains a separate open scientific task |
| Saved force/Hessian reconstruction and mode comparison | E2 in `research/evidence-audit/`; N1 in `research/mode-comparison/`; arithmetic and convergence support only |
| Site and reduced-handle observations | A1 in `research/site-selectivity/`; A2 in `research/candidate-feasibility/` |
| Resumable characterization | H1 research prototype and H2 review; not an integrated production-resume guarantee |
| Cost / possible GPU work | C2 planning and G2 offline protocol; no verified GPU adapter, parity or speedup |
| Integration and assignment | C1; current P2 author and P3 reviewer handle this protocol only |
| IRC implementation; continuous state tracking; finite mount/bath; cycle rates; experiments | No dedicated completed implementation/evidence lane identified; C1 must scope these explicitly |

### Ordered decision path

1. **Freeze the claim and existing evidence.** Fill the decision inputs for the next
   intended use, retain source hashes and failed attempts, and review S1/A1/A2
   terminal records when their owners publish them. Do not launch duplicate work.
2. **Resolve the small calibration surface.** Use S1's candidate structures and
   state findings to decide whether to follow separate branches, refine the method,
   or obtain a more suitable reference. If G1 is unresolved, state-dependent G3/G4
   conclusions remain conditional; do not select a barrier solely by convergence.
3. **Converge and connect the calibration reaction.** Complete G2/G3, then assess
   G4 against justified references. Missing IRC/state evidence requires a scoped
   implementation or an independently documented external calculation.
4. **Test transfer before expensive ranking.** Complete site/reduced-handle
   observations, define mechanically meaningful candidate boundaries and competing
   channels, then obtain state-assessed relaxed supported endpoints. Stop or revise
   the candidate if the intended products/basins do not survive.
5. **Calculate supported pathways under a bounded plan.** At explicitly chosen
   poses, converge target and competitor paths, Hessians and connectivity. The
   141-coordinate central-difference Hessian routine needs `1 + 2 × 141 = 283`
   force requests. The current workflow wrapper makes 284 force-array requests
   with a cached baseline, giving 283 new quantum evaluations for a complete first
   attempt; early rejection and retries change the count. See the
   [cost assumptions](../research/compute-planning/README.md). Reuse compatible
   evidence only with provenance checks.
6. **Extend to a physical operation.** Close G6/G7 and then G8 for a declared
   operating envelope. A failed criterion triggers redesign or a narrower claim;
   missing evidence remains inconclusive, not a favorable score.

### Commands available now

Run from the repository with the documented environment installed. These first
commands inspect or prepare evidence; they do not execute electronic calculations.
Use a new output directory for each plan or bundle.

```sh
nanodesign check examples/h-abstraction/design.json
nanodesign campaign-report examples/pose-campaign
nanodesign state-scan-create examples/h-abstraction/initial.xyz \
  --settings examples/h-abstraction/design.json \
  --guesses minao atom 1e huckel --out runs/validation-tip-guesses
nanodesign state-scan-report runs/validation-tip-guesses
nanodesign bundle-create data/validation/paired-ccpvdz-atom \
  --out runs/validation-reference-bundle
nanodesign bundle-verify runs/validation-reference-bundle
```

Only after C1/owner scheduling and a declared resource bound, this separate command
starts **one real energy-plus-gradient job**:

```sh
nanodesign state-scan-run runs/validation-tip-guesses --max-jobs 1
```

It does not optimize a geometry, evaluate orbital stability or prove state identity.
Do not run it merely to follow this document during the current host-overload
advisory. [State-scan behavior](../research/state-scan/README.md),
[bundle integrity limits](EVIDENCE_BUNDLES.md), [calculation accuracy](ACCURACY.md)
and [campaign behavior](CAMPAIGNS.md) describe the implemented contracts.
There is presently no single command that executes and certifies all eight gates.
Software tests, evidence transport and numerical replay remain necessary support;
chemical and operating validation require the additional observations above.
