# A2 — Is the 53-atom candidate reaction path computable on this hardware?

Lane `andresarriaga-a8`. Draft of 2026-09-16; sections marked **pending** are
not yet measured and are named rather than estimated.

**Verdict: not at default settings on this host as it stands — but the margin
is contention, not chemistry, and that changes what to do about it.** Using the
archived timings at face value, a path on *one* pose costs between 4.3 and 19.5
days of uninterrupted single-core compute. However, those timings are wall clock
taken under unrecorded contention, and this host's measured contention factor
reaches 10.4x. If the archived runs carried anything like that factor, the same
central scenario is **20 hours, not 8.6 days**. The verdict genuinely hangs on a
number nobody has measured yet.

Note what is *not* blocked: the nine-pose campaign in `examples/pose-campaign`
is configured for `stage: "singlepoint"`, one evaluation per pose, roughly
1.8 hours in total. It is affordable today. The cliff is the path stage, which
costs 530-2400 evaluations where the campaign costs one.

So the deliverable splits in two. The evaluation *count* is solid and is the
real structural problem. The per-evaluation *cost* is not yet established to
better than an order of magnitude, and establishing it is cheap. A reduced
model brings one pose into range under either reading; the reduction is
defensible on electronics and undemonstrated on mechanics.

**What this verdict is not.** It says nothing about whether the reaction works.
An affordable path would not be evidence that the tool abstracts the right
hydrogen, and an unaffordable one is not evidence against it. Site selectivity
belongs to A1 and locating a genuine saddle to S1.

**Two findings this lane produced by asking what the cost was *for*.** Both are
written up separately in `PHYSICS_COMPLETENESS.md` and cost no quantum time.

*Positional control has a real constraint with a shape.* The apex may move
2.495 Å in any direction before a different hydrogen is nearest, which needs
4.8 N/m of mount stiffness at 298 K for a 10⁻¹⁵ error rate and 1.25 N/m at
77 K. Set against the buildable range, a long handle worked in bending (~2 N/m)
misses by nine orders of magnitude at room temperature and passes cold. The
design rule is **mount stiffly and short, or operate cold** — conditional
throughout on the nearest hydrogen being the one that reacts, which is
unestablished.

*Whether a classical barrier can become a rate is unresolved, and one number
decides it.* The crossover temperature is 0.229 × ω*[cm⁻¹]. On the only
measured imaginary frequency available (259i) tunneling is a 7% correction; on
the value implied by the experimental activation energy (~1648i) it is the
mechanism. S1's saddle produces the first ω* on a functional's own surface, and
it should be reported as a headline beside the barrier.

---

## 1. What one evaluation costs

From the archived 53-atom records, reused rather than recomputed:

| | direct | density-fitted |
|---|---|---|
| Wall clock, one energy+gradient | 2044.8 s | 701.4 s |
| SCF / gradient split | not instrumented | 477.2 s / 224.1 s |
| Basis functions | 463 | 463 |
| SCF cycles | 14 | 14 |
| Effective threads | not recorded | 1 (`threads_honored` false) |

System: C22H31, 53 atoms, neutral doublet, UKS, PBE0-D3(BJ)/def2-SVP.

Three caveats that travel with these numbers.

The direct run predates the thread diagnostics added in `121c30f`, so its
effective thread count was never recorded. One thread is the expected value
because this PySCF build has no OpenMP, but that is an inference from the
build, not a measurement of that run.

Both runs were concurrent with other jobs, and the load average at the time was
not recorded. Each figure is therefore an upper bound on its uncontended cost
by an unknown factor, and the 2.9x direct-over-density-fitted ratio is an
observation, not a speedup measurement.

A peer lane has since measured this host's contention factor directly: an
identical five-species benchmark took 6.1 s on a quiet host and 63.4 s at load
201, a factor of 10.4. That is a property of the box, not of the method. It
means the archived figures cannot be corrected without the load at the time,
and it is the reason this lane measures CPU time rather than wall clock for
everything it runs itself.

**Pending:** a CPU-time measurement of the same 53-atom evaluation, which
would replace the uncalibrated 701.4 s with a transferable number. This is the
single most valuable outstanding measurement in the lane.

## 2. What the path costs

`workflow.run(stage="path")` is four stages, read from the code rather than
assumed: two endpoint relaxations, each up to `steps` evaluations; a
preliminary NEB and then a climbing-image NEB, each up to `steps` optimizer
steps costing `images - 2` evaluations apiece. Images are evaluated serially in
one process and no SCF density is carried between changed geometries, so every
evaluation pays full price.

At the default 7 images and 200 steps the ceiling is **2400 evaluations**.

The ceiling is what the code permits, not what it should need. Three scenarios
with declared step counts, using the archived per-evaluation costs:

| Scenario | Endpoint / pre-NEB / CI-NEB steps | Evaluations | Density-fitted | Direct |
|---|---|---|---|---|
| Optimistic | 40 / 30 / 60 | 530 | **4.3 d** | 12.5 d |
| Central | 80 / 60 / 120 | 1060 | **8.6 d** | 25.1 d |
| Ceiling | 200 / 200 / 200 | 2400 | **19.5 d** | 56.8 d |

The step counts are declared assumptions. Nothing in this repository has ever
run a NEB on this system, so there is no calibration for them; that is a real
weakness of the projection and the reason the answer is a range.

Even the optimistic scenario exceeds four days for a single pose on a machine
that is shared.

**A distinction worth getting right, because it changes what is actually
blocked.** The campaign in `examples/pose-campaign` has nine poses, but its
controls specify `stage: "singlepoint"`, not `"path"`. As configured it is one
evaluation per pose — about 1.8 hours in total at the archived density-fitted
cost, and affordable today. Nothing about the campaign is blocked by this
report.

What is infeasible is switching that campaign, or any single pose, to
`stage: "path"`. Nine poses at the central scenario would be roughly 77 days
at the pessimistic end of the contention range and about 7.5 days at the
optimistic end. The gap between "run the campaign" and "run a path" is three
orders of magnitude in evaluations, and it is the whole of the problem.

**How much of that is the machine rather than the method.** The archived
per-evaluation cost is wall clock under contention of unrecorded magnitude, so
the projection inherits that uncertainty. Bounding it at both ends, for the
central scenario, density-fitted:

| Assumed contention factor on the archived runs | Implied cost per evaluation | Central scenario |
|---|---|---|
| 1.0 (the archived number is the true cost) | 701.4 s | 8.6 days |
| 10.4 (this host's measured factor at load 201) | 67.5 s | **20 hours** |

The true factor lies somewhere between, and cannot be recovered because the
load during those runs was not recorded. This is a bound, not a measurement.
It is also the whole ballgame: one end says wait for a cluster, the other says
run it overnight. One CPU-time measurement of a single 53-atom energy and
gradient would collapse the range to a number, which is why it is the top
outstanding item in this lane.

Note what does *not* change across that range: the evaluation count. Nine poses
at 530-2400 serial evaluations each is a structural cost that a faster machine
rescales but does not remove.

**Process parallelism does not rescue it.** The five interior images are
independent within an optimizer step and could run in five processes, which
would cut the band term by about fivefold and leave the endpoint relaxations
untouched — roughly 8.6 days to 2.5 in the central scenario. Three problems:
ASE's NEB here evaluates images serially and the code would have to be changed;
this host's eight cores are already shared by other lanes; and the machine is
currently memory-bound rather than core-bound, so adding processes would not
buy the ideal factor. It is a real lever, but it is a code change plus a quiet
machine, not a setting.

## 3. Does a smaller tool handle preserve the chemistry?

The handle is 27 of the 53 atoms and most of the cost. Shrinking it changes two
separable things, and conflating them is the trap.

### Electronic substituent effect — measured

Tip H-affinity `A(R) = E(R-CC-H) - E(R-CC*) - E(H)` at PBE0-D3(BJ)/def2-SVP
with density fitting, rigid geometries truncated from the candidate's own 3.6 Å
pose, four-guess SCF scan on every open-shell species:

| Handle | Tool fragment | A(R), kcal/mol |
|---|---|---|
| Hydrogen | HC≡C· | -136.72 |
| Methyl | CH₃-C≡C· | -136.52 |
| Adamantyl | pending | pending |

**Sanity check against known chemistry.** The hydrogen-handle case is just
acetylene: `A(H) = -136.72 kcal/mol` is the negative of the acetylene C-H bond
dissociation energy, so this calculation puts that bond at 136.7 kcal/mol.
Acetylene's C-H bond is experimentally about 133 kcal/mol as an enthalpy at
298 K. These numbers are bare electronic energies with no zero-point or
thermal correction, and zero-point energy lowers a C-H dissociation energy by
roughly 4-5 kcal/mol, so an electronic value a few kcal/mol above the
experimental enthalpy is what agreement looks like here. It is a coarse check
— rigid unrelaxed geometries, a small basis, no ZPE computed — but it is the
first number in this lane anchored to something outside the repository, and it
lands where it should rather than somewhere absurd.

Hydrogen versus methyl differ by **0.19 kcal/mol**. For scale, this project's
own CCSD(T) calibration puts tertiary abstraction 7.4 kcal/mol more exothermic
than primary, and published kinetic selectivity in adamantane corresponds to
roughly 1.0-1.3 kcal/mol. A 0.19 kcal/mol handle effect is well below both.

### Mechanical and steric fidelity — not measured, and not implied

The handle also sets the mount's stiffness, its mass, its sterics against the
substrate, and where the anchors sit. The reduced models fix a single terminal
handle atom where the candidate fixes three cage carbons: that is a *different*
boundary condition, not a scaled one. Nothing here bears on it, and a small
H-affinity spread must not be read as licence to treat the reduced model as the
same mounted tool. A relayed finding from another lane sharpens this: the
nominal pose sits at only +0.200 Å van der Waals clearance, so steric detail
near the tip is not a free parameter.

The honest status of the reduction is therefore: **defensible for screening the
electronics of the bond being broken, undemonstrated as a mechanical
substitute.**

## 4. The cheap knobs

**Density fitting.** Archived, at one geometry: maximum force-component
deviation 4.7e-4 eV/Å and RMS 1.8e-4 against the 0.03 eV/Å optimization
tolerance, with S² agreeing to 1e-6. That is roughly 60x inside tolerance, and
being a ratio it is immune to the contention that spoils the timings. Measured
at the unrelaxed pose only; the surface is flatter near a saddle and the check
does not automatically transfer there.

**Pending:** a controlled same-process direct-versus-fitted CPU-time
comparison, and the basis and force-tolerance knobs.

## 5. Traps this lane hit

The default `minao` SCF guess converged 11.2 kcal/mol **above** the solution
the other three guesses agree on, for the propynyl radical at
PBE0-D3(BJ)/def2-SVP. Two things make this worth recording beyond this lane.

It is DFT, not Hartree-Fock. The repository's trap narrative had treated DFT as
the safe side; two independent lanes have now falsified that, here and at the
published transition structure.

The misleading solution has the *cleaner* spin expectation value: minao gives
S² = 0.7521 against 0.7876 for the correct one, and 0.75 is the ideal doublet
value. Across all three cases this project has scanned, the lower-energy
solution is the more spin-contaminated one every time. Picking the cleanest S²
— the shortcut a reader reaches for once they learn that scanning is expensive
— would have chosen wrong in every case tested. That is a warning against a
heuristic, not a selection rule, and it does not relieve anyone of scanning.

## 6. Recommended protocol

Runnable this week, with its compromises stated.

1. **Run the nine-pose campaign as configured — it is affordable.** Its
   controls say `stage: "singlepoint"`, about 1.8 hours in total density-fitted.
   Do not switch it to `stage: "path"`: that is 530-2400 evaluations per pose
   instead of one.
2. **Use density fitting throughout.** Measured 4.7e-4 eV/Å force deviation
   against a 0.03 tolerance is the best-supported approximation available, and
   it is the one number here that contention cannot touch.
3. **Establish the path on the reduced model first**, at the same pose, and
   treat the result as electronics screening rather than as the mounted tool's
   barrier. Carry the boundary-condition difference explicitly into any
   comparison.
4. **Scan SCF guesses for every open-shell species at every new geometry**, and
   record the selected guess in the evidence. Do not select on S².
5. **Do not cost or run poses tighter than 3.6 Å** without revisiting the van
   der Waals clearance finding first.
6. **Re-measure per-evaluation cost on a quiet host** before quoting any
   wall-clock figure as a hardware verdict, or quote CPU time instead.

On "does this need a GPU cluster": not established, and the numbers do not yet
support saying so. The dominant term is 530 to 2400 serial evaluations, and the
first lever to pull is reducing that count — a smaller model, a better starting
band, parallel images — not a faster evaluation. A GPU would help the
per-evaluation cost and change none of the evaluation count.
