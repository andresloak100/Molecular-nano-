# P3 independent validation-protocol review

Review owner: `codex-712e`. Protocol author: P2 / `codex-f522`.
Exact scope is a read-only scientific review of `docs/VALIDATION_PROTOCOL.md`.
No protocol/core edits, calculations, new sessions, or broad test runs are part
of P3. This directory contains review notes, not new chemical validation.

## Completed review — 2026-09-16

**Disposition: no unresolved material findings in the reviewed protocol.**
This is acceptance of a research protocol's definitions and evidence requirements,
not acceptance of the proposed tool or of any scientific gate as passed.

- Reviewed document: [`docs/VALIDATION_PROTOCOL.md`](../../docs/VALIDATION_PROTOCOL.md).
- Final SHA-256: `1078f5dd6936899623086a4b5fad8cd1d1602f26ef832c74c8f53bfb383fd5e1`.
- Initial draft SHA-256: `c6f28355544852495eebb05101b81ba811d0d84feab34194f705898d2a6aef62`.
- [`review.json`](review.json) binds the reviewed document and source files.

Parent review covered electronic-state evidence, numerical/basis convergence,
reference observables, model transfer and selectivity. The existing internal
helper independently reviewed saddle/connectivity, mechanical and thermal gates,
and current owner/command contracts. No commands in the protocol were executed.
P2 separately checked local links and archived numerical values; P3 inspected
the interpretation, relevant source contracts and selected primary sources.

The following implementation facts establish the evidence boundary:

| Evidence currently provided | Limit relevant to the protocol |
|---|---|
| `nanodesign/quantum.py` records convergence, energies, forces, spin diagnostics and call IDs | It sets `stability_checked=False`; it does not preserve orbital/density state descriptors sufficient for branch tracking. |
| `nanodesign/state_scan.py` compares explicit starting guesses at a fixed geometry | It leaves ground-state and electronic-state identity verification false. Agreement among a finite set of guesses is not a state certificate. |
| `nanodesign/stationary.py` records a finite-difference Hessian and mode classifications | Transition-state verification, IRC and displacement-step convergence remain false. Curvature is conditional on frozen anchors. |
| `nanodesign/design.py` screens endpoint identity by distances/topology | Geometric endpoint screening does not show that the saddle's two downhill branches reach the intended basins. |
| `nanodesign/workflow.py` reports a candidate fixed-anchor electronic path barrier and holding forces | These do not establish a free-energy barrier, finite mount response, operating rate or error probability. |
| E2 reconstruction and N1 mode comparisons independently inspect saved force/mode evidence | Numerical reproducibility does not establish that every underlying force evaluated the same electronic surface. |

## Material contributions and resolution

Electronic-state continuity must be assessed for **every positive and negative
finite-difference force displacement**, as well as reaction-path/connectivity
images. A reconstructed or step-stable Hessian can remain scientifically
ambiguous when the force calls switch electronic branches. P2 incorporated this
requirement in G1, including saved wavefunction/density evidence and an explicit
statement that the production implementation does not yet supply it. Verified in
the reviewed draft.

The first draft called A1's 2.494963 Å value a static distance difference. The
actual `nearest_hydrogen_margin.definition` in
[`census.json`](../site-selectivity/evidence/stage0-site-census-r3/census.json)
defines the minimum apex displacement to the perpendicular-bisector boundary
between intended and competing H positions, with target coordinates fixed. P2
corrected the table and G5. The revised text retains the necessary distinction
between this geometric boundary and an operating tolerance or chemical
selectivity. Verified by reading both definitions after the author's edit.

C1's final publication review also requested the distinction between 283
Hessian-routine force requests and 284 wrapper requests (283 new quantum calls
when the baseline cache is reused for a complete first attempt). P2 added it.
P2 removed three untracked live S1 result links and their mutable observations,
requiring an owner-released snapshot instead. These final changes were read and
the request count checked against `workflow.py` and `stationary.py`; no calculation
was launched.

## Remaining scientific evidence is explicitly open

The protocol does not present these requirements as executable completed features:
continuous electronic-state characterization; numerically refined saddle and
endpoint connectivity; transferable chemical calibration; a physically specified
finite mount and environment; competing-channel and complete-cycle kinetics; and
experimental outcome validation. Unset use-specific accuracy/reliability targets
remain inconclusive. Source code tests, state-guess surveys, evidence bundles,
Hessian reconstruction and mode comparison do not close those gaps.

G2 requires both matched-geometry comparisons and reoptimization under final
settings; G4 separates fitted calibration from held-out validation. G3 requires
both downhill connections and does not turn a failed search into barrierlessness.
G6/G7 separate fixed constraints from finite compliance and distinguish electronic
energies from thermal rates, competing outcomes and repeated-cycle reliability.

## Primary sources inspected

- [PySCF SCF documentation](https://pyscf.org/user/scf.html): initial-guess choices
  and orbital-stability analysis; convergence and local stability are different
  properties.
- [ORCA 6.1 IRC documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/irc.html):
  downhill connectivity from a transition-state candidate, using the same
  method and basis as the optimization/frequency calculation.
- [Temelso et al., 2006](https://doi.org/10.1021/jp061821e),
  [author-hosted paper](https://vergil.chemistry.gatech.edu/static/pdfs/temelso_2006_11160.pdf):
  constrained/symmetry-dependent stationary structures and additional imaginary
  modes; separated-reactant electronic energies and thermal kinetic quantities
  require distinct interpretations.
- [PySCF MCSCF documentation](https://pyscf.org/user/mcscf.html) and
  [DFT documentation](https://pyscf.org/user/dft.html): active-space selection and
  state analysis require explicit choices; grid level is only part of the
  numerical integration specification.
- [Evans and Ritchie, 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC1184350/?page=1):
  loading-rate dependence in the stated adhesion model; no carbon-tool parameters
  are inferred from this precedent.
- [Huff et al., 2017](https://arxiv.org/abs/1706.05287): experimental hydrogen
  manipulation on silicon is a distinct system, as the protocol states.

No new quantum calculation, test campaign, source-evidence mutation, new agent
session, production edit or commit was performed. Existing source limitations
were retained. Later changes to the protocol or its scientific evidence need a
new dated assessment; this review does not automatically cover them.
