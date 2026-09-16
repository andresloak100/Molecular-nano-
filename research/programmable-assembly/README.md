# Programmable assembly: a second, possibly easier route to custom structures

Owner: `molecular-nano-6c` (cloud session). Lane opened 2026-09-16 ~23:40 UTC
on the user's strategic direction: *"Might be easier to build a programmable
assembler/enzyme and direct it to build the custom structure — like directing a
ribosome to build novel proteins."* New files only; no core module or other
lane touched. This is a **proposed second research track**, offered to the lead
for the roster — not a replacement for the positional-mechanosynthesis work the
project is built around.

## The idea, stated honestly

The project's current path is *positional mechanosynthesis*: a rigid tool
placed with sub-angstrom precision makes one bond at a chosen location. The
user's proposal is the other well-known paradigm: a *programmable assembler*
that reads a sequence and builds the corresponding structure, the way the
ribosome reads mRNA and builds an arbitrary protein.

This is not speculative in the way a diamondoid assembler is. Template-directed,
sequence-programmed molecular assembly **already exists and works**, in several
forms with published, reproduced results:

- **The ribosome** — the existence proof. A programmable assembler that builds
  arbitrary sequences from 20 monomers at ~10⁻³–10⁻⁴ error per residue.
- **Engineered / orthogonal ribosomes and flexizymes** — incorporate
  non-canonical monomers, extending the alphabet beyond biology.
- **DNA-templated synthesis** — a nucleic-acid sequence programs the order of
  small-molecule/polymer bond-forming reactions.
- **Sequence-defined synthetic polymers** — artificial molecular machines that
  add monomers in a programmed order along a track.
- **DNA origami / structural DNA nanotechnology** — a sequence programs the
  self-assembly of near-arbitrary 3D shapes at nanometer resolution.

So the strategic question is real and worth putting to the team: is a
programmable-assembler route to the project's target structures more tractable
than positional mechanosynthesis?

## Where the analogy holds, and where it breaks (the honest limit)

The ribosome makes a **linear polymer** that then folds; DNA origami builds 3D
shapes but out of **one structural material** (duplex DNA) held by a uniform,
weak, reversible interaction. Neither does arbitrary **3D covalent** assembly of
a stiff solid. The reason the ribosome is programmable is that every step is the
*same* reaction (peptide-bond formation) between mutually compatible chemistries;
the program only chooses which side chain rides along. A general assembler for
arbitrary covalent 3D structures would need a uniform, programmable, high-yield
bond-forming chemistry that does not yet exist for diamondoid-class targets.

That is the crux, and it is the same crux the mechanosynthesis track keeps
hitting: **selective, high-fidelity single-step chemistry.** The two paradigms
are not as different as they look. Both live or die on how cleanly one intended
reaction can be favored over its competitors. Which means the validated
quantum-chemistry engine this repository already has is the right tool for
*either* path.

## The one quantitative fact that governs feasibility: error propagation

An assembler that adds N blocks, each correct with probability f, produces
correct full-length product with probability **f^N**. This compounds brutally,
and it is why the ribosome spends so much energy on accuracy. Per-step fidelity
is set by the free-energy discrimination ΔΔG‡ between the intended reaction and
its nearest competitor — *exactly* the quantity the site-selectivity lane (A1)
and the quantum backend compute. So the feasibility of a programmable assembler
reduces to chemistry this project already measures.

`assembler_model.py` implements this (transition-state theory + Hopfield/Ninio
kinetic proofreading; 12 passing tests). Room-temperature results:

**Per-step discrimination → yield of a 100-unit structure (no proofreading):**

| ΔΔG‡ (kcal/mol) | per-step fidelity | 100-mer yield |
|---:|---:|---:|
| 2 | 0.967 | 3.5% |
| 3 | 0.994 | 53% |
| 4 | 0.9988 | 89% |
| 5 | 0.99978 | 98% |
| 6 | 0.99996 | 99.6% |

**Discrimination required for 90% yield of an N-unit structure:**

| N | no proofreading | one proofreading stage |
|---:|---:|---:|
| 10 | 2.7 kcal/mol | 1.3 kcal/mol |
| 100 | 4.1 kcal/mol | 2.0 kcal/mol |
| 1000 | 5.4 kcal/mol | 2.7 kcal/mol |

**Model is calibrated against reality:** ~2.7 kcal/mol intrinsic discrimination
with one proofreading stage reproduces the ribosome's ~10⁻⁴ per-residue error
and ~97% yield for a 300-residue protein — the measured biological regime.

### Why this matters to the whole project, not just this lane

Three consequences fall straight out, and each ties to work already underway:

1. **A programmable assembler needs ~4 kcal/mol of clean discrimination per
   step (or ~2 with proofreading) for structures of interesting size.** That is
   the target A1's site-selectivity numbers should be read against. A1 is
   probing the ~2 kcal/mol scale — enough for one step, not for a 100-step
   program without a proofreading mechanism.
2. **The ~5 kcal/mol error that DFT makes on this reaction is *larger* than the
   discrimination window that decides success.** A method that cannot rank
   competing pathways to ~1 kcal/mol cannot design either an assembler or a
   mechanosynthetic sequence. This is the strongest argument yet for the
   reference-reliability and saddle work: accuracy is not academic, it is the
   whole ballgame.
3. **Proofreading is the highest-leverage design choice.** One kinetic-
   proofreading stage roughly halves the required per-step discrimination. Any
   serious programmable-assembler proposal should build in proofreading rather
   than demand heroic single-step selectivity — exactly what evolution did.

## Proposed next steps (for lead/team discussion, not unilateral)

1. **Pick a target structure class and ask which paradigm fits it.** Linear or
   graftable → programmable assembly is plausibly easier. Stiff 3D covalent
   solid → mechanosynthesis, and the assembler analogy mostly fails.
2. **Feed A1's real ΔΔG‡ into this model** to state a concrete yield for a
   concrete target once A1's numbers land. (Hook is ready; it takes one number.)
3. **Model a proofreading step chemically:** what reversible discard pathway,
   driven by what free-energy source, could a synthetic assembler use? This is
   where the quantum engine and this layer would genuinely combine.
4. **Survey the real precedents** (DNA-templated synthesis, sequence-defined
   polymers, orthogonal ribosomes) for the closest existing chemistry to the
   project's targets, and cost a build path from there.

## Scope and safety boundary

This lane models **abstract building blocks and reaction thermodynamics/kinetics
for materials/molecular-structure assembly.** It is deliberately
chemistry-agnostic: `AssemblyProgram` holds opaque labels and carries no
biological-sequence semantics. It is **not** a tool for designing biological
agents, pathways, or sequences, and must not be extended into one. Engineering
real assemblers (especially ribosome/enzyme engineering) carries biosecurity
considerations that belong with the user and appropriate review, not with an
autonomous session; this lane stops at the physics of fidelity and yield.

## Files

- `assembler_model.py` — fidelity / proofreading / yield model. Stdlib only.
- `test_assembler_model.py` — 12 tests, all passing.
