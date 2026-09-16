# Molecular nano

Quantum-chemistry research software for **positional diamondoid mechanosynthesis**: the atomically precise assembly of stiff molecular structures envisioned in Drexler's molecular machinery proposals.

This repository starts with one elementary operation: a supported ethynyl radical removing a selected hydrogen from a diamondoid site. It generates atom-resolved candidate geometries, computes electronic energies and forces, relaxes mechanically constrained structures, and searches reaction paths using climbing-image nudged elastic band (CI-NEB).

**Status: an executable research foundation, not a validated assembler-design system.** It does not design molecular nanomachines perfectly. Numerical convergence establishes that a chosen calculation finished; it does not establish that its chemistry, manufacturing process, or device will work. The software never labels an uncalibrated candidate a validated design.

## What is implemented

- A 53-atom, neutral-doublet starting candidate: an adamantane target and an adamantane-supported ethynyl tip. The same hydrogen and atom identities are preserved in both reaction endpoints.
- Real PySCF restricted/unrestricted DFT and analytical forces, with explicit D3(BJ) dispersion through the official `dftd3` library. Default: PBE0/def2-SVP with D3(BJ), a starting method requiring calibration for this chemistry.
- Fixed distal carbon anchors, adjustable tool separation and lateral offset, endpoint optimization, and CI-NEB with independent electronic calculations for each image.
- Checks for clashes, inconsistent charge/spin, moving anchors, failed SCF/geometry/path convergence, and endpoints that collapse to the same reaction state.
- Reproducible outputs: input geometries, settings, input hashes, forces, spin diagnostics, trajectories, timings, failure records and outstanding validation work.
- Published small-molecule geometry seeds for method calibration, with provenance, original units, and identified source inconsistencies.

The first structure is a **finite diamondoid cluster**, not a converged diamond surface. The initial coordinates are constructed geometry, not optimized coordinates. A reaction path describes hydrogen transfer at one fixed tool pose; approach, withdrawal, regeneration and entire assembly sequences require additional calculations.

## Install and run

Python 3.11 or newer is required. Quantum calculations run locally on the CPU; cost grows rapidly with atom count, basis size, optimization steps and path images.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'

# Build a candidate and inspect its inputs, without running quantum chemistry.
nanodesign candidate --out designs/h-abstraction --separation 3.6 --offset 0.0
nanodesign check designs/h-abstraction/design.json

# Calculate energy and forces at the initial geometry.
nanodesign calculate designs/h-abstraction/design.json \
  --stage singlepoint --out runs/h-abstraction-energy

# Relax one endpoint with the specified mechanical anchors.
nanodesign calculate designs/h-abstraction/design.json \
  --stage relax --state initial --out runs/h-abstraction-relax

# Relax BOTH endpoints, check their identities, then search the reaction path.
# This is substantially more expensive than a single-point calculation.
nanodesign calculate designs/h-abstraction/design.json \
  --stage path --images 7 --steps 200 --fmax 0.03 --out runs/h-abstraction-path

nanodesign audit runs/h-abstraction-path/result.json
pytest -q
```

Each output directory must be new, so a new calculation cannot silently overwrite an earlier result. A stopped or unconverged calculation returns a nonzero exit code. Post-start failures are recorded in `result.json`. `electronic.jsonl` records electronic convergence iterations and gradient stages for inspecting long calculations. Geometry trajectories can be opened with `ase gui runs/.../initial.traj` or any viewer supporting extended XYZ.

Use the supplied `examples/h-abstraction/design.json` directly to inspect the default candidate. Change `quantum`, input coordinates, or `fixed_indices` in a copied design file to investigate another explicitly defined model. All atom indices are zero-based. Charge is in elementary-charge units; `spin` means **2S = Nα − Nβ**, not multiplicity. Coordinates use Å, output energies eV and forces eV/Å.

## What would make this accurate enough to design a tool?

Accuracy needs a target: material and surface, elementary reaction, environment, operating temperature, placement tolerance, competing reactions and acceptable failure probability. A perfect-design guarantee is not scientifically available.

Before accepting a design, establish reaction-specific quantum benchmarks; converge the basis, grid, cluster size and mechanical boundaries; verify electronic states; confirm transition states and connectivity; quantify competing pathways and finite-temperature effects; then validate predicted behavior experimentally. Those are outstanding research tasks, not boxes that this implementation silently checks off.

The program does **not** currently predict assembly error rates, synthesis accessibility, arbitrary mechanosynthesis reactions, tool lifetime or a whole nanofactory. It does not replace electronic structure with a Lennard-Jones animation or use a nonreactive force field to infer bond-making chemistry.

Read [the physical model](docs/MODEL.md), [source notes](SOURCE_NOTES.md), and [reference geometry provenance](data/reference/provenance.json). The reference data have their own stated licensing terms; see [data/reference/README.md](data/reference/README.md).

For the implemented software layers and the CPU/GPU extension boundary, see [architecture](docs/ARCHITECTURE.md).
