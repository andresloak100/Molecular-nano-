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
- A fixed-geometry comparison with a published methane/ethynyl reaction reference, preserving results species by species. Its discrepancy includes method and geometry differences; it is not a measured assembler error rate.
- Small-system coupled-cluster energy references, with an explicit same-geometry, same-basis comparison against DFT. Open-shell UHF/UCCSD(T) is identified separately from the paper's ROHF-based reference.
- Central-difference vibrational characterization of the free coordinates, reporting negative-curvature and unresolved soft modes, numerical Hessian asymmetry and the remaining transition-state checks.
- Bounded, resumable campaigns over explicit tool poses, with preserved input snapshots, attempt history, worker locking and evidence reports. They do not infer a scientifically validated best tool.
- A local 3D workbench for inspecting actual atomic coordinates, measuring distances, switching poses, and reading calculation evidence and its limitations.

The first structure is a **finite diamondoid cluster**, not a converged diamond surface. The initial coordinates are constructed geometry, not optimized coordinates. A reaction path describes hydrogen transfer at one fixed tool pose; approach, withdrawal, regeneration and entire assembly sequences require additional calculations.

## Install and run

Python 3.11 or newer is required. Quantum calculations run locally on the CPU; cost grows rapidly with atom count, basis size, optimization steps and path images.

To inspect the recorded candidate and results without starting calculations:

```bash
python3 workbench/server.py --port 8765
```

Open [the local workbench](http://127.0.0.1:8765). The viewer needs only Python's standard library and a browser. It shows the real 53-atom candidate, nine uncomputed poses, direct/density-fitting results and both coupled-cluster starting-guess comparisons. See [workbench instructions](workbench/README.md) for importing a local campaign.

To generate inputs and run quantum chemistry:

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

# Compare the default method at published small-molecule geometries.
nanodesign benchmark --out runs/methane-reference

# Compare DFT and CCSD(T) on the same coordinates and cc-pVDZ basis.
# Compare explicit HF starting guesses: this reference has multiple solutions.
nanodesign compare-methods --out runs/paired-reference-minao --cc-initial-guess minao
nanodesign compare-methods --out runs/paired-reference-atom --cc-initial-guess atom

# Optional: characterize a converged structure. This is expensive:
# two force evaluations for each free Cartesian coordinate, plus the original.
# This 53-atom candidate has 141 free coordinates (283 evaluations).
nanodesign characterize designs/h-abstraction/design.json \
  --structure runs/h-abstraction-relax/structure.extxyz \
  --max-free-coordinates 141 --out runs/h-abstraction-modes

pytest -q
```

Each output directory must be new, so a new calculation cannot silently overwrite an earlier result. A stopped or unconverged calculation returns a nonzero exit code. `calculate` and `characterize` preserve post-start failures in `result.json`; `electronic.jsonl` records electronic convergence iterations and gradient stages. `benchmark` writes `benchmark.json` and per-species logs under `electronic/`; `compare-methods` writes `method_comparison.json`, with DFT records and logs under `dft/`. Geometry trajectories can be opened with `ase gui runs/.../initial.traj` or any viewer supporting extended XYZ.

Use the supplied `examples/h-abstraction/design.json` directly to inspect the default candidate. Change `quantum`, input coordinates, or `fixed_indices` in a copied design file to investigate another explicitly defined model. All atom indices are zero-based. Charge is in elementary-charge units; `spin` means **2S = Nα − Nβ**, not multiplicity. Coordinates use Å, output energies eV and forces eV/Å.

Set `quantum.scf_initial_guess` to `minao` (the default), `atom`, `1e`, or `huckel` to select one DFT starting guess explicitly. The choice is recorded with every new calculation. It does not scan alternatives or verify an electronic ground state. The paired comparison's `--cc-initial-guess` controls its separate HF/CC reference only.

## What would make this accurate enough to design a tool?

Accuracy needs a target: material and surface, elementary reaction, environment, operating temperature, placement tolerance, competing reactions and acceptable failure probability. A perfect-design guarantee is not scientifically available.

Before accepting a design, establish reaction-specific quantum benchmarks; converge the basis, grid, cluster size and mechanical boundaries; verify electronic states; confirm transition states and connectivity; quantify competing pathways and finite-temperature effects; then validate predicted behavior experimentally. Those are outstanding research tasks, not boxes that this implementation silently checks off.

The program does **not** currently predict assembly error rates, synthesis accessibility, arbitrary mechanosynthesis reactions, tool lifetime or a whole nanofactory. It does not replace electronic structure with a Lennard-Jones animation or use a nonreactive force field to infer bond-making chemistry.

The reference comparison uses the attributed data in this source checkout. When running a separately installed package, supply `--reference-dir /path/to/data/reference`. A reference-geometry energy difference is not a newly optimized DFT activation barrier. Likewise, one negative vibrational mode alone does not verify a reaction transition state: its direction and connectivity still need checking. See [accuracy checks](docs/ACCURACY.md).

Read [the physical model](docs/MODEL.md), [source notes](SOURCE_NOTES.md), and [reference geometry provenance](data/reference/provenance.json). The reference data have their own stated licensing terms; see [data/reference/README.md](data/reference/README.md).

For the implemented software layers and the CPU/GPU extension boundary, see [architecture](docs/ARCHITECTURE.md).

## Evaluate a family of tool poses

Create a finite grid without launching calculations, then run an explicit number of jobs. The default is one serial job per invocation. A geometry optimization or reaction-path job can itself require many quantum evaluations; `--max-jobs` is not a time limit.

```bash
nanodesign campaign-create --out campaigns/pose-study \
  --separations 3.4 3.6 3.8 --offsets -0.2 0.0 0.2 --stage singlepoint
nanodesign campaign-report campaigns/pose-study
nanodesign campaign-run campaigns/pose-study --max-jobs 1
# Continue the remaining jobs; completed work is not repeated.
nanodesign campaign-run campaigns/pose-study --max-jobs 1
```

Alternatively, pass existing design files to `campaign-create`. They must have matching element order, quantum settings, anchor indices and reaction identities. The campaign snapshots their inputs and records hashes. Changing those snapshots or the plan invalidates the campaign. Each retry gets a new directory; `--retry-incomplete` explicitly restarts incomplete work from the original inputs while retaining all earlier evidence. It does not resume an interrupted optimizer trajectory.

Reports preserve plan order and expose energies, force residuals, anchor loads, reaction-basin checks and validation status when available. They do not silently rank an uncalibrated barrier or a failed reaction basin as a useful design. Campaign execution currently supports macOS and Linux. See [campaign interfaces](docs/CAMPAIGNS.md) for integration details.
