# GPU readiness and numerical comparison protocol

G1 support for the molecular-design project, checked 2026-09-16. This directory
adds a working environment inspector, an exact-input benchmark-plan builder, and
an offline CPU/GPU record checker. It contains **no GPU solver adapter**. None of
the preparation or local checks executes quantum chemistry, installs software,
starts cloud machines, or spends money.

**G2 update:** the default plan/compare commands now use version 2. They verify
the actual snapshot files and require resolved basis, grid, pruning and float64
metadata before accepting a recorded comparison. See [the v2 contract](PROTOCOL_V2.md).
The original v1 plan, preflight and verification records are preserved unchanged.

## Do we need a GPU?

We can continue calibration and small-system development on the existing CPU.
Larger pose and reaction-path studies justify testing acceleration. The archived
53-atom evaluation took about 701.5 seconds with density fitting and 2,044.8
seconds with direct DFT. They overlapped other workloads; their ratio is not a
controlled speedup. The independent [C2 planner](../compute-planning/README.md)
projects about 19.54 DF days or 56.97 direct days for one explicit scenario with
200 updates in each of four path stages. That is an assumption-dependent
projection, not an expected completion date or convergence guarantee.

[Local preflight](local-preflight.json) measured Darwin/arm64, eight logical
CPUs and 24 GiB memory. The optional subprocess probe requested two PySCF threads
and observed one. CuPy and GPU4PySCF were absent. The current production adapter
still runs CPU PySCF.

The supported GPU4PySCF route is **NVIDIA CUDA**; the Mac's Apple/Metal GPU is a
different platform. Its compiled packages document compute capability 7.0 or
later. Use a compatible NVIDIA system for an initial measured trial before
choosing or buying hardware. Neither minimum VRAM nor the speedup of this
particular workload has been measured. [GPU4PySCF project](https://github.com/pyscf/gpu4pyscf)

The official support table includes density-fitted hybrid and unrestricted DFT
energies and gradients, relevant to this project. Dispersion remains a CPU
contribution. These advertised capabilities need to be checked on the exact
settings and derivative options used here. [PySCF GPU documentation](https://pyscf.org/user/gpu.html)

Use a separate environment with matching GPU4PySCF, CUDA, CuPy and cuTENSOR
versions; the project warns against arbitrary combinations. The current
upstream CUDA-12 requirements include PySCF 2.14.0, CuPy 13.4.1 and cuTENSOR
2.2.0. Those are a dated reference, not an installation lock or proof that a
particular Python/platform wheel is available. This tool records the versions
actually installed and does not install them. [Upstream environment requirements](https://github.com/pyscf/gpu4pyscf/blob/master/requirements.txt)

## Run the inspector

From the repository root:

```sh
.venv/bin/python research/gpu-readiness/gpu_readiness.py inspect
.venv/bin/python research/gpu-readiness/gpu_readiness.py inspect --probe-runtime
```

The default inventories package metadata and available NVIDIA management output.
`--probe-runtime` additionally imports libraries and queries device visibility in
a child process with a 20-second timeout. It may initialize a CUDA context, but
does not solve a molecule. The PySCF thread probe runs in that separate process
and restores its setting; it does not affect other agents' running jobs.

`prerequisites_observed` means only that the required package metadata, imports
and a supported visible device were observed. It does not establish numerical
parity, method support, enough device memory, useful speed, or production
integration. A missing `nvidia-smi` alone is not treated as proof of GPU absence.

Add `--output NEW_FILE.json` to preserve a report; existing files are rejected.

## Prepare a bounded benchmark

```sh
.venv/bin/python research/gpu-readiness/gpu_readiness.py plan --output /tmp/new-gpu-protocol
```

The shipped [v2 protocol](protocol-v2/plan.json) contains five cases using the exact
bytes of existing geometries, with per-file hashes, original settings records,
recorded changes and reference attribution. It is **planned, not executed**:

| Case | Atoms | Purpose |
|---|---:|---|
| H2, restricted direct DFT without D3 | 2 | Basic energy/force conversion and direct route |
| Methane, restricted DF + D3 | 5 | Closed-shell density fitting and dispersion |
| Ethynyl, unrestricted DF + D3, `minao` | 3 | Open-shell energy, force and spin diagnostics |
| Ethynyl, unrestricted DF + D3, `atom` | 3 | Explicit starting-guess propagation |
| Nominal reference transition geometry, unrestricted DF + D3 | 8 | A difficult supplied geometry; not a verified saddle |

One CPU and one GPU evaluation per case gives **ten planned energy-plus-force
calls**, serially, with no implicit retries or optimization. The two ethynyl
guesses are separate checks, not automatic state selection. Some settings differ
from the archived source calculations (for example density fitting); every change
is recorded. The old calculated energies are not substituted for new matched CPU
controls.

Only after the small checks pass should a separate bundle be prepared with
`--include-candidate`. It adds one pair on the actual 53-atom starting structure,
retaining all atom identities and the six anchors as context. It increases the
planned count to twelve but still launches nothing. Compare **raw forces on all
atoms**, including anchors, because anchor loads matter to mounting the tool.

## Capture matched CPU and GPU results

A future reviewed runner must produce one JSON record per backend and case.
The production adapter currently has no backend selector; creating a plan does
not change that. Required record fields are:

| Field | Contract |
|---|---|
| `schema_version`, `case_id`, `geometry_sha256`, `plan_sha256` | Version 2, exact planned case/coordinate identity, and digest of the complete plan bytes |
| `settings` | Complete exact case settings, including charge, spin, basis, XC, dispersion, grid, tolerances, DF and guess |
| `backend`, `execution_device` | `pyscf_cpu` / `cpu` or `gpu4pyscf` / `nvidia_cuda`; CPU fallback is not a GPU result |
| `scf_initial_guess` | Matches planned settings; no implicit fallback |
| `scf_converged`, `gradient_completed` | Both explicitly true |
| `energy_hartree`, `forces_ev_per_angstrom`, `s2` | Finite energy, N×3 raw forces in original atom order, determinant spin diagnostic |
| `grid_response`, `auxiliary_basis_response` | Grid response true; auxiliary response agrees with DF use |
| `dispersion_evaluations` | Exactly one if D3 is requested, otherwise zero |
| `versions` | Nonempty installed-version strings; GPU also identifies GPU4PySCF and CuPy |
| `gpu_synchronized` | Explicitly true for GPU; required operations completed before timing/output |
| `resolved_numerics` | Actual precision, expanded orbital/auxiliary basis identities and canonical paired-grid identity; see [v2 schema](PROTOCOL_V2.md) |

The CPU implementation explicitly disables native dispersion and adds
simple-dftd3 once. Its gradient includes grid response and, for DF, auxiliary
basis response. The GPU implementation must demonstrate those operations are
honored, not merely assign unused attributes. GPU arrays must be transferred to
host arrays and converted with the same force sign and unit constants as the
current adapter. Check finite differences on a small displaced molecule before
extending to geometry optimization.

Record synchronized total energy-plus-force wall time, initialization separately
when available, actual CPU threads, device name, memory usage, software versions
and concurrent load alongside these fields. A single cold pair is a correctness
screen, not a trustworthy benchmark of sustained throughput. The checker does
not infer a speedup or extrapolate memory from atom count.

## Compare records

```sh
.venv/bin/python research/gpu-readiness/equivalence.py compare \
  --plan research/gpu-readiness/protocol-v2/plan.json \
  --cpu /path/to/cpu-case.json --gpu /path/to/gpu-case.json
```

Provisional comparison thresholds are 1e-6 Hartree for total energy, 1e-4 eV/Å
for the maximum force-component difference, and 1e-5 for determinant S². These
are engineering gates for the same input and approximation, not calibrated
chemical-error bounds. Differences require investigation of grids, auxiliary
bases, derivative terms, library behavior and electronic solutions. Agreement in
energy alone is insufficient, and matching S² does not prove state identity.

The checker validates the supplied records and numerical differences. It cannot
independently attest that their claimed hardware was used; the future runner
must capture that provenance. Missing data never becomes a passing score.
`numerical_parity_passed` applies only to that supplied pair; scientific
validation remains false.

The default CLI additionally reads and verifies every declared coordinate,
source-record and attribution file, checks geometry/order and requested-setting
changes, and binds each execution record to the complete plan hash. Historical
v1 records can be inspected explicitly with `--records-only`; that mode always
leaves `input_files_verified`, `resolved_numerics_verified` and
`parity_evidence_accepted` false, even if their supplied numbers agree.

## What is still needed to make the design model useful?

1. **Calibrate one elementary operation.** Locate stationary reactants and a
   first-order saddle candidate, verify the path connects the intended states,
   compare alternative electronic solutions, and check competing reaction sites.
   Existing S1/A1/source studies address parts of this work; they do not yet
   validate a complete assembly operation.
2. **Establish the model's accuracy range.** Quantify basis, functional,
   handle/substrate-size and mechanical-boundary sensitivity against suitable
   reference calculations or measurements. A2 owns the reduced-handle fidelity
   experiments. A smaller model cannot simply inherit a larger model's accuracy.
3. **Build a usable computational search.** Integrate a tested accelerated
   adapter or bounded parallel CPU execution, then use calibrated reaction
   objectives to compare poses. The current campaign runner is serial and
   preserves attempts; it does not discover arbitrary tools automatically.
4. **Close the physical operating cycle.** Approach, selective transfer,
   withdrawal and tool regeneration need compatible structures, forces and
   chemistry. Fabrication, mounting and operation require experimental evidence;
   a successful simulation is not a manufactured device.

The near-term milestone is a defensible model of a particular reaction with
stated accuracy and failure conditions. GPU access helps us test more cases;
it does not supply the missing physical validation.

## Verification

```sh
.venv/bin/python -m pytest research/gpu-readiness -q
```

The G1 checkpoint passed **98 offline tests** locally in 3.82 seconds. Tests use real archived
inputs for snapshot checks and explicitly synthetic records for comparison/error
cases. They run no quantum calculations. Original archives and other agents'
source files remain unchanged. `verification.json` records the command, output
and source hashes. Actual GPU execution, numerical parity and acceleration remain
untested because this host lacks the required CUDA environment.
G2 passed **229 offline tests** in 2.96 seconds, including snapshot, integration
and resolved-numerics checks. Its command, output, source hashes and verification
of 12 bundled files are recorded separately in `verification-g2.json`; the
original verification record stays intact.
