# Electronic snapshots and determinant comparisons

Electronic energies and forces alone cannot show whether successive SCF calls
followed the same electronic solution. This optional evidence layer preserves the
occupied orbitals and their atomic-orbital (AO) metric, so different converged
solutions can be inspected. **It does not identify the physical state, select a
ground state, perform stability analysis or validate a reaction barrier.**

## Capture scope and lifecycle

`nanodesign.electronic_state` accepts an already converged, real, all-electron
molecular RKS or UKS calculation with integer occupations. It does not run SCF,
compute nuclear gradients, choose occupations, or deserialize executable molecular
objects. Calling `get_ovlp()` may evaluate one-electron AO overlap integrals; capture
also copies arrays, evaluates matrix diagnostics and writes data. Those operations
have computational and storage cost even though they do not repeat SCF.

Every snapshot has `phase: converged_scf_only`. The optional backend hook captures
after SCF/spin checks and before gradients. A saved snapshot can therefore survive
a later gradient, dispersion, logging or output failure. Its presence says nothing
about force availability or success of the complete energy-and-force request.
The surrounding calculation's per-call diagnostics and terminal result determine
that separate outcome. SCF energy in the snapshot excludes the separately added
D3 energy; the orbitals describe the Kohn–Sham determinant.

Export is off by default. An explicitly requested export error stops the requested
calculation through the backend's failure handling. Cached energy/force reads reuse
their evidence reference rather than generating another snapshot. Files are
created exclusively: an existing destination, including a symlink, is rejected.
This API does not make a file immune to modification by other programs. An I/O
failure can leave an incomplete file; it is not silently repaired or overwritten.

## API

```python
from nanodesign.electronic_state import (
    capture_snapshot, save_snapshot, read_snapshot, read_snapshot_bytes,
    validate_snapshot, compare_snapshots,
)

# mf is the live, already converged RKS/UKS object from the current calculation.
# settings is the complete dictionary from dataclasses.asdict(QuantumSettings).
snapshot = capture_snapshot(mf, settings=settings, call_id=call_id)
reference = save_snapshot(existing_output_directory / "new-call.json", snapshot)

left = read_snapshot(reference["path"])
right = read_snapshot(other_snapshot_path)
report = compare_snapshots(left, right)
```

`capture_snapshot` additionally accepts `coordinate_frame_id=None`. Supply an
explicit shared identifier only when the caller knows the frames agree; a generic
`solver_cartesian` label does not establish a common laboratory frame for unrelated
calculations. The initial comparator requires identical saved coordinates and
context. It performs no alignment, atom permutation or cross-geometry comparison.

`save_snapshot` returns `path`, `sha256` and `size_bytes`. The parent directory must
already exist. File reads/writes reject symlinks in the file and its parent path,
and readers reject nonregular files. `read_snapshot_bytes(raw)` validates the exact
immutable bytes already captured by a caller, allowing a checksum and parsing to
refer to the same read. Duplicate JSON keys and nonfinite data are rejected. There
is no pickle, HDF5 checkpoint or `gto.loads` importer in this layer.

The default and maximum file size is **16 MiB**; `max_bytes` may lower it. The AO
dimension is limited to **1,024**. These bounds accommodate the candidate's 463 AO
dimension, but a particular serialized snapshot must still fit the byte limit.
Occupied coefficients and the overlap are saved; full virtual-orbital arrays and
redundant dense density matrices are not. Each explicitly captured call creates
one file, so a long path or Hessian campaign can consume substantial disk space.
No automatic deletion, retention policy or campaign-wide byte budget is supplied.

## Schema 1 and the information it preserves

| Field | Meaning |
|---|---|
| `call_id`, `phase`, `producer` | Exact source-call identifier, SCF-only scope and PySCF version/class |
| `settings` | Complete explicit `QuantumSettings` dictionary; missing fields are rejected rather than filled with today's defaults |
| `geometry` | Ordered symbols/nuclear charges, actual solver coordinates in Bohr, recorded Bohr-to-Å conversion, charge/spin, frame and nonperiodic convention |
| `basis` | Cartesian/spherical convention; ordered AO labels; expanded construction basis; Gaussian normalization flag; per-shell centers, angular momenta, kappa, exponents and contraction coefficients |
| `ao_overlap` | Actual same-geometry AO overlap evaluated at capture |
| `source_occupations` | Original RKS spatial `0/2` or UKS separate alpha/beta `0/1` occupation arrays |
| `channels` | Explicit alpha/beta electron counts, occupied source-column indices and occupied coefficient matrices |
| `source_array_dtypes` | Source MO, occupation and overlap dtypes observed before conversion to JSON numbers |
| `scf_energy_hartree` | Converged SCF energy before the external dispersion addition |
| Scientific flags | State identity and ground-state verification remain false; capture itself does not assess orbital stability |

RKS occupied spatial orbitals are used once in each spin channel. UKS channels
remain distinct, including when their summed density is identical to that of
another determinant. Fractional occupations, complex orbitals, ECP/pseudopotential,
periodic, spinor and other reference forms are unsupported and must not be rounded
or silently cast into this contract.

Two coefficient conventions are deliberately recorded: public `bas_ctr_coeff`
removes primitive normalization, while `_libcint_ctr_coeff` describes the actual
integral representation. Expanded `_basis` alone may be insufficient to reproduce
an object modified after construction; ordered actual shell data must also match.
The named private accessors are version-bound PySCF details, not a portable
cross-engine basis format. See the [capture-field mapping](../research/gpu-readiness/CAPTURE_FIELDS.md)
and [PySCF molecular implementation](https://pyscf.org/_modules/pyscf/gto/mole.html).

Requested settings are not a complete measured numerical-environment record.
This first snapshot does not save actual SCF quadrature arrays, transient force
response grids, auxiliary-basis tensors or a GPU equivalence record. In particular,
the final `mf.grids` arrays cannot be relabeled as the transient full-response
force grid. The comparison therefore reports determinant evidence under the saved
context, not proof that all numerical settings or electronic surfaces are equivalent.

## Arithmetic checks and comparison

For each spin channel, let `C` contain occupied AO coefficients and `S` the AO
overlap. The reader requires metric orthonormality `Cᵀ S C ≈ I` and checks the
electron trace `Tr(P S) ≈ N`, where the derived density is `P = C Cᵀ`. The overlap
must be symmetric and positive definite. The arithmetic tolerance is `1e-8`, and
the smallest overlap eigenvalue must exceed `1e-12` times the largest. These are
implementation guards, not chemical-accuracy or physical-state thresholds;
ill-conditioned evidence is rejected instead of silently reducing the basis.

At identical geometry and AO context, the comparison forms
`O = C_leftᵀ S C_right`. Its singular values are the cosines of occupied-subspace
principal angles. Rotating occupied orbitals, changing their signs or permuting
their order preserves these values. The reported spin-channel density distance is
`Tr[(P_left − P_right) S (P_left − P_right) S]`. It uses the nonorthogonal AO metric,
not a raw coefficient or density-array difference. [PySCF orbital-overlap code](https://pyscf.org/_modules/pyscf/tools/mo_mapping.html)

Both spin channels are reported separately. There is no aggregate score that hides
a spin-density change behind agreement of the total density. Small numerical
excursions of singular values above one are reported and clipped only for the
arccos calculation, within the declared allowance. The original singular values
remain available. Empty spin channels have empty angle arrays and zero density
distance.

The result is one of `compared`, `insufficient_evidence`, or
`incompatible_context`. Only `compared` has metrics. Geometry/order/frame,
expanded basis/AO representation, reference, producer and physical/numerical
settings must match. Starting guess, requested threads and memory may differ;
those differences are preserved in `comparison_axis`. Initial guesses are starting
conditions and must be allowed to differ when comparing a guess survey.

The two saved AO overlap matrices must be **exactly equal** in this first version.
An absolute elementwise tolerance can hide a large relative change along a small
metric eigenvalue and produce misleading, direction-dependent overlaps. Evidence
with even a small metric difference is therefore declined rather than silently
choosing one matrix; no tolerance-based cross-metric bridge is implemented.

The output binds both normalized input records with SHA-256 digests and includes
their call IDs, numerical diagnostics and limits. File reference hashes identify
the exact written bytes; comparison digests identify canonical JSON content.
Hashes establish consistency, not provenance authenticity or a correct physical
interpretation. No similarity cutoff, automatic state selection, ground-state
decision or scientific pass/fail is inferred.

## Boundaries and the next scientific use

Occupied-subspace agreement can persist for an incorrect or unstable solution.
Orbital stability, spin/density character, suitability of the electronic method
and reference calculations remain separate evidence needs. The maximum-overlap
literature supports following selected SCF solutions but does not make overlap a
physical-state certificate. This module performs post-calculation comparison, not
maximum-overlap SCF optimization. [Gilbert, Besley and Gill, 2008](https://rsc.anu.edu.au/~pgill/papers/118MOM.pdf)

For different nuclear geometries, `C_leftᵀ C_right` or reuse of either same-geometry
overlap is wrong: one needs the cross-AO overlap for the two explicitly represented
molecules. PySCF documents such cross-geometry overlap construction. The production
comparator currently declines those comparisons; the separate X1 research bridge
is not automatically consumed. [PySCF cross-geometry wavefunction example](https://pyscf.org/user/ci.html#wavefunction-overlap)

The next scientific exercise is an explicitly bounded capture of the already
planned guess/geometry evaluations, followed by expert assessment of observed
differences. This document launches none. Synthetic tests verify arithmetic,
parsing, binding and output behavior; they do not establish accuracy of the methane
calibration, the 53-atom candidate, or a working molecular machine. See the
[validation protocol](VALIDATION_PROTOCOL.md) for the remaining scientific gates.
