# Cross-geometry AO overlap bridge (X1)

Research adapter for the electronic snapshots owned by E3. This is a basis-
reconstruction and integral-consistency tool, not electronic-state identification,
branch tracking, a ground-state search, or chemical validation. No SCF, gradient,
optimization, or new quantum campaign is performed here.

`cross_overlap.py` consumes the concrete E3 snapshot schema. The bridge
computes `S_AB[mu,nu] = <AO_mu(left) | AO_nu(right)>` in the declared common
Cartesian frame. Nuclear coordinates are never independently aligned/recentered.
E3 owns the occupied-subspace mathematics; this adapter does not duplicate it.
Its current comparator still rejects different geometries. A computed bridge is
a separate research artifact, not an automatically integrated continuity result.

## Use

From the repository root, with two E3 snapshot JSON files:

```sh
.venv/bin/python research/ao-overlap-bridge/cross_overlap.py left.json right.json \
  --common-frame-id study-lab-frame --atom-mapping 0 1 --output new-bridge.json
```

The example mapping is for two atoms; give every zero-based atom index in order.
Only identity correspondence is supported. The common-frame identifier is an
explicit caller declaration; any recorded non-null frame ID must agree. Neither
the mapping nor the frame label proves physical correspondence. No alignment,
translation or rotation is applied. A whole molecule's displacement can change
orbital overlap without changing its chemical state.

Python API: `build_cross_overlap(left_path, right_path, *, common_frame_id,
atom_mapping)`. It returns `status="computed"`, the dimensionless `cross_overlap`
matrix with left bra rows/right ket columns, exact source-file hashes and separate
canonical-content hashes, call IDs, geometry/basis digests, reconstruction
residuals/tolerances, and non-certification fields. Unsupported input raises
`BridgeError`. The CLI records `status="not_computed"`, no matrix, and returns 2
for these failures; it never overwrites an output file.

The runtime PySCF version, normalization convention and Bohr conversion must
match capture. Ordered elements, full expanded/actual basis, charge/spin,
reference, producer, and physical/numerical method settings must agree across the
pair. Initial guess, requested threads and memory may differ and are reported.
Reusing one call ID for contradictory content is rejected. The original bytes
are never modified; no old energy-only record is upgraded into an orbital record.

Bounds: at most 16 MiB per snapshot, 256 atoms, 512 AOs, 256 primitives per shell,
and angular momentum 0 through 6. E3 also validates the integer real RKS/UKS
determinants and AO metric. This adapter supports untagged real elements,
nonperiodic all-electron Gaussian bases, and zero-kappa shells. Basis string
aliases, executable molecule dumps, and unsupported representations are rejected.
These implementation bounds are not a promise of runtime/memory performance.

## Required reconstruction contract

1. Bind both exact snapshot artifacts, calculation IDs, ordered atom identities,
   charge/spin, full-precision coordinates and their units/common frame.
2. Reconstruct from explicit numeric expanded basis shells and the captured
   Cartesian/spherical and normalization conventions. A basis name is inadequate.
   Reject unsupported ECP/pseudo/ghost/spinor representations rather than infer.
3. Check the actual AO dimensions/order and shell definitions, then compare each
   recomputed same-geometry overlap with its captured matrix. Preserve negative
   contraction coefficients; reject nonfinite coefficients and invalid exponents.
4. Only after those checks, compute the directional cross-overlap and bind it to
   left/right hashes, actual geometry/basis definitions, operator, software and
   declared numerical tolerance. Matching inputs is not proof of a physical state.

Use numeric atom tuples and allowlisted fields. Do not deserialize an arbitrary
`Mole.dumps()` payload: installed PySCF's `loads()` evaluates stored strings.
Do not mutate PySCF's process-global normalization or thread controls in a bridge.
Expanded constructor input, `bas_ctr_coeff` and actual `_libcint_ctr_coeff` are
checked separately: the latter includes primitive normalization. Pointer offsets
in an internal environment array are not treated as basis function identity.
The declared shell tolerance is absolute/relative `1e-12`; the overlap absolute
tolerance is `1e-10`. The reconstruction must additionally satisfy a spectral
norm of at most `1e-8` after whitening its error by the captured AO metric:
`||S^(-1/2) (S_rebuilt - S) S^(-1/2)||₂`. This refuses an apparently small
absolute error that is large along a nearly linearly dependent AO direction.
Near-singular metrics may therefore fail closed even when entrywise agreement
looks good. These are reconstruction checks, not chemical error bounds.

## Preliminary bounded checks

An integral-only two-atom/two-AO check with PySCF 2.14.0 reconstructed H2/STO-3G
from expanded basis shells and Bohr coordinates with zero maximum absolute
self-overlap discrepancy. Changing the H-H separation from 1.4 to 1.5 Bohr
produced a cross-overlap whose reverse-pair transpose agreed exactly. These are
integral checks with specified geometries, not solved electronic structures.

Independent review also checked normalized single-atom 1s Gaussians: changing
the exponent can leave AO labels and self-overlap unchanged while changing the
cross-overlap. Thus self-overlap matching alone cannot authenticate the basis;
the expanded function definitions must also be bound and checked.

Actual verification is limited by assignment to at most four atoms and ten AOs.
The existing scientific calculations and their original evidence are unchanged.

Tests use actual E3 capture/save APIs over clearly synthetic, already populated
determinants and real tiny Gaussian integral matrices. Fixture `kernel` and
gradient methods raise if invoked. They check analytical Gaussian overlaps,
Angstrom-to-Bohr conversion, swapped-pair transposition, spherical/Cartesian
angular functions, exact file binding and adversarial metadata/basis cases.
Run:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q research/ao-overlap-bridge
```

All **37 tests pass** on the completed handoff: seven real integral/producer
checks plus 30 contract cases. The initial integration test caught E3 importing
an unavailable `pyscf.lib.param` module; its owner fixed that import to use the
actual `lib.param` alias, and the full suite passed after repair. No producer
files were patched by this lane. These results verify software and tiny integral
arithmetic, not solved electronic states or an operating molecular tool.
The two final regressions address E4's reported metric-conditioning boundary:
small entrywise errors must also be small relative to the captured AO metric.

## Primary references

- [PySCF cross-molecule AO integral implementation](https://pyscf.org/_modules/pyscf/gto/mole.html#intor_cross)
  defines the bra-left/ket-right overlap and Cartesian/spherical handling.
- [PySCF determinant-overlap implementation](https://pyscf.org/_modules/pyscf/scf/uhf.html#det_ovlp)
  uses cross-AO overlaps between occupied orbitals. The subsequent electronic
  comparison belongs to E3, with independent E4 review.
- [Gilbert, Besley and Gill, 2008](https://doi.org/10.1021/jp801738f)
  describes overlap-guided SCF solution following. This bridge does not implement
  that algorithm or establish continuity from an isolated overlap value.
