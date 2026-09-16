# E4 independent electronic-state evidence review

Owner: `codex-712e`; implementation owner E3 / `codex-f522`.
Scope: `nanodesign/electronic_state.py` mathematics, saved-state contracts and
interpretation. C1 owns solver integration; S3 separately reviews its lifecycle.
This review uses synthetic exact fixtures and launches no quantum calculations.

## Completed review — 2026-09-16

**All 35 independent synthetic checks pass. No unresolved material finding remains
within this review's scope.** The reviewed electronic-state module SHA-256 is
`e67b2a072ab415422dd4907ebf8198e70ef82313d96e2e52f9525772ba3bbaf6`.
[`verification-after.json`](verification-after.json) records exact output,
test hashes and unchanged before/after source hashes. E4 did not edit the
implementation or its production tests.

The original overlap-compatibility defect is preserved in
[`verification-before-overlap.json`](verification-before-overlap.json): 33 passed,
two failed manifestations of one issue. These checks establish synthetic
arithmetic/contracts only, not physical state identity, a calibrated electronic
method or operation of an assembly tool.

Independent checks cover:

- AO-metric electron counts and occupied-orbital normalization.
- RKS double occupations versus UKS separate alpha/beta single occupations.
- Invariance under occupied-orbital sign, permutation and real rotation.
- Known occupied-subspace changes in a nonorthogonal AO basis.
- Distinct spin densities concealed by an unchanged total density.
- Exact geometry, atom order, method, resolved basis and AO representation binding.
- Explicit rejection of cross-geometry comparison without cross-AO overlaps.
- Saved-file integrity and unsupported/malformed mathematical inputs.
- Numeric diagnostics that never claim physical state identity or a ground state.

## Findings and repair status

1. **Matching invalid or missing settings were accepted.** Two records could both
   omit `grid_level`, `density_fit`, `conv_tol` or `dispersion`, or contain the same
   invalid values, and still receive comparison metrics. E3 independently added
   complete explicit `QuantumSettings` validation during review. All eight E4
   reproductions now pass without reconstructing missing defaults.
2. **The capture function could not import its Bohr conversion.**
   `from pyscf.lib.param import BOHR` fails in the installed package because
   `param` is an exported alias for `pyscf.lib.parameters`. E4's fake-mean-field
   capture tests and X1's actual tiny producer tests independently caught it.
   E3 changed capture to `lib.param.BOHR`; all five E4 capture-convention tests
   pass without running a solver or evaluating an AO integral.
3. **Absolute overlap agreement hid a large relative metric change.** Synthetic
   `S_left = diag(1,1e-9,1)` and `S_right = diag(1,2e-9,1)` each meet the current
   conditioning guard. Independently normalized occupied coefficients pass each
   record's validation. Yet the `1e-8` absolute matrix tolerance accepts them as
   a common metric. Forward comparison reports an apparent `1/sqrt(2)` occupied
   overlap; reverse comparison raises on `sqrt(2)`. The records have inconsistent
   metrics and should receive no comparison metrics in either direction. E3 now
   requires exactly equal saved same-context metrics, returning
   `incompatible_context` in both directions. It also turns numerical allowance
   failures into structured `insufficient_evidence`. Both reproductions pass.

The third finding is preserved against electronic-state module SHA-256
`7ba2de99179d3998a162803421e75f75cf13f6d1665c0005ba9bab52b943cdff`.
These are synthetic arithmetic/context records, not Gaussian integrals for their
placeholder geometry metadata. Their purpose is to test the comparator's declared
matrix domain. Actual basis reconstruction is the separate X1 lane.

X1 was informed of the same conditioning concern and independently hardened its
reconstructed-overlap check using an absolute bound plus a metric-whitened
spectral residual. Its owner reports 37 passing checks, including a small-eigenvalue
case. E4 did not rerun or certify the separately owned X1 suite.

## Run the bounded review

```sh
.venv/bin/python research/electronic-state-review/run_review.py \
  --output research/electronic-state-review/NEW_RECEIPT.json
```

The runner refuses an existing receipt. It records test output, reviewed source
and test hashes, and whether source changed during execution. It launches only
the E4 synthetic tests; no SCF, gradient or AO-integral work is performed.

Different initial guesses are intentional in the fixed-geometry survey and must
remain visible without being mistaken for a Hamiltonian change. A saved converged
SCF snapshot is not proof that gradients, the complete call or a scientific gate
succeeded; backend lifecycle testing belongs to the separate S3 assignment.

## Independent mathematical references

With AO overlap matrix S and spin density P, electron count is Tr(P S). An
occupied coefficient matrix C must obey C-transpose S C = I in the supported real
representation. For two valid states in the same AO representation, the squared
Hilbert-Schmidt density distance is Tr[(P1-P2) S (P1-P2) S].

The occupied cross-overlap is C1-transpose S C2. Singular values describe overlap
between the occupied subspaces and are unchanged by rotations within either
occupied basis. Different geometries require a cross-AO matrix S12, as shown by
the official [PySCF orbital mapping source](https://pyscf.org/_modules/pyscf/tools/mo_mapping.html)
and [PySCF determinant-overlap implementation](https://pyscf.org/_modules/pyscf/scf/uhf.html).
Equal array sizes and matching basis names alone do not justify subtraction.

These metrics quantify representation overlap. Neither a universal threshold nor
energy ordering is treated as a certificate of physical electronic-state identity.
The complete protocol still requires reference suitability and branch evidence.

The final suite consists of six mathematical-invariance cases, 24 context/record
cases and five fake-mean-field capture-convention cases. The latter use the actual
installed PySCF import path but supply synthetic arrays and overlap; the fake
solver's `kernel` raises if called. A captured snapshot remains SCF-only evidence,
and raw-file integrity is not independent authentication of the computation.

Module/support fixtures are uniquely prefixed to avoid collision with B2's
separate fixtures when collected together. The earlier failure receipt retains
the original fixture filenames/hashes from that stage. All E4 files were left
uncommitted for C1 integration; no other owner's files or calculations were
changed.
