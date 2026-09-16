# G3: capture actual CPU numerical evidence for E3

Status: source-based collection map, 2026-09-16. This document supports E3's
electronic snapshots; it implements no capture hook, GPU adapter or calculation.
E3 owns its snapshot schema and C1 owns production integration. G2 code and its
historical records remain unchanged. No runtime probe or solver import was used
for this mapping.

The immediate useful capture is the **live molecule, ordered AO representation
and final orbital arrays**. Store what exists at the declared capture phase.
Missing fields stay unavailable; do not rebuild an object or run a new integral,
grid, gradient, SCF or stability calculation to make a snapshot appear complete.
The initial E3 scope is AO/subspace evidence, not extending G2 parity acceptance.

## 1. Capture boundary and three distinct kinds of information

In the inspected baseline [calculator](../../nanodesign/quantum.py), `molecule`,
`mean_field` and `gradients` are local to `calculate()`. Capture after the relevant
solver stage while the needed objects and call ID remain available, before reset
or discard. The wrapper converts the SCF result with `float(...)` and the
gradient with `np.asarray(..., dtype=float)`. A capture made afterward cannot
recover the original return dtype. C1/E3 may preserve the original return value
alongside the conversion when integrating their hook; this document makes no
production change.

During G3 review, C1's ongoing integration added an E3 hook **after converged SCF
and before gradients**, labeled `converged_scf_only`. That phase is appropriate
for the state/SCF fields below. The gradient object and force results do not yet
exist at that boundary. A later force failure does not erase the earlier SCF
snapshot or turn it into a completed force evaluation. Fields described below
for a post-gradient boundary apply only if that later boundary is separately
supplied; this map does not ask E3 to move or extend its current hook.

Keep three explicit categories in an E3 record:

- **Observed object data:** existing arrays, shapes/dtypes, populated molecule
  and grid objects, actual callable identities, and solver output flags.
- **Derived from captured data:** a density matrix reconstructed from saved
  coefficients/occupations, or a canonical hash of captured arrays. Retain the
  derivation and inputs; these are not independently observed arrays.
- **Requested or unavailable:** input basis names, grid level/override recipes,
  uncaptured transient arrays and internal arithmetic precision. Request labels
  must not be promoted to measured resolved values.

Bind the capture to the exact call, geometry/atom order, charge/spin, constraints,
settings, versions and result. Retain failure/nonconvergence state. A complete
snapshot is neither electronic-state identity nor a ground-state certificate.

## 2. Fields safely available from an existing molecule

Use `mf.mol` as the identity source; confirm it is the object associated with the
calculation rather than reconstructing a new molecule from the design file.
Accessors below inspect stored tables or perform simple layout arithmetic, not
electronic-structure calculations.

| Evidence | Actual source | Interpretation / limit |
|---|---|---|
| Atom identity and coordinates | Existing `mol._atm` / coordinate slices, or `mol.atom_coords(unit='Bohr')`, ordered atom symbols/charges | Coordinates in internal Bohr; preserve geometry binding and recorded conversion. Do not infer units from the original input label. |
| Basis representation | `mol.cart`, `mol.nbas`, `mol.nao_nr()`, `mol.ao_loc_nr()` | Cartesian versus real spherical representation; shell offsets and dimensions in actual order. |
| AO labels | `mol.ao_labels(fmt=False, base=0)` | Preserve tuple content, zero-based atom indices and ordering. Labels alone do not establish phases or normalization. |
| Ordered expanded shells | Each shell's `bas_atom`, `bas_angular`, `bas_nprim`, `bas_nctr`, `bas_kappa`, `bas_exp`, and `_libcint_ctr_coeff` | Bind each shell to its center; keep primitive and contraction order. Copy views before object mutation. |
| Orbital/state arrays | Existing `mf.mo_coeff`, `mf.mo_occ`, `mf.mo_energy` | RKS arrays versus UKS alpha/beta arrays must be explicit. Preserve native dtype, shape and ordering; never flatten spin channels together. |
| Final density, if desired | Existing retained density, or algebra from captured MO coefficients and occupations | Label reconstruction separately. PySCF's RKS/UHF `make_rdm1` contracts occupied columns and occupations; no extra SCF is needed, but memory allocation must be bounded. |
| Overlap evidence | An overlap matrix retained from the existing call, if provided to the E3 boundary | It is not generally an already stored `mf` attribute. Calling `get_ovlp()`/`intor()` computes integrals. If E3's assigned capture explicitly includes that evaluation, label it capture-stage computation with its cost; it is not a passive field read or recovery of the original array. |

The stored coefficient convention matters. In PySCF 2.14.0,
`_libcint_ctr_coeff(shell)` returns the coefficients stored in `_env`, with the
normalization used by the integral representation. `bas_ctr_coeff(shell)` divides
out primitive normalization. They must not share a hash namespace without an
explicit conversion. Shell construction can also normalize contractions.
[Source: shell construction and accessors](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/gto/mole.py#L984).

For a future G2 `orbital_basis.sha256`, use an agreed versioned serialization of
the **ordered semantic shell payload**, its coefficient convention and AO
representation. Whole `_env` hashes are useful raw-source fingerprints but are
not a normalized basis identity: `_env` also contains coordinates/operator data,
and `_bas` contains layout-dependent pointers. Do not sort shells or AO labels.
Retain atom-center identity/order and bind coordinates separately. G3 does not
implement or declare a completed cross-backend basis serializer.

AO labels and matching shell hashes still do not establish a CPU/GPU permutation,
phase convention or occupied-subspace equivalence. E3/F1 must preserve the basis
and overlap context for their own comparisons; direct entrywise MO comparison is
not supplied by this mapping.

## 3. Density-fitting evidence: resolved object versus name

With density fitting enabled, inspect the existing `mf.with_df`. Its `auxbasis`
is a request/resolution label; capture the actual **existing `with_df.auxmol`**
with the same ordered shell/AO representation above for expanded-basis evidence.
Normal `DF.build()` assigns that object. Automatic selection can choose predefined
or even-tempered bases. A label such as `None` does not describe the basis that
was ultimately used. [DF lifecycle](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/df/df.py#L103),
[basis resolution](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/df/addons.py#L230).

Preserve these exceptions:

- With no DF requested/used, G2 permits `auxiliary_basis.status: not_applicable`
  and no hash. With DF enabled but no retained auxiliary molecule, the evidence
  is unavailable, not `not_applicable`.
- A supplied `_cderi` can bypass auxiliary-molecule construction. Some DF gradient
  paths construct an auxiliary molecule locally without retaining it on
  `with_df`. A local object not passed to the capture boundary remains missing.
- Do not call `build()`, `make_auxmol()` or `get_naoaux()` to fill the gap;
  these can construct objects or integrals. Do not open external `_cderi` paths
  automatically. Capture a storage-kind/status description instead.
- Auxiliary AO count is not necessarily the retained rank of a metric
  factorization. Matching the auxiliary basis does not establish matching
  three-center tensors, metric cutoffs, J-only versus J/K fitting, or auxiliary
  response derivatives. Preserve available `only_dfj`/DF-class and response
  settings separately; do not infer them from the basis name.

## 4. Critical distinction: SCF grid and force-response grid

The current calculator enables `gradients.grid_response = True`. That is an
**algorithm setting**, not evidence that force-grid arrays were captured.

The inspected RKS/UKS gradient implementation calls `grids_response_cc`, which
dispatches to a generator such as `grids_response_becke`. It regenerates atomic
grids and yields local `(coords, weights, weight_derivatives)` blocks. Those
blocks are not retained as `mf.grids.coords/weights`. Meanwhile SCF initialization
can density-prune `mf.grids`. Consequently, post-call `mf.grids` arrays establish
the final **SCF grid only**, not the actual force-response grid.
[Gradient generator](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/grad/rks.py#L483),
[SCF initialization](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/dft/rks.py#L499).

For initial E3 capture, explicitly record a scope such as `final_scf_grid` and
`force_response_grid: unavailable_not_retained`. Capturing the transient force
blocks would require instrumentation during the existing gradient call. Do not
invoke the generator again and label the result as the arrays used by that call.
E3 is not assigned that instrumentation by this document.

| SCF grid field | Safe collection / interpretation |
|---|---|
| Actual points and partitioned weights | Copy existing `mf.grids.coords` `(N,3)` and `mf.grids.weights` `(N,)`, including dtype, units Bohr/Bohr³ and finite-value checks. `N` is an observed row count. |
| Per-row atom and unpartitioned weights | Existing `atm_idx` and `quadrature_weights`, if populated. They are distinct from partitioned `weights`; do not substitute one for the other. |
| Padding and ordering | Record `alignment`, raw row count and existing `atm_idx == -1` padding markers. Preserve duplicate and zero-weight rows. A zero weight alone is not proof of padding. |
| Actual callable identities | Read `prune`, `radi_method`, `radii_adjust`, `becke_scheme`: type/module/qualified name where identifiable, with installed source/version. Explicit disabled/none is a state, not a placeholder default. For custom closures/objects, name alone cannot establish behavior; retain unresolved identity. |
| Radius and screening configuration | Copy available `atomic_radii`; record `cutoff`, `mf.small_rho_cutoff`, and existing `non0tab`/`screen_index` array metadata within the capture budget. These supplement, not replace, actual grid arrays. |
| Recipe request fields | Preserve `level` and `atom_grid` exactly as observed configuration. An empty `atom_grid` is an unresolved per-element recipe, not an empty integration grid. |
| Transient generation details | Actual resolved radial counts, angular counts per radial shell and pruning-mask history are generally not stored on the completed `Grids` object. Leave unavailable unless captured during that original execution. |

`gen_atomic_grids` resolves absent recipes from tables and can reinterpret an
angular-order request as a point count. Callable pruning then varies angular
counts by radius. The final points cannot generally reconstruct that history.
Reading today's tables from a level is a source-derived reconstruction, not an
observed generation trace. Do not populate G2's required per-element
`atom_grid.{element}.{radial,angular}` with such a reconstruction under a measured
label. This deliberate missing field may prevent a partial E3 snapshot from
passing G2. [Grid implementation](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/dft/gen_grid.py#L254).

Custom per-atom recipes can also exceed G2's current per-element schema. Preserve
the richer E3 evidence and mark the mapping unsupported; do not collapse distinct
atom recipes just to satisfy the older schema. Capture a separate NLC grid only
if one actually exists and was used; it cannot stand in for the ordinary XC grid.

## 5. Exact existing G2 interface and its limits

The frozen [resolved_numerics.py](resolved_numerics.py) exposes:

```python
validate_resolved_numerics(cpu_record, gpu_record, case)  # list[str]
canonical_grid_signature(rows)  # SHA-256 hexadecimal string
```

`rows` must be a nonempty list/tuple of finite ordinary Python-number quadruples
`(x, y, z, partitioned_weight)`, in Bohr and Bohr³. The function converts to
binary64, normalizes signed zero, sorts **complete paired rows**, and hashes a
little-endian uint64 row count followed by little-endian float64 quadruples.
The identifier is `paired-xyz-weight-f64le-sorted-v1`. It does not discard
duplicates or padding, round coordinates, or sort points and weights separately.
This is a derived hash of the saved SCF grid if the input is `mf.grids`.

Check E3's declared capture size bound **before** array-to-list conversion:
canonicalization allocates rows and sorts them. If too large, keep a bounded
unavailable/omitted status with the reason; do not truncate the grid and hash it
as complete. Independent raw coordinate and weight hashes are optional and do
not replace the paired signature.

| Required G2 field | E3 mapping readiness |
|---|---|
| `precision.energy/forces/density == float64` | Capture native return/array dtypes where available; post-cast output storage does not establish internal arithmetic precision. A reconstructed density must be labeled derived. Missing pre-cast evidence stays missing; this map does not authorize a blanket float64 claim. |
| `orbital_basis.sha256` | Ordered expanded-shell evidence is collectable; the agreed versioned semantic serializer is still needed. An arbitrary raw archive hash is not automatically this identity. |
| `auxiliary_basis.status/sha256` | Depends on DF applicability and the retained actual auxiliary molecule, with the same serializer requirement. |
| `quadrature.level/pruning/radial_method/becke_scheme` | Live configuration identities can be recorded, distinguishing explicit none/custom behavior from missing data. |
| `quadrature.atom_grid` | Actual resolved per-element counts are not generally retained. Do not fill from defaults or confuse them with final pruned row counts. |
| `quadrature.point_count/signature_sha256/canonicalization` | Derivable from complete saved SCF points and partitioned weights, with explicit SCF-only scope. |

`validate_resolved_numerics` returns `unavailable:`, `malformed:` or `mismatch:`
messages; an empty list means only its supplied-field checks passed. It does not
authenticate arrays, verify execution or infer scientific state identity.
G2's single quadrature object does **not** distinguish SCF and force-response
streams or require all the richer recipe/screening fields above. Its accepted
record comparison must not be interpreted as complete force-grid provenance.
C1 has accepted this scope clarification; no G2 contract change is made here.

## 6. Source binding and handoff

Initial inspected checkout HEAD: `026d3fc69a03c67ef91e3f0d2252093a84840320`.
Installed package metadata: PySCF `2.14.0`; local source root is
`.venv/lib/python3.14/site-packages/pyscf/`. Line references refer to these
inspected bytes. The calculator row is bound to that committed baseline, not the
concurrently changing working file. The subsequent before-gradient E3 hook was
read as a scope update above; this is not a review of that implementation.

| Source | Relevant lines | SHA-256 of inspected file |
|---|---|---|
| `nanodesign/quantum.py` at the baseline commit above | 208–315: actual objects, conversions, response flags and results | `9e70a81c5658d23e64568f30bf455e1de9875d76a174d0a687ed675ee1231b2b` |
| `research/gpu-readiness/resolved_numerics.py` | 22–62, 179–242: existing signature and validation API | `0516754a5443b0b33b49148d3da7934acafdfc735eab54d61d44b8915f7d3875` |
| PySCF `gto/mole.py` | 984–1018, 1456–1473, 1658–1674, 3390–3440 | `005bd63c5fa100a20cc264a995c38f44f2e05c7c82d37bd36eba2355dc651392` |
| PySCF `df/df.py` | 103–169, 248–266: auxmol lifecycle and implicit build | `e88c889248ddfd211bbd58a825e63d441c47040bbbd49a54b2994d49f138ff3d` |
| PySCF `df/addons.py` | 230 onward: actual auxiliary-molecule resolution | `0041fedea02f200aae8ea199856a5f3fdabebfd9dffb1af0a9e21eb87b2497c7` |
| PySCF `df/grad/rhf.py` | 85–92: local auxiliary-molecule fallback | `a52525f269974884b28846b55d021c169d430048e0755b55557480f733128735` |
| PySCF `dft/gen_grid.py` | 254–338, 565–606, 649–741 | `44aa665fb25efc544a42cb80fdad737767c787c20479d4df0ae4b1649028ff14` |
| PySCF `dft/rks.py` | 499–520: initialization and density pruning | `fcc8d3dbd03df368f3b70f3ce5999be6a19c608a6ff0bdc1890074616f0befa8` |
| PySCF `grad/rks.py` | 111–118, 483–488, 605–614, 732–739 | `cffc3e9efa10cf1b547d757e9d165c1000377e9159c794f87464a79c8c93ab46` |
| PySCF `grad/uks.py` | 30–52, 160 onward: shared full-response generator | `57ba4d9528be5d5dbd1eadb7f91e539e03008b74b28ebfbbb9a33f50b841484a` |
| PySCF `scf/hf.py` | 855–869: RKS/RHF density construction | `34ea1a9057da3a2d5aa37fe51d4e7305580bcb8e0118a5f17ba95ba5bb36c22f` |
| PySCF `scf/uhf.py` | 191–206: spin-resolved density construction | `b3ee3e679ab64d51ad610f04d02f4ac332c3f9238aeede7ba6bd80ad962a8fee` |

The existing internal helper independently mapped the orbital/AO and auxiliary
objects. G3 checked source fields and document links only; no chemical result,
runtime capture, solver import, quantum job or repeated test suite was produced.
Both C1 and E3 received the force-grid distinction before this handoff. Unknown
capture fields remain unknown; the observed numerical representation supports
later investigation and does not validate the proposed molecular operation.
