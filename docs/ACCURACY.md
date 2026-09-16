# Accuracy checks and their limits

The implementation now distinguishes three different questions: whether the energy/force code is numerically consistent, whether a structure has the expected local curvature, and how a selected electronic method compares with published calculations for a specific reaction. None establishes a reliable complete assembly tool.

## Numerical verification

The test suite checks analytical quantum forces against independent central differences of electronic energies, including dispersion and density-fitting variants. It checks coordinate translation and rotation, charge/spin validity, missing convergence, topology screening, data provenance and known analytic reaction-path surfaces. Analytic test potentials test optimizer software only; production calculations use the quantum backend.

## Fixed-geometry reaction comparison

`nanodesign benchmark --out runs/methane-reference` evaluates methane, the ethynyl radical, a published transition-structure geometry, acetylene and the methyl radical using consistent functional, basis, grid and dispersion settings. Electronic states differ as required by each species. All atom identities, coordinate conversions and source hashes are checked before calculation.

It reports the energy at the supplied transition-structure geometry relative to separated reactants, and the reaction energy between the supplied product/reactant geometries. Comparison with a published high-level bare electronic barrier includes differences in electronic method, basis and potentially geometry. It is **not** the activation barrier of a newly optimized DFT pathway or a calibrated estimate of the error for a diamondoid tool.

The published 0 K enthalpic barrier includes a vibrational correction; it is not compared to the uncorrected electronic energies. Source supporting-information absolute energies remain excluded as targets because of the inconsistency documented in `SOURCE_NOTES.md`.

The supplied methane transition geometry also has unresolved source limitations: the main paper reports three imaginary modes at UCCSD(T)/cc-pVDZ, and its tabulated distances differ from the SI coordinates. It is treated as a nominal collinear transition-structure geometry, not a verified first-order saddle. Coordinates are not silently changed to fit the comparison.

Inputs and completed species are saved incrementally, with per-species electronic diagnostics. A failure never produces a complete comparison from an incomplete set of species. This is a reference comparison, not an automatically passed physical-validation gate.

## Local curvature and vibrational modes

`nanodesign characterize` uses the structure and electronic method selected in the design, with the same fixed anchors. It first checks the free-coordinate force residual. If that exceeds the requested tolerance, it stops before undertaking the much more expensive Hessian calculation.

For free coordinates q, the Hessian is estimated as H_ij = −[F_i(q+h e_j)−F_i(q−h e_j)]/(2h). Both the raw asymmetry and symmetric Hessian are recorded. Dividing by the square roots of the relevant atomic masses gives the mass-weighted curvature matrix. Its eigenvalues determine the squared harmonic frequencies. Negative eigenvalues are reported as negative signed frequencies (imaginary normal modes).

Frozen coordinates are excluded. No global translations or rotations are projected away. Thus an unanchored molecule retains external modes among its 3N coordinates, and an anchored structure describes the conditional fixed-boundary dynamics only. Frequencies near zero are classified as unresolved under the explicitly stated threshold. Repeat with different displacements, tighter quantum settings and threshold sensitivity before treating a small mode as physical.

Every free Cartesian coordinate requires two displaced force calculations, plus a force evaluation at the original structure. A 53-atom candidate with six fixed atoms therefore requires 283 force evaluations. The command has an explicit free-coordinate cost guard. It rejects unsupported constraint types instead of silently changing the Hessian being computed.

No negative modes at low force residual suggests a local minimum under these constraints. One resolved negative mode suggests a first-order saddle candidate. A reliable transition-state assignment also needs the correct reaction direction, forward/backward connectivity, electronic-state verification and numerical convergence. Zero-point/free-energy corrections, tunnelling and reaction rates remain unimplemented.

## Outstanding work before design selection

- Reoptimize and characterize reference stationary points at the selected method, then establish consistent high-level comparison targets.
- Converge basis, grid, molecular size, frozen boundaries and handle stiffness.
- Explore competing chemistry, different tool poses and the approach/withdrawal/regeneration cycle.
- Estimate uncertainty and operating failure probabilities using appropriately validated dynamics and temperature-dependent methods.
- Compare with experiments before claiming fabrication or operating reliability.
