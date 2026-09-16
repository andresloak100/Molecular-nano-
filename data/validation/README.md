# Executed calculation records

These are actual calculations run on 2026-09-16, not illustrative numbers.

| Record | Calculation | Interpretation |
|---|---|---|
| `pbe0-svp/benchmark.json` | PBE0-D3(BJ)/def2-SVP on five published reference geometries | Fixed-geometry comparison; neither calibrated accuracy nor a verified barrier |
| `pbe0-tzvp/benchmark.json` | Same with def2-TZVP | Basis sensitivity at identical coordinates |
| `h2-integration/relax/result.json` | Actual DFT H2 relaxation | Small-system integration check, not a mechanosynthesis model |
| `h2-integration/modes/result.json` | DFT force-difference Hessian after H2 relaxation | Checks the complete geometry-to-vibration workflow |
| `h2-integration/modes-half-step/result.json` | Same H2 structure with displacement halved to 0.0015 Å | Numerical step sensitivity; stretch changes by −0.0896 cm⁻¹ |
| `paired-ccpvdz/method_comparison.json` | PBE0-D3(BJ) and CCSD(T) at identical cc-pVDZ reference geometries | Electronic-state identity is **unverified**; these are method differences, not accuracy estimates |
| `paired-ccpvdz-atom/method_comparison.json` | Same paired comparison with explicit HF `atom` guess | Reproduces the nominal +2.3986 kcal/mol UCCSD(T)/cc-pVDZ relative energy; state/accuracy still unverified |
| `h-abstraction-direct-initial/result.json` | Direct PBE0-D3(BJ)/def2-SVP energy and forces on the actual 53-atom candidate | Completed single point; large residual forces, unrelaxed geometry |
| `h-abstraction-df-initial/result.json` | Density-fitting variant at identical coordinates | Numerical approximation comparison at one geometry, not chemical calibration |

The reference-geometry electronic energy relative to reactants is −3.62484 kcal/mol with SVP and −2.89058 kcal/mol with TZVP. The basis change shifts this quantity by +0.73426 kcal/mol. The source's reported 2.2 kcal/mol bare high-level barrier is a different method/geometry calculation; discrepancies of −5.82484 and −5.09058 kcal/mol are **combined differences**, not isolated functional errors. The nominal transition geometries retain substantial forces under these DFT methods.

The supplied transition geometry has documented source limitations, including three reported imaginary modes and an SI/main-table geometry discrepancy. The earlier SVP record predates the additional source review; it is retained as written, and the updated interpretation above and in `docs/ACCURACY.md` applies. No archived calculation is retroactively described as a verified transition state.

Each reference comparison preserves input snapshots, input hashes, per-species energies/forces and convergence events. The source data retain their CC BY-NC attribution in each `inputs` directory. Portable exports omit the machine-specific `source_directory` field from the input manifest and record the original manifest hash; the remaining input and calculation files are copied unchanged.

The H2 calculation yields one resolved positive frequency, approximately 4383.90 cm⁻¹, and five unresolved external/near-zero modes at the stated 20 cm⁻¹ threshold. It has no resolved negative mode and a residual maximum force of 0.0002874 eV/Å. These are computed results at the chosen method; no experimental-accuracy claim is made.

Halving the H2 displacement from 0.003 to 0.0015 Å gives 4383.8070 cm⁻¹ for the stretch, a change of −0.0896 cm⁻¹. Both runs retain five unresolved near-zero modes and no resolved negative mode. The exact comparison is in `h2-integration/step-comparison.json`; this only checks the software and this small system's displacement sensitivity.

The original paired cc-pVDZ run with the default `minao` HF guess gives nominal transition-structure relative energies of −3.5574 kcal/mol (DFT) and −11.9239 kcal/mol (CCSD(T)). Its UHF transition-structure reference has S²=1.2143 versus the intended doublet value 0.75. The ethynyl reference converged to a different solution from the source calculation. This original evidence is retained, including its schema version 1.

[A retrospective interpretation sidecar](paired-ccpvdz/interpretation.json) binds these cautions to the original result's SHA-256. It labels the inferred MINAO choice separately from the absent original field and leaves the original record unchanged.

The later `atom`-guess run gives **+2.3986156 kcal/mol** for the CCSD(T) nominal transition-structure relative energy, with the same DFT value of −3.5573880. The DFT−CC difference is −5.9560036 kcal/mol. Ethynyl and transition-structure HF spin squares are 1.22118 and 1.21426. Source forensics separately reproduce four SI energies and isolate the inconsistent methane entry; see `si-energy-reproduction/`. The paired run uses computed values for all five species, without substituting published energies. An explicit guess is a reproducible control, not proof of the ground state. Current outputs record `electronic_state_identity_verified: false`; no actual first-order saddle or chemically accurate tool prediction is established.

The 53-atom direct and density-fitting runs both converged SCF and completed analytical forces. Their total energies differ by −0.00652168 eV (DF minus direct), and their maximum per-atom force-vector difference is 0.000473513 eV/Å. The density-fitting structure still has a maximum free-atom force of 0.95806 eV/Å, so it is not relaxed. Observed elapsed times were 2044.84 and 701.53 seconds under overlapping, different workloads: the ratio is not a controlled hardware benchmark. Exact inputs, raw results and the independently checked comparison are archived. Agreement at this one structure does not bound errors along a reaction path or validate the electronic method.
