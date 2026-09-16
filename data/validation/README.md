# Executed calculation records

These are actual calculations run on 2026-09-16, not illustrative numbers.

| Record | Calculation | Interpretation |
|---|---|---|
| `pbe0-svp/benchmark.json` | PBE0-D3(BJ)/def2-SVP on five published reference geometries | Fixed-geometry comparison; neither calibrated accuracy nor a verified barrier |
| `pbe0-tzvp/benchmark.json` | Same with def2-TZVP | Basis sensitivity at identical coordinates |
| `h2-integration/relax/result.json` | Actual DFT H2 relaxation | Small-system integration check, not a mechanosynthesis model |
| `h2-integration/modes/result.json` | DFT force-difference Hessian after H2 relaxation | Checks the complete geometry-to-vibration workflow |

The reference-geometry electronic energy relative to reactants is −3.62484 kcal/mol with SVP and −2.89058 kcal/mol with TZVP. The basis change shifts this quantity by +0.73426 kcal/mol. The source's reported 2.2 kcal/mol bare high-level barrier is a different method/geometry calculation; discrepancies of −5.82484 and −5.09058 kcal/mol are **combined differences**, not isolated functional errors. The nominal transition geometries retain substantial forces under these DFT methods.

The supplied transition geometry has documented source limitations, including three reported imaginary modes and an SI/main-table geometry discrepancy. The earlier SVP record predates the additional source review; it is retained as written, and the updated interpretation above and in `docs/ACCURACY.md` applies. No archived calculation is retroactively described as a verified transition state.

Each reference comparison preserves input snapshots, input hashes, per-species energies/forces and convergence events. The source data retain their CC BY-NC attribution in each `inputs` directory. Portable exports omit the machine-specific `source_directory` field from the input manifest and record the original manifest hash; the remaining input and calculation files are copied unchanged.

The H2 calculation yields one resolved positive frequency, approximately 4383.90 cm⁻¹, and five unresolved external/near-zero modes at the stated 20 cm⁻¹ threshold. It has no resolved negative mode and a residual maximum force of 0.0002874 eV/Å. These are computed results at the chosen method; no experimental-accuracy claim is made.
