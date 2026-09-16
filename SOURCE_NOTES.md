# Mechanosynthesis model: evidence and validation plan

Prepared 2026-09-16. Scope: an auditable quantum-chemistry workbench for a proposed hydrogen-abstraction tool. Numerical convergence, chemical accuracy, and experimentally demonstrated tool operation are different claims. None establishes perfection.

## Direct primary sources

1. **Temelso, Sherrill, Merkle & Freitas (2006), “High-Level ab Initio Studies of Hydrogen Abstraction from Prototype Hydrocarbon Systems.”** [Article DOI](https://doi.org/10.1021/jp061821e), [author-hosted published paper](https://vergil.chemistry.gatech.edu/static/pdfs/temelso_2006_11160.pdf), [supporting data DOI](https://doi.org/10.1021/jp061821e.s001). This is a computational study, including ethynyl-radical abstraction and small diamond-site proxies. It identifies substantial method dependence: the tested DFT functionals underestimate barriers, while MP2 can be compromised by spin contamination. The reported isobutane result is a small-molecule calculation, not validation of an attached tool operating on a real diamond surface. Table 5 reports a methane-abstraction bare electronic barrier of 2.2 kcal/mol for RCCSD(T)/cc-pVTZ; its 0 K enthalpic barrier is 1.7 kcal/mol with the stated vibrational correction. These are distinct quantities. Table 1 also flags extra imaginary modes in several small-basis stationary structures; supplied geometries must be rechecked.

2. **Freitas & Merkle (2008), “A Minimal Toolset for Positional Diamond Mechanosynthesis.”** [Published paper](https://www.molecularassembler.com/Papers/MinToolset.pdf), DOI 10.1166/jctn.2008.002. A theoretical toolset and reaction-sequence proposal, useful for selecting structures, intended chemistry, and competing reactions. Its proposed closure is a property of the modeled reaction scheme; it is not a demonstrated fabrication system. [Authors’ overview](https://www.molecularassembler.com/Nanofactory/DMS.htm) explicitly distinguishes computationally studied sequences from experimental realization.

3. **Drexler (1999), “Building Molecular Machine Systems.”** [Author’s preprint](https://www.imm.org/reports/rep008/). Establishes the connection between the requested vision and rigid graphitic/diamondoid structures in vacuum. This historical framing supplies no modern performance guarantee.

## Published geometry package delivered

`data/reference/temelso_2006_si.txt` is the downloaded original supporting data. `data/reference/provenance.json` records source URL, SHA-256, attribution, per-structure charge/multiplicity, and transformations. The publisher repository labels this supplement **CC BY-NC 4.0**: retain attribution and account for its noncommercial condition when redistributing derived files.

The source coordinates are **Bohr**, not Å. The eight extracted `.xyz` files use Å, converted with 1 Bohr = 0.529177210544 Å. Only the nonphysical `X` dummy in the tert-butyl structure was omitted; real coordinates were otherwise unchanged. These are published starting geometries, not newly optimized or independently validated structures.

Files: `methane.xyz`, `ethynyl_radical.xyz`, `methyl_radical.xyz`, `acetylene.xyz`, `methane_ethynyl_ts.xyz`, `isobutane.xyz`, `tert_butyl_radical.xyz`, `isobutane_ethynyl_ts.xyz`.

**Data-quality hold:** directly subtracting the SI's labeled absolute methane and ethynyl energies from its methane transition-structure energy gives approximately −14.76 kcal/mol, inconsistent with the main paper's positive barrier. No cause is assigned here. Do not use these SI absolute energies as acceptance targets until reconciled with the original calculation settings or an author correction. Metadata marks them unverified. Geometry provenance remains useful independently.

Additional source checks: the methane UCCSD(T)/cc-pVDZ collinear structure has three imaginary frequencies in main-paper Table 1 (259i, 50i, 50i cm⁻¹), so it is a nominal transition-structure candidate, not a verified first-order saddle. Its converted SI donor–H and acceptor–H distances are 1.149220 and 1.672376 Å, whereas Table 2 lists 1.148 and 1.678 Å; this discrepancy is unresolved. Original coordinates are retained. Table 5 footnote d uses RCCSD(T)/cc-pVDZ vibrational corrections for the reported 1.7 kcal/mol 0 K barrier; the 2.2 kcal/mol bare comparator must remain separate.

## Proposed first executable chemical operation

Start with the gas-phase calibration reaction **C₂H• + CH₄ → C₂H₂ + CH₃•**, all species neutral. Methane and acetylene are singlets; ethynyl and methyl are doublets; the combined reaction is a neutral doublet. In PySCF conventions, the combined system has `charge=0, spin=1` (spin means Nα−Nβ, not multiplicity). There are eight atoms and 23 electrons. Use the ground-state ethynyl solution, inspect spin density and electronic stability, and watch for state switching along the path.

In `methane_ethynyl_ts.xyz` with **zero-based indices**: acceptor carbon = 2, transferred hydrogen = 3, donor carbon = 4. Atoms 0–2 are the original ethynyl group; atoms 3–7 are the methane group. The same 2/3/4 reactive indices apply to `isobutane_ethynyl_ts.xyz`.

Then extend to ethynyl + constrained isobutane, followed by an ethynyl group mounted on a diamondoid handle approaching a hydrogen-terminated diamond cluster/slab. Label each level explicitly. A free ethynyl radical is a chemical proxy, not a fabricated positional tool. A relaxed isolated isobutane is not a constrained diamond surface.

## Proposed equations and units

Use a Born–Oppenheimer electronic energy E(R), with forces Fᵢ = −∇ᵢE. For positional control, optimize free coordinates at fixed anchor coordinates or with documented restraints. A useful reaction coordinate is q = r(C_donor,H) − r(C_acceptor,H); the full multidimensional path must still be relaxed. Finite restraints add a separate mechanical energy, ½k|R_anchor−R_target|², which must be reported separately from molecular electronic energy.

Report separately: ΔE = E_products−E_reactants; ΔE‡ = E_saddle−E_reactants; zero-point correction; temperature-dependent free-energy correction; and any restraint work. A maximum along an unrelaxed interpolation is only a path-sampling maximum, never a transition-state barrier. A maximum along a converged NEB still needs saddle refinement and mode/connectivity checks. Submerged saddles must not be silently clamped to zero.

Useful conversions: 1 eV/Å = 1.602176634 nN; 1 Hartree ≈ 27.21138625 eV ≈ 627.50947 kcal/mol. Record exact software versions, functional, basis, integration grid, convergence thresholds, charge, spin, constraints, geometries and random seeds.

## Required validation gates (proposed engineering acceptance process)

1. **Input integrity:** atom identity/count/units and charge/electron parity checks; no clashes; correct transfer atom mapping; explicit anchor set. No invented bonds treated as evidence of chemistry.
2. **Electronic calculation:** converged SCF and energy gradients; stable intended electronic solution; spin contamination and density tracked across all images. Numerical tightening and grid/basis sensitivity checked where results affect decisions.
3. **Stationary structures:** optimize endpoints at the same method. Verify minima within permitted degrees of freedom. Refine the candidate saddle and establish exactly one chemically relevant negative-curvature mode; follow both directions to the intended endpoints. Extra modes cannot be dismissed automatically.
4. **Independent chemical calibration:** compare with a consistent high-level calculation and explicitly matching published quantities. Recompute the small ethynyl/methane benchmark rather than trusting the inconsistent SI absolute energies. Lower-cost DFT is a screening method until calibrated.
5. **Transferability:** grow cluster/slab and handle models; vary surface termination, frozen boundaries, restraint stiffness, functional and basis. Confirm reaction preference remains stable. Investigate undesired attachment, hydrogen exchange, handle self-abstraction and surface reconstruction.
6. **Operation reliability:** sample positional and angular offsets, finite-temperature fluctuations, approach/retraction speed and competing channels. Assess hydrogen tunneling where relevant. Gas-phase entropy and harmonic rates cannot be copied unchanged onto an anchored tool. No universal success probability follows from a barrier alone.
7. **Experiment:** validate fabrication, attachment, registration, force control, product identity and tool recharge. Computational checks alone do not establish a buildable machine.

Unpassed gates should remain visible in results. Report “calculated under these assumptions,” “numerically checked,” and “chemically benchmarked” as distinct statuses; reserve claims of experimental validation for measured evidence.
