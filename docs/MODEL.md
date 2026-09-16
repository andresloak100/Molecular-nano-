# Physical model and interpretation

## Defined system

Initial implementation: neutral, isolated finite clusters in the Born–Oppenheimer approximation. Nuclei are classical coordinates. Electronic energies and gradients come from Kohn–Sham DFT with an explicitly selected functional and basis. The default supported ethynyl radical plus closed-shell substrate has an odd number of electrons and a doublet electronic state.

For an atom position R, the computed force is F = −∂E/∂R. The underlying PySCF energy in hartree and gradient in hartree/bohr are converted to eV and eV/Å. Dispersion energy and gradient must use the same method/damping settings and must be included exactly once.

Fixed anchors define a conditional mechanical model: selected distal carbon coordinates remain fixed while mobile atoms relax. The reported sum of anchor forces is the **sum of internal forces on all anchors**, not automatically the tip actuator force. The external reaction force is its negative. The candidate additionally records separate tool and substrate anchor groups and their holding forces. Finite anchor stiffness and thermal motion are not represented.

## H-abstraction candidate

The ideal carbon cage has diamond tetrahedral connectivity. Hydrogens terminate its missing valences. The substrate presents a bridgehead C–H bond to an ethynyl radical on a second cage. The candidate is C22H31, neutral, spin 1 in the PySCF 2S convention. The product guess transfers the same H atom to the tip apex; no atoms are created or deleted.

This is a designed starting configuration. Its bonds, approach distance and anchors require relaxation and sensitivity studies. A small hydrogenated cage cannot establish the behavior of an extended diamond facet or its elastic surroundings. Neither endpoint establishes a synthesis route for the tip.

## Reaction path

The workflow first optimizes both endpoints with identical anchors. It rejects unconverged endpoints, H-transfer endpoints that no longer match the intended H-transfer basins, and identical relaxed endpoints. The H-basin predicate checks donor–H and acceptor–H distances only. When nominal endpoint bond lists are supplied, a separate geometry screen flags stretched expected bonds and unexpected close contacts using declared covalent-radius thresholds. This screening is a heuristic, not electronic bond-order analysis. A band with independently evaluated images is initialized with IDPP, relaxed without climbing, and then refined with a climbing image.

Only a converged band with an interior energy maximum yields a `candidate_electronic_barrier_ev`: max(Eimage) − Einitial. A missing interior maximum is not proof of barrierless chemistry. The endpoints are states at the specified tool pose; this energy difference is not automatically the gas-phase separated-reactant activation barrier.

NEB convergence alone does not verify a first-order saddle. The optional characterization command computes a constrained finite-difference Hessian and signed vibrational modes. Transition-state refinement, displacement-step convergence and forward/backward connectivity remain required. There is no zero-point, vibrational entropy, tunnelling, solvent, field, thermal bath or alternative-spin correction in these energies. No rate or assembly reliability is inferred from them.

## Validation layers

1. **Numerical verification:** force finite differences, unit conversion, atom conservation, fixed-coordinate consistency, SCF convergence, gradient convergence and path convergence. Tests exercise the implementation, not physical truth.
2. **Electronic-method calibration:** compare consistently defined reaction energies and barriers with reliable high-level calculations or measurements, including open-shell state diagnostics. DFT convergence is insufficient to establish the correct electronic state.
3. **Physical representation:** enlarge cluster and handle, vary boundary conditions and pose, investigate mechanical compliance, surface orientation, defects and environmental effects.
4. **Chemical selectivity:** competing abstraction/deposition reactions, accidental bonding, handle damage, approach, withdrawal and regeneration. A successful elementary path does not show a complete operational cycle.
5. **Device validation:** free-energy and dynamical sampling, positioning errors, temperature, actuation, fabrication and experimental comparison.

The result audit deliberately never emits a validated design. Physical validation is not yet implemented or supplied. A numerical pass is kept separate from the missing evidence.

## Primary implementation references

- [PySCF DFT](https://pyscf.org/user/dft.html) and [SCF methods](https://pyscf.org/user/scf.html): electronic method and state conventions.
- [Simple DFT-D3 Python interface](https://dftd3.readthedocs.io/en/latest/api/python.html): dispersion energies and gradients.
- [ASE NEB](https://docs.ase-lib.org/ase/neb.html): numerical path search and image constraints.
- [ASE optimizers](https://docs.ase-lib.org/ase/optimize.html): force-based geometry optimization.

Scientific reaction sources and their limitations are recorded in the repository's source notes, separately from software documentation.
