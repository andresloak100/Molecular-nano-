# F2 independent energy/force review

Owner: `codex-9395`, assigned by C1. F1/4280 owns
`research/force-energy-consistency/`; F2 writes only this review directory and
its own status. This review uses exact analytic software fixtures, not quantum
calculations. E2's Hessian reconstruction and N1's mode comparisons are outside
this review.

**Result: 39 independent tests passed in 2.33 seconds against the frozen F1
implementation, with no unresolved finding in this bounded review.** F1 repaired
the issues below. No molecular or electronic-state validation is claimed.

The final run started on 2026-09-16 at 22:04:32 UTC. Its exact command, output
and before/after SHA256 values are in `test-output.txt` and `verification.json`.
All seven implementation/test files stayed unchanged during execution. The F1
implementation SHA256 is
`7c8962708b83ce3c0a97e27d9805c31cd926915db5312e2e108a87e25de56b77`.
The 188 warnings concern ASE's existing NumPy shape-assignment deprecation.

To reproduce from the repository root:

```sh
.venv/bin/python -m pytest -q research/force-energy-review
```

## Independent numerical oracle

`analytic_fixtures.py` evaluates a polynomial as rational arithmetic at the
actual supplied coordinates, then rounds once when writing a simulated saved
energy. Its force is known independently by differentiating the polynomial.
The fixture retains a nonzero force on a declared fixed atom.

Let the actual displacement lengths be `a = x0 - x_minus > 0` and
`b = x_plus - x0 > 0`. Differentiating the three-point interpolating polynomial
at the reference gives

```text
D = b/(a+b) * (E0-E_minus)/a + a/(a+b) * (E_plus-E0)/b
F_difference = -D
```

For the fixture `E = C + g*x + k*x²/2 + c*x³ + q*x⁴`, with `x=0` at the
reference, this evaluates to `F_difference = -g - c*a*b - q*a*b*(b-a)`.
The quadratic term cancels for arbitrary positive `a,b`. In contrast,
`-(E_plus-E_minus)/(a+b)` has quadratic bias when `a != b`: changing only the
denominator does not estimate the derivative at the reference correctly.
This oracle is derived from the polynomial and independent of F1's code.
[NIST DLMF §3.4](https://dlmf.nist.gov/3.4) describes numerical differentiation
through differentiated interpolation and its remainder.

For symmetric steps, the cubic residual scales as `h²`; three halved steps
therefore provide an exact manufactured trend against which to test reporting.
Observing such a trend in real saved data would not by itself prove a bound on
the derivative error or chemical accuracy. Subtracting a shared energy reference
can improve arithmetic but cannot recover energy digits already lost in saved
absolute energies. The large-offset fixture deliberately tests that distinction.

## Scientific evidence boundaries

The comparison must use total energy for the same model as the raw force.
Current production DFT diagnostics distinguish `dft_energy_hartree` from
`total_energy_hartree`; separately added dispersion belongs in both the total
energy and force. H1 records can be used only where the total and its conversion
are actually recorded. SCF-only energies cannot silently replace them.
Analytic gradient response terms are documented in the
[PySCF gradient API](https://pyscf.org/pyscf_api_docs/pyscf.grad.html).

Matching settings, starting guesses, convergence flags or an opaque declared
state ID do not establish that a single electronic branch was followed.
Numerical residuals and electronic-state verification must remain separate.
Likewise, source digests inside supplied JSON are declarations unless bound to
actual bytes by the adapter; matching strings do not authenticate a calculation.
SCF convergence tolerance is not a certified absolute energy-error bound.

Read-only archive inspection found that the older
`data/validation/h2-integration/modes/` and `modes-half-step/` results have energy
logs but lack explicit displaced-geometry/call-ID mappings. Log order will not
be promoted into a validated join. E2 saved evidence is intentionally force-only;
missing energies stay missing. No original evidence will be rewritten.

## Bounded checks

- Correct force sign and units; symmetric and actual asymmetric geometry.
- Polynomial exactness and three-step truncation sensitivity without a claimed
  universal bound; large-energy subtraction and unresolved saved precision.
- Raw versus constraint-masked forces, including nonzero retained anchor loads.
- Missing, failed or inconsistent records; displaced-coordinate/atom/axis errors;
  incompatible settings or energy scope and contradictory calculation identities.
- Electronic-state evidence remains separate from numerical agreement.

`test_independent_consistency.py` exercises these boundaries independently of
F1's test fixtures. `test_checkpoint_bridge.py` creates actual H1-format
checkpoints using an explicitly substituted analytic calculator, then runs the
F1 reader and analyzer. It verifies two displacement sizes, total-energy
conversion, missing totals despite available SCF-only values, preserved failed
attempt details, unchanged source bytes and zero new calculation calls during
analysis.

## Findings addressed during development

1. **Checkpoint software metadata shape.** H1 stores a nested package-version
   object; F1 initially required flat version values. The adapter now flattens
   the comparison fields while retaining the original manifest. The real-format
   fixture produces six comparisons over two steps successfully.
2. **Failed-attempt fields.** The adapter initially looked for a field named
   `failure`, while H1 saves `error` and `failed_calculation`. Both are now
   preserved, and an interrupted fixture retains the actual failure reason
   without filling its missing displacement slots.
3. **Finite-input span overflow.** For a linear extreme-range fixture with
   distances `a=b=1e308`, their sum overflowed and silently zeroed the derivative
   weights. The owner added a finite-positive-span guard. The original failed
   observation is retained in `development-overflow-failure.txt`; the new
   regression passes. This is numerical input validation, not an observed
   molecular-data problem.
4. **Uninformative agreement after energy rounding.** A constant energy offset
   of `1e16 eV` can erase the small changes in the synthetic cubic fixture.
   The saved-data residual is then zero, although its binary64 input-rounding
   sensitivity is `8 eV/angstrom`, above the explicit `0.001 eV/angstrom`
   residual criterion. F1 now reports
   `input_rounding_exceeds_declared_tolerance` and
   `numerical_resolution_findings`. The literal residual comparison stays
   separate from that warning; no true-error bound is invented.

F1 owns and applied all repairs. This directory has not altered producer or
production code. The passing checks establish the observed recovery from these
software defects and consistency with the manufactured numerical oracles.
They do not validate forces from a real electronic calculation or close the
validation protocol's full G2 gate.

## Handoff

Released to C1 uncommitted: this report, `analytic_fixtures.py`,
`test_independent_consistency.py`, `test_checkpoint_bridge.py`,
`development-overflow-failure.txt`, `test-output.txt` and `verification.json`.
F1 received all concrete findings and the final frozen-source result directly.
No original calculation evidence, producer code or other lane's files changed.
