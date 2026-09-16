# F1: saved energy and force consistency

This research tool asks whether recorded raw forces agree with the change in
recorded **total** energy when one free Cartesian coordinate is displaced. It
does not acquire energies or forces, run a quantum solver, or change an input
archive. It complements force-to-Hessian reconstruction (E2) and mode comparison
(N1), which test different relationships.

The current scientific motivation is gate G2 in
[the validation protocol](../../docs/VALIDATION_PROTOCOL.md). A wrong force sign,
missing energy contribution, mismatched method or unsuitable displacement can
spoil an optimization or curvature calculation. Agreement at a few saved points
is useful numerical evidence, but does not prove electronic branch continuity,
chemical accuracy, reaction connectivity or an operating molecular machine.

## Use existing checkpoint evidence

With the repository environment installed, run:

```sh
.venv/bin/python research/force-energy-consistency/consistency.py \
  --checkpoints /path/to/H1/coarse /path/to/H1/fine \
  --output /path/to/new-force-energy-report.json
```

This reads the existing H1 `force_checkpoint.json` format under its shared
execution lock and verifies its native envelope, manifest and record bindings.
It uses each accepted record's `total_energy_hartree` and recorded `hartree_eV`
conversion. That total includes the configured dispersion contribution; SCF-only
`e_tot` or `dft_energy_hartree` is never substituted for it. The raw force is the
one saved for that exact geometry and calculation call. New files are written
only when an explicit, previously absent output path is supplied.

There is deliberately no acquisition command here. H1 remains a separate research
prototype with its own ownership and resource limits. Existing force-only records
produce an unavailable comparison. Missing native checkpoint files are unsupported
input; this tool never invents a stencil from electronic-log order. In particular,
the historical `data/validation/h2-integration/modes*` archives lack the required
per-displacement call/geometry binding and cannot be used to infer this check.

For another producer, supply the explicit [normalized schema](SCHEMA.md):

```sh
.venv/bin/python research/force-energy-consistency/consistency.py \
  --input research/force-energy-consistency/synthetic-example.json
```

`synthetic-example.json` and its report demonstrate mathematics using an exact
polynomial fixture. They are **not molecular results**. An optional
`--force-tolerance VALUE` in eV/angstrom compares observed residuals against a
caller-declared numerical threshold; no default scientific tolerance is chosen.
Exit codes: 0 means every submitted stencil produced a comparison; 1 means some
evidence was unavailable; 2 means contradictory/malformed input. A code of 0
does not assert that the residual met a threshold or that any scientific gate
passed. Read the report's individual fields.

## What is calculated

Let the exact saved positions have positive distances
`a = x0 - x_minus` and `b = x_plus - x0`. The reported energy derivative is
the derivative at `x0` of the three-point quadratic interpolant:

```text
D = [b/(a+b)] * (E0 - E_minus)/a
  + [a/(a+b)] * (E_plus - E0)/b
F_FD = -D
residual = F_FD - F_raw(x0)
```

When `a=b=h`, this becomes `-(E_plus-E_minus)/(2h)`. A simple secant divided by
`a+b` is biased at the baseline for unequal steps, even for a quadratic energy;
it is not silently used. Requested and actual offsets, asymmetry and the three
energy weights are retained. A displacement must change exactly the requested
free coordinate; no alignment, atom reordering, off-axis motion or displacement
of a fixed atom is allowed. The one-part-per-million representation guard is a
coordinate-format check, not a force-accuracy budget. Overflowing arithmetic is
rejected rather than producing an artificial zero.

For each coordinate, the report gives the spread of finite-difference forces
and separately the spread of baseline analytical forces across distinct saved
stencils. One available stencil has null spreads. These are observed sensitivities;
no asymptotic convergence, Richardson extrapolation or error bound is claimed.

The report also propagates half a binary64 unit-in-the-last-place for each stored
energy through the interpolation weights. This is only an **input-rounding
sensitivity indicator**, excluding subtraction/division rounding, unknown solver
noise, SCF error and finite-step truncation. When it exceeds the declared
tolerance, `input_rounding_exceeds_declared_tolerance` and a finding explicitly
warn that a small observed residual cannot provide an informative precision
check. `within_declared_tolerance` still describes the literal residual only.
The SCF convergence setting is never used as an energy-error bound.

These design choices follow the standard finite-difference tradeoff between
step size and numerical precision described in
[SciPy's differentiation documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.differentiate.derivative.html).
PySCF also supplies an
[analytical-versus-finite-difference gradient example](https://github.com/pyscf/pyscf/blob/master/examples/grad/16-scan_force.py).
This tool uses the formula above directly; it does not invoke either package's
derivative or electronic calculation routines.

## Interpretation and provenance

`comparison_produced`, `partial`, `unavailable` and `invalid` distinguish usable
numerical comparisons, incomplete evidence and contradictory data. A failed
calculation remains unusable even when an energy was recorded before failure.
Native H1 errors and failed-calculation details are retained in the normalized
document's source records. All source checkpoint bytes remain unchanged.

The normalized pure API checks the **declared** source hashes, call IDs, exact
geometry, units and repeated context. It does not open arbitrary referenced
files or authenticate a physical calculation. The H1 adapter additionally checks
the native bytes and embedded integrity digests. Either record can be fabricated
by its author; hashes are not an electronic-state or experimental certificate.

Opaque state-evidence IDs are carried beside their calculation call and artifact
digest. They are not interpreted as branch assignments. Reusing an ID across
different geometries is flagged. Matching starting guesses or spin labels is
also insufficient. The report always leaves electronic-state identity, branch
continuity, numerical convergence and scientific model validation unverified.
Future state-snapshot consumers need their own reviewed comparison.

The current producer does not record the complete resolved quadrature and
auxiliary-basis specification. `resolved_numerics: null` exposes that gap;
recorded response flags and available software versions are retained separately.
Agreement can therefore support only a conditional numerical observation, not
closure of the full G2 gate.

## Software verification and ownership

```sh
.venv/bin/python -m pytest -q research/force-energy-consistency --disable-warnings
```

Tests use conservative polynomial energies, inconsistent forces, unequal steps,
rounding/overflow cases and the actual H1 producer with an explicitly synthetic
calculator. Producer/consumer checks include source preservation, full-precision
coordinates, raw anchor forces, partial failures and missing total energies.
No quantum calculation has been launched for F1. Independent F2 review is in
`research/force-energy-review/`; its receipts describe its tested source revision.

Final owned selection: **42 tests passed** in 3.71 s. Independent F2 selection:
**39 tests passed** in 2.33 s, with unchanged implementation/test hashes and no
unresolved finding in its bounded review. These were separate runs. F2 preserves
the original overflowing-span reproduction alongside its passing verification.
`verification.json` records the owned check and source hashes; the independent
receipt is under `research/force-energy-review/`. Existing ASE/NumPy deprecation
warnings in the native producer path are not scientific findings.

F1 handed this directory to C1 on 2026-09-16 for review/integration, leaving files
uncommitted as requested. C1 owns subsequent changes and any production wiring.
Existing source evidence and other agents' implementations are unchanged.
