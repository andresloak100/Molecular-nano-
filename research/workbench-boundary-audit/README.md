# V2 workbench boundary and evidence review

Assigned by C1 to support task `01a0ac0c-712e-71a0-8160-cab84676bacf`.
This review uses synthetic files in temporary directories, imports no solver,
and leaves production code and recorded scientific evidence untouched.
Implementation belongs to V1; C1 integrates the corrections. No commits are made
from this lane.

## Verification update — 2026-09-16 21:17 UTC

V1 repaired all four findings below. All six independent regression tests now
pass, including the valid singlepoint/path positive control. The final receipt
is `verification-after.json`; it records exact output and unchanged production
source hashes across the run. `audit-results.json` preserves the original
failures. The detailed findings below remain as the historical explanation of
why these checks exist. No chemical-accuracy claim follows from these passes.

## Reproduce

From the repository root:

```sh
.venv/bin/python -m unittest discover -s research/workbench-boundary-audit -p 'test_*.py' -v
```

These checks assert the desired behavior. Failures recorded in `audit-results.json`
are reproducible gaps on the recorded source hash, not intentional application
behavior or failed chemistry calculations. All fixtures are explicitly synthetic.
The ordinary workbench suite separately passed 19 tests at 2026-09-16 21:11 UTC.

## Findings handed to C1/V1

### V2-1: a previously selected campaign directory can escape its read boundary

`SafeRoot.directory()` validates a selected directory once, then returns a new
boundary containing its resolved pathname. `read_bytes()` opens that pathname
without protecting or pinning the root directory itself. The `O_NOFOLLOW`
protection applies to children only.

A deterministic fixture selects `allowed/selected`, renames that directory, and
replaces its original name with a symbolic link to `outside`. A read through the
already-created boundary returns `outside/probe.json` instead of rejecting the
replacement or retaining the original directory. No race or private file is
needed to reproduce it. A running store retains these selected boundaries across
Refresh, so this is relevant after startup too.

Expected: keep reads bound to the originally selected directory or fail closed.
Suggested repair: retain a directory descriptor and traverse children relative
to it, or use an equivalent safe root-identity strategy that also protects
ancestor replacement. Merely adding `O_NOFOLLOW` to the final root open is not
sufficient for symlink replacement of an ancestor.

Initial locations: `workbench/server.py` lines 61–83.
Regression: `test_selected_directory_replacement_cannot_escape_original_boundary`.

### V2-2: a refresh silently changes what earlier artifact IDs mean

`WorkbenchStore.catalog()` replaces one global snapshot, while structure and
source IDs are reused across snapshots. A first reader can receive a catalog
showing a 0.74 Å synthetic pose. After a second reader refreshes a legitimately
updated 1.48 Å campaign, the first reader's unchanged structure ID and source
link return the new coordinates. The two campaigns have internally consistent,
matching input and plan hashes; this is a cross-request consistency failure,
not tampering.

Expected: catalog links retain their snapshot identity, or stale requests fail
explicitly and require Refresh. Suggested repair: snapshot-specific opaque IDs
or required snapshot tokens. Keep retained snapshots bounded; alternatively
reject stale tokens. The browser currently requests a structure using only its
ID, so per-request locking alone cannot fix this case.

Initial locations: `workbench/server.py` lines 508–525;
`workbench/static/app.js` `loadStructure` and `refresh`.
Regression: `test_catalog_link_remains_consistent_after_another_reader_refreshes`.

### V2-3: equal invalid quantum dictionaries pass a completed-status check

The importer tests that quantum settings are a nonempty dictionary and that
recorded dictionaries agree. It does not validate their basic field types.
Matching `{"basis": 7, "spin": "banana"}` in the plan, design and synthetic
result is accepted as `status="completed"`, `status_verified=true` with no
integrity errors, although the production settings contract rejects both
values. All relevant hashes in this fixture match.

Expected: unsupported setting types make evidence invalid, even when the same
invalid value is repeated everywhere. Suggested repair: validate the documented
JSON settings contract without importing or running a solver, and then compare
settings with documented legacy-default handling. This does not require proving
that a chosen functional or basis is scientifically accurate.

Initial locations: `workbench/server.py` lines 382–384, 408–410, 467–468.

### V2-4: malformed path endpoints abort catalog construction

For an otherwise well-formed synthetic completed path, `endpoints[0]=null` is
recognized by the structure check but is immediately dereferenced by
`endpoint.get(...)`. A valid endpoint with `topology_screen=null` similarly
causes an uncaught `AttributeError`. The per-row exception handler does not
isolate it. The result is an aborted catalog build rather than an invalid
campaign row; another valid campaign cannot finish loading in that refresh.

Expected: an invalid row and readable error while other evidence remains
available. Suggested repair: type-check endpoint and topology objects before
accessing fields. Do not hide arbitrary application errors under an overly broad
catch.

Initial locations: `workbench/server.py` lines 234–239 and 473.

## Scope limits and checks that already work

The ordinary suite covers explicit endpoint units, changed input/plan/result
hashes, wrong result settings, missing results, fixed-anchor geometry, initial
symlink/traversal checks, denied write routes and foreign Host headers. Those
checks passed; this review does not re-report them. Earlier malformed controls
were corrected by V1 during this review and are not listed as an open finding.

Historical attempts are exposed through the raw ledger, but only the latest
attempt is numerically checked. That is a limitation, not proof that a latest
valid result is invalid: the current UI does not claim every historical attempt
was checked. A future full-history viewer should label or audit each attempt.

These are software evidence-handling findings. They do not validate a molecular
machine, reaction barrier, electronic state, or synthesis route. No quantum
calculation was launched or repeated.
