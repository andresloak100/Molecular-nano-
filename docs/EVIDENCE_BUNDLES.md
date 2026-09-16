# Portable evidence bundles

`nanodesign.bundle` copies one explicitly selected directory to a new directory
and records the exact bytes copied. This supports handoffs between researchers,
machines and review tools without rerunning a calculation or altering its record.

The exporter and verifier use only Python's standard library. They perform no
network requests, electronic calculations, campaign recovery or archive unpacking.
The API currently requires POSIX no-follow directory access (macOS or Linux).

## Create and verify

From the installed repository environment:

```sh
nanodesign bundle-create data/validation/paired-ccpvdz --out paired-reference-bundle
nanodesign bundle-verify paired-reference-bundle
```

Both commands print a JSON integrity report and return 0 for verified contents
or 2 for an invalid bundle or operational failure. An export verifies the written
bundle before reporting success. The destination's parent must already exist.
Optional positive integer bounds are `--max-files`, `--max-file-bytes`,
`--max-total-bytes`, `--max-entries`, and `--max-depth`; their defaults appear
below. An interrupted export returns 130 and may leave an incomplete directory.

The same operations are available in Python:

```python
from nanodesign.bundle import export_bundle, verify_bundle

manifest = export_bundle(
    "data/validation/paired-ccpvdz",
    "paired-reference-bundle",  # Must not exist; its parent must exist.
)
report = verify_bundle("paired-reference-bundle")
assert report["integrity_verified"], report["issues"]
```

The result has two parts:

```text
paired-reference-bundle/
  manifest.json
  payload/
    method_comparison.json
    interpretation.json
    dft/...
    ccsd_t/...
```

Every regular file and directory in the selected source tree is included,
including dotfiles, empty files/directories, failed attempts, malformed JSON,
source attribution and original filenames. Files named `manifest.json` in the
source are ordinary evidence under `payload/`. No source JSON is rewritten or
interpreted. No file is silently excluded. Select an appropriately scoped
evidence directory instead of an entire checkout or environment.

The exported tree can be moved, and verification does not need the original
source. Existing relative paths within a campaign remain relative to `payload/`.
Absolute paths and references to files outside the selected tree remain exactly
as recorded; they are not silently rewritten or followed. Reference closure and
full reproducibility are **not** established by this transport operation.

## What verification means

The version-1 manifest records each relative file path, byte length and SHA-256,
the included directory list, capture time, count/size totals and applied limits.
It is written once after capture and checks finish; the API has no manifest-update
operation. This does not make the operating-system files write-protected.

`verify_bundle` returns `integrity_verified: true` only when the supplied manifest
and current bundle inventory match, including all file lengths and hashes. It
detects changed, missing and extra files or directories, unsafe paths, symlinks,
malformed records and limits exceeded. It never repairs or deletes anything.
Invalid or incomplete bundles return `integrity_verified: false` with `issues`;
invalid API limit arguments raise `BundleError`.

Both APIs label `scientific_validation: "not_assessed"` and
`reference_closure_verified: false`. A faithful export can contain a failed job,
an incomplete campaign, contradictory records or scientific claims that need
review. Their original bytes remain evidence. Matching export hashes is neither
scientific approval nor authentication of the source: somebody who changes both
the payload and manifest can produce a matching new bundle. Save the returned
`manifest_sha256` from verification separately when an external digest is needed.

An exported `.campaign.lock` is only a copied file. It does not transfer a running
process's lock, prove job liveness or authorize simultaneous workers.

## Bounds and failure handling

```python
from nanodesign.bundle import BundleLimits

limits = BundleLimits(
    max_files=10_000,
    max_file_bytes=32 * 1024**2,
    max_total_bytes=256 * 1024**2,
    max_entries=20_000,  # Files plus directories in the source/payload.
    max_depth=32,
)
report = verify_bundle("paired-reference-bundle", limits=limits)
```

These are the defaults. Bounds must be positive integers; file size zero is still
allowed. Verification applies its caller's limits independently of limits written
in an untrusted manifest. The manifest has a separate 16 MiB read/write bound.
Paths have a 4,096-byte bound; traversal components, absolute manifest paths and
backslashes are rejected. All source entries must be regular files or directories;
symlinks and devices/sockets/FIFOs are rejected rather than followed or skipped.

The destination must be new and outside the source. Existing destinations,
including empty directories, are never reused or overwritten. A failure before
manifest publication can leave an incomplete new destination without
`manifest.json`; it cannot verify as a complete bundle. An interruption during
publication may leave an extra `.manifest.pending`, which verification also
rejects. Inspect or preserve such a directory and choose a new destination for
a new attempt. The exporter never cleans up an arbitrary existing directory.

Use completed, quiescent evidence where possible. The exporter pins directory
access, rejects detected file/root replacement and rechecks source inventory and
bytes. Nevertheless, copying a filesystem tree is **not a globally simultaneous
snapshot**. The manifest explicitly records `snapshot.atomic: false` and
`snapshot.source_rechecked: true`; undetected changes outside those observation
windows remain possible. Verification likewise describes the bytes it checked,
not an ongoing monitor of future writes.

## Focused verification

```sh
.venv/bin/python -m pytest -q tests/test_bundle.py
```

These tests use synthetic files only. Independent adversarial review lives under
`research/bundle-audit/`; it has a separate owner from the implementation.
