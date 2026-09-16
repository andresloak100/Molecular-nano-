# Molecular design workbench

A read-only, offline workspace for inspecting real atom coordinates, pose studies,
and saved calculation evidence. It does not run calculations, change designs,
control hardware, or infer that a molecular machine can be manufactured.

From the repository root on macOS or Linux, start it with Python 3.11 or newer:

```sh
.venv/bin/python workbench/server.py --port 8765
```

Open [the local workbench](http://127.0.0.1:8765). Stop it with Ctrl-C. The server
binds only to `127.0.0.1`; there is no external hosting or cloud service. Plain
`python3` also works when its version is recent enough. No package installation,
API key, GPU, browser extension, or internet connection is required.

## What is loaded

- Exact archived initial coordinates for the 53-atom C22H31 candidate, with six
  fixed carbon anchors and the three identified hydrogen-transfer atoms.
- Both recorded direct and density-fitting single-point calculations, their
  method settings, force/spin diagnostics, interpretation sidecars and raw JSON.
- The nine pending poses in `examples/pose-campaign`, including their initial and
  proposed final coordinates. Selecting a pose loads its actual saved file.
- Both paired `minao` and `atom` CC reference calculations. These found different
  HF solutions; neither record verifies ground-state or reaction-saddle identity.

Drag the structure to rotate it, scroll to zoom, and select atoms to inspect
their supplied coordinates and pair distances in ångströms. Display changes do
not alter saved coordinates. Bond lines use nominal design connectivity when it
exists, otherwise a labeled covalent-radius distance heuristic; neither is a
computed electronic bond-order measurement.

The candidate's path and vibration results remain **not computed**. A completed
single point is neither a relaxed geometry nor a reaction barrier. Elapsed times
come from overlapping jobs under different workloads and are not a controlled
speed comparison. Missing and failed records remain visible instead of becoming
zeroes or favorable scores. Reload the page to read a new evidence snapshot; no
campaign recovery or status-changing code is invoked.

## Select a local campaign

Campaigns are selected at startup, through an explicit directory boundary:

```sh
.venv/bin/python workbench/server.py --campaign-root /absolute/path/to/studies --campaign study-one --campaign study-two
```

Each `--campaign` must be a relative directory inside `--campaign-root` with the
existing `plan.json`, `campaign.json`, and saved design files. Repeated flags load
multiple campaigns. Without flags, the bundled nine-pose study is loaded. To
select a campaign within this repository, use `--campaign relative/directory`.
`--project-root` can select another checkout for the fixed archive locations.

Absolute imported artifact paths, `..` traversal, symlinks inside a selected
boundary, nonregular files, files larger than 8 MiB, more than 100 campaign poses,
and structures larger than 1,000 atoms are rejected. All HTTP reads use registered
artifact IDs; the server has no arbitrary file browser or write/job API. Source
links return snapshot bytes, with structure SHA-256 hashes available in the API.
Links carry a snapshot token and expire when another page loads a new catalog;
reload the page if a source or structure link reports that it has expired. A
selected campaign directory that is replaced after startup is rejected instead
of following a replacement symlink or silently changing the read boundary.
Source JSON is data, never executable markup. Imported status claims are displayed
as recorded evidence; a matching plan hash is not scientific validation.

Imported designs must explicitly declare `schema_version: 1` and
`length_unit: "angstrom"`. The viewer verifies plan, design, coordinate and
recorded result hashes, matching method/optimization settings, and basic completed
energy/force evidence. Missing or changed completed evidence is labeled **invalid
evidence**, with the original ledger status preserved for inspection. These checks
establish consistency of the saved records, not their chemical accuracy. The
curated original direct archive predates the design unit field; its extended-XYZ
coordinates use the documented ASE ångström convention, with that historical
context preserved in its interpretation sidecar.

The historical paired reference has a separate, source-hash-bound retrospective
annotation. Its inferred `minao` CC guess is labeled as an inference, while the
original missing guess remains “Not recorded.” An invalid source binding cannot
supply inferred settings. DFT and CC initial guesses are displayed independently;
the `atom` reference label describes only the explicitly recorded CC guess.

## Verify

The server and tests use only Python's standard library; the browser UI uses
local HTML, CSS, JavaScript, and a canvas projection without a remote dependency.

```sh
.venv/bin/python -m unittest workbench.test_server -v
```

Tests exercise the actual archived counts, coordinates, anchors, 3.6 Å default
donor–apex separation, nine pending poses, source values/hashes, missing and
failed results, campaign import bounds, symlink rejection, read-only HTTP routes,
and loopback host checks. They do not run electronic-structure calculations.

The canvas viewer is for inspection, not a force-field simulation or an editable
molecular CAD system. Coordinates are never relaxed or interpolated for display.
The first XYZ frame is used; initial and final endpoint files are separate views.
