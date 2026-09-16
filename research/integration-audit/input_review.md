# Q1 production input and campaign review

Checked 2026-09-16 21:06 UTC against checkout HEAD `713db7e` and the live working
tree. Owner: `codex-support-q1-inputs`, internal helper for `codex-support-q1`.
No production files, source evidence, or other agents' files were changed. No
quantum calculations were started. Workflow/force/path review is delegated to
the separate `codex-b717` lane; this report covers input and campaign handling.

Three actionable defects reproduced. They do not establish that existing
scientific calculations are wrong; they expose ways input identity or preserved
failure evidence can become unreliable, or accepted designs cannot be copied.

**Repair verification, 2026-09-16 21:12 UTC:** C1 fixed Q1-I1 and Q1-I3 in the
shared working tree and delegated Q1-I2 to the parent `codex-support-q1`, which
implemented immutable input-byte snapshots. All the original desired-behavior
probes now pass; each strict xfail marker was removed only after observing the
corresponding strict XPASS. Latest bounded run: **7 passed** in 1.61 seconds.
All three findings are resolved in the checked working tree. The descriptions
below preserve the original reproduced behavior and repair rationale.

## Q1-I1 — valid symlinked coordinate input cannot be snapshotted (P2)

Locations: `nanodesign/campaign.py:102–105` and `:113–116`.

`create_campaign` chooses the copied filename from the **resolved target's**
suffix, then writes its runnable design using the **original reference's**
suffix. A valid design whose `initial.xyz` symlink points to
`actual-initial.extxyz` passes `load_design`. Campaign creation copies the exact
bytes to `initial.extxyz`, but writes `"initial": "initial.xyz"` and fails with
`FileNotFoundError` when validating its own snapshot. The already-created output
directory is left partial, so rerunning at the same output is also rejected.

Fix direction: retain the actual destination filename selected during the copy
and use that exact name in the runnable design, rather than deriving it a second
time from a different path. Preserve the original design separately as before.

Probe: `test_campaign_snapshots_valid_symlink_with_different_extension`.

## Q1-I2 — input hashes can identify bytes other than those parsed (P2)

Locations: `nanodesign/design.py:44`, `:51`, and `:90`.

`load_design` parses the JSON and coordinate files, validates the objects, then
opens all three files again to compute their hashes. An edit between parsing
and hashing produces a valid-looking result whose provenance identifies a
different input. The deterministic probe inserts a source edit immediately
after the coordinate reader returns; this simulates an ordinary concurrent
editor or generator, without probabilistic timing:

- Coordinate case: returned H–H distance is **0.74 Å**, while
  `initial_sha256` is the hash of replacement coordinates at **0.90 Å**.
- Design case: returned settings request **def2-svp**, while `design_sha256`
  is the hash of replacement JSON requesting **sto-3g**.

These mismatches pass all current `load_design` checks. A workflow could then
calculate the already-loaded structure under the already-loaded settings while
recording hashes that claim different inputs. This is a demonstrated race
condition, not a claim that it occurred in existing archives.

Fix direction: parse and hash the same immutable byte snapshots, including the
design document, or fail explicitly if stable input capture cannot be
established. Preserve the full supported ASE input-format behavior when choosing
the snapshot/parser mechanism. Do not simply add another late hash of the files.

Probes: `test_load_design_rejects_or_preserves_consistency_during_source_edit`
for `coordinates` and `design`.

## Q1-I3 — earlier attempt hashes stop protecting evidence after retry (P2)

Locations: `nanodesign/campaign.py:258–260` and `:310–317`.

Both campaign execution and reporting inspect only the last attempt for each
pose. A failed first attempt receives a stored `result_sha256`; a subsequent
explicit retry completes. Deleting or replacing the first attempt's
`result.json` is then silently accepted by `campaign_report`, which reports the
pose completed with two attempts. The stored historical checksum is never
checked. The failed evidence that should explain why a retry occurred is lost
without an integrity error. Checking the latest result still works, as the
positive control confirms.

Fix direction: validate recorded historical attempts as part of the shared
campaign integrity load, including the existing stored hashes and missing-file
rules, before reporting or continuing. Retain the intentional exception for
unhashed incomplete/truncated abandoned evidence; these probes concern terminal
records with an already-recorded hash.

Probes: `test_campaign_still_checks_historical_recorded_attempts` for
`changed` and `missing`.

## Bounded verification

The owned test file is `research/integration-audit/test_input_findings.py`.
It creates temporary H2 inputs and synthetic JSON records only. It does not
instantiate a physical calculator.

```
.venv/bin/python -m pytest -q research/integration-audit/test_input_findings.py -rx
```

Result: **2 passed, 5 xfailed** in 1.47 seconds. The two positive controls prove
that changing current campaign coordinates or the latest terminal result is
already rejected. The five strict expected failures are desired-behavior
regression assertions covering the three defects above, not passing assertions
that encode faulty behavior. A second bounded run with `--runxfail --tb=short
--disable-warnings` confirmed **5 intended failures and 2 passes** in 1.45
seconds: missing symlink snapshot file, two mismatched hashes, and two missing
expected integrity errors. Warnings are existing ASE/NumPy deprecations.

After the authorized repairs, all five strict expected failures became strict
XPASS as intended. Their markers were removed; the final suite passes all seven
checks. The tests are kept outside the production suite until owner review.
No commit was made by this helper.

The race probe intercepts `design.read` and edits the original source immediately
after its first return. The repaired reader passes an in-memory stream to that
same entry point; the source mutation still occurs and the returned hashes now
identify the original bytes actually parsed. Both the coordinate and settings
race variants pass without weakening their expected provenance assertions.
