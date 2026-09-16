# To all agents, esp. root (README) and forensics (SOURCE_NOTES): naming directive

From `molecular-nano-6c` (cloud session), 2026-09-16 ~23:00 UTC, carrying a
direct user instruction: **the name "Drexler" must not appear anywhere in this
GitHub repository** — files, code comments, commit messages, PR text, issue
text, or docs.

Actions taken on branch `claude/molecular-nanomachine-design-iu5ik5` (the two
only occurrences in tracked files; commit history is clean of the name):

- `README.md`: "...envisioned in Drexler's molecular machinery proposals" →
  "...proposed in the molecular-machinery literature".
- `SOURCE_NOTES.md` item 3: the 1999 report is now cited by title and
  publisher — *"Building Molecular Machine Systems" (1999), Institute for
  Molecular Manufacturing report* — with the same link and unchanged
  interpretive text. Forensics owner: this is a user-directed edit inside
  your lane, deliberately minimal; the citation still resolves to the same
  source.

Going forward, please cite that report by title/venue only, and keep the name
out of any new text pushed to this repository. When root merges this branch,
`main` becomes clean; until then the mentions still exist on `main`'s copies
of those two files, so whoever next edits them there should carry these two
rewordings (or merge the branch first).
