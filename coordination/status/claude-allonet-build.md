# AlloNet build session (Claude Code remote, 2026-09-17)

New top-level directory `allonet/` — user-directed, outside all science
lanes, new files only. No existing lane's files were touched and no quantum
jobs were started.

AlloNet is a standalone, zero-dependency Python library ("the trust layer
for AI x biology"): bio identifier validation with check digits, sequence
utilities and database checksums (CRC64-ISO as UniProt publishes, SEGUID),
provenance-recording fetch clients for UniProt/RCSB/NCBI/Ensembl/doi.org,
and claim verification returning CONFIRMED/REFUTED/UNVERIFIABLE with an
evidence chain. 560 offline tests pass without network; a live contract
script (`allonet/scripts/live_check.py`) passes 14/14 against the real
APIs.

It lives here only until it has its own repository; it imports nothing
from `nanodesign` and nothing here imports it. Agents in this repository
are welcome to use it — the Evidence pattern (URL + UTC time + sha256 of
exact bytes next to every published number) is the same discipline
`FINDINGS.md` already asks for.
