# AlloNet

**The trust layer for AI × biology.** One file, zero dependencies.

Every AI agent that touches biological data eventually states something
confidently wrong: an accession number that was never issued, a sequence that
doesn't match the database it's attributed to, a binding energy with no source
attached. These failures are silent, they propagate into papers, orders and
downstream models, and every team re-implements the guards against them badly,
or not at all.

AlloNet is that guard, done once, properly:

```
pip install allonet        # or just copy allonet.py into your project
```

```python
import allonet

# 1. Is this identifier even real? (the hallucination check)
allonet.identify("P69905")        # [IdMatch(database='uniprot', ...)]
allonet.identify("XP123NOTREAL")  # [] — never issued by anyone. Alarm.
allonet.is_valid("50-78-2", "cas")   # True  (check digit verified)
allonet.is_valid("50-78-3", "cas")   # False (format ok, checksum wrong)

# 2. Checksum sequences the way the databases do
allonet.crc64_iso(seq)   # the checksum UniProt publishes for every entry
allonet.seguid(seq)      # SEGUID, Biopython-compatible

# 3. Fetch with provenance — every byte gets an Evidence record
seq, evidence = allonet.uniprot_sequence("P69905")
evidence.as_dict()
# {'url': 'https://rest.uniprot.org/uniprotkb/P69905.fasta',
#  'status': 200, 'sha256': 'aecc3cb7…', 'size': 269,
#  'retrieved_at': '2026-09-17T01:12:03Z', 'from_cache': False}

# 4. Verify claims against the primary source
v = allonet.verify_uniprot_sequence("P69905", seq)   # verdict: CONFIRMED
v = allonet.verify_pdb_uniprot_link("1MBN", "P02185")
v = allonet.verify_doi("10.1093/nar/gkae1010")
v = allonet.verify_gene_symbol("TP53")
v.verdict   # CONFIRMED / REFUTED / UNVERIFIABLE — never a silent guess
v.evidence  # the fetches that back the verdict, hashes and all
```

Or from the command line:

```console
$ allonet id Q9Y6K9
$ allonet seq stats hemoglobin.fasta
$ allonet checksum hemoglobin.fasta
$ allonet verify uniprot-seq P69905 hemoglobin.fasta
$ allonet verify gene TP53
$ allonet units kd2dg 1e-9        # → -12.28 kcal/mol at 298.15 K
```

## What's inside

| Layer | What it does |
|---|---|
| **Identifiers** | Recognize/validate/normalize UniProt, PDB (incl. extended `pdb_`), RefSeq, GenBank, Ensembl, GO, EC, DOI, PMID, dbSNP, ChEMBL, InChIKey, CAS, Pfam, InterPro, HGNC, NCBI Taxonomy, ORCID — with check-digit verification where the scheme defines one. |
| **Sequences** | Clean copy-paste noise, guess/validate alphabets (full IUPAC), reverse-complement, transcribe, translate (NCBI table 1), GC content, protein MW, strict FASTA round-trip. |
| **Checksums** | CRC64-ISO (what UniProt publishes), SEGUID, SHA-256, MD5 — compare sequences to database entries without re-downloading them. |
| **Fetching** | Polite client for UniProt, RCSB PDB, NCBI E-utilities, Ensembl, doi.org: on-disk cache, per-host rate limits, retries with backoff, HTTPS-only host allowlist, offline mode for reproducible runs. Every response wrapped in an `Evidence` record (URL, UTC time, SHA-256 of the exact bytes). |
| **Verification** | Claims turned into checks against primary sources, returning `CONFIRMED` / `REFUTED` / `UNVERIFIABLE` plus the evidence chain — never a plausible-sounding guess. |
| **Units** | The conversions people get wrong under deadline: kcal↔kJ (thermochemical), Kd↔ΔG at temperature. |

## Design principles

- **Zero dependencies.** Standard library only. Vendor the single file if you
  don't want a dependency at all.
- **Loud failure beats silent damage.** Partial codons raise instead of being
  dropped; ambiguous residues refuse to average into a molecular weight;
  unknown checksum formats say `UNVERIFIABLE` instead of guessing.
- **Evidence or it didn't happen.** Anything fetched carries the URL,
  timestamp and hash of the exact bytes. Put `evidence.as_dict()` next to
  every number you publish.
- **Offline mode is a feature.** `ALLONET_OFFLINE=1` makes reruns byte-stable
  and network-free — cache misses raise instead of silently refetching.
- **Read-only.** AlloNet validates data; it never modifies anything remote.

## For AI agents

If you are an AI agent (or you build them), adopt the three-rule protocol in
[`examples/AGENT_PROTOCOL.md`](examples/AGENT_PROTOCOL.md):

1. Never emit a database identifier that fails `identify()`.
2. Never attribute a sequence to an entry without `verify_uniprot_sequence`
   (or a checksum comparison) passing.
3. Never publish a fetched fact without its `Evidence` record.

## Environment

| Variable | Effect |
|---|---|
| `ALLONET_CACHE` | cache directory (default `~/.allonet/cache`) |
| `ALLONET_OFFLINE` | `1` = cache only, never touch the network |
| `ALLONET_NCBI_EMAIL` | contact email forwarded to NCBI (their policy recommends it) |
| `ALLONET_NCBI_KEY` | NCBI API key (raises your rate limit) |

## Testing

```console
python -m pytest        # offline suite; no network required
python scripts/live_check.py   # optional live contract check against the real APIs
```

## License

MIT.
