# The AlloNet protocol for AI agents

Three rules. They cost milliseconds and prevent the three ways AI agents
poison biological work.

## Rule 1 — Never emit an identifier that fails `identify()`

LLMs invent plausible-looking accessions. Format-validate every identifier
*you are about to output*, not just ones you receive:

```python
import allonet

matches = allonet.identify(candidate_id)
if not matches:
    # This string was never issued by any database you know.
    # Do not output it. Say you could not verify the identifier.
    ...
```

`identify()` returning a match means the *format* is real. For anything that
matters, escalate to Rule 3 and fetch the record — a well-formed accession can
still be unissued (fetching it 404s, which `FetchError` surfaces loudly).

## Rule 2 — Never attribute a sequence without a passing check

"The sequence of human hemoglobin alpha is …" must be backed by either a live
comparison:

```python
v = allonet.verify_uniprot_sequence("P69905", seq)
assert v.ok, v.detail
```

or, offline, a checksum comparison against the entry's published CRC64:

```python
v = allonet.verify_sequence_checksum(seq, "15E13666573BBBAE")
```

A `REFUTED` verdict comes with the differing-position detail — show it to the
user instead of papering over it.

## Rule 3 — Never publish a fetched fact without its Evidence

Every AlloNet fetch returns an `Evidence` record: URL, UTC timestamp, HTTP
status, SHA-256 of the exact bytes received. Attach it to the claim:

```python
seq, ev = allonet.uniprot_sequence("P69905")
report["claims"].append({
    "claim": "P69905 canonical sequence, 142 aa",
    "evidence": ev.as_dict(),
})
```

Six months later, anyone can re-fetch the URL and diff the hash. A claim with
evidence is auditable; a claim without is an anecdote.

## Reproducible runs

Pin the cache and go offline for reruns:

```console
ALLONET_CACHE=./evidence_cache ALLONET_OFFLINE=1 python your_pipeline.py
```

Cache misses now raise instead of silently refetching, so a "reproduction"
cannot quietly become a different dataset. Commit the cache directory next to
your results if your evidence policy allows it.

## Numbers

Report energies with the conversion done by the library, not by mental
arithmetic under deadline:

```python
allonet.kd_to_delta_g(1e-9)          # -12.28 kcal/mol at 298.15 K
allonet.kcal_to_kj(-12.28)           # -51.4 kJ/mol
```

And state the temperature — `kd_to_delta_g(kd, temp_k=310.15)` for
physiological claims.
