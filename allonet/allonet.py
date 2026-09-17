"""AlloNet — the trust layer for AI x biology.

Every AI agent that touches biological data eventually states something
confidently wrong: a hallucinated accession number, a sequence that doesn't
match the database, a number with no source attached. AlloNet is the
correctness layer that catches those failures before they propagate.

It is a single file with zero dependencies beyond the Python standard
library, so any agent or pipeline can vendor it. It does four things:

1.  **Identifiers** — recognize, validate and normalize the identifiers of
    the major biological databases (UniProt, PDB, RefSeq, Ensembl, GO, EC,
    CAS, InChIKey, DOI, PubMed, dbSNP, Pfam, InterPro, ChEMBL, ORCID, ...),
    including check-digit verification where the scheme has one. An ID that
    fails `identify()` was never issued by anyone — the canonical
    hallucination check.

2.  **Sequences** — clean copy-pasted DNA/RNA/protein sequences, guess and
    validate their alphabet, reverse-complement, transcribe, translate
    (NCBI table 1), and compute the checksums the databases themselves use
    (CRC64-ISO as in UniProt, SEGUID, SHA-256, MD5) so a sequence can be
    compared to a database entry without shipping the sequence.

3.  **Provenance-first fetching** — a polite HTTP client (caching, per-host
    rate limits, retries with backoff) for UniProt, RCSB PDB, NCBI
    E-utilities and Ensembl. Every response is wrapped in an `Evidence`
    record: URL, UTC timestamp, HTTP status, SHA-256 of the exact bytes.
    A claim with an Evidence record is auditable; one without is an
    anecdote.

4.  **Verification** — turn common factual claims into checks against the
    primary source: "this is the sequence of UniProt P69905", "PDB 1MBN is
    a structure of P02185", "this DOI exists". Each returns a
    `Verification` with a verdict (CONFIRMED / REFUTED / UNVERIFIABLE) and
    the evidence chain.

Command line: ``python allonet.py --help`` (or the ``allonet`` entry point).

Environment variables:
    ALLONET_CACHE       cache directory (default: ~/.allonet/cache)
    ALLONET_OFFLINE     "1" = never touch the network, serve cache only
    ALLONET_NCBI_EMAIL  contact email forwarded to NCBI E-utilities (optional,
                        recommended by NCBI; never defaulted for you)
    ALLONET_NCBI_KEY    NCBI API key (optional, raises rate limits)

This module validates data; it never modifies anything remote.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import http.client
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

__version__ = "0.1.0"
__all__ = [
    # identifiers
    "IdMatch", "identify", "is_valid", "normalize",
    # sequences
    "clean_sequence", "guess_seq_type", "validate_sequence",
    "reverse_complement", "transcribe", "translate", "gc_content",
    "protein_mw", "parse_fasta", "write_fasta",
    # checksums
    "crc64_iso", "seguid", "seq_sha256", "seq_md5",
    # fetching
    "Evidence", "Record", "Fetcher", "get_fetcher",
    "uniprot_fasta", "uniprot_json", "uniprot_sequence",
    "rcsb_entry", "rcsb_polymer_entities", "pdb_structure",
    "ncbi_efetch", "ncbi_esearch", "ensembl_get",
    # verification
    "Verification", "verify_sequence_checksum", "verify_uniprot_sequence",
    "verify_pdb_uniprot_link", "verify_doi", "verify_gene_symbol",
    # units
    "kcal_to_kj", "kj_to_kcal", "kd_to_delta_g", "delta_g_to_kd",
    "R_KCAL_PER_MOL_K",
]


# ---------------------------------------------------------------------------
# 1. Identifiers
# ---------------------------------------------------------------------------

# Anchored patterns for identifiers whose format is published by the issuing
# database. Sources: UniProt help (accession format), wwPDB (entry ID and
# extended entry ID), NCBI RefSeq accession prefixes, Ensembl stable IDs,
# GO/EC/CAS/InChIKey/DOI/ORCID public specifications.
_ID_PATTERNS: dict[str, tuple[re.Pattern[str], str]] = {
    "uniprot": (
        re.compile(
            r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
            r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
            r"(?:-(?P<isoform>[1-9][0-9]*))?",
            re.IGNORECASE,
        ),
        "UniProtKB accession (optionally with isoform suffix)",
    ),
    "pdb": (
        re.compile(r"[1-9][A-Za-z0-9]{3}"),
        "PDB entry ID (4-character)",
    ),
    "pdb_extended": (
        # wwPDB's official regex: pdb_ + 8 alphanumerics.
        re.compile(r"pdb_[A-Za-z0-9]{8}", re.IGNORECASE),
        "PDB extended entry ID",
    ),
    "refseq": (
        # NZ_ (genomes/WGS) embeds an INSDC-style accession with letters.
        re.compile(
            r"(?:(?:AC|AP|NC|NG|NM|NP|NR|NT|NW|WP|XM|XP|XR|YP)_[0-9]+"
            r"|NZ_[A-Z]{0,6}[0-9]+)"
            r"(?:\.[0-9]+)?"
        ),
        "NCBI RefSeq accession",
    ),
    "genbank_protein": (
        # NCBI issues 3 letters + 5 digits (legacy) or 3 + 7, never 3 + 6.
        re.compile(r"[A-Z]{3}(?:[0-9]{5}|[0-9]{7})(?:\.[0-9]+)?"),
        "GenBank/DDBJ/ENA protein accession",
    ),
    "ensembl": (
        re.compile(
            r"ENS(?:[A-Z]{0,5}(?:G|T|P|E|R)[0-9]{11}"
            r"|(?:FM|GT)[0-9]{14})(?:\.[0-9]+)?"
        ),
        "Ensembl stable ID (gene/transcript/protein/exon/regulatory/"
        "gene-tree/family)",
    ),
    "go": (
        re.compile(r"GO:[0-9]{7}"),
        "Gene Ontology term",
    ),
    "ec": (
        # Dash placeholders are only valid from the tail inward (1.2.-.- is
        # real, 1.-.3.4 is hierarchically impossible).
        re.compile(
            r"(?:EC[ :])?[1-7]"
            r"\.(?:[0-9]{1,3}\.(?:[0-9]{1,3}\.(?:n?[0-9]{1,3}|-)|-\.-)"
            r"|-\.-\.-)"
        ),
        "Enzyme Commission number",
    ),
    "doi": (
        # Angle brackets are legal in legacy SICI-style DOIs, and the prefix
        # may carry dot-separated sub-elements.
        re.compile(
            r"(?:doi:|https?://(?:dx\.)?doi\.org/)?"
            r"10\.[0-9]{4,9}(?:\.[0-9]+)*/[^\s\"]+",
            re.IGNORECASE,
        ),
        "Digital Object Identifier",
    ),
    "pmid": (
        re.compile(r"PMID:? ?[0-9]{1,9}", re.IGNORECASE),
        "PubMed ID (only recognized with PMID prefix)",
    ),
    "dbsnp": (
        re.compile(r"rs[0-9]{1,12}"),
        "dbSNP reference SNP ID",
    ),
    "chembl": (
        re.compile(r"CHEMBL[0-9]+"),
        "ChEMBL compound/target ID",
    ),
    "inchikey": (
        # 14 chars, hyphen, 8 chars + flag (S/N) + version (A), hyphen,
        # protonation char.
        re.compile(r"[A-Z]{14}-[A-Z]{8}[SN]A-[A-Z]"),
        "InChIKey (structure hash; format-valid does not mean the compound exists)",
    ),
    "cas": (
        re.compile(r"[1-9][0-9]{1,6}-[0-9]{2}-[0-9]"),
        "CAS Registry Number (check digit verified)",
    ),
    "pfam": (
        re.compile(r"PF[0-9]{5}(?:\.[0-9]+)?"),
        "Pfam family",
    ),
    "interpro": (
        re.compile(r"IPR[0-9]{6}"),
        "InterPro entry",
    ),
    "hgnc": (
        re.compile(r"HGNC:[0-9]+"),
        "HGNC gene ID",
    ),
    "ncbi_taxon": (
        re.compile(r"(?:taxid:|NCBITaxon:)[0-9]{1,7}"),
        "NCBI Taxonomy ID (only recognized with a taxid:/NCBITaxon: prefix)",
    ),
    "orcid": (
        # ORCID's allocated ISNI blocks start 0000- or 0009-.
        re.compile(r"(?:(?:https?://)?orcid\.org/)?"
                   r"000[09]-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]"),
        "ORCID researcher ID (MOD 11-2 check digit verified)",
    ),
}


@dataclass(frozen=True)
class IdMatch:
    """One plausible interpretation of an identifier string."""

    database: str          # key into _ID_PATTERNS
    normalized: str        # canonical form (case, prefixes fixed)
    description: str
    checksum_ok: Optional[bool] = None   # None = scheme has no check digit
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _cas_checksum_ok(cas: str) -> bool:
    digits = cas.replace("-", "")
    body, check = digits[:-1], int(digits[-1])
    total = sum(int(d) * i for i, d in enumerate(reversed(body), start=1))
    return total % 10 == check


def _orcid_checksum_ok(orcid: str) -> bool:
    digits = orcid.replace("-", "")
    total = 0
    for ch in digits[:-1]:
        total = (total + int(ch)) * 2
    remainder = total % 11
    result = (12 - remainder) % 11
    expected = "X" if result == 10 else str(result)
    return digits[-1] == expected


def identify(text: str) -> list[IdMatch]:
    """Return every database identifier interpretation of ``text``.

    An empty list means the string matches no known identifier scheme —
    for an LLM-produced "accession", that is the hallucination alarm.
    Ambiguity is real (e.g. ``1ABC`` could be a PDB ID and looks nothing
    like anything else, but many 6-char strings match several schemes), so
    the result is a list, not a single verdict.
    """
    s = text.strip()
    matches: list[IdMatch] = []
    if not s:
        return matches
    for db, (pattern, desc) in _ID_PATTERNS.items():
        m = pattern.fullmatch(s)
        if not m:
            continue
        checksum_ok: Optional[bool] = None
        note = ""
        if db == "cas":
            checksum_ok = _cas_checksum_ok(s)
            if not checksum_ok:
                note = "matches CAS format but the check digit is WRONG"
        elif db == "orcid":
            checksum_ok = _orcid_checksum_ok(s.rsplit("/", 1)[-1])
            if not checksum_ok:
                note = "matches ORCID format but the check digit is WRONG"
        normalized = _normalize_for(db, s)
        matches.append(IdMatch(db, normalized, desc, checksum_ok, note))
    # A long pure-sequence string is not an identifier; say so.
    if not matches:
        cleaned = clean_sequence(s)
        if len(cleaned) >= 10:
            kind, conf, _ = guess_seq_type(cleaned)
            if kind != "unknown" and conf >= 0.9:
                matches.append(IdMatch(
                    "sequence", cleaned[:30] + ("..." if len(cleaned) > 30 else ""),
                    f"not an identifier — looks like a {kind} sequence "
                    f"({len(cleaned)} residues)",
                ))
    return matches


def _normalize_for(db: str, s: str) -> str:
    if db == "uniprot":
        return s.upper()
    if db == "pdb":
        return s.upper()
    if db == "pdb_extended":
        return s.lower()
    if db in ("refseq", "genbank_protein", "ensembl", "chembl", "inchikey",
              "pfam", "interpro"):
        return s.upper()
    if db == "doi":
        s2 = re.sub(r"^(doi:|https?://(dx\.)?doi\.org/)", "", s, flags=re.IGNORECASE)
        # DOI names are case-insensitive; Crossref's convention is lowercase.
        return s2.lower()
    if db == "pmid":
        return re.sub(r"^PMID:? ?", "", s, flags=re.IGNORECASE)
    if db == "ec":
        return re.sub(r"^EC[ :]", "", s)
    if db == "dbsnp":
        return s.lower()
    if db == "orcid":
        return s.rsplit("/", 1)[-1].upper()
    if db == "ncbi_taxon":
        return re.sub(r"^(taxid:|NCBITaxon:)", "", s)
    return s


def is_valid(text: str, database: str) -> bool:
    """True if ``text`` is a well-formed identifier of ``database``
    (including check digit where the scheme has one)."""
    return any(
        m.database == database and m.checksum_ok is not False
        for m in identify(text)
    )


def normalize(text: str, database: str) -> str:
    """Canonical form of ``text`` as a ``database`` identifier.

    Raises ValueError if it is not a well-formed identifier of that scheme.
    """
    for m in identify(text):
        if m.database == database:
            if m.checksum_ok is False:
                raise ValueError(f"{text!r}: {m.note}")
            return m.normalized
    raise ValueError(f"{text!r} is not a well-formed {database} identifier")


# ---------------------------------------------------------------------------
# 2. Sequences
# ---------------------------------------------------------------------------

DNA_UNAMBIGUOUS = set("ACGT")
DNA_ALPHABET = set("ACGTRYSWKMBDHVN")          # IUPAC ambiguity codes
RNA_ALPHABET = set("ACGURYSWKMBDHVN")
PROTEIN_UNAMBIGUOUS = set("ACDEFGHIKLMNPQRSTVWY")
PROTEIN_ALPHABET = PROTEIN_UNAMBIGUOUS | set("BJZXUO")  # ambiguity + Sec/Pyl

_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)

# NCBI translation table 1 (standard code), TCAG ordering.
_T1_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_BASES = "TCAG"
CODON_TABLE_1: dict[str, str] = {
    b1 + b2 + b3: _T1_AAS[16 * i + 4 * j + k]
    for i, b1 in enumerate(_BASES)
    for j, b2 in enumerate(_BASES)
    for k, b3 in enumerate(_BASES)
}

# Average (isotopically weighted) residue masses in Da, per Expasy; the
# residue mass is the amino-acid mass minus one water.
_RESIDUE_MASS_AVG = {
    "G": 57.0519, "A": 71.0788, "S": 87.0782, "P": 97.1167, "V": 99.1326,
    "T": 101.1051, "C": 103.1388, "L": 113.1594, "I": 113.1594,
    "N": 114.1038, "D": 115.0886, "Q": 128.1307, "K": 128.1741,
    "E": 129.1155, "M": 131.1926, "H": 137.1411, "F": 147.1766,
    "R": 156.1875, "Y": 163.1760, "W": 186.2132,
    "U": 150.0388, "O": 237.3018,
}
_WATER_MASS_AVG = 18.01524


def clean_sequence(raw: str) -> str:
    """Strip the noise a copy-pasted sequence carries.

    Removes whitespace, digits (GenBank-style position numbering) and
    trailing ``*`` stop markers; uppercases. Does NOT remove letters, so a
    FASTA header pasted along with the sequence still fails validation
    loudly instead of being silently mangled — strip headers with
    :func:`parse_fasta` first.
    """
    s = re.sub(r"[\s0-9]+", "", raw)
    s = s.rstrip("*")
    return s.upper()


def guess_seq_type(seq: str) -> tuple[str, float, str]:
    """Guess whether ``seq`` is dna / rna / protein.

    Returns ``(kind, confidence, note)`` where kind is one of
    ``"dna" | "rna" | "protein" | "unknown"``. Short ACGT-only strings are
    genuinely ambiguous (every base is also an amino-acid letter); the
    note says so rather than hiding it.
    """
    s = clean_sequence(seq)
    if not s:
        return "unknown", 0.0, "empty sequence"
    letters = set(s)
    nuc_core = sum(s.count(c) for c in "ACGTUN") / len(s)
    if "U" in letters and "T" not in letters and letters <= RNA_ALPHABET:
        return "rna", 0.95, ""
    if letters <= DNA_ALPHABET and nuc_core >= 0.9:
        note = ""
        conf = 0.95
        if len(s) < 20:
            conf = 0.6
            note = ("short ACGT-only string: could equally be a peptide of "
                    "Ala/Cys/Gly/Thr")
        return "dna", conf, note
    if letters <= PROTEIN_ALPHABET:
        return "protein", 0.95, ""
    bad = sorted(letters - PROTEIN_ALPHABET - DNA_ALPHABET - RNA_ALPHABET)
    return "unknown", 0.0, f"characters outside every alphabet: {''.join(bad)}"


def validate_sequence(seq: str, kind: str) -> tuple[bool, str]:
    """Strictly validate ``seq`` against alphabet ``kind``
    (``dna``/``rna``/``protein``, IUPAC ambiguity codes allowed).

    Returns ``(ok, message)``; on failure the message names the offending
    characters and their positions (1-based, in the cleaned sequence).
    """
    alphabet = {"dna": DNA_ALPHABET, "rna": RNA_ALPHABET,
                "protein": PROTEIN_ALPHABET}.get(kind)
    if alphabet is None:
        raise ValueError(f"kind must be dna, rna or protein, not {kind!r}")
    s = clean_sequence(seq)
    if not s:
        return False, "empty sequence after cleaning"
    bad = [(i + 1, c) for i, c in enumerate(s) if c not in alphabet]
    if bad:
        shown = ", ".join(f"{c!r}@{i}" for i, c in bad[:5])
        more = f" (+{len(bad) - 5} more)" if len(bad) > 5 else ""
        return False, f"invalid {kind} characters: {shown}{more}"
    return True, f"valid {kind}, {len(s)} residues"


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA/RNA sequence (IUPAC codes handled).
    U complements to A; output is DNA-style (A pairs to T) unless the
    input contains U, in which case T->A mappings still apply per IUPAC."""
    s = clean_sequence(seq)
    ok, msg = validate_sequence(s, "dna")
    if not ok:
        ok_rna, _ = validate_sequence(s, "rna")
        if not ok_rna:
            raise ValueError(f"not a nucleotide sequence: {msg}")
    return s.translate(_COMPLEMENT)[::-1]


def transcribe(dna: str) -> str:
    """DNA coding strand -> mRNA (T becomes U)."""
    s = clean_sequence(dna)
    ok, msg = validate_sequence(s, "dna")
    if not ok:
        raise ValueError(msg)
    return s.replace("T", "U")


def translate(seq: str, to_stop: bool = False) -> str:
    """Translate a DNA or RNA sequence with NCBI table 1 (standard code).

    Translation starts at position 0 of the given string — pick your
    reading frame before calling. Trailing partial codons raise, because
    silently dropping them is how off-by-one frame bugs survive.
    Ambiguity codes raise unless every expansion agrees (only ``N`` in a
    third position that is fourfold degenerate is resolved).
    """
    s = clean_sequence(seq).replace("U", "T")
    if len(s) % 3 != 0:
        raise ValueError(
            f"length {len(s)} is not a multiple of 3; refusing to guess "
            "the reading frame"
        )
    out = []
    for i in range(0, len(s), 3):
        codon = s[i:i + 3]
        aa = CODON_TABLE_1.get(codon)
        if aa is None:
            if codon[2] == "N" and codon[:2].isalpha():
                fam = {CODON_TABLE_1.get(codon[:2] + b) for b in "TCAG"}
                if len(fam) == 1 and None not in fam:
                    aa = fam.pop()
            if aa is None:
                raise ValueError(f"cannot translate codon {codon!r} at {i}")
        if aa == "*" and to_stop:
            break
        out.append(aa)
    return "".join(out)


def gc_content(seq: str) -> float:
    """GC fraction of a nucleotide sequence. S (G/C) counts toward GC and
    W (A/T) counts against it; all other ambiguity codes are excluded from
    both numerator and denominator."""
    s = clean_sequence(seq)
    unambiguous = [c for c in s if c in "ACGTUSW"]
    if not unambiguous:
        raise ValueError("no unambiguous bases")
    gc = sum(1 for c in unambiguous if c in "GCS")
    return gc / len(unambiguous)


def protein_mw(seq: str) -> float:
    """Average molecular weight of a protein in Da (Expasy residue masses).

    Raises on B/J/Z/X ambiguity codes — an averaged guess presented as a
    mass is exactly the kind of silent error this library exists to stop.
    """
    s = clean_sequence(seq)
    if not s:
        raise ValueError("empty sequence after cleaning")
    total = _WATER_MASS_AVG
    for i, c in enumerate(s):
        m = _RESIDUE_MASS_AVG.get(c)
        if m is None:
            raise ValueError(
                f"residue {c!r} at position {i + 1} has no defined mass "
                "(ambiguity codes B/J/Z/X are refused, not averaged)"
            )
        total += m
    return total


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into ``[(header, sequence), ...]``.

    Strict: content before the first ``>``, records with empty sequences,
    and duplicate IDs (first whitespace-token of the header) all raise
    ValueError. Sequences are cleaned and uppercased.
    """
    entries: list[tuple[str, str]] = []
    header: Optional[str] = None
    chunks: list[str] = []
    seen_ids: set[str] = set()

    def flush() -> None:
        if header is None:
            return
        seq = clean_sequence("".join(chunks))
        if not seq:
            raise ValueError(f"record {header!r} has an empty sequence")
        fid = header.split()[0] if header.split() else ""
        if fid in seen_ids:
            raise ValueError(f"duplicate FASTA id {fid!r}")
        seen_ids.add(fid)
        entries.append((header, seq))

    for lineno, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            flush()
            header = stripped[1:].strip()
            chunks = []
            if not header:
                raise ValueError(f"line {lineno}: empty FASTA header")
        else:
            if header is None:
                raise ValueError(
                    f"line {lineno}: sequence data before any '>' header"
                )
            chunks.append(stripped)
    flush()
    if not entries:
        raise ValueError("no FASTA records found")
    return entries


def write_fasta(records: Iterable[tuple[str, str]], width: int = 60) -> str:
    """Serialize ``[(header, seq), ...]`` to FASTA text (round-trips with
    :func:`parse_fasta`). Headers containing newlines are rejected — they
    would forge extra records that round-trip as genuine."""
    lines: list[str] = []
    for header, seq in records:
        if not header.strip():
            raise ValueError("empty FASTA header")
        if "\n" in header or "\r" in header:
            raise ValueError(f"newline in FASTA header {header!r}")
        s = clean_sequence(seq)
        lines.append(f">{header}")
        lines.extend(s[i:i + width] for i in range(0, len(s), width))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# 2b. Checksums the databases themselves use
# ---------------------------------------------------------------------------

def _crc64_table() -> list[int]:
    poly = 0xD800000000000000  # ISO 3309 polynomial, reversed
    table = []
    for i in range(256):
        crc = i
        for _ in range(8):
            crc = (crc >> 1) ^ poly if crc & 1 else crc >> 1
        table.append(crc)
    return table


_CRC64_TABLE = _crc64_table()


def crc64_iso(seq: str) -> str:
    """CRC64 (ISO 3309) of the cleaned, uppercased sequence — the checksum
    UniProt publishes for every entry, so equality against the ``crc64``
    field of a UniProt record verifies a sequence without downloading it
    twice. Returns 16 uppercase hex characters."""
    crc = 0
    for b in clean_sequence(seq).encode("ascii"):
        crc = _CRC64_TABLE[(crc ^ b) & 0xFF] ^ (crc >> 8)
    return f"{crc:016X}"


def seguid(seq: str) -> str:
    """SEGUID (SEquence Globally Unique IDentifier): base64 of the SHA-1
    of the cleaned, uppercased sequence, padding stripped. 27 characters,
    compatible with Biopython's ``seguid``."""
    digest = hashlib.sha1(clean_sequence(seq).encode("ascii")).digest()
    return base64.b64encode(digest).decode("ascii").rstrip("=")


def seq_sha256(seq: str) -> str:
    """SHA-256 hex digest of the cleaned, uppercased sequence."""
    return hashlib.sha256(clean_sequence(seq).encode("ascii")).hexdigest()


def seq_md5(seq: str) -> str:
    """MD5 hex digest of the cleaned, uppercased sequence (as used by
    Ensembl and some archives). Not collision-resistant; fine as a
    cross-reference key, not as an integrity guarantee."""
    return hashlib.md5(clean_sequence(seq).encode("ascii")).hexdigest()


# ---------------------------------------------------------------------------
# 3. Provenance-first fetching
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Evidence:
    """Where a piece of data came from, verifiably.

    Attach this (via ``as_dict()``) to any number or sequence you publish;
    the sha256 is of the exact bytes received, so anyone can re-fetch and
    diff."""

    url: str
    status: int
    sha256: str
    size: int
    retrieved_at: str          # ISO-8601 UTC
    from_cache: bool
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class Record:
    """Response bytes plus their Evidence."""

    content: bytes
    evidence: Evidence

    @property
    def text(self) -> str:
        return self.content.decode("utf-8", errors="replace")

    def json(self) -> Any:
        return json.loads(self.content)


class FetchError(RuntimeError):
    """Network/HTTP failure after retries, or offline-mode cache miss."""


# A transport takes (url, headers, timeout) and returns
# (status_code, body_bytes). Injectable so tests never touch the network.
Transport = Callable[[str, dict[str, str], float], tuple[int, bytes]]

_MAX_RESPONSE_BYTES = 64 * 1024 * 1024


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Refuse redirects. urlopen otherwise follows them silently — even
    cross-host and https->http — which would let an allowlisted host hand
    us bytes from anywhere while Evidence.url still named the original
    URL. A 3xx therefore surfaces as its own status and the fetch fails
    loudly instead of lying about provenance."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: D102
        return None


_OPENER = urllib.request.build_opener(_NoRedirect)


def _default_transport(url: str, headers: dict[str, str],
                       timeout: float) -> tuple[int, bytes]:
    req = urllib.request.Request(url, headers=headers)
    try:
        with _OPENER.open(req, timeout=timeout) as resp:
            body = resp.read(_MAX_RESPONSE_BYTES + 1)
            if len(body) > _MAX_RESPONSE_BYTES:
                raise FetchError(
                    f"response body exceeds {_MAX_RESPONSE_BYTES} bytes: {url}"
                )
            return resp.status, body
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read()


_REDACT_PARAMS = {"api_key", "email"}


def _redact_url(url: str) -> str:
    """Mask credential-bearing query parameters (NCBI api_key/email) so
    they never land in Evidence records or cache metadata."""
    parts = urllib.parse.urlsplit(url)
    pairs = urllib.parse.parse_qsl(parts.query, keep_blank_values=True)
    if not any(k in _REDACT_PARAMS for k, _ in pairs):
        return url
    query = urllib.parse.urlencode(
        [(k, "REDACTED" if k in _REDACT_PARAMS else v) for k, v in pairs]
    )
    return urllib.parse.urlunsplit(
        (parts.scheme, parts.netloc, parts.path, query, parts.fragment)
    )


def _utc_now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


class Fetcher:
    """Polite, cached, provenance-recording HTTP GET.

    - On-disk cache keyed by SHA-256 of URL + Accept header (no
      path-traversal risk), written atomically, integrity-checked on every
      read (a corrupt entry is refetched, never served as authentic).
    - Per-host minimum interval between requests (NCBI gets 0.34 s to
      respect their 3 req/s unkeyed limit; everyone else 0.25 s).
    - Retries with fixed exponential backoff (1/2/4 s) on 429/5xx.
      Retry-After is NOT consulted — the transport contract deliberately
      carries no response headers.
    - HTTPS on the default port only, hosts on the allowlist only (unless
      ``allow_any_host=True``), and redirects refused — a 3xx fails
      loudly rather than fetching bytes Evidence would misattribute.
    - NCBI credentials (api_key/email) are redacted from Evidence.url and
      cache metadata; the cache directory is created mode 0700.
    - ``offline=True`` (or ALLONET_OFFLINE=1) serves cache only and raises
      FetchError on a miss — reproducible runs cannot silently refetch.
    """

    ALLOWED_HOSTS = {
        "rest.uniprot.org", "www.uniprot.org",
        "data.rcsb.org", "files.rcsb.org", "www.rcsb.org",
        "eutils.ncbi.nlm.nih.gov",
        "rest.ensembl.org",
        "doi.org", "www.ebi.ac.uk",
    }

    def __init__(self,
                 cache_dir: Optional[str] = None,
                 offline: Optional[bool] = None,
                 transport: Optional[Transport] = None,
                 timeout: float = 30.0,
                 max_retries: int = 3,
                 allow_any_host: bool = False,
                 user_agent: Optional[str] = None) -> None:
        if cache_dir is None:
            cache_dir = os.environ.get(
                "ALLONET_CACHE",
                str(Path.home() / ".allonet" / "cache"),
            )
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(self.cache_dir, 0o700)
        except OSError:
            pass
        if offline is None:
            offline = os.environ.get("ALLONET_OFFLINE", "") == "1"
        self.offline = offline
        self.transport: Transport = transport or _default_transport
        self.timeout = timeout
        if max_retries < 0:
            raise ValueError("max_retries must be >= 0")
        self.max_retries = max_retries
        self.allow_any_host = allow_any_host
        self.user_agent = user_agent or f"allonet/{__version__}"
        self._last_request: dict[str, float] = {}
        self._sleep = time.sleep  # injectable for tests

    # -- internals ---------------------------------------------------------

    def _cache_paths(self, url: str,
                     accept: Optional[str] = None) -> tuple[Path, Path]:
        key = hashlib.sha256(f"{url}\n{accept or ''}".encode("utf-8")).hexdigest()
        return (self.cache_dir / f"{key}.body",
                self.cache_dir / f"{key}.meta.json")

    def _min_interval(self, host: str) -> float:
        if host == "eutils.ncbi.nlm.nih.gov" or host.endswith(".ncbi.nlm.nih.gov"):
            return 0.34
        return 0.25

    def _throttle(self, host: str) -> None:
        last = self._last_request.get(host)
        if last is not None:
            wait = self._min_interval(host) - (time.monotonic() - last)
            if wait > 0:
                self._sleep(wait)
        self._last_request[host] = time.monotonic()

    def _check_url(self, url: str) -> str:
        parsed = urllib.parse.urlsplit(url)
        if parsed.scheme != "https":
            raise ValueError(f"only https URLs are fetched, got {url!r}")
        if parsed.port not in (None, 443):
            raise ValueError(f"only port 443 is fetched, got {url!r}")
        host = parsed.hostname or ""
        if not self.allow_any_host and host not in self.ALLOWED_HOSTS:
            raise ValueError(
                f"host {host!r} is not on the allowlist; pass "
                "allow_any_host=True to a Fetcher you construct yourself"
            )
        return host

    # -- public ------------------------------------------------------------

    def _read_cache(self, body_path: Path,
                    meta_path: Path) -> Optional[Record]:
        """Return the cached Record, or None if the entry is absent,
        corrupt, or fails its integrity check."""
        try:
            meta = json.loads(meta_path.read_text())
            content = body_path.read_bytes()
            digest = hashlib.sha256(content).hexdigest()
            if digest != meta["sha256"] or len(content) != int(meta["size"]):
                return None
            return Record(content, Evidence(
                url=meta["url"], status=int(meta["status"]),
                sha256=digest, size=len(content),
                retrieved_at=meta["retrieved_at"], from_cache=True,
            ))
        except (OSError, ValueError, KeyError, TypeError):
            return None

    @staticmethod
    def _atomic_write(path: Path, data: bytes) -> None:
        tmp = path.parent / (path.name + ".tmp")
        tmp.write_bytes(data)
        os.replace(tmp, path)

    def get(self, url: str, accept: Optional[str] = None,
            force_refresh: bool = False) -> Record:
        host = self._check_url(url)
        body_path, meta_path = self._cache_paths(url, accept)
        safe_url = _redact_url(url)

        if not force_refresh and body_path.exists() and meta_path.exists():
            cached = self._read_cache(body_path, meta_path)
            if cached is not None:
                return cached
            if self.offline:
                raise FetchError(
                    f"offline mode and the cache entry for {safe_url} is "
                    "corrupt or unreadable"
                )

        if self.offline:
            raise FetchError(f"offline mode and no cache entry for {safe_url}")

        headers = {"User-Agent": self.user_agent}
        if accept:
            headers["Accept"] = accept

        delay = 1.0
        last_status = -1
        last_body = b""
        for attempt in range(self.max_retries + 1):
            self._throttle(host)
            try:
                status, body = self.transport(url, headers, self.timeout)
            except (OSError, http.client.HTTPException) as exc:
                if attempt == self.max_retries:
                    raise FetchError(f"GET {safe_url} failed: {exc}") from exc
                self._sleep(delay)
                delay *= 2
                continue
            last_status, last_body = status, body
            if status in (429, 500, 502, 503, 504) and attempt < self.max_retries:
                self._sleep(delay)
                delay *= 2
                continue
            break

        if last_status != 200:
            snippet = last_body[:200].decode("utf-8", errors="replace")
            raise FetchError(
                f"GET {safe_url} returned HTTP {last_status}: {snippet}"
            )

        retrieved_at = _utc_now_iso()
        digest = hashlib.sha256(last_body).hexdigest()
        self._atomic_write(body_path, last_body)
        self._atomic_write(meta_path, json.dumps({
            "url": safe_url, "status": last_status,
            "retrieved_at": retrieved_at, "sha256": digest,
            "size": len(last_body),
        }).encode("utf-8"))
        return Record(last_body, Evidence(
            url=safe_url, status=last_status, sha256=digest,
            size=len(last_body), retrieved_at=retrieved_at,
            from_cache=False,
        ))


_DEFAULT_FETCHER: Optional[Fetcher] = None


def get_fetcher() -> Fetcher:
    """The process-wide default Fetcher (created lazily)."""
    global _DEFAULT_FETCHER
    if _DEFAULT_FETCHER is None:
        _DEFAULT_FETCHER = Fetcher()
    return _DEFAULT_FETCHER


# -- database clients -------------------------------------------------------

def uniprot_fasta(accession: str, fetcher: Optional[Fetcher] = None) -> Record:
    acc = normalize(accession, "uniprot")
    f = fetcher or get_fetcher()
    return f.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")


def uniprot_json(accession: str, fetcher: Optional[Fetcher] = None) -> Record:
    acc = normalize(accession, "uniprot")
    f = fetcher or get_fetcher()
    return f.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                 accept="application/json")


def uniprot_sequence(accession: str,
                     fetcher: Optional[Fetcher] = None) -> tuple[str, Evidence]:
    """The canonical sequence of a UniProt entry, plus its Evidence."""
    rec = uniprot_fasta(accession, fetcher)
    entries = parse_fasta(rec.text)
    return entries[0][1], rec.evidence


def rcsb_entry(pdb_id: str, fetcher: Optional[Fetcher] = None) -> Record:
    pid = normalize(pdb_id, "pdb")
    f = fetcher or get_fetcher()
    return f.get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}",
                 accept="application/json")


def rcsb_polymer_entities(pdb_id: str,
                          fetcher: Optional[Fetcher] = None) -> list[Record]:
    """All polymer-entity records of a PDB entry (they carry the UniProt
    cross-references)."""
    pid = normalize(pdb_id, "pdb")
    f = fetcher or get_fetcher()
    entry = rcsb_entry(pid, f).json()
    entity_ids = (entry.get("rcsb_entry_container_identifiers", {})
                  .get("polymer_entity_ids", []))
    return [
        f.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{eid}",
              accept="application/json")
        for eid in entity_ids
    ]


def pdb_structure(pdb_id: str, fmt: str = "cif",
                  fetcher: Optional[Fetcher] = None) -> Record:
    """Download a structure file ('cif' or 'pdb') from files.rcsb.org."""
    pid = normalize(pdb_id, "pdb")
    if fmt not in ("cif", "pdb"):
        raise ValueError("fmt must be 'cif' or 'pdb'")
    f = fetcher or get_fetcher()
    return f.get(f"https://files.rcsb.org/download/{pid}.{fmt}")


def _ncbi_params(extra: dict[str, str]) -> str:
    params = dict(extra)
    params.setdefault("tool", f"allonet-{__version__}")
    email = os.environ.get("ALLONET_NCBI_EMAIL", "")
    if email:
        params["email"] = email
    key = os.environ.get("ALLONET_NCBI_KEY", "")
    if key:
        params["api_key"] = key
    return urllib.parse.urlencode(params)


def ncbi_efetch(db: str, id_: str, rettype: str = "fasta",
                retmode: str = "text",
                fetcher: Optional[Fetcher] = None) -> Record:
    """NCBI E-utilities efetch. Set ALLONET_NCBI_EMAIL / ALLONET_NCBI_KEY
    to identify yourself to NCBI (recommended by their usage policy)."""
    qs = _ncbi_params({"db": db, "id": id_, "rettype": rettype,
                       "retmode": retmode})
    f = fetcher or get_fetcher()
    return f.get(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{qs}"
    )


def ncbi_esearch(db: str, term: str,
                 fetcher: Optional[Fetcher] = None) -> Record:
    qs = _ncbi_params({"db": db, "term": term, "retmode": "json"})
    f = fetcher or get_fetcher()
    return f.get(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{qs}"
    )


def ensembl_get(path: str, fetcher: Optional[Fetcher] = None) -> Record:
    """GET an Ensembl REST path, e.g. ``lookup/id/ENSG00000139618``.

    The path must be a plain endpoint path: each segment is
    percent-encoded, and query strings, fragments, empty and dot segments
    are rejected so an untrusted ID cannot rewrite the request.
    """
    f = fetcher or get_fetcher()
    clean = path.strip("/")
    if not clean or "?" in clean or "#" in clean:
        raise ValueError(
            f"ensembl path must not be empty or contain '?'/'#': {path!r}"
        )
    segments = clean.split("/")
    if any(seg in ("", ".", "..") for seg in segments):
        raise ValueError(f"bad ensembl path segment in {path!r}")
    quoted = "/".join(urllib.parse.quote(seg, safe="") for seg in segments)
    return f.get(
        f"https://rest.ensembl.org/{quoted}?content-type=application/json",
        accept="application/json",
    )


# ---------------------------------------------------------------------------
# 4. Verification — claims checked against primary sources
# ---------------------------------------------------------------------------

CONFIRMED = "CONFIRMED"
REFUTED = "REFUTED"
UNVERIFIABLE = "UNVERIFIABLE"


@dataclass
class Verification:
    """The outcome of checking one claim against a primary source."""

    claim: str
    verdict: str               # CONFIRMED / REFUTED / UNVERIFIABLE
    detail: str
    evidence: list[Evidence] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.verdict == CONFIRMED

    def as_dict(self) -> dict[str, Any]:
        return {
            "claim": self.claim,
            "verdict": self.verdict,
            "detail": self.detail,
            "evidence": [e.as_dict() for e in self.evidence],
        }


def verify_sequence_checksum(seq: str, claimed: str) -> Verification:
    """Check a sequence against a claimed checksum, auto-detecting the
    scheme from the checksum's format (CRC64 = 16 hex, MD5 = 32 hex,
    SHA-256 = 64 hex, SEGUID = 27 base64 chars)."""
    c = claimed.strip()
    claim = f"sequence has checksum {c}"
    if re.fullmatch(r"[0-9A-Fa-f]{16}", c):
        actual, scheme = crc64_iso(seq), "CRC64-ISO"
        match = actual.upper() == c.upper()
    elif re.fullmatch(r"[0-9A-Fa-f]{32}", c):
        actual, scheme = seq_md5(seq), "MD5"
        match = actual.lower() == c.lower()
    elif re.fullmatch(r"[0-9A-Fa-f]{64}", c):
        actual, scheme = seq_sha256(seq), "SHA-256"
        match = actual.lower() == c.lower()
    elif re.fullmatch(r"[A-Za-z0-9+/]{27}", c):
        actual, scheme = seguid(seq), "SEGUID"
        match = actual == c
    else:
        return Verification(claim, UNVERIFIABLE,
                            f"unrecognized checksum format: {c!r}")
    verdict = CONFIRMED if match else REFUTED
    return Verification(claim, verdict,
                        f"{scheme} of the given sequence is {actual}")


def verify_uniprot_sequence(accession: str, seq: str,
                            fetcher: Optional[Fetcher] = None) -> Verification:
    """Does ``seq`` match the current canonical sequence of ``accession``?"""
    acc = normalize(accession, "uniprot")
    claim = f"the given sequence is the canonical sequence of UniProt {acc}"
    try:
        db_seq, ev = uniprot_sequence(acc, fetcher)
    except (FetchError, ValueError) as exc:
        return Verification(claim, UNVERIFIABLE, str(exc))
    given = clean_sequence(seq)
    if given == db_seq:
        return Verification(claim, CONFIRMED,
                            f"exact match, {len(db_seq)} residues, "
                            f"CRC64 {crc64_iso(db_seq)}", [ev])
    detail = (f"mismatch: given {len(given)} residues "
              f"(CRC64 {crc64_iso(given)}), database has {len(db_seq)} "
              f"(CRC64 {crc64_iso(db_seq)})")
    if len(given) == len(db_seq):
        diffs = [i + 1 for i, (a, b) in enumerate(zip(given, db_seq)) if a != b]
        detail += (f"; same length, {len(diffs)} differing positions, "
                   f"first at {diffs[0]}")
    return Verification(claim, REFUTED, detail, [ev])


def verify_pdb_uniprot_link(pdb_id: str, accession: str,
                            fetcher: Optional[Fetcher] = None) -> Verification:
    """Is PDB entry ``pdb_id`` actually a structure of UniProt
    ``accession``, per RCSB's own cross-references?"""
    pid = normalize(pdb_id, "pdb")
    acc = normalize(accession, "uniprot").split("-")[0]
    claim = f"PDB {pid} is a structure of UniProt {acc}"
    try:
        entities = rcsb_polymer_entities(pid, fetcher)
    except (FetchError, ValueError) as exc:
        return Verification(claim, UNVERIFIABLE, str(exc))
    found: set[str] = set()
    evidence = [e.evidence for e in entities]
    for rec in entities:
        data = rec.json()
        refs = (data.get("rcsb_polymer_entity_container_identifiers", {})
                .get("reference_sequence_identifiers") or [])
        for ref in refs:
            if ref.get("database_name", "").lower() == "uniprot":
                found.add(str(ref.get("database_accession", "")).upper())
    if acc in found:
        return Verification(claim, CONFIRMED,
                            f"RCSB lists UniProt refs {sorted(found)}",
                            evidence)
    if not found:
        return Verification(claim, UNVERIFIABLE,
                            "RCSB lists no UniProt reference for any "
                            "polymer entity of this entry", evidence)
    return Verification(claim, REFUTED,
                        f"RCSB lists {sorted(found)}, not {acc}", evidence)


def verify_doi(doi: str, fetcher: Optional[Fetcher] = None) -> Verification:
    """Does this DOI resolve, per the doi.org handle registry?"""
    d = normalize(doi, "doi")
    claim = f"DOI {d} exists"
    f = fetcher or get_fetcher()
    url = "https://doi.org/api/handles/" + urllib.parse.quote(d, safe="/")
    try:
        rec = f.get(url, accept="application/json")
    except FetchError as exc:
        # The handle API answers 404 with a JSON body for unknown handles.
        msg = str(exc)
        if "HTTP 404" in msg:
            return Verification(claim, REFUTED,
                                "doi.org handle registry: not found")
        return Verification(claim, UNVERIFIABLE, msg)
    data = rec.json()
    code = data.get("responseCode")
    if code == 1:
        return Verification(claim, CONFIRMED,
                            "doi.org handle registry resolves it",
                            [rec.evidence])
    return Verification(claim, REFUTED,
                        f"doi.org responseCode={code}", [rec.evidence])


def verify_gene_symbol(symbol: str, organism: str = "Homo sapiens",
                       fetcher: Optional[Fetcher] = None) -> Verification:
    """Is ``symbol`` a current gene symbol in ``organism``, per NCBI Gene?

    Matches current symbols only (``[sym]`` field); an alias or a retired
    symbol REFUTES the claim as stated, which is what you want when the
    symbol is about to go into a figure or an order form.
    """
    sym = symbol.strip()
    claim = f"{sym!r} is a current gene symbol in {organism}"
    if not re.fullmatch(r"[A-Za-z0-9_.@-]{1,64}", sym):
        return Verification(claim, REFUTED,
                            "not even a well-formed gene symbol")
    term = f"{sym}[sym] AND {organism}[orgn]"
    try:
        rec = ncbi_esearch("gene", term, fetcher)
        count = int(rec.json()["esearchresult"]["count"])
    except (FetchError, KeyError, ValueError) as exc:
        return Verification(claim, UNVERIFIABLE, str(exc))
    if count > 0:
        return Verification(claim, CONFIRMED,
                            f"NCBI Gene has {count} matching record(s)",
                            [rec.evidence])
    return Verification(claim, REFUTED,
                        "NCBI Gene has no record with this current symbol "
                        "in this organism (it may be an alias — check "
                        "esearch without [sym])", [rec.evidence])


# ---------------------------------------------------------------------------
# 5. Units — the conversions people get wrong under deadline
# ---------------------------------------------------------------------------

# CODATA R = 8.31446261815324 J/(mol K), thermochemical calorie = 4.184 J.
R_KCAL_PER_MOL_K = 8.31446261815324 / 4184   # gas constant, kcal / (mol K)
_KJ_PER_KCAL = 4.184                         # thermochemical calorie, exact


def kcal_to_kj(kcal: float) -> float:
    """kcal/mol -> kJ/mol (thermochemical calorie, factor exactly 4.184)."""
    return kcal * _KJ_PER_KCAL


def kj_to_kcal(kj: float) -> float:
    """kJ/mol -> kcal/mol (thermochemical calorie)."""
    return kj / _KJ_PER_KCAL


def kd_to_delta_g(kd_molar: float, temp_k: float = 298.15) -> float:
    """Dissociation constant (molar) -> binding free energy, kcal/mol.

    ΔG = RT ln(Kd) with the 1 M standard state; a 1 nM binder at 298.15 K
    gives about -12.3 kcal/mol. Negative = favorable binding.
    """
    if kd_molar <= 0:
        raise ValueError("Kd must be positive (molar units)")
    if temp_k <= 0:
        raise ValueError("temperature must be positive kelvin")
    return R_KCAL_PER_MOL_K * temp_k * math.log(kd_molar)


def delta_g_to_kd(dg_kcal: float, temp_k: float = 298.15) -> float:
    """Binding free energy (kcal/mol) -> dissociation constant (molar)."""
    if temp_k <= 0:
        raise ValueError("temperature must be positive kelvin")
    return math.exp(dg_kcal / (R_KCAL_PER_MOL_K * temp_k))


# ---------------------------------------------------------------------------
# 6. Command line
# ---------------------------------------------------------------------------

def _print(obj: Any, as_json: bool) -> None:
    if as_json:
        print(json.dumps(obj, indent=2, sort_keys=True))
    elif isinstance(obj, dict):
        for k, v in obj.items():
            print(f"{k}: {v}")
    else:
        print(obj)


def _seq_from_file(p: Path) -> str:
    text = p.read_text()
    if text.lstrip().startswith(">"):
        records = parse_fasta(text)
        if len(records) > 1:
            print(f"note: {p} has {len(records)} FASTA records; using the "
                  "first", file=sys.stderr)
        return records[0][1]
    return clean_sequence(text)


def _read_seq_arg(value: str) -> str:
    """A sequence argument: ``@path`` forces file input; otherwise a bare
    value naming an existing file is read as one (with a stderr notice, so
    a literal shadowed by a file name is never reinterpreted silently)."""
    if value.startswith("@"):
        p = Path(value[1:])
        if not p.is_file():
            raise ValueError(f"no such file: {p}")
        return _seq_from_file(p)
    p = Path(value)
    if p.is_file():
        print(f"note: reading {value!r} as a file (prefix with @ to make "
              "this explicit)", file=sys.stderr)
        return _seq_from_file(p)
    return clean_sequence(value)


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="allonet",
        description="AlloNet: validate identifiers, checksum sequences, "
                    "fetch with provenance, verify claims against primary "
                    "biological databases.",
    )
    parser.add_argument("--version", action="version",
                        version=f"allonet {__version__}")
    parser.add_argument("--json", action="store_true",
                        help="emit JSON instead of text")
    sub = parser.add_subparsers(dest="command", required=True)

    p_id = sub.add_parser("id", help="identify / validate an identifier")
    p_id.add_argument("value")

    p_seq = sub.add_parser("seq", help="sequence utilities")
    p_seq.add_argument("op", choices=["stats", "validate", "revcomp",
                                      "translate", "clean"])
    p_seq.add_argument("value", help="sequence literal, or @path / existing file path")
    p_seq.add_argument("--kind", choices=["dna", "rna", "protein"],
                       help="required for validate")

    p_ck = sub.add_parser("checksum", help="all checksums of a sequence")
    p_ck.add_argument("value", help="sequence literal, or @path / existing file path")

    p_fetch = sub.add_parser("fetch", help="fetch a record with provenance")
    p_fetch.add_argument("source", choices=["uniprot", "pdb", "ncbi",
                                            "ensembl"])
    p_fetch.add_argument("id")
    p_fetch.add_argument("--db", default="nuccore",
                         help="NCBI database for source=ncbi")

    p_ver = sub.add_parser("verify", help="verify a claim against a "
                                          "primary source")
    ver_sub = p_ver.add_subparsers(dest="what", required=True)
    v1 = ver_sub.add_parser("uniprot-seq")
    v1.add_argument("accession")
    v1.add_argument("value", help="sequence literal, or @path / existing file path")
    v2 = ver_sub.add_parser("pdb-uniprot")
    v2.add_argument("pdb_id")
    v2.add_argument("accession")
    v3 = ver_sub.add_parser("checksum")
    v3.add_argument("value", help="sequence literal, or @path / existing file path")
    v3.add_argument("claimed")
    v4 = ver_sub.add_parser("doi")
    v4.add_argument("doi")
    v5 = ver_sub.add_parser("gene")
    v5.add_argument("symbol")
    v5.add_argument("--organism", default="Homo sapiens")

    p_units = sub.add_parser("units", help="unit conversions")
    p_units.add_argument("op", choices=["kd2dg", "dg2kd", "kcal2kj",
                                        "kj2kcal"])
    p_units.add_argument("value", type=float)
    p_units.add_argument("--temp", type=float, default=298.15)

    args = parser.parse_args(argv)

    try:
        if args.command == "id":
            matches = identify(args.value)
            if not matches:
                _print({"input": args.value, "matches": [],
                        "verdict": "no known identifier scheme matches — "
                                   "treat as unverified or hallucinated"},
                       args.json)
                return 1
            _print({"input": args.value,
                    "matches": [m.as_dict() for m in matches]}, args.json)
            return 0

        if args.command == "seq":
            seq = _read_seq_arg(args.value)
            if args.op == "clean":
                _print(seq, args.json)
            elif args.op == "stats":
                kind, conf, note = guess_seq_type(seq)
                stats: dict[str, Any] = {
                    "length": len(seq), "guessed_type": kind,
                    "confidence": conf, "note": note,
                    "crc64": crc64_iso(seq), "seguid": seguid(seq),
                }
                if kind in ("dna", "rna"):
                    stats["gc_content"] = round(gc_content(seq), 4)
                if kind == "protein":
                    try:
                        stats["mw_avg_da"] = round(protein_mw(seq), 2)
                    except ValueError as exc:
                        stats["mw_avg_da"] = f"undefined ({exc})"
                _print(stats, args.json)
            elif args.op == "validate":
                if not args.kind:
                    parser.error("seq validate requires --kind")
                ok, msg = validate_sequence(seq, args.kind)
                _print({"ok": ok, "detail": msg}, args.json)
                return 0 if ok else 1
            elif args.op == "revcomp":
                _print(reverse_complement(seq), args.json)
            elif args.op == "translate":
                _print(translate(seq), args.json)
            return 0

        if args.command == "checksum":
            seq = _read_seq_arg(args.value)
            _print({"length": len(seq), "crc64_iso": crc64_iso(seq),
                    "seguid": seguid(seq), "sha256": seq_sha256(seq),
                    "md5": seq_md5(seq)}, args.json)
            return 0

        if args.command == "fetch":
            if args.source == "uniprot":
                rec = uniprot_fasta(args.id)
            elif args.source == "pdb":
                rec = rcsb_entry(args.id)
            elif args.source == "ncbi":
                rec = ncbi_efetch(args.db, args.id)
            else:
                rec = ensembl_get(f"lookup/id/{normalize(args.id, 'ensembl')}")
            out = {"evidence": rec.evidence.as_dict(),
                   "content": rec.text}
            if args.json:
                _print(out, True)
            else:
                print(rec.text)
                ev = rec.evidence
                print(f"--- evidence: {ev.url}\n--- sha256 {ev.sha256}"
                      f"  retrieved {ev.retrieved_at}"
                      f"  cache={ev.from_cache}", file=sys.stderr)
            return 0

        if args.command == "verify":
            if args.what == "uniprot-seq":
                v = verify_uniprot_sequence(args.accession,
                                            _read_seq_arg(args.value))
            elif args.what == "pdb-uniprot":
                v = verify_pdb_uniprot_link(args.pdb_id, args.accession)
            elif args.what == "checksum":
                v = verify_sequence_checksum(_read_seq_arg(args.value),
                                             args.claimed)
            elif args.what == "doi":
                v = verify_doi(args.doi)
            else:
                v = verify_gene_symbol(args.symbol, args.organism)
            _print(v.as_dict(), args.json)
            return 0 if v.ok else 1

        if args.command == "units":
            if args.op == "kd2dg":
                out = kd_to_delta_g(args.value, args.temp)
            elif args.op == "dg2kd":
                out = delta_g_to_kd(args.value, args.temp)
            elif args.op == "kcal2kj":
                out = kcal_to_kj(args.value)
            else:
                out = kj_to_kcal(args.value)
            _print(out, args.json)
            return 0

    except (ValueError, FetchError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    parser.error("unhandled command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
