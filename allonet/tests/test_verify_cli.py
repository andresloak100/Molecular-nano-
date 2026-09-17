"""Offline tests for AlloNet's verification layer and CLI.

Everything here is deterministic and network-free: every Fetcher is built
with an injected fake transport and a per-test cache directory, and only
CLI paths that never touch the network are exercised.

Checksum ground truths were computed independently of allonet's own
implementations: CRC64 values come from a separate bit-wise ISO-3309
engine validated against the published CRC-64/GO-ISO check value
(0xB90956C775A41001 for b"123456789"); MD5/SHA-256/SEGUID come straight
from hashlib per each scheme's definition. The HBA_HUMAN CRC64 below
(15E13666573BBBAE) also matches the checksum UniProt publishes for
P69905, tying the implementation to the database convention it claims.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
import sys
import urllib.parse
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import allonet
from allonet import (
    CONFIRMED,
    REFUTED,
    UNVERIFIABLE,
    Evidence,
    Fetcher,
    Verification,
    verify_doi,
    verify_gene_symbol,
    verify_pdb_uniprot_link,
    verify_sequence_checksum,
    verify_uniprot_sequence,
)

# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

# Human hemoglobin alpha, UniProt P69905, 142 residues.
HBA_SEQ = (
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
    "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP"
    "AVHASLDKFLASVSTVLTSKYR"
)
HBA_FASTA = (
    ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens "
    "OX=9606 GN=HBA1 PE=1 SV=2\n"
    + "\n".join(HBA_SEQ[i:i + 60] for i in range(0, len(HBA_SEQ), 60))
    + "\n"
)
HBA_CRC64 = "15E13666573BBBAE"          # UniProt's published checksum

ELVIS = "ELVISLIVES"                     # ten valid amino acids
ELVIS_CRC64 = "C59213C05735B042"         # independent bit-wise engine
ACGT_CRC64 = "71A87EBDB0000000"          # independent bit-wise engine

KD2DG_1NM = -12.278223117308158          # R * 298.15 K * ln(1e-9)

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/P69905.fasta"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/"
RCSB_ENTITY = "https://data.rcsb.org/rest/v1/core/polymer_entity/"
DOI_HANDLES = "https://doi.org/api/handles/"
ESEARCH_PREFIX = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"


def _entry_json(pdb_id: str, entity_ids: list[str]) -> str:
    return json.dumps({
        "rcsb_id": pdb_id,
        "rcsb_entry_container_identifiers": {
            "entry_id": pdb_id,
            "polymer_entity_ids": entity_ids,
        },
        "struct": {"title": f"structure {pdb_id}"},
    })


def _entity_json(pdb_id: str, entity_id: str, refs) -> str:
    return json.dumps({
        "rcsb_id": f"{pdb_id}_{entity_id}",
        "rcsb_polymer_entity_container_identifiers": {
            "entry_id": pdb_id,
            "entity_id": entity_id,
            "reference_sequence_identifiers": refs,
        },
    })


def _esearch_json(count: int, ids: list[str]) -> str:
    return json.dumps({
        "header": {"type": "esearch", "version": "0.3"},
        "esearchresult": {
            "count": str(count),
            "retmax": str(len(ids)),
            "retstart": "0",
            "idlist": ids,
        },
    })


# ---------------------------------------------------------------------------
# Fake transports and fetcher factory
# ---------------------------------------------------------------------------

class CannedTransport:
    """Transport serving canned (status, body) pairs by URL prefix.

    An unrouted URL raises AssertionError so a test cannot silently hit
    an unexpected endpoint (or the real network).
    """

    def __init__(self, routes: list[tuple[str, int, str]]) -> None:
        self.routes = routes
        self.calls: list[str] = []

    def __call__(self, url: str, headers: dict[str, str],
                 timeout: float) -> tuple[int, bytes]:
        self.calls.append(url)
        for prefix, status, body in self.routes:
            if url.startswith(prefix):
                return status, body.encode("utf-8")
        raise AssertionError(f"test transport got unexpected URL: {url}")


class FailingTransport:
    """Transport that always raises OSError (network down)."""

    def __init__(self) -> None:
        self.calls = 0

    def __call__(self, url, headers, timeout):
        self.calls += 1
        raise OSError("connection refused (fake)")


def make_fetcher(tmp_path, transport, max_retries: int = 0) -> Fetcher:
    f = Fetcher(cache_dir=str(tmp_path / "cache"), offline=False,
                transport=transport, max_retries=max_retries)
    f._sleep = lambda s: None   # documented injection point; no real waiting
    return f


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Keep fetch URLs and offline behavior deterministic."""
    for var in ("ALLONET_NCBI_EMAIL", "ALLONET_NCBI_KEY", "ALLONET_OFFLINE"):
        monkeypatch.delenv(var, raising=False)


# ---------------------------------------------------------------------------
# verify_sequence_checksum — format auto-detection
# ---------------------------------------------------------------------------

def test_checksum_crc64_confirmed():
    v = verify_sequence_checksum(ELVIS, ELVIS_CRC64)
    assert v.verdict == CONFIRMED
    assert v.ok is True
    assert "CRC64-ISO" in v.detail
    assert ELVIS_CRC64 in v.detail
    assert v.evidence == []          # local computation, no fetch


def test_checksum_crc64_case_and_whitespace_insensitive():
    # lowercase hex with surrounding whitespace still confirms
    v = verify_sequence_checksum(ELVIS, f"  {ELVIS_CRC64.lower()} ")
    assert v.verdict == CONFIRMED


def test_checksum_crc64_refuted():
    wrong = "0000000000000000"
    v = verify_sequence_checksum(ELVIS, wrong)
    assert v.verdict == REFUTED
    assert v.ok is False
    assert ELVIS_CRC64 in v.detail   # tells you what it actually is


def test_checksum_md5_confirmed():
    md5 = hashlib.md5(ELVIS.encode()).hexdigest()
    v = verify_sequence_checksum(ELVIS, md5)
    assert v.verdict == CONFIRMED
    assert "MD5" in v.detail


def test_checksum_md5_refuted():
    wrong = "d" * 32
    v = verify_sequence_checksum(ELVIS, wrong)
    assert v.verdict == REFUTED
    assert "MD5" in v.detail


def test_checksum_sha256_confirmed():
    sha = hashlib.sha256(ELVIS.encode()).hexdigest()
    v = verify_sequence_checksum(ELVIS, sha.upper())   # case-insensitive
    assert v.verdict == CONFIRMED
    assert "SHA-256" in v.detail


def test_checksum_sha256_refuted():
    v = verify_sequence_checksum(ELVIS, "e" * 64)
    assert v.verdict == REFUTED
    assert "SHA-256" in v.detail


def test_checksum_seguid_confirmed():
    seg = (hashlib.sha1(ELVIS.encode()).digest())
    seg_b64 = __import__("base64").b64encode(seg).decode().rstrip("=")
    assert len(seg_b64) == 27
    v = verify_sequence_checksum(ELVIS, seg_b64)
    assert v.verdict == CONFIRMED
    assert "SEGUID" in v.detail


def test_checksum_seguid_refuted():
    # well-formed 27-char base64, but not this sequence's SEGUID
    v = verify_sequence_checksum(ELVIS, "A" * 27)
    assert v.verdict == REFUTED
    assert "SEGUID" in v.detail


@pytest.mark.parametrize("claimed", [
    "",                                              # empty
    "xyz",                                           # junk
    "0123456789abcdef0123456789abcdef01234567",      # 40 hex = SHA-1, unsupported
    "A" * 26 + "-",                                  # 27 chars, not base64
    "CRC-C59213C05735B042",                          # Biopython-prefixed CRC64
    "0123456789ABCDE",                               # 15 hex, wrong length
])
def test_checksum_unrecognized_format_is_unverifiable(claimed):
    v = verify_sequence_checksum(ELVIS, claimed)
    assert v.verdict == UNVERIFIABLE
    assert v.ok is False
    assert v.evidence == []


def test_checksum_sequence_is_cleaned_before_hashing():
    messy = "  elvis\n  lives 123 *"
    md5 = hashlib.md5(ELVIS.encode()).hexdigest()
    v = verify_sequence_checksum(messy, md5)
    assert v.verdict == CONFIRMED


# ---------------------------------------------------------------------------
# verify_uniprot_sequence
# ---------------------------------------------------------------------------

def _uniprot_fetcher(tmp_path, status=200, body=HBA_FASTA):
    transport = CannedTransport([(UNIPROT_FASTA_URL, status, body)])
    return make_fetcher(tmp_path, transport), transport


def test_uniprot_sequence_confirmed(tmp_path):
    fetcher, transport = _uniprot_fetcher(tmp_path)
    v = verify_uniprot_sequence("P69905", HBA_SEQ, fetcher)
    assert v.verdict == CONFIRMED
    assert v.ok is True
    assert "142 residues" in v.detail
    assert HBA_CRC64 in v.detail
    assert len(v.evidence) == 1
    ev = v.evidence[0]
    assert isinstance(ev, Evidence)
    assert ev.url == UNIPROT_FASTA_URL
    assert ev.status == 200
    assert ev.sha256 == hashlib.sha256(HBA_FASTA.encode()).hexdigest()
    assert ev.size == len(HBA_FASTA.encode())
    assert ev.from_cache is False
    assert transport.calls == [UNIPROT_FASTA_URL]


def test_uniprot_sequence_confirmed_from_messy_pasted_sequence(tmp_path):
    fetcher, _ = _uniprot_fetcher(tmp_path)
    # lowercase, GenBank-style numbering and line breaks must not matter
    messy = "\n".join(
        f"{i + 1:>4} {HBA_SEQ[i:i + 60].lower()}"
        for i in range(0, len(HBA_SEQ), 60)
    )
    v = verify_uniprot_sequence("P69905", messy, fetcher)
    assert v.verdict == CONFIRMED


def test_uniprot_sequence_refuted_different_length(tmp_path):
    fetcher, _ = _uniprot_fetcher(tmp_path)
    v = verify_uniprot_sequence("P69905", HBA_SEQ[:100], fetcher)
    assert v.verdict == REFUTED
    assert v.ok is False
    assert "given 100 residues" in v.detail
    assert "database has 142" in v.detail
    assert len(v.evidence) == 1      # the refutation is still evidence-backed


def test_uniprot_sequence_refuted_same_length_reports_positions(tmp_path):
    fetcher, _ = _uniprot_fetcher(tmp_path)
    mutant = "G" + HBA_SEQ[1:]                  # M1G point difference
    v = verify_uniprot_sequence("P69905", mutant, fetcher)
    assert v.verdict == REFUTED
    assert "1 differing positions" in v.detail
    assert "first at 1" in v.detail


def test_uniprot_sequence_unverifiable_on_http_404(tmp_path):
    fetcher, _ = _uniprot_fetcher(tmp_path, status=404, body="Error")
    v = verify_uniprot_sequence("P69905", HBA_SEQ, fetcher)
    assert v.verdict == UNVERIFIABLE
    assert "404" in v.detail
    assert v.evidence == []


def test_uniprot_sequence_unverifiable_on_non_fasta_body(tmp_path):
    fetcher, _ = _uniprot_fetcher(tmp_path, body="<html>maintenance</html>")
    v = verify_uniprot_sequence("P69905", HBA_SEQ, fetcher)
    assert v.verdict == UNVERIFIABLE


def test_uniprot_sequence_second_check_served_from_cache(tmp_path):
    fetcher, transport = _uniprot_fetcher(tmp_path)
    first = verify_uniprot_sequence("P69905", HBA_SEQ, fetcher)
    second = verify_uniprot_sequence("P69905", HBA_SEQ, fetcher)
    assert first.evidence[0].from_cache is False
    assert second.evidence[0].from_cache is True
    assert second.verdict == CONFIRMED
    assert len(transport.calls) == 1             # one real fetch only
    assert first.evidence[0].sha256 == second.evidence[0].sha256


# ---------------------------------------------------------------------------
# verify_pdb_uniprot_link
# ---------------------------------------------------------------------------

def _mbn_fetcher(tmp_path):
    routes = [
        (RCSB_ENTRY + "1MBN", 200, _entry_json("1MBN", ["1"])),
        (RCSB_ENTITY + "1MBN/1", 200, _entity_json(
            "1MBN", "1",
            [{"database_name": "UniProt", "database_accession": "P02185"}])),
    ]
    transport = CannedTransport(routes)
    return make_fetcher(tmp_path, transport), transport


def test_pdb_uniprot_link_confirmed(tmp_path):
    fetcher, transport = _mbn_fetcher(tmp_path)
    v = verify_pdb_uniprot_link("1mbn", "P02185", fetcher)   # lowercase PDB id ok
    assert v.verdict == CONFIRMED
    assert v.ok is True
    assert "P02185" in v.detail
    assert "PDB 1MBN" in v.claim
    assert len(v.evidence) == 1                  # one polymer entity
    assert v.evidence[0].url == RCSB_ENTITY + "1MBN/1"
    # both the entry and the entity endpoints were consulted
    assert transport.calls[0] == RCSB_ENTRY + "1MBN"


def test_pdb_uniprot_link_isoform_suffix_is_stripped(tmp_path):
    fetcher, _ = _mbn_fetcher(tmp_path)
    v = verify_pdb_uniprot_link("1MBN", "P02185-2", fetcher)
    assert v.verdict == CONFIRMED
    assert "UniProt P02185" in v.claim


def test_pdb_uniprot_link_refuted_lists_actual_references(tmp_path):
    routes = [
        (RCSB_ENTRY + "4HHB", 200, _entry_json("4HHB", ["1", "2"])),
        (RCSB_ENTITY + "4HHB/1", 200, _entity_json(
            "4HHB", "1",
            [{"database_name": "UniProt", "database_accession": "P69905"}])),
        (RCSB_ENTITY + "4HHB/2", 200, _entity_json(
            "4HHB", "2",
            [{"database_name": "UniProt", "database_accession": "P68871"}])),
    ]
    fetcher = make_fetcher(tmp_path, CannedTransport(routes))
    v = verify_pdb_uniprot_link("4HHB", "P02185", fetcher)
    assert v.verdict == REFUTED
    assert v.ok is False
    assert "P69905" in v.detail and "P68871" in v.detail
    assert "not P02185" in v.detail
    assert len(v.evidence) == 2                  # one per polymer entity


def test_pdb_uniprot_link_unverifiable_without_uniprot_refs(tmp_path):
    # a DNA-only entry: entities carry no UniProt cross-reference
    routes = [
        (RCSB_ENTRY + "1BNA", 200, _entry_json("1BNA", ["1", "2"])),
        (RCSB_ENTITY + "1BNA/1", 200, _entity_json("1BNA", "1", None)),
        (RCSB_ENTITY + "1BNA/2", 200, _entity_json(
            "1BNA", "2",
            [{"database_name": "GenBank", "database_accession": "X00001"}])),
    ]
    fetcher = make_fetcher(tmp_path, CannedTransport(routes))
    v = verify_pdb_uniprot_link("1BNA", "P02185", fetcher)
    assert v.verdict == UNVERIFIABLE
    assert "no UniProt reference" in v.detail
    assert len(v.evidence) == 2


def test_pdb_uniprot_link_unverifiable_on_fetch_error(tmp_path):
    routes = [(RCSB_ENTRY + "9XYZ", 404, '{"status": 404}')]
    fetcher = make_fetcher(tmp_path, CannedTransport(routes))
    v = verify_pdb_uniprot_link("9XYZ", "P02185", fetcher)
    assert v.verdict == UNVERIFIABLE
    assert "404" in v.detail
    assert v.evidence == []


def test_pdb_uniprot_link_malformed_accession_raises(tmp_path):
    fetcher, _ = _mbn_fetcher(tmp_path)
    with pytest.raises(ValueError):
        verify_pdb_uniprot_link("1MBN", "NOTANACC", fetcher)


# ---------------------------------------------------------------------------
# verify_doi
# ---------------------------------------------------------------------------

GOOD_DOI = "10.1038/nature12373"
GOOD_DOI_BODY = json.dumps({
    "responseCode": 1,
    "handle": GOOD_DOI,
    "values": [{"index": 1, "type": "URL",
                "data": {"format": "string",
                         "value": "https://www.nature.com/articles/nature12373"}}],
})


def test_doi_confirmed(tmp_path):
    transport = CannedTransport([(DOI_HANDLES + GOOD_DOI, 200, GOOD_DOI_BODY)])
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_doi(GOOD_DOI, fetcher)
    assert v.verdict == CONFIRMED
    assert v.ok is True
    assert len(v.evidence) == 1
    assert v.evidence[0].url == DOI_HANDLES + GOOD_DOI


def test_doi_prefix_forms_are_normalized(tmp_path):
    transport = CannedTransport([(DOI_HANDLES + GOOD_DOI, 200, GOOD_DOI_BODY)])
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_doi(f"https://doi.org/{GOOD_DOI}", fetcher)
    assert v.verdict == CONFIRMED
    assert v.claim == f"DOI {GOOD_DOI} exists"
    assert transport.calls == [DOI_HANDLES + GOOD_DOI]


def test_doi_refuted_on_handle_404(tmp_path):
    fake = "10.99999/does-not-exist"
    body_404 = json.dumps({"responseCode": 100, "handle": fake})
    transport = CannedTransport([(DOI_HANDLES + fake, 404, body_404)])
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_doi(fake, fetcher)
    assert v.verdict == REFUTED
    assert v.ok is False
    assert "not found" in v.detail


def test_doi_refuted_on_nonresolving_response_code(tmp_path):
    body = json.dumps({"responseCode": 2, "handle": GOOD_DOI})
    transport = CannedTransport([(DOI_HANDLES + GOOD_DOI, 200, body)])
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_doi(GOOD_DOI, fetcher)
    assert v.verdict == REFUTED
    assert "responseCode=2" in v.detail
    assert len(v.evidence) == 1


def test_doi_unverifiable_on_network_failure_after_retries(tmp_path):
    transport = FailingTransport()
    fetcher = make_fetcher(tmp_path, transport, max_retries=2)
    v = verify_doi(GOOD_DOI, fetcher)
    assert v.verdict == UNVERIFIABLE
    assert v.evidence == []
    assert transport.calls == 3                  # initial try + 2 retries


# ---------------------------------------------------------------------------
# verify_gene_symbol
# ---------------------------------------------------------------------------

def test_gene_symbol_confirmed(tmp_path):
    transport = CannedTransport([(ESEARCH_PREFIX, 200,
                                  _esearch_json(1, ["7157"]))])
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_gene_symbol("TP53", fetcher=fetcher)
    assert v.verdict == CONFIRMED
    assert v.ok is True
    assert "1 matching record" in v.detail
    assert len(v.evidence) == 1
    # the query really was a current-symbol search in the right organism
    qs = urllib.parse.parse_qs(urllib.parse.urlsplit(transport.calls[0]).query)
    assert qs["db"] == ["gene"]
    assert qs["term"] == ["TP53[sym] AND Homo sapiens[orgn]"]
    assert qs["retmode"] == ["json"]


def test_gene_symbol_organism_is_forwarded(tmp_path):
    transport = CannedTransport([(ESEARCH_PREFIX, 200,
                                  _esearch_json(1, ["22059"]))])
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_gene_symbol("Trp53", organism="Mus musculus", fetcher=fetcher)
    assert v.verdict == CONFIRMED
    qs = urllib.parse.parse_qs(urllib.parse.urlsplit(transport.calls[0]).query)
    assert qs["term"] == ["Trp53[sym] AND Mus musculus[orgn]"]


def test_gene_symbol_refuted_on_zero_count(tmp_path):
    transport = CannedTransport([(ESEARCH_PREFIX, 200,
                                  _esearch_json(0, []))])
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_gene_symbol("NOTAGENE1", fetcher=fetcher)
    assert v.verdict == REFUTED
    assert v.ok is False
    assert "alias" in v.detail                   # points at the likely cause
    assert len(v.evidence) == 1


@pytest.mark.parametrize("symbol", ["TP53 human", "", "TP53;DROP"])
def test_gene_symbol_malformed_is_refuted_without_fetching(tmp_path, symbol):
    transport = CannedTransport([])              # any fetch would assert
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_gene_symbol(symbol, fetcher=fetcher)
    assert v.verdict == REFUTED
    assert "well-formed" in v.detail
    assert v.evidence == []
    assert transport.calls == []


def test_gene_symbol_unverifiable_on_malformed_esearch_body(tmp_path):
    transport = CannedTransport([(ESEARCH_PREFIX, 200, '{"esearchresult": {}}')])
    fetcher = make_fetcher(tmp_path, transport)
    v = verify_gene_symbol("TP53", fetcher=fetcher)
    assert v.verdict == UNVERIFIABLE


def test_gene_symbol_unverifiable_on_http_500(tmp_path):
    transport = CannedTransport([(ESEARCH_PREFIX, 500, "oops")])
    fetcher = make_fetcher(tmp_path, transport)  # max_retries=0
    v = verify_gene_symbol("TP53", fetcher=fetcher)
    assert v.verdict == UNVERIFIABLE
    assert "500" in v.detail


# ---------------------------------------------------------------------------
# Verification dataclass
# ---------------------------------------------------------------------------

def test_verification_ok_only_for_confirmed():
    assert Verification("c", CONFIRMED, "d").ok is True
    assert Verification("c", REFUTED, "d").ok is False
    assert Verification("c", UNVERIFIABLE, "d").ok is False


def test_verification_as_dict_round_trips_evidence(tmp_path):
    fetcher, _ = _uniprot_fetcher(tmp_path)
    v = verify_uniprot_sequence("P69905", HBA_SEQ, fetcher)
    d = v.as_dict()
    assert set(d) == {"claim", "verdict", "detail", "evidence"}
    assert d["verdict"] == CONFIRMED
    assert isinstance(d["evidence"], list) and len(d["evidence"]) == 1
    ev = d["evidence"][0]
    assert set(ev) == {"url", "status", "sha256", "size", "retrieved_at",
                       "from_cache", "note"}
    json.dumps(d)                                # fully JSON-serializable


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def run_cli(argv):
    """Run allonet.main(argv) capturing stdout/stderr; return (code, out, err)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = allonet.main(argv)
    return code, out.getvalue(), err.getvalue()


def test_cli_id_valid_exit_0():
    code, out, _ = run_cli(["id", "P69905"])
    assert code == 0
    assert "uniprot" in out


def test_cli_id_valid_json():
    code, out, _ = run_cli(["--json", "id", "P69905"])
    assert code == 0
    data = json.loads(out)
    assert data["input"] == "P69905"
    assert any(m["database"] == "uniprot" and m["normalized"] == "P69905"
               for m in data["matches"])


def test_cli_id_hallucinated_exit_1():
    code, out, _ = run_cli(["id", "NOT_A_REAL_ID"])
    assert code == 1
    assert "hallucinated" in out


def test_cli_id_hallucinated_json_empty_matches():
    code, out, _ = run_cli(["--json", "id", "NOT_A_REAL_ID"])
    assert code == 1
    data = json.loads(out)
    assert data["matches"] == []


def test_cli_seq_stats_json():
    seq = "ACGTACGTACGTACGTACGTACGT"
    code, out, _ = run_cli(["--json", "seq", "stats", seq])
    assert code == 0
    stats = json.loads(out)
    assert stats["length"] == 24
    assert stats["guessed_type"] == "dna"
    assert stats["gc_content"] == 0.5
    assert stats["crc64"] == "68AFDFF47FFEDFF4"   # independent bit-wise engine
    assert len(stats["seguid"]) == 27


def test_cli_seq_validate_ok_exit_0():
    code, out, _ = run_cli(["seq", "validate", "ACGT", "--kind", "dna"])
    assert code == 0
    assert "valid dna" in out


def test_cli_seq_validate_bad_exit_1():
    code, out, _ = run_cli(["seq", "validate", "ACGTJ", "--kind", "dna"])
    assert code == 1
    assert "'J'@5" in out


def test_cli_seq_validate_missing_kind_is_usage_error():
    with pytest.raises(SystemExit) as ei:
        run_cli(["seq", "validate", "ACGT"])
    assert ei.value.code == 2


def test_cli_seq_revcomp():
    code, out, _ = run_cli(["seq", "revcomp", "ATGC"])
    assert code == 0
    assert out.strip() == "GCAT"


def test_cli_seq_translate():
    code, out, _ = run_cli(["seq", "translate", "ATGGCC"])
    assert code == 0
    assert out.strip() == "MA"


def test_cli_seq_clean():
    code, out, _ = run_cli(["seq", "clean", " atg 123 gca *"])
    assert code == 0
    assert out.strip() == "ATGGCA"


def test_cli_seq_reads_fasta_file(tmp_path):
    fasta = tmp_path / "hba.fasta"
    fasta.write_text(HBA_FASTA)
    code, out, _ = run_cli(["seq", "clean", str(fasta)])
    assert code == 0
    assert out.strip() == HBA_SEQ


def test_cli_checksum_json():
    code, out, _ = run_cli(["--json", "checksum", ELVIS])
    assert code == 0
    data = json.loads(out)
    assert data["length"] == 10
    assert data["crc64_iso"] == ELVIS_CRC64
    assert data["md5"] == hashlib.md5(ELVIS.encode()).hexdigest()
    assert data["sha256"] == hashlib.sha256(ELVIS.encode()).hexdigest()
    seg = __import__("base64").b64encode(
        hashlib.sha1(ELVIS.encode()).digest()).decode().rstrip("=")
    assert data["seguid"] == seg


def test_cli_checksum_of_fasta_file(tmp_path):
    fasta = tmp_path / "hba.fasta"
    fasta.write_text(HBA_FASTA)
    code, out, _ = run_cli(["--json", "checksum", str(fasta)])
    assert code == 0
    data = json.loads(out)
    assert data["length"] == 142
    assert data["crc64_iso"] == HBA_CRC64
    assert data["md5"] == hashlib.md5(HBA_SEQ.encode()).hexdigest()


def test_cli_units_kd2dg():
    code, out, _ = run_cli(["units", "kd2dg", "1e-9"])
    assert code == 0
    assert math.isclose(float(out.strip()), KD2DG_1NM, rel_tol=1e-12)


def test_cli_units_kd2dg_json_parses_as_number():
    code, out, _ = run_cli(["--json", "units", "kd2dg", "1e-9"])
    assert code == 0
    assert math.isclose(json.loads(out), KD2DG_1NM, rel_tol=1e-12)


def test_cli_units_kd2dg_respects_temperature():
    code, out, _ = run_cli(["units", "kd2dg", "1e-9", "--temp", "310.15"])
    assert code == 0
    expected = allonet.R_KCAL_PER_MOL_K * 310.15 * math.log(1e-9)
    assert math.isclose(float(out.strip()), expected, rel_tol=1e-12)


def test_cli_units_nonpositive_kd_is_error_exit_2():
    code, out, err = run_cli(["units", "kd2dg", "0"])
    assert code == 2
    assert out == ""
    assert "error:" in err and "positive" in err


def test_cli_verify_checksum_confirmed_exit_0():
    code, out, _ = run_cli(["verify", "checksum", "ACGT", ACGT_CRC64])
    assert code == 0
    assert "CONFIRMED" in out


def test_cli_verify_checksum_refuted_exit_1():
    code, out, _ = run_cli(["verify", "checksum", "ACGT",
                            "FFFFFFFFFFFFFFFF"])
    assert code == 1
    assert "REFUTED" in out
    assert ACGT_CRC64 in out                     # the correct value is shown


def test_cli_verify_checksum_unrecognized_exit_1():
    code, out, _ = run_cli(["verify", "checksum", "ACGT", "not-a-checksum"])
    assert code == 1
    assert "UNVERIFIABLE" in out


def test_cli_verify_checksum_json():
    code, out, _ = run_cli(["--json", "verify", "checksum", "ACGT",
                            ACGT_CRC64])
    assert code == 0
    data = json.loads(out)
    assert data["verdict"] == "CONFIRMED"
    assert data["claim"] == f"sequence has checksum {ACGT_CRC64}"
    assert data["evidence"] == []


def test_cli_error_exit_2_on_untranslatable_length():
    code, _, err = run_cli(["seq", "translate", "ATGG"])
    assert code == 2
    assert "error:" in err and "multiple of 3" in err


def test_cli_error_exit_2_on_non_nucleotide_revcomp():
    code, _, err = run_cli(["seq", "revcomp", "ELVISLIVES"])
    assert code == 2
    assert "error:" in err


def test_cli_missing_command_is_usage_error():
    with pytest.raises(SystemExit) as ei:
        run_cli([])
    assert ei.value.code == 2
