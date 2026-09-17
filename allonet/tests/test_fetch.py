"""Offline tests for AlloNet's Fetcher / Evidence layer.

Every test injects a fake transport — a callable
``(url, headers, timeout) -> (status, bytes)`` — so the suite never
touches the network, regardless of environment or cache state. Cache
directories are pytest ``tmp_path`` dirs, and ``fetcher._sleep`` is
patched to record requested delays instead of sleeping, so the suite is
fast and deterministic. Live checks against the real databases are
deliberately out of scope here.
"""

import hashlib
import json
import os
import re
import sys
import time
from unittest import mock
from urllib.parse import parse_qs, urlsplit

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import allonet
from allonet import (
    Evidence,
    FetchError,
    Fetcher,
    Record,
    ensembl_get,
    ncbi_efetch,
    ncbi_esearch,
    pdb_structure,
    rcsb_entry,
    rcsb_polymer_entities,
    uniprot_fasta,
    uniprot_json,
    uniprot_sequence,
)

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/P69905.fasta"
UNIPROT_URL_2 = "https://rest.uniprot.org/uniprotkb/P02185.fasta"
RCSB_URL = "https://data.rcsb.org/rest/v1/core/entry/1MBN"


# ---------------------------------------------------------------------------
# Test doubles
# ---------------------------------------------------------------------------

class FakeTransport:
    """Scripted stand-in for the (url, headers, timeout) -> (status, bytes)
    transport. Records every call; never touches the network.

    ``script`` entries are consumed one per call: a ``(status, bytes)``
    tuple is returned, an exception instance is raised. When the script is
    exhausted, ``by_url`` (exact-URL routing) is consulted, then
    ``default``.
    """

    def __init__(self, *script, default=(200, b"ok"), by_url=None):
        self.script = list(script)
        self.default = default
        self.by_url = dict(by_url or {})
        self.calls = []

    def __call__(self, url, headers, timeout):
        self.calls.append((url, dict(headers), timeout))
        if self.script:
            step = self.script.pop(0)
            if isinstance(step, BaseException):
                raise step
            return step
        if url in self.by_url:
            return self.by_url[url]
        return self.default

    @property
    def urls(self):
        return [c[0] for c in self.calls]


def make_fetcher(tmp_path, transport=None, *, throttle=False, **kw):
    """A Fetcher wired to a tmp cache dir and a fake transport.

    ``_sleep`` is replaced by a recorder (``fetcher.sleeps``). Unless
    ``throttle=True``, the per-host minimum interval is zeroed so that
    recorded sleeps are exactly the retry-backoff sleeps.
    """
    kw.setdefault("offline", False)  # never depend on ALLONET_OFFLINE
    f = Fetcher(cache_dir=str(tmp_path / "cache"),
                transport=transport if transport is not None else FakeTransport(),
                **kw)
    f.sleeps = []
    f._sleep = f.sleeps.append
    if not throttle:
        f._min_interval = lambda host: 0.0
    return f


# ---------------------------------------------------------------------------
# Cache behavior + Evidence contents
# ---------------------------------------------------------------------------

def test_first_get_is_a_cache_miss_and_records_evidence(tmp_path):
    t = FakeTransport((200, b"payload"))
    f = make_fetcher(tmp_path, t)
    rec = f.get(UNIPROT_URL)
    assert isinstance(rec, Record)
    assert rec.content == b"payload"
    ev = rec.evidence
    assert isinstance(ev, Evidence)
    assert ev.url == UNIPROT_URL
    assert ev.status == 200
    assert ev.from_cache is False
    assert ev.size == len(b"payload")
    assert ev.sha256 == hashlib.sha256(b"payload").hexdigest()
    assert len(t.calls) == 1


def test_second_get_is_served_from_cache_without_transport(tmp_path):
    t = FakeTransport((200, b"payload"))
    f = make_fetcher(tmp_path, t)
    rec1 = f.get(UNIPROT_URL)
    rec2 = f.get(UNIPROT_URL)
    assert len(t.calls) == 1, "cache hit must not call the transport"
    assert rec2.evidence.from_cache is True
    assert rec2.content == rec1.content
    assert rec2.evidence.sha256 == rec1.evidence.sha256
    assert rec2.evidence.status == 200
    # the timestamp of the original retrieval is preserved, not re-stamped
    assert rec2.evidence.retrieved_at == rec1.evidence.retrieved_at


def test_cache_is_keyed_by_full_url(tmp_path):
    t = FakeTransport()
    f = make_fetcher(tmp_path, t)
    f.get(UNIPROT_URL)
    f.get(UNIPROT_URL_2)          # different URL: its own miss
    assert len(t.calls) == 2
    f.get(UNIPROT_URL)            # first URL again: hit
    assert len(t.calls) == 2


def test_force_refresh_bypasses_cache_and_updates_it(tmp_path):
    t = FakeTransport((200, b"version-1"), (200, b"version-2"))
    f = make_fetcher(tmp_path, t)
    assert f.get(UNIPROT_URL).content == b"version-1"
    rec = f.get(UNIPROT_URL, force_refresh=True)
    assert len(t.calls) == 2, "force_refresh must call the transport"
    assert rec.content == b"version-2"
    assert rec.evidence.from_cache is False
    # the refreshed bytes replaced the cached copy
    rec3 = f.get(UNIPROT_URL)
    assert rec3.evidence.from_cache is True
    assert rec3.content == b"version-2"


def test_evidence_sha256_is_of_exact_bytes(tmp_path):
    payload = bytes(range(256))   # binary, not valid UTF-8
    f = make_fetcher(tmp_path, FakeTransport((200, payload)))
    rec = f.get(UNIPROT_URL)
    assert rec.content == payload
    assert rec.evidence.sha256 == hashlib.sha256(payload).hexdigest()
    assert rec.evidence.size == 256
    # byte-exact through the cache as well
    rec2 = f.get(UNIPROT_URL)
    assert rec2.content == payload
    assert rec2.evidence.sha256 == rec.evidence.sha256


def test_evidence_retrieved_at_is_utc_iso8601(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport())
    ev = f.get(UNIPROT_URL).evidence
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", ev.retrieved_at)


def test_evidence_as_dict_carries_all_fields(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport((200, b"x")))
    d = f.get(UNIPROT_URL).evidence.as_dict()
    assert set(d) == {"url", "status", "sha256", "size", "retrieved_at",
                      "from_cache", "note"}
    assert d["url"] == UNIPROT_URL
    assert d["status"] == 200
    assert d["from_cache"] is False


# ---------------------------------------------------------------------------
# Record.text / Record.json
# ---------------------------------------------------------------------------

def test_record_text_decodes_utf8(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport((200, "α-helix".encode("utf-8"))))
    assert f.get(UNIPROT_URL).text == "α-helix"


def test_record_text_replaces_undecodable_bytes(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport((200, b"\xff\xfeok")))
    assert f.get(UNIPROT_URL).text == "��ok"


def test_record_json_parses_content(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport((200, b'{"a": [1, 2], "b": null}')))
    assert f.get(UNIPROT_URL).json() == {"a": [1, 2], "b": None}


def test_record_json_raises_on_non_json(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport((200, b"not json")))
    with pytest.raises(ValueError):
        f.get(UNIPROT_URL).json()


# ---------------------------------------------------------------------------
# Retries with backoff
# ---------------------------------------------------------------------------

def test_retry_on_429_then_success_backs_off(tmp_path):
    t = FakeTransport((429, b"slow down"), (200, b"done"))
    f = make_fetcher(tmp_path, t)
    rec = f.get(UNIPROT_URL)
    assert rec.content == b"done"
    assert rec.evidence.status == 200
    assert rec.evidence.from_cache is False
    assert len(t.calls) == 2
    assert f.sleeps == [1.0]


def test_retry_on_500_twice_then_success_doubles_backoff(tmp_path):
    t = FakeTransport((500, b""), (500, b""), (200, b"recovered"))
    f = make_fetcher(tmp_path, t)
    rec = f.get(UNIPROT_URL)
    assert rec.content == b"recovered"
    assert len(t.calls) == 3
    assert f.sleeps == [1.0, 2.0], "exponential backoff: 1 s then 2 s"


def test_gives_up_after_max_retries_with_status_in_error(tmp_path):
    t = FakeTransport(default=(503, b"unavailable"))
    f = make_fetcher(tmp_path, t, max_retries=2)
    with pytest.raises(FetchError) as exc:
        f.get(UNIPROT_URL)
    msg = str(exc.value)
    assert "503" in msg
    assert "unavailable" in msg          # body snippet included
    assert UNIPROT_URL in msg
    assert len(t.calls) == 3, "max_retries=2 means 3 attempts total"
    assert f.sleeps == [1.0, 2.0]


def test_error_reports_the_final_attempts_status(tmp_path):
    t = FakeTransport((500, b"first"), (429, b"rate limited"))
    f = make_fetcher(tmp_path, t, max_retries=1)
    with pytest.raises(FetchError) as exc:
        f.get(UNIPROT_URL)
    assert "429" in str(exc.value)
    assert "rate limited" in str(exc.value)
    assert len(t.calls) == 2
    assert f.sleeps == [1.0]


def test_non_retryable_status_raises_immediately(tmp_path):
    t = FakeTransport((404, b"no such entry"))
    f = make_fetcher(tmp_path, t)
    with pytest.raises(FetchError) as exc:
        f.get(UNIPROT_URL)
    assert "404" in str(exc.value)
    assert "no such entry" in str(exc.value)
    assert len(t.calls) == 1, "404 is not retryable"
    assert f.sleeps == []


def test_failed_fetch_is_not_cached(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport(default=(404, b"gone")))
    with pytest.raises(FetchError):
        f.get(UNIPROT_URL)
    assert list(f.cache_dir.glob("*.body")) == []
    assert list(f.cache_dir.glob("*.meta.json")) == []


def test_transport_oserror_is_retried_then_succeeds(tmp_path):
    t = FakeTransport(OSError("connection reset"), (200, b"recovered"))
    f = make_fetcher(tmp_path, t)
    rec = f.get(UNIPROT_URL)
    assert rec.content == b"recovered"
    assert len(t.calls) == 2
    assert f.sleeps == [1.0]


def test_transport_oserror_exhausts_retries(tmp_path):
    t = FakeTransport(OSError("boom 1"), OSError("boom 2"))
    f = make_fetcher(tmp_path, t, max_retries=1)
    with pytest.raises(FetchError) as exc:
        f.get(UNIPROT_URL)
    assert "boom 2" in str(exc.value)
    assert isinstance(exc.value.__cause__, OSError)
    assert len(t.calls) == 2
    assert f.sleeps == [1.0]


# ---------------------------------------------------------------------------
# Offline mode
# ---------------------------------------------------------------------------

def poison_transport(url, headers, timeout):
    raise AssertionError(f"offline fetcher touched the transport: {url}")


def test_offline_cache_miss_raises_without_touching_transport(tmp_path):
    f = make_fetcher(tmp_path, poison_transport, offline=True)
    with pytest.raises(FetchError) as exc:
        f.get(UNIPROT_URL)
    assert "offline" in str(exc.value)
    assert UNIPROT_URL in str(exc.value)


def test_offline_cache_hit_is_served(tmp_path):
    # prime the cache with an online fetcher on the same directory
    online = make_fetcher(tmp_path, FakeTransport((200, b"cached bytes")))
    online.get(UNIPROT_URL)
    off = make_fetcher(tmp_path, poison_transport, offline=True)
    rec = off.get(UNIPROT_URL)
    assert rec.content == b"cached bytes"
    assert rec.evidence.from_cache is True


def test_offline_force_refresh_refuses_network_even_with_cache(tmp_path):
    online = make_fetcher(tmp_path, FakeTransport((200, b"cached bytes")))
    online.get(UNIPROT_URL)
    off = make_fetcher(tmp_path, poison_transport, offline=True)
    with pytest.raises(FetchError):
        off.get(UNIPROT_URL, force_refresh=True)


def test_offline_default_comes_from_env(tmp_path):
    with mock.patch.dict(os.environ, {"ALLONET_OFFLINE": "1"}):
        f = Fetcher(cache_dir=str(tmp_path / "c1"), transport=poison_transport)
        assert f.offline is True
        with pytest.raises(FetchError):
            f.get(UNIPROT_URL)
    with mock.patch.dict(os.environ):
        os.environ.pop("ALLONET_OFFLINE", None)
        f2 = Fetcher(cache_dir=str(tmp_path / "c2"), transport=FakeTransport())
        assert f2.offline is False


def test_cache_dir_default_comes_from_env(tmp_path):
    envdir = tmp_path / "envcache"
    with mock.patch.dict(os.environ, {"ALLONET_CACHE": str(envdir)}):
        f = Fetcher(transport=FakeTransport(), offline=False)
    assert f.cache_dir == envdir
    assert envdir.is_dir()


# ---------------------------------------------------------------------------
# URL policy: https only, host allowlist
# ---------------------------------------------------------------------------

def test_http_scheme_rejected(tmp_path):
    t = FakeTransport()
    f = make_fetcher(tmp_path, t)
    with pytest.raises(ValueError) as exc:
        f.get("http://rest.uniprot.org/uniprotkb/P69905.fasta")
    assert "https" in str(exc.value)
    assert t.calls == []


def test_unknown_host_rejected(tmp_path):
    t = FakeTransport()
    f = make_fetcher(tmp_path, t)
    with pytest.raises(ValueError) as exc:
        f.get("https://evil.example.com/steal")
    assert "allowlist" in str(exc.value)
    assert "evil.example.com" in str(exc.value)
    assert t.calls == []


def test_allow_any_host_permits_unknown_hosts(tmp_path):
    t = FakeTransport((200, b"fine"))
    f = make_fetcher(tmp_path, t, allow_any_host=True)
    rec = f.get("https://evil.example.com/ok")
    assert rec.content == b"fine"
    assert t.urls == ["https://evil.example.com/ok"]
    # https-only still enforced even with allow_any_host
    with pytest.raises(ValueError):
        f.get("http://evil.example.com/ok")


def test_every_allowlisted_host_is_accepted(tmp_path):
    t = FakeTransport()
    f = make_fetcher(tmp_path, t)
    for host in sorted(Fetcher.ALLOWED_HOSTS):
        f.get(f"https://{host}/ping")
    assert len(t.calls) == len(Fetcher.ALLOWED_HOSTS)


# ---------------------------------------------------------------------------
# Per-host throttling
# ---------------------------------------------------------------------------

def test_throttle_sleeps_between_rapid_requests_to_same_host(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport(), throttle=True)
    f.get(UNIPROT_URL)
    assert f.sleeps == [], "first contact with a host never sleeps"
    # simulate the previous request having finished just now, so the
    # second request is unambiguously "rapid" even on a slow test host
    f._last_request["rest.uniprot.org"] = time.monotonic()
    f.get(UNIPROT_URL_2)
    assert len(f.sleeps) == 1
    assert 0 < f.sleeps[0] <= 0.25


def test_no_throttle_across_different_hosts(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport(), throttle=True)
    f.get(UNIPROT_URL)
    f._last_request["rest.uniprot.org"] = time.monotonic()
    f.get(RCSB_URL)   # different host: independent budget
    assert f.sleeps == []


def test_cache_hit_never_throttles(tmp_path):
    t = FakeTransport()
    f = make_fetcher(tmp_path, t, throttle=True)
    f.get(UNIPROT_URL)
    f._last_request["rest.uniprot.org"] = time.monotonic()
    rec = f.get(UNIPROT_URL)          # served from cache
    assert rec.evidence.from_cache is True
    assert len(t.calls) == 1
    assert f.sleeps == []


def test_ncbi_gets_a_slower_min_interval(tmp_path):
    f = make_fetcher(tmp_path, FakeTransport(), throttle=True)
    assert f._min_interval("eutils.ncbi.nlm.nih.gov") == 0.34
    assert f._min_interval("rest.uniprot.org") == 0.25


# ---------------------------------------------------------------------------
# Headers and timeout reach the transport
# ---------------------------------------------------------------------------

def test_default_headers_and_timeout_reach_transport(tmp_path):
    t = FakeTransport()
    f = make_fetcher(tmp_path, t, timeout=7.5)
    f.get(UNIPROT_URL, accept="application/json")
    url, headers, timeout = t.calls[0]
    assert headers["User-Agent"] == f"allonet/{allonet.__version__}"
    assert headers["Accept"] == "application/json"
    assert timeout == 7.5
    f.get(UNIPROT_URL_2)              # no accept given
    assert "Accept" not in t.calls[1][1]


def test_custom_user_agent(tmp_path):
    t = FakeTransport()
    f = make_fetcher(tmp_path, t, user_agent="acme-pipeline/2.0")
    f.get(UNIPROT_URL)
    assert t.calls[0][1]["User-Agent"] == "acme-pipeline/2.0"


# ---------------------------------------------------------------------------
# Client URL builders (through the fake transport)
# ---------------------------------------------------------------------------

def test_uniprot_fasta_builds_expected_url(tmp_path):
    t = FakeTransport((200, b">sp|P69905|\nMV\n"))
    f = make_fetcher(tmp_path, t)
    rec = uniprot_fasta("P69905", fetcher=f)
    assert t.urls == ["https://rest.uniprot.org/uniprotkb/P69905.fasta"]
    assert rec.evidence.url == "https://rest.uniprot.org/uniprotkb/P69905.fasta"


def test_uniprot_json_url_and_accept_header(tmp_path):
    t = FakeTransport((200, b"{}"))
    f = make_fetcher(tmp_path, t)
    uniprot_json("P69905", fetcher=f)
    url, headers, _ = t.calls[0]
    assert url == "https://rest.uniprot.org/uniprotkb/P69905.json"
    assert headers["Accept"] == "application/json"


def test_uniprot_fasta_rejects_malformed_accession(tmp_path):
    t = FakeTransport()
    f = make_fetcher(tmp_path, t)
    with pytest.raises(ValueError):
        uniprot_fasta("XYZ123", fetcher=f)
    assert t.calls == [], "a hallucinated accession must never reach the network"


def test_uniprot_sequence_parses_fasta_and_returns_evidence(tmp_path):
    fasta = b">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\nMVLSPADKTN\nVKAAWGKVGA\n"
    f = make_fetcher(tmp_path, FakeTransport((200, fasta)))
    seq, ev = uniprot_sequence("P69905", fetcher=f)
    assert seq == "MVLSPADKTNVKAAWGKVGA"
    assert isinstance(ev, Evidence)
    assert ev.url == "https://rest.uniprot.org/uniprotkb/P69905.fasta"
    assert ev.sha256 == hashlib.sha256(fasta).hexdigest()


def test_rcsb_entry_builds_expected_url_and_normalizes_case(tmp_path):
    t = FakeTransport((200, b"{}"))
    f = make_fetcher(tmp_path, t)
    rcsb_entry("1mbn", fetcher=f)
    url, headers, _ = t.calls[0]
    assert url == "https://data.rcsb.org/rest/v1/core/entry/1MBN"
    assert headers["Accept"] == "application/json"


def test_rcsb_polymer_entities_follows_entity_ids(tmp_path):
    entry = {"rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1", "2"]}}
    base = "https://data.rcsb.org/rest/v1/core"
    t = FakeTransport(by_url={
        f"{base}/entry/1MBN": (200, json.dumps(entry).encode()),
        f"{base}/polymer_entity/1MBN/1": (200, b'{"entity": 1}'),
        f"{base}/polymer_entity/1MBN/2": (200, b'{"entity": 2}'),
    })
    f = make_fetcher(tmp_path, t)
    records = rcsb_polymer_entities("1mbn", fetcher=f)
    assert t.urls == [f"{base}/entry/1MBN",
                      f"{base}/polymer_entity/1MBN/1",
                      f"{base}/polymer_entity/1MBN/2"]
    assert [r.json()["entity"] for r in records] == [1, 2]


def test_pdb_structure_url_and_fmt_validation(tmp_path):
    t = FakeTransport((200, b"data_1MBN"))
    f = make_fetcher(tmp_path, t)
    pdb_structure("1mbn", fmt="cif", fetcher=f)
    assert t.urls == ["https://files.rcsb.org/download/1MBN.cif"]
    with pytest.raises(ValueError):
        pdb_structure("1mbn", fmt="exe", fetcher=f)
    assert len(t.calls) == 1, "invalid fmt must be rejected before any request"


def test_ncbi_efetch_url_without_env_credentials(tmp_path):
    t = FakeTransport((200, b">seq\nACGT\n"))
    f = make_fetcher(tmp_path, t)
    with mock.patch.dict(os.environ):
        os.environ.pop("ALLONET_NCBI_EMAIL", None)
        os.environ.pop("ALLONET_NCBI_KEY", None)
        ncbi_efetch("nuccore", "NM_000518.5", fetcher=f)
    parts = urlsplit(t.urls[0])
    assert parts.scheme == "https"
    assert parts.hostname == "eutils.ncbi.nlm.nih.gov"
    assert parts.path == "/entrez/eutils/efetch.fcgi"
    q = parse_qs(parts.query)
    assert q == {"db": ["nuccore"], "id": ["NM_000518.5"],
                 "rettype": ["fasta"], "retmode": ["text"],
                 "tool": [f"allonet-{allonet.__version__}"]}
    assert "email" not in q and "api_key" not in q


def test_ncbi_efetch_url_with_env_credentials(tmp_path):
    t = FakeTransport((200, b"data"))
    f = make_fetcher(tmp_path, t)
    with mock.patch.dict(os.environ, {"ALLONET_NCBI_EMAIL": "agent@example.org",
                                      "ALLONET_NCBI_KEY": "TESTKEY123"}):
        ncbi_efetch("protein", "NP_000509.1", rettype="gb", retmode="xml",
                    fetcher=f)
    q = parse_qs(urlsplit(t.urls[0]).query)
    assert q["email"] == ["agent@example.org"]
    assert q["api_key"] == ["TESTKEY123"]
    assert q["db"] == ["protein"]
    assert q["id"] == ["NP_000509.1"]
    assert q["rettype"] == ["gb"]
    assert q["retmode"] == ["xml"]


def test_ncbi_esearch_url(tmp_path):
    t = FakeTransport((200, b'{"esearchresult": {"count": "1"}}'))
    f = make_fetcher(tmp_path, t)
    with mock.patch.dict(os.environ):
        os.environ.pop("ALLONET_NCBI_EMAIL", None)
        os.environ.pop("ALLONET_NCBI_KEY", None)
        ncbi_esearch("gene", "HBB[sym] AND Homo sapiens[orgn]", fetcher=f)
    parts = urlsplit(t.urls[0])
    assert parts.hostname == "eutils.ncbi.nlm.nih.gov"
    assert parts.path == "/entrez/eutils/esearch.fcgi"
    q = parse_qs(parts.query)
    assert q == {"db": ["gene"], "term": ["HBB[sym] AND Homo sapiens[orgn]"],
                 "retmode": ["json"],
                 "tool": [f"allonet-{allonet.__version__}"]}


def test_ensembl_get_appends_content_type(tmp_path):
    t = FakeTransport((200, b"{}"))
    f = make_fetcher(tmp_path, t)
    ensembl_get("lookup/id/ENSG00000139618", fetcher=f)
    assert t.urls == [
        "https://rest.ensembl.org/lookup/id/ENSG00000139618"
        "?content-type=application/json"
    ]
    assert t.calls[0][1]["Accept"] == "application/json"


def test_ensembl_get_rejects_query_strings_and_strips_leading_slash(tmp_path):
    t = FakeTransport((200, b"{}"), (200, b"{}"))
    f = make_fetcher(tmp_path, t)
    # Query strings, fragments and dot segments in the path are request
    # rewriting, not IDs — refused since the injection audit.
    with pytest.raises(ValueError):
        ensembl_get("/overlap/region/human/7?feature=gene", fetcher=f)
    with pytest.raises(ValueError):
        ensembl_get("lookup/../secrets", fetcher=f)
    with pytest.raises(ValueError):
        ensembl_get("lookup//id", fetcher=f)
    ensembl_get("/sequence/id/ENST00000288602", fetcher=f)
    assert t.urls[0] == (
        "https://rest.ensembl.org/sequence/id/ENST00000288602"
        "?content-type=application/json"
    )
