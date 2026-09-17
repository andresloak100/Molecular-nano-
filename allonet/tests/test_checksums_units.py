"""Offline tests for AlloNet's checksum and unit-conversion layers.

Covers section 2b (crc64_iso, seguid, seq_sha256, seq_md5, and the offline
verify_sequence_checksum built on them) and section 5 (kcal/kJ conversions,
Kd <-> delta-G, R_KCAL_PER_MOL_K). Every expectation is computed either from
a documented external vector or independently with hashlib/math inside the
test, never by calling the library twice. No network access anywhere.
"""

import base64
import hashlib
import math
import string

import pytest

import allonet
from allonet import (
    R_KCAL_PER_MOL_K,
    crc64_iso,
    delta_g_to_kd,
    kcal_to_kj,
    kd_to_delta_g,
    kj_to_kcal,
    seguid,
    seq_md5,
    seq_sha256,
    verify_sequence_checksum,
)

# A realistic protein fragment (start of human hemoglobin alpha, P69905)
# and a nucleotide sequence, used across the checksum tests.
PROT = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF"
DNA = "ACGTACGTACGT"


@pytest.fixture(params=[PROT, DNA, "IHATEMATH", "A"])
def any_seq(request):
    """A cleaned, uppercase sequence with no noise to strip."""
    return request.param


def messy(seq: str) -> str:
    """The same sequence as a copy-paste would carry it: lowercased,
    GenBank-style numbering, line breaks, internal whitespace."""
    low = seq.lower()
    mid = len(low) // 2
    return f"  1 {low[:mid]}\n 61\t{low[mid:]}  \n"


# ---------------------------------------------------------------------------
# CRC64-ISO (the checksum UniProt publishes)
# ---------------------------------------------------------------------------

class TestCrc64Iso:
    def test_documented_spcrc_vector(self):
        # The canonical test vector shipped with the SWISS-PROT SPcrc tool.
        assert crc64_iso("IHATEMATH") == "E3DCADD69B01ADD1"

    def test_returns_16_uppercase_hex(self, any_seq):
        out = crc64_iso(any_seq)
        assert len(out) == 16
        assert set(out) <= set(string.digits + "ABCDEF")

    def test_case_invariance(self, any_seq):
        assert crc64_iso(any_seq.lower()) == crc64_iso(any_seq)
        # Mixed case too.
        mixed = "".join(
            c.lower() if i % 2 else c for i, c in enumerate(any_seq)
        )
        assert crc64_iso(mixed) == crc64_iso(any_seq)

    def test_whitespace_and_numbering_invariance(self, any_seq):
        assert crc64_iso(messy(any_seq)) == crc64_iso(any_seq)

    def test_trailing_stop_marker_stripped(self):
        # UniProt sequences carry no trailing '*'; cleaning removes it so
        # a translated ORF checksums identically to the database entry.
        assert crc64_iso(PROT + "*") == crc64_iso(PROT)

    def test_empty_sequence_is_zero(self):
        # SPcrc initializes the register to 0, so no input leaves it 0.
        assert crc64_iso("") == "0" * 16

    def test_distinguishes_near_identical_sequences(self):
        assert crc64_iso(PROT) != crc64_iso(PROT[:-1] + "Y")
        assert crc64_iso("ACGT") != crc64_iso("ACGA")

    def test_matches_independent_bitwise_implementation(self, any_seq):
        # Bit-at-a-time CRC64-ISO (reflected polynomial 0xD800000000000000,
        # init 0, no final xor) written here from the ISO 3309 definition —
        # an oracle independent of the library's table-driven code.
        crc = 0
        for b in any_seq.encode("ascii"):
            crc ^= b
            for _ in range(8):
                crc = (crc >> 1) ^ 0xD800000000000000 if crc & 1 else crc >> 1
        assert crc64_iso(any_seq) == f"{crc:016X}"


# ---------------------------------------------------------------------------
# SEGUID
# ---------------------------------------------------------------------------

class TestSeguid:
    def test_is_27_chars_of_base64(self, any_seq):
        out = seguid(any_seq)
        assert len(out) == 27
        assert set(out) <= set(
            string.ascii_letters + string.digits + "+/"
        )
        assert not out.endswith("=")  # padding must be stripped

    def test_matches_manual_base64_sha1(self, any_seq):
        digest = hashlib.sha1(any_seq.encode("ascii")).digest()
        expected = base64.b64encode(digest).decode("ascii").rstrip("=")
        assert seguid(any_seq) == expected

    def test_biopython_documented_vector(self):
        # Vector from Biopython's Bio.SeqUtils.CheckSum.seguid docs.
        assert seguid("ACGTACGTACGT") == "If6HIvcnRSQDVNiAoefAzySc6i4"

    def test_cleaning_invariance(self, any_seq):
        assert seguid(messy(any_seq)) == seguid(any_seq)

    def test_distinguishes_sequences(self):
        assert seguid(PROT) != seguid(DNA)


# ---------------------------------------------------------------------------
# SHA-256 / MD5
# ---------------------------------------------------------------------------

class TestPlainDigests:
    def test_sha256_matches_hashlib(self, any_seq):
        expected = hashlib.sha256(any_seq.encode("ascii")).hexdigest()
        assert seq_sha256(any_seq) == expected

    def test_md5_matches_hashlib(self, any_seq):
        expected = hashlib.md5(any_seq.encode("ascii")).hexdigest()
        assert seq_md5(any_seq) == expected

    def test_digest_lengths(self):
        assert len(seq_sha256(PROT)) == 64
        assert len(seq_md5(PROT)) == 32

    def test_cleaning_invariance(self, any_seq):
        assert seq_sha256(messy(any_seq)) == seq_sha256(any_seq)
        assert seq_md5(messy(any_seq)) == seq_md5(any_seq)


# ---------------------------------------------------------------------------
# verify_sequence_checksum — the offline verification built on the above
# ---------------------------------------------------------------------------

class TestVerifySequenceChecksum:
    def test_crc64_confirmed_and_case_insensitive_claim(self):
        v = verify_sequence_checksum("IHATEMATH", "e3dcadd69b01add1")
        assert v.verdict == allonet.CONFIRMED
        assert v.ok

    def test_crc64_refuted(self):
        v = verify_sequence_checksum("IHATEMATH", "0" * 16)
        assert v.verdict == allonet.REFUTED
        assert not v.ok

    def test_md5_and_sha256_detected_by_length(self):
        ok_md5 = verify_sequence_checksum(PROT, seq_md5(PROT).upper())
        ok_sha = verify_sequence_checksum(PROT, seq_sha256(PROT))
        assert ok_md5.verdict == allonet.CONFIRMED
        assert ok_sha.verdict == allonet.CONFIRMED
        bad_sha = verify_sequence_checksum(PROT, "f" * 64)
        assert bad_sha.verdict == allonet.REFUTED

    def test_seguid_detected(self):
        v = verify_sequence_checksum(DNA, "If6HIvcnRSQDVNiAoefAzySc6i4")
        assert v.verdict == allonet.CONFIRMED

    def test_unrecognized_format_is_unverifiable(self):
        v = verify_sequence_checksum(PROT, "not-a-checksum!")
        assert v.verdict == allonet.UNVERIFIABLE


# ---------------------------------------------------------------------------
# Units: kcal <-> kJ
# ---------------------------------------------------------------------------

class TestCalorieJoule:
    def test_thermochemical_factor_is_exactly_4_184(self):
        assert kcal_to_kj(1.0) == 4.184
        assert kcal_to_kj(0.0) == 0.0
        assert kcal_to_kj(-2.0) == -8.368
        assert kj_to_kcal(4.184) == 1.0

    def test_known_conversion(self):
        # ATP hydrolysis, textbook value: -7.3 kcal/mol = -30.5 kJ/mol.
        assert kcal_to_kj(-7.3) == pytest.approx(-30.5432, rel=1e-12)

    @pytest.mark.parametrize("x", [1.0, -12.278, 0.5924, 1e6, -1e-6])
    def test_roundtrip(self, x):
        assert kj_to_kcal(kcal_to_kj(x)) == pytest.approx(x, rel=1e-14)
        assert kcal_to_kj(kj_to_kcal(x)) == pytest.approx(x, rel=1e-14)


# ---------------------------------------------------------------------------
# Units: Kd <-> delta-G and the gas constant
# ---------------------------------------------------------------------------

class TestBindingThermodynamics:
    def test_gas_constant_matches_codata(self):
        # R = 8.31446261815324 J/(mol K) exactly (2019 SI); one
        # thermochemical kcal = 4184 J exactly.
        assert R_KCAL_PER_MOL_K == pytest.approx(
            8.31446261815324 / 4184, rel=1e-15
        )

    def test_one_nanomolar_at_room_temperature(self):
        dg = kd_to_delta_g(1e-9, 298.15)
        assert dg == pytest.approx(-12.28, abs=0.01)
        # And exactly RT ln(Kd), computed independently here.
        assert dg == pytest.approx(
            R_KCAL_PER_MOL_K * 298.15 * math.log(1e-9), rel=1e-15
        )

    def test_default_temperature_is_298_15(self):
        assert kd_to_delta_g(1e-9) == kd_to_delta_g(1e-9, 298.15)
        assert delta_g_to_kd(-12.28) == delta_g_to_kd(-12.28, 298.15)

    def test_standard_state_anchor(self):
        # Kd = 1 M is the standard state: delta-G exactly zero.
        assert kd_to_delta_g(1.0) == 0.0
        assert delta_g_to_kd(0.0) == 1.0

    def test_signs(self):
        assert kd_to_delta_g(1e-6) < 0        # binder: favorable
        assert kd_to_delta_g(10.0) > 0        # sub-standard-state affinity
        # Tighter binding is more negative.
        assert kd_to_delta_g(1e-12) < kd_to_delta_g(1e-9)

    @pytest.mark.parametrize("kd", [1e-12, 1e-9, 1e-6, 1e-3, 1.0, 5.0])
    @pytest.mark.parametrize("temp", [277.0, 298.15, 310.0])
    def test_kd_roundtrip(self, kd, temp):
        assert delta_g_to_kd(kd_to_delta_g(kd, temp), temp) == pytest.approx(
            kd, rel=1e-12
        )

    @pytest.mark.parametrize("dg", [-18.0, -12.28, -5.0, 0.0, 3.0])
    def test_delta_g_roundtrip(self, dg):
        assert kd_to_delta_g(delta_g_to_kd(dg, 310.0), 310.0) == pytest.approx(
            dg, abs=1e-12
        )

    @pytest.mark.parametrize("bad_kd", [0.0, -1e-9, -1.0])
    def test_nonpositive_kd_raises(self, bad_kd):
        with pytest.raises(ValueError):
            kd_to_delta_g(bad_kd, 298.15)

    @pytest.mark.parametrize("bad_temp", [0.0, -273.15, -1.0])
    def test_nonpositive_temperature_raises(self, bad_temp):
        with pytest.raises(ValueError):
            kd_to_delta_g(1e-9, bad_temp)
        with pytest.raises(ValueError):
            delta_g_to_kd(-12.28, bad_temp)
