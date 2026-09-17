"""Offline tests for AlloNet's sequence utilities.

Covers clean_sequence, guess_seq_type, validate_sequence,
reverse_complement, transcribe, translate, gc_content, protein_mw,
parse_fasta and write_fasta. Deterministic; never touches the network.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from allonet import (
    clean_sequence,
    gc_content,
    guess_seq_type,
    parse_fasta,
    protein_mw,
    reverse_complement,
    transcribe,
    translate,
    validate_sequence,
    write_fasta,
)

# Expasy average mass of one water; free peptide = sum(residues) + water.
WATER = 18.01524
GLYCINE_RESIDUE = 57.0519


# ---------------------------------------------------------------------------
# clean_sequence
# ---------------------------------------------------------------------------

def test_clean_removes_whitespace_and_uppercases():
    assert clean_sequence(" ac g\tt\r\nacgt ") == "ACGTACGT"


def test_clean_removes_digits():
    assert clean_sequence("A1C22G333T") == "ACGT"


def test_clean_genbank_numbered_lines():
    raw = (
        "        1 atgcgtaaag ctggtctcat\n"
        "       21 gcaatgctaa\n"
    )
    assert clean_sequence(raw) == "ATGCGTAAAGCTGGTCTCATGCAATGCTAA"


def test_clean_strips_trailing_stop_markers():
    assert clean_sequence("MKVL*") == "MKVL"
    assert clean_sequence("MKVL**") == "MKVL"
    # trailing whitespace/digits after the * are removed first, so the
    # stop marker is still recognized as trailing
    assert clean_sequence("mkvl* \n") == "MKVL"


def test_clean_keeps_internal_stop_marker():
    assert clean_sequence("MK*VL") == "MK*VL"


def test_clean_does_not_remove_letters_or_header_punctuation():
    # a pasted FASTA header must survive so validation fails loudly
    out = clean_sequence(">sp|P69905|HBA_HUMAN")
    assert out.startswith(">")
    assert "|" in out


def test_clean_empty_and_noise_only():
    assert clean_sequence("") == ""
    assert clean_sequence(" \n\t 123 456 ") == ""
    assert clean_sequence("***") == ""


# ---------------------------------------------------------------------------
# guess_seq_type
# ---------------------------------------------------------------------------

def test_guess_dna_long():
    kind, conf, note = guess_seq_type("ACGT" * 10)
    assert kind == "dna"
    assert conf == pytest.approx(0.95)
    assert note == ""


def test_guess_dna_short_is_flagged_ambiguous():
    kind, conf, note = guess_seq_type("ACGTACGTACGTACGT")  # 16 nt
    assert kind == "dna"
    assert conf == pytest.approx(0.6)
    assert "short" in note.lower()


def test_guess_dna_length_20_boundary_full_confidence():
    kind, conf, note = guess_seq_type("ACGT" * 5)  # exactly 20
    assert kind == "dna"
    assert conf == pytest.approx(0.95)
    assert note == ""


def test_guess_dna_with_n_ambiguity():
    kind, conf, _ = guess_seq_type("ACGTN" * 5)
    assert kind == "dna"
    assert conf == pytest.approx(0.95)


def test_guess_rna():
    kind, conf, note = guess_seq_type("AUGGCUUACGGAUCCUAA")
    assert kind == "rna"
    assert conf == pytest.approx(0.95)
    assert note == ""


def test_guess_protein():
    kind, conf, note = guess_seq_type("MKVLANQEFDHIRSTWYCGP")
    assert kind == "protein"
    assert conf == pytest.approx(0.95)
    assert note == ""


def test_guess_unknown_names_offending_characters():
    kind, conf, note = guess_seq_type("ACGT!?")
    assert kind == "unknown"
    assert conf == 0.0
    assert "!" in note and "?" in note


def test_guess_empty():
    assert guess_seq_type("") == ("unknown", 0.0, "empty sequence")
    # digits/whitespace-only cleans to empty too
    assert guess_seq_type(" 123 \n") == ("unknown", 0.0, "empty sequence")


def test_guess_cleans_input_first():
    kind, conf, _ = guess_seq_type("  acg t\n acgt acgt acgt acgt 99")
    assert kind == "dna"
    assert conf == pytest.approx(0.95)


# ---------------------------------------------------------------------------
# validate_sequence
# ---------------------------------------------------------------------------

def test_validate_dna_accepts_all_iupac_codes():
    ok, msg = validate_sequence("ACGTRYSWKMBDHVN", "dna")
    assert ok
    assert "15" in msg  # residue count reported


def test_validate_rna_accepts_u_rejects_t():
    ok, _ = validate_sequence("ACGURYSWKMBDHVN", "rna")
    assert ok
    ok, msg = validate_sequence("ACGT", "rna")
    assert not ok
    assert "'T'" in msg


def test_validate_dna_rejects_u():
    ok, msg = validate_sequence("ACGU", "dna")
    assert not ok
    assert "'U'@4" in msg


def test_validate_reports_1based_positions_in_cleaned_sequence():
    # raw has whitespace and a digit; positions refer to "ACGTE"
    ok, msg = validate_sequence("AC 1gtE", "dna")
    assert not ok
    assert "'E'@5" in msg


def test_validate_reports_multiple_positions_and_truncates():
    ok, msg = validate_sequence("ACGTEEEEEEE", "dna")  # 7 bad chars
    assert not ok
    assert "'E'@5" in msg
    assert "+2 more" in msg


def test_validate_protein_accepts_ambiguity_and_rare_residues():
    ok, _ = validate_sequence("ACDEFGHIKLMNPQRSTVWYBJZXUO", "protein")
    assert ok


def test_validate_protein_rejects_non_residue():
    ok, msg = validate_sequence("MKV@LF", "protein")
    assert not ok
    assert "'@'@4" in msg


def test_validate_bad_kind_raises():
    with pytest.raises(ValueError, match="kind"):
        validate_sequence("ACGT", "peptide")


def test_validate_empty_after_cleaning():
    ok, msg = validate_sequence("  1234 ", "dna")
    assert not ok
    assert "empty" in msg


# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------

def test_revcomp_basic_dna():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAACCC") == "GGGTTTT"


def test_revcomp_is_an_involution_on_unambiguous_dna():
    s = "ATGCGTTAGCCGGAT"
    assert reverse_complement(reverse_complement(s)) == s


def test_revcomp_iupac_ambiguity_codes():
    # per IUPAC: R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, N<->N
    assert reverse_complement("RYSWKMBDHVN") == "NBDHVKMWSRY"


def test_revcomp_rna_input_gives_dna_style_output():
    # U complements to A; A complements to T
    assert reverse_complement("AUG") == "CAT"
    assert reverse_complement("AUGC") == "GCAT"


def test_revcomp_cleans_input():
    assert reverse_complement(" at\n1gc ") == "GCAT"


def test_revcomp_non_nucleotide_raises():
    with pytest.raises(ValueError, match="not a nucleotide"):
        reverse_complement("EFILP")


def test_revcomp_empty_raises():
    with pytest.raises(ValueError):
        reverse_complement("")


# ---------------------------------------------------------------------------
# transcribe
# ---------------------------------------------------------------------------

def test_transcribe_basic():
    assert transcribe("ATGC") == "AUGC"
    assert transcribe("TTTT") == "UUUU"


def test_transcribe_preserves_iupac_codes():
    assert transcribe("ATNRYS") == "AUNRYS"


def test_transcribe_cleans_and_uppercases():
    assert transcribe("  at gc 12") == "AUGC"


def test_transcribe_rejects_rna_input():
    with pytest.raises(ValueError):
        transcribe("AUGC")


def test_transcribe_rejects_protein():
    with pytest.raises(ValueError):
        transcribe("EFILP")


# ---------------------------------------------------------------------------
# translate
# ---------------------------------------------------------------------------

def test_translate_standard_code_spot_checks():
    assert translate("ATG") == "M"
    assert translate("TTT") == "F"
    assert translate("TGG") == "W"
    assert translate("ATGTTT") == "MF"


def test_translate_accepts_rna():
    assert translate("AUGUUU") == "MF"


def test_translate_all_three_stops():
    assert translate("TAA") == "*"
    assert translate("TAG") == "*"
    assert translate("TGA") == "*"


def test_translate_keeps_stop_and_continues_by_default():
    assert translate("ATGTAATTT") == "M*F"


def test_translate_to_stop():
    assert translate("ATGTAATTT", to_stop=True) == "M"
    # stop in the first codon yields an empty protein
    assert translate("TAAATG", to_stop=True) == ""
    # no stop present: identical either way
    assert translate("ATGTTT", to_stop=True) == "MF"


def test_translate_non_multiple_of_three_raises():
    with pytest.raises(ValueError, match="multiple of 3"):
        translate("ATGT")


def test_translate_fully_ambiguous_codon_raises():
    with pytest.raises(ValueError, match="NNN"):
        translate("NNN")


def test_translate_third_position_n_fourfold_degenerate_resolves():
    assert translate("GGN") == "G"   # glycine: GGT/GGC/GGA/GGG
    assert translate("GCN") == "A"   # alanine
    assert translate("CTN") == "L"   # leucine (CTx block)
    assert translate("GTN") == "V"   # valine
    assert translate("TCN") == "S"   # serine (TCx block)
    assert translate("CCN") == "P"   # proline
    assert translate("ACN") == "T"   # threonine
    assert translate("CGN") == "R"   # arginine (CGx block)
    # and inside a longer frame
    assert translate("ATGGGNTTT") == "MGF"


def test_translate_third_position_n_non_degenerate_raises():
    with pytest.raises(ValueError):
        translate("AAN")  # AAT/AAC = N, AAA/AAG = K
    with pytest.raises(ValueError):
        translate("ATN")  # ATT/ATC/ATA = I, ATG = M
    with pytest.raises(ValueError):
        translate("TGN")  # C/C/*/W


def test_translate_other_ambiguity_positions_raise():
    with pytest.raises(ValueError):
        translate("NGG")
    with pytest.raises(ValueError):
        translate("GRG")


def test_translate_empty_is_empty():
    assert translate("") == ""


def test_translate_cleans_input():
    assert translate(" atg\n ttt 1") == "MF"


# ---------------------------------------------------------------------------
# gc_content
# ---------------------------------------------------------------------------

def test_gc_basic():
    assert gc_content("ACGT") == pytest.approx(0.5)
    assert gc_content("GGCC") == pytest.approx(1.0)
    assert gc_content("AATT") == pytest.approx(0.0)


def test_gc_s_counts_as_gc():
    assert gc_content("SSSS") == pytest.approx(1.0)
    assert gc_content("ATSS") == pytest.approx(0.5)


def test_gc_u_counts_in_denominator_as_at():
    assert gc_content("AUGC") == pytest.approx(0.5)


def test_gc_other_ambiguity_excluded_from_both_sides():
    # N and R drop out of numerator and denominator alike
    assert gc_content("ACGTNNNN") == pytest.approx(0.5)
    assert gc_content("GCR") == pytest.approx(1.0)


def test_gc_empty_raises():
    with pytest.raises(ValueError, match="no unambiguous bases"):
        gc_content("")


def test_gc_all_ambiguous_raises():
    with pytest.raises(ValueError, match="no unambiguous bases"):
        gc_content("NNRYKM")


def test_gc_w_counts_against_gc_symmetrically_with_s():
    # S is unambiguously G/C, W unambiguously A/T; counting S in the
    # numerator but dropping W from the denominator would bias GC upward.
    assert gc_content("SSWW") == pytest.approx(0.5)
    assert gc_content("GWWW") == pytest.approx(0.25)


def test_gc_cleans_input():
    assert gc_content(" ac gt 12") == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# protein_mw
# ---------------------------------------------------------------------------

def test_mw_free_glycine():
    # residue (57.0519) + water (18.01524) = free glycine, ~75.07 Da
    assert protein_mw("G") == pytest.approx(GLYCINE_RESIDUE + WATER, abs=1e-3)
    assert protein_mw("G") == pytest.approx(75.07, abs=0.01)


def test_mw_additivity_one_water_per_chain():
    # extending a chain adds exactly one residue mass, no extra water
    assert protein_mw("GG") - protein_mw("G") == pytest.approx(
        GLYCINE_RESIDUE, abs=1e-6)
    # joining two peptides into one costs one water
    assert protein_mw("GAVL") == pytest.approx(
        protein_mw("GA") + protein_mw("VL") - WATER, abs=1e-6)


def test_mw_selenocysteine_and_pyrrolysine_have_masses():
    assert protein_mw("U") == pytest.approx(150.0388 + WATER, abs=1e-3)
    assert protein_mw("O") == pytest.approx(237.3018 + WATER, abs=1e-3)


@pytest.mark.parametrize("bad", ["B", "J", "Z", "X"])
def test_mw_ambiguity_codes_refused(bad):
    with pytest.raises(ValueError, match="no defined mass"):
        protein_mw(f"GG{bad}G")


def test_mw_error_reports_position():
    with pytest.raises(ValueError, match="position 3"):
        protein_mw("GGXG")


def test_mw_cleans_input():
    assert protein_mw(" g 1\n") == pytest.approx(protein_mw("G"), abs=1e-9)


# ---------------------------------------------------------------------------
# parse_fasta
# ---------------------------------------------------------------------------

MULTI_FASTA = (
    ">seq1 first record\n"
    "ACGT ACGT\n"
    "acgtacgt\n"
    "\n"
    ">seq2 second record\n"
    "  10 mkvl anqe\n"
    "\n"
    ">seq3\n"
    "TTTT*\n"
)


def test_parse_multi_record_cleans_and_uppercases():
    entries = parse_fasta(MULTI_FASTA)
    assert entries == [
        ("seq1 first record", "ACGTACGTACGTACGT"),
        ("seq2 second record", "MKVLANQE"),
        ("seq3", "TTTT"),
    ]


def test_parse_tolerates_blank_lines_and_indentation():
    text = "\n\n>only one\n   ACGT\n\nACGT\n\n"
    assert parse_fasta(text) == [("only one", "ACGTACGT")]


def test_parse_data_before_header_raises():
    with pytest.raises(ValueError, match="before any"):
        parse_fasta("ACGT\n>late header\nACGT\n")


def test_parse_empty_sequence_raises_middle_record():
    with pytest.raises(ValueError, match="empty sequence"):
        parse_fasta(">a\n>b\nACGT\n")


def test_parse_empty_sequence_raises_last_record():
    with pytest.raises(ValueError, match="empty sequence"):
        parse_fasta(">a\nACGT\n>b\n")


def test_parse_sequence_of_only_noise_is_empty():
    with pytest.raises(ValueError, match="empty sequence"):
        parse_fasta(">a\n123 456\n")


def test_parse_duplicate_ids_raise_even_with_different_descriptions():
    text = ">P1 human\nACGT\n>P1 mouse\nTTTT\n"
    with pytest.raises(ValueError, match="duplicate"):
        parse_fasta(text)


def test_parse_same_description_different_ids_ok():
    text = ">P1 desc\nACGT\n>P2 desc\nTTTT\n"
    assert len(parse_fasta(text)) == 2


def test_parse_empty_header_raises():
    with pytest.raises(ValueError, match="empty FASTA header"):
        parse_fasta(">\nACGT\n")
    with pytest.raises(ValueError, match="empty FASTA header"):
        parse_fasta(">   \nACGT\n")


def test_parse_no_records_raises():
    with pytest.raises(ValueError, match="no FASTA records"):
        parse_fasta("")
    with pytest.raises(ValueError, match="no FASTA records"):
        parse_fasta("\n  \n\n")


# ---------------------------------------------------------------------------
# write_fasta and round-trip
# ---------------------------------------------------------------------------

def test_write_wraps_at_default_width_60():
    seq = "A" * 130
    out = write_fasta([("long", seq)])
    lines = out.splitlines()
    assert lines[0] == ">long"
    assert [len(l) for l in lines[1:]] == [60, 60, 10]
    assert out.endswith("\n")


def test_write_custom_width_and_exact_multiple():
    out = write_fasta([("x", "ACGTACGT")], width=4)
    assert out == ">x\nACGT\nACGT\n"


def test_write_cleans_sequence():
    out = write_fasta([("y", " ac 1gt* ")])
    assert out == ">y\nACGT\n"


def test_write_parse_round_trip():
    records = [
        ("seq1 first record", "ACGTACGTACGTACGT"),
        ("seq2 second record", "MKVLANQE"),
        ("seq3", "TTTT" * 40),  # long enough to wrap
    ]
    assert parse_fasta(write_fasta(records)) == records


def test_round_trip_preserves_full_headers_with_narrow_wrap():
    records = parse_fasta(MULTI_FASTA)
    assert parse_fasta(write_fasta(records, width=5)) == records
