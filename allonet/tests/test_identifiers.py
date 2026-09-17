"""Offline tests for the AlloNet identifier layer.

Covers ``identify()``, ``is_valid()``, ``normalize()``, ``IdMatch`` and the
CAS / ORCID check-digit verification, exercising every scheme in
``allonet._ID_PATTERNS`` with both genuine and corrupted examples.

No test touches the network: the identifier layer is pure string logic.
Live database checks are deliberately out of scope here.
"""

import dataclasses
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import allonet
from allonet import IdMatch, identify, is_valid, normalize


# ---------------------------------------------------------------------------
# Example tables: at least one valid and one invalid example per scheme.
# Valid examples are real issued identifiers wherever the database is public
# (P69905 = human hemoglobin alpha, 1MBN = myoglobin, 50-78-2 = aspirin,
# 0000-0002-1825-0097 = the ORCID documentation example, ...).
# ---------------------------------------------------------------------------

VALID = [
    ("uniprot", "P69905"),
    ("uniprot", "Q9Y6K9"),
    ("uniprot", "O43526"),
    ("uniprot", "A0A024R161"),          # 10-character accession
    ("uniprot", "P69905-2"),            # isoform suffix
    ("pdb", "1MBN"),
    ("pdb", "4HHB"),
    ("pdb", "6vxx"),                    # lowercase accepted, normalized up
    ("pdb_extended", "pdb_00001abc"),
    ("pdb_extended", "PDB_00006VXX"),   # case-insensitive, normalized down
    ("refseq", "NC_000001"),
    ("refseq", "NM_000546"),            # unversioned
    ("refseq", "NM_000546.6"),          # versioned
    ("refseq", "NP_000537.3"),
    ("refseq", "WP_003131952"),
    ("refseq", "XP_011541469.1"),
    ("genbank_protein", "AAA98665"),
    ("genbank_protein", "CAA12345.1"),
    ("genbank_protein", "EAW79549"),
    ("ensembl", "ENSG00000139618"),     # human BRCA2 gene
    ("ensembl", "ENSG00000139618.15"),  # versioned
    ("ensembl", "ENST00000380152"),
    ("ensembl", "ENSP00000369497"),
    ("ensembl", "ENSMUSG00000017167"),  # mouse: species-prefixed
    ("go", "GO:0008150"),
    ("go", "GO:0003677"),
    ("ec", "1.1.1.1"),
    ("ec", "EC 1.1.1.1"),
    ("ec", "EC:2.7.11.1"),
    ("ec", "3.4.-.-"),                  # partial classification
    ("ec", "6.3.2.n1"),                 # preliminary n-number
    ("doi", "10.1038/nature12373"),
    ("doi", "doi:10.1038/nature12373"),
    ("doi", "https://doi.org/10.1093/nar/gkaa1100"),
    ("doi", "https://dx.doi.org/10.1000/182"),
    ("doi", "10.2210/pdb1MBN/pdb"),     # RCSB DOIs contain slashes
    ("pmid", "PMID:12345"),
    ("pmid", "PMID 22997153"),
    ("pmid", "PMID12345678"),
    ("pmid", "pmid:9606"),
    ("dbsnp", "rs334"),                 # sickle-cell SNP
    ("dbsnp", "rs123456789012"),        # 12 digits, upper bound
    ("chembl", "CHEMBL25"),             # aspirin
    ("chembl", "CHEMBL1201585"),
    ("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"),   # aspirin
    ("inchikey", "XLYOFNOQVPJJNP-UHFFFAOYSA-N"),   # water
    ("cas", "50-78-2"),                 # aspirin
    ("cas", "7732-18-5"),               # water
    ("cas", "64-17-5"),                 # ethanol
    ("cas", "50-00-0"),                 # formaldehyde, check digit 0
    ("pfam", "PF00001"),
    ("pfam", "PF13649.6"),              # versioned
    ("interpro", "IPR000001"),
    ("interpro", "IPR036291"),
    ("hgnc", "HGNC:5"),
    ("hgnc", "HGNC:1100"),              # BRCA1
    ("ncbi_taxon", "taxid:9606"),
    ("ncbi_taxon", "NCBITaxon:10090"),
    ("orcid", "0000-0002-1825-0097"),
    ("orcid", "0000-0002-1694-233X"),   # X check digit
    ("orcid", "https://orcid.org/0000-0002-1825-0097"),
    ("orcid", "http://orcid.org/0000-0002-1694-233X"),
]

INVALID = [
    ("uniprot", "P6990"),               # truncated
    ("uniprot", "P699050"),             # 7 characters
    ("uniprot", "B12345"),              # shape wrong for a B accession
    ("uniprot", "AOA024R161"),          # letter O where a zero belongs
    ("pdb", "0ABC"),                    # cannot start with 0
    ("pdb", "1AB"),                     # too short
    ("pdb", "1ABCD"),                   # too long
    ("pdb_extended", "pdb_0001abc"),    # tail too short
    ("refseq", "AB_123456"),            # AB is not a RefSeq prefix
    ("refseq", "NM000546"),             # missing underscore
    ("refseq", "nm_000546"),            # lowercase
    ("genbank_protein", "AA12345"),     # two letters
    ("genbank_protein", "AAAA12345"),   # four letters
    ("genbank_protein", "AAA1234"),     # only four digits
    ("ensembl", "ENSG0000013961"),      # 10 digits
    ("ensembl", "ENSG000001396181"),    # 12 digits
    ("ensembl", "ENSMUSXYZG00000017167"),  # 6-letter species prefix
    ("go", "GO:123"),                   # too few digits
    ("go", "GO:00081501"),              # 8 digits
    ("go", "go:0008150"),               # lowercase prefix
    ("ec", "8.1.1.1"),                  # EC classes stop at 7
    ("ec", "1.1.1"),                    # three fields
    ("ec", "1.1.1.1.1"),                # five fields
    ("ec", "EC1.1.1.1"),                # prefix needs a space or colon
    ("doi", "10.123/abc"),              # registrant must be 4-9 digits
    ("doi", "11.1038/nature12373"),     # DOIs start with 10.
    ("doi", "10.1038/"),                # empty suffix
    ("pmid", "12345"),                  # bare number: prefix required
    ("pmid", "PMID:"),                  # no digits
    ("pmid", "PMID:1234567890"),        # 10 digits
    ("dbsnp", "rs"),                    # no digits
    ("dbsnp", "RS334"),                 # rs prefix is lowercase
    ("dbsnp", "rs1234567890123"),       # 13 digits
    ("chembl", "CHEMBL"),               # no digits
    ("chembl", "chembl25"),             # lowercase
    ("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYS-N"),    # 9-char middle block
    ("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYSA-NN"),  # 2-char last block
    ("inchikey", "bsynrymutxbxsq-uhfffaoysa-n"),   # lowercase
    ("cas", "50-78-3"),                 # corrupted check digit
    ("cas", "5-78-2"),                  # first group needs >= 2 digits
    ("cas", "50-782"),                  # two groups
    ("cas", "50-78-22"),                # two-digit check group
    ("pfam", "PF001"),                  # three digits
    ("pfam", "PF000001"),               # six digits
    ("pfam", "pf00001"),                # lowercase
    ("interpro", "IPR0001"),            # four digits
    ("interpro", "IPR00001"),           # five digits (GenBank-shaped)
    ("interpro", "ipr000001"),          # lowercase
    ("hgnc", "HGNC:"),                  # no number
    ("hgnc", "HGNC5"),                  # missing colon
    ("hgnc", "hgnc:5"),                 # lowercase
    ("ncbi_taxon", "9606"),             # bare number: prefix required
    ("ncbi_taxon", "taxid:"),           # no number
    ("ncbi_taxon", "taxid:12345678"),   # 8 digits
    ("ncbi_taxon", "Taxid:9606"),       # prefix is case-sensitive
    ("orcid", "0000-0002-1825-0098"),   # corrupted check digit
    ("orcid", "0000-0002-1825-009X"),   # X where 7 belongs
    ("orcid", "0000-0002-1825-97"),     # short last group
    ("orcid", "0000-0002-1694-233x"),   # lowercase x not permitted
]

# Strings an LLM might emit as an "accession" that match no scheme at all.
# (Long purely-alphabetic inventions are excluded here on purpose: those
# fall through to the sequence detector, pinned separately below.)
HALLUCINATED = [
    "P6990",              # truncated UniProt
    "Q9Y6K99",            # UniProt with an extra character
    "AOA024R161",         # UniProt with O-for-0 transcription error
    "B12345",             # plausible-looking six characters
    "0ABC",               # PDB-like but starts with zero
    "ZZ_12345",           # RefSeq-like with a prefix nobody issues
    "ENSG1234",           # Ensembl-like, wrong digit count
    "GO:123",             # GO with too few digits
    "10.42/nope",         # DOI with a 2-digit registrant
    "CHEMBL_25",          # ChEMBL with an invented underscore
    "not-an-id",
]


def _ids(pairs):
    return [f"{scheme}:{value}" for scheme, value in pairs]


# ---------------------------------------------------------------------------
# Coverage guard: the tables above must cover every published pattern.
# ---------------------------------------------------------------------------

def test_example_tables_cover_every_pattern():
    schemes = set(allonet._ID_PATTERNS)
    assert {s for s, _ in VALID} == schemes
    assert {s for s, _ in INVALID} == schemes


# ---------------------------------------------------------------------------
# The three core behaviors, across all schemes.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("scheme,value", VALID, ids=_ids(VALID))
def test_valid_example_is_recognized(scheme, value):
    assert scheme in {m.database for m in identify(value)}
    assert is_valid(value, scheme)


@pytest.mark.parametrize("scheme,value", INVALID, ids=_ids(INVALID))
def test_invalid_example_is_rejected(scheme, value):
    assert not is_valid(value, scheme)


@pytest.mark.parametrize("scheme,value", INVALID, ids=_ids(INVALID))
def test_invalid_example_fails_normalize(scheme, value):
    with pytest.raises(ValueError):
        normalize(value, scheme)


@pytest.mark.parametrize("value", HALLUCINATED)
def test_hallucinated_id_matches_nothing(value):
    assert identify(value) == []


# ---------------------------------------------------------------------------
# identify(): basics
# ---------------------------------------------------------------------------

def test_identify_empty_string_returns_empty_list():
    assert identify("") == []


def test_identify_whitespace_only_returns_empty_list():
    assert identify("   \t\n") == []


def test_identify_strips_surrounding_whitespace():
    matches = identify("  P69905\n")
    assert len(matches) == 1
    assert matches[0].database == "uniprot"
    assert matches[0].normalized == "P69905"


def test_uniprot_accession_is_a_single_unambiguous_match():
    matches = identify("P69905")
    assert len(matches) == 1
    m = matches[0]
    assert isinstance(m, IdMatch)
    assert m.database == "uniprot"
    assert m.normalized == "P69905"
    assert m.checksum_ok is None       # scheme has no check digit
    assert m.note == ""


def test_uniprot_isoform_suffix_is_kept():
    (m,) = identify("P69905-2")
    assert m.database == "uniprot"
    assert m.normalized == "P69905-2"


def test_uniprot_lowercase_is_accepted_and_normalized_upper():
    # The scheme is case-insensitive on input; canonical form is uppercase.
    assert normalize("p69905", "uniprot") == "P69905"


def test_uniprot_isoform_zero_is_rejected():
    # UniProt never issues isoform -0.
    assert not is_valid("P69905-0", "uniprot")


def test_pdb_extended_allows_eight_alphanumerics():
    # wwPDB's official regex is pdb_ + 8 alphanumerics; new-style IDs
    # (pdb_1xxxxxxx) must not trip the hallucination alarm.
    assert is_valid("pdb_abcd1234", "pdb_extended")
    assert normalize("PDB_00001ABC", "pdb_extended") == "pdb_00001abc"


def test_refseq_nz_with_embedded_insdc_letters_is_valid():
    # Most prokaryotic RefSeq genome records look like this; rejecting
    # them was the audit's critical finding.
    assert is_valid("NZ_CP007542.1", "refseq")
    assert is_valid("NZ_AAAA02000001.1", "refseq")


def test_sici_doi_with_angle_brackets_is_valid():
    doi = "10.1002/(sici)1099-1085(199805)12:6<821::aid-hyp668>3.0.co;2-y"
    assert is_valid(doi, "doi")
    # DOI names are case-insensitive; normalization case-folds.
    assert normalize(doi.upper(), "doi") == normalize(doi, "doi")


def test_ensembl_gene_tree_and_family_ids_are_valid():
    assert is_valid("ENSGT00390000003602", "ensembl")
    assert is_valid("ENSFM00250000000731", "ensembl")


def test_cas_leading_zero_garbage_is_rejected():
    # 00-00-0 passes the check-digit arithmetic but is not an issuable
    # CAS number.
    assert identify("00-00-0") == []


def test_ec_dash_only_valid_from_tail_inward():
    assert is_valid("1.2.-.-", "ec")
    assert not is_valid("1.-.3.4", "ec")


def test_orcid_outside_allocated_blocks_is_rejected():
    # Valid MOD 11-2 checksum, but not in ORCID's 0000-/0009- ISNI blocks.
    assert not is_valid("1234-5678-9012-3456", "orcid")


# ---------------------------------------------------------------------------
# Ambiguity: one string, several schemes — identify() must return them all.
# ---------------------------------------------------------------------------

def test_interpro_id_is_not_a_genbank_accession():
    # NCBI issues 3 letters + 5 or 7 digits, never 3 + 6, so IPR + 6
    # digits is unambiguously InterPro.
    dbs = {m.database for m in identify("IPR123456")}
    assert dbs == {"interpro"}
    assert not is_valid("IPR123456", "genbank_protein")


def test_five_digit_interpro_lookalike_is_only_genbank():
    dbs = {m.database for m in identify("IPR00001")}
    assert dbs == {"genbank_protein"}
    assert not is_valid("IPR00001", "interpro")


def test_versioned_six_digit_accession_matches_nothing():
    # 3 letters + 6 digits is a format NCBI never issued, versioned or not.
    assert identify("IPR123456.1") == []
    # The genuinely issued widths pass with a version.
    assert is_valid("AAA12345.1", "genbank_protein")
    assert is_valid("AXY1234567.2", "genbank_protein")


def test_bare_number_reads_as_pdb_not_taxon_or_pmid():
    # 9606 is the human taxid, but without a taxid:/NCBITaxon: prefix the
    # only scheme it satisfies is a 4-character PDB entry ID.
    matches = identify("9606")
    assert {m.database for m in matches} == {"pdb"}
    assert not is_valid("9606", "ncbi_taxon")
    assert not is_valid("9606", "pmid")


# ---------------------------------------------------------------------------
# CAS check digits
# ---------------------------------------------------------------------------

def test_cas_valid_check_digit_is_flagged_ok():
    (m,) = identify("50-78-2")
    assert m.database == "cas"
    assert m.checksum_ok is True
    assert m.note == ""


def test_cas_corrupted_check_digit_still_matches_but_is_flagged():
    (m,) = identify("50-78-3")
    assert m.database == "cas"
    assert m.checksum_ok is False
    assert "WRONG" in m.note


def test_cas_corrupted_check_digit_fails_is_valid_and_normalize():
    assert not is_valid("50-78-3", "cas")
    with pytest.raises(ValueError, match="WRONG"):
        normalize("50-78-3", "cas")


@pytest.mark.parametrize("digit", [d for d in range(10) if d != 2])
def test_cas_every_wrong_check_digit_is_rejected(digit):
    # Aspirin's true check digit is 2; the other nine must all fail.
    assert not is_valid(f"50-78-{digit}", "cas")


def test_cas_normalize_returns_the_id_unchanged():
    assert normalize("7732-18-5", "cas") == "7732-18-5"


# ---------------------------------------------------------------------------
# ORCID check digits (ISO 7064 MOD 11-2)
# ---------------------------------------------------------------------------

def test_orcid_valid_check_digit():
    (m,) = identify("0000-0002-1825-0097")
    assert m.database == "orcid"
    assert m.checksum_ok is True
    assert m.note == ""


def test_orcid_x_check_digit_is_valid():
    (m,) = identify("0000-0002-1694-233X")
    assert m.checksum_ok is True
    assert m.normalized == "0000-0002-1694-233X"


def test_orcid_url_form_normalizes_to_bare_id():
    assert (normalize("https://orcid.org/0000-0002-1825-0097", "orcid")
            == "0000-0002-1825-0097")


def test_orcid_url_form_checksum_is_checked_on_the_bare_id():
    (m,) = identify("https://orcid.org/0000-0002-1825-0098")
    assert m.database == "orcid"
    assert m.checksum_ok is False
    assert "WRONG" in m.note


def test_orcid_corrupted_check_digit_fails_is_valid_and_normalize():
    assert not is_valid("0000-0002-1825-0098", "orcid")
    with pytest.raises(ValueError, match="WRONG"):
        normalize("0000-0002-1825-0098", "orcid")


@pytest.mark.parametrize("check",
                         [str(d) for d in range(10) if d != 7] + ["X"])
def test_orcid_every_wrong_check_digit_is_rejected(check):
    # 0000-0002-1825-0097's true check digit is 7; every substitute,
    # including a wrongly placed X, must fail.
    assert not is_valid(f"0000-0002-1825-009{check}", "orcid")


# ---------------------------------------------------------------------------
# normalize(): canonical forms per scheme
# ---------------------------------------------------------------------------

def test_pdb_normalizes_to_uppercase():
    assert normalize("1mbn", "pdb") == "1MBN"
    assert normalize("6vxx", "pdb") == "6VXX"


def test_pdb_extended_normalizes_to_lowercase():
    assert normalize("PDB_00006VXX", "pdb_extended") == "pdb_00006vxx"


def test_refseq_version_is_part_of_the_normalized_id():
    assert normalize("NM_000546.6", "refseq") == "NM_000546.6"
    assert normalize("NM_000546", "refseq") == "NM_000546"


@pytest.mark.parametrize("raw", [
    "10.1038/nature12373",
    "doi:10.1038/nature12373",
    "DOI:10.1038/nature12373",
    "https://doi.org/10.1038/nature12373",
    "http://dx.doi.org/10.1038/nature12373",
])
def test_doi_prefix_forms_all_normalize_to_the_bare_doi(raw):
    assert normalize(raw, "doi") == "10.1038/nature12373"


@pytest.mark.parametrize("raw", [
    "PMID:22997153",
    "PMID 22997153",
    "PMID22997153",
    "pmid:22997153",
])
def test_pmid_prefix_forms_all_normalize_to_digits(raw):
    assert normalize(raw, "pmid") == "22997153"


def test_ec_prefix_is_stripped():
    assert normalize("EC 1.1.1.1", "ec") == "1.1.1.1"
    assert normalize("EC:2.7.11.1", "ec") == "2.7.11.1"
    assert normalize("1.1.1.1", "ec") == "1.1.1.1"


def test_taxon_prefixes_normalize_to_the_bare_number():
    assert normalize("taxid:9606", "ncbi_taxon") == "9606"
    assert normalize("NCBITaxon:9606", "ncbi_taxon") == "9606"


def test_normalize_strips_surrounding_whitespace():
    assert normalize("  Q9Y6K9 ", "uniprot") == "Q9Y6K9"


def test_normalize_raises_for_the_wrong_scheme():
    with pytest.raises(ValueError, match="not a well-formed"):
        normalize("P69905", "pdb")


def test_normalize_raises_for_garbage():
    with pytest.raises(ValueError):
        normalize("total garbage!!!", "uniprot")


# ---------------------------------------------------------------------------
# is_valid(): scheme discrimination
# ---------------------------------------------------------------------------

def test_is_valid_requires_the_matching_database():
    assert is_valid("P69905", "uniprot")
    assert not is_valid("P69905", "pdb")
    assert not is_valid("1MBN", "uniprot")


def test_is_valid_empty_and_whitespace_input():
    assert not is_valid("", "uniprot")
    assert not is_valid("   ", "pdb")


def test_is_valid_tolerates_surrounding_whitespace():
    assert is_valid("  P69905 ", "uniprot")


def test_unknown_database_name_is_false_not_an_error():
    assert not is_valid("P69905", "frobnicate")
    with pytest.raises(ValueError):
        normalize("P69905", "frobnicate")


# ---------------------------------------------------------------------------
# The sequence-not-an-identifier fallback
# ---------------------------------------------------------------------------

HBA_PEPTIDE = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF"   # 37 aa, HBA N-term


def test_protein_sequence_is_reported_as_sequence_not_identifier():
    (m,) = identify(HBA_PEPTIDE)
    assert m.database == "sequence"
    assert "not an identifier" in m.description
    assert "protein" in m.description
    assert f"{len(HBA_PEPTIDE)} residues" in m.description
    assert m.checksum_ok is None


def test_sequence_preview_is_truncated_to_30_residues():
    (m,) = identify(HBA_PEPTIDE)
    assert m.normalized == HBA_PEPTIDE[:30] + "..."


def test_short_sequence_preview_is_not_truncated():
    (m,) = identify("ACGTACGTACGTACGTACGTACGT")
    assert m.normalized == "ACGTACGTACGTACGTACGTACGT"


def test_dna_sequence_fallback():
    (m,) = identify("ACGTACGTACGTACGTACGTACGT")
    assert m.database == "sequence"
    assert "dna" in m.description
    assert "24 residues" in m.description


def test_rna_sequence_fallback():
    (m,) = identify("AUGGCCAUUGUAAUGGGCCGCUGA")
    assert m.database == "sequence"
    assert "rna" in m.description


def test_numbered_genbank_style_paste_is_still_detected_as_sequence():
    # Position numbers and internal whitespace are cleaned away first.
    (m,) = identify("1 MVLSPADKTN VKAAWGKVGA")
    assert m.database == "sequence"
    assert "protein" in m.description
    assert "20 residues" in m.description


def test_short_ambiguous_dna_is_not_claimed_as_sequence():
    # 12-mer ACGT-only: guess confidence is below the 0.9 bar, so
    # identify() stays silent instead of guessing.
    assert identify("ACGTACGTACGT") == []


def test_sub_10_residue_string_is_not_claimed_as_sequence():
    assert identify("MVLSPADKT") == []


def test_long_alphabetic_gibberish_reads_as_probable_sequence():
    # Every A-Z letter is in the extended protein alphabet, so a long
    # purely alphabetic invented "accession" is reported as a probable
    # sequence rather than an identifier — still not a validated ID for
    # any database, which is what matters for the hallucination check.
    (m,) = identify("NOTANACCESSIONATALL")
    assert m.database == "sequence"
    assert "not an identifier" in m.description
    assert not is_valid("NOTANACCESSIONATALL", "uniprot")


# ---------------------------------------------------------------------------
# IdMatch itself
# ---------------------------------------------------------------------------

def test_idmatch_is_frozen():
    m = identify("P69905")[0]
    with pytest.raises(dataclasses.FrozenInstanceError):
        m.database = "pdb"


def test_idmatch_as_dict_has_the_full_shape():
    (m,) = identify("50-78-2")
    assert m.as_dict() == {
        "database": "cas",
        "normalized": "50-78-2",
        "description": m.description,
        "checksum_ok": True,
        "note": "",
    }


def test_idmatch_as_dict_checksum_none_for_schemes_without_one():
    (m,) = identify("GO:0008150")
    d = m.as_dict()
    assert d["checksum_ok"] is None
    assert d["normalized"] == "GO:0008150"
