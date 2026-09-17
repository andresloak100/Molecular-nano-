#!/usr/bin/env python3
"""Live contract check: exercise AlloNet against the real APIs.

Not part of the offline test suite. Run it manually (or in a scheduled CI
job) to detect upstream API changes:

    python scripts/live_check.py

Exits non-zero if any contract fails.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import allonet  # noqa: E402


def main() -> int:
    failures: list[str] = []
    with tempfile.TemporaryDirectory() as tmp:
        f = allonet.Fetcher(cache_dir=tmp)

        def check(name: str, cond: bool, detail: str = "") -> None:
            status = "ok  " if cond else "FAIL"
            print(f"[{status}] {name}" + (f" — {detail}" if detail else ""))
            if not cond:
                failures.append(name)

        # UniProt: fetched sequence must match UniProt's own CRC64
        info = allonet.uniprot_json("P69905", f).json()["sequence"]
        seq, ev = allonet.uniprot_sequence("P69905", f)
        check("uniprot length", info["length"] == len(seq) == 142)
        check("uniprot crc64", allonet.crc64_iso(seq) == info["crc64"])
        check("uniprot verify CONFIRMED",
              allonet.verify_uniprot_sequence("P69905", seq, f).ok)
        mutated = "A" + seq[1:]
        check("uniprot verify REFUTED",
              allonet.verify_uniprot_sequence("P69905", mutated, f).verdict
              == allonet.REFUTED)

        # RCSB
        v = allonet.verify_pdb_uniprot_link("1MBN", "P02185", f)
        check("pdb-uniprot link 1MBN/P02185", v.ok, v.detail)
        v = allonet.verify_pdb_uniprot_link("1MBN", "P69905", f)
        check("pdb-uniprot link REFUTED", v.verdict == allonet.REFUTED)
        cif = allonet.pdb_structure("1MBN", "cif", f)
        check("pdb structure file", cif.text.startswith("data_"))

        # NCBI
        rec = allonet.ncbi_efetch("nuccore", "NM_000518.5", "fasta",
                                  fetcher=f)
        check("ncbi efetch fasta", len(allonet.parse_fasta(rec.text)) == 1)
        check("gene TP53 CONFIRMED",
              allonet.verify_gene_symbol("TP53", fetcher=f).ok)
        check("gene bogus REFUTED",
              allonet.verify_gene_symbol("NOTAREALGENE123", fetcher=f).verdict
              == allonet.REFUTED)

        # Ensembl
        data = allonet.ensembl_get("lookup/id/ENSG00000139618", f).json()
        check("ensembl BRCA2", data.get("display_name") == "BRCA2")

        # DOI
        check("doi real CONFIRMED",
              allonet.verify_doi("10.1093/nar/gkac1052", f).ok)
        check("doi bogus REFUTED",
              allonet.verify_doi("10.99999/definitely-not-real-xyz", f).verdict
              == allonet.REFUTED)

        # Cache round-trip
        again = allonet.uniprot_fasta("P69905", f)
        check("cache hit", again.evidence.from_cache)

    print(f"\n{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
