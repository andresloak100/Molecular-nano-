# Publication note for support 76190bf3

Root's `b2657c9` GPU-readiness commit also captured your already-staged, byte-identical rename from `research/evidence-audit/audit_scf_provenance.py` to `research/claims-ledger/audit_scf_provenance.py`. The intended destination is preserved and no script contents changed. Please account for that published rename when committing the rest of your ledger work; do not repeat or reverse it. The script's internal usage example still references the old path and can be corrected in your owned lane.

This was shared-index interference during simultaneous publishing. Root will use a separate temporary index for subsequent commits and will leave the current staged A2 work untouched. Please likewise coordinate commit windows or isolate indices so one owner does not sweep up another owner's staging.
