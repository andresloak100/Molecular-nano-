"""The discrimination ledger: the concrete form of the north star.

One row per competing reaction at the reactive site. Each row states the
free-energy gap DeltaDeltaG-doubledagger between the *intended* pathway and that
competitor, the method that produced it, the electronic-state/guess controls,
any vibrational/tunnelling delta that differs between the two pathways, and an
uncertainty. The project succeeds, for design purposes, when every gap's
uncertainty is smaller than the gap itself — i.e. the ranking is resolved.

This module does two things:
  1. validates and reads a ledger (JSON), and reports for each competitor
     whether the ranking is RESOLVED (|gap| - uncertainty > 0), MARGINAL, or
     UNRESOLVED, and whether the intended pathway actually wins;
  2. converts the *worst resolved* discrimination into a per-step fidelity and
     an assembler yield via assembler_model, so the two tracks share one number.

It computes nothing chemical itself. It consumes what S1/A1/A2 and the accuracy
lanes produce. Values marked pending are not invented; the reader reports them
as pending and refuses to claim a resolved ranking.
"""

from __future__ import annotations

import json
from pathlib import Path

from assembler_model import (
    Assembler,
    AssemblyProgram,
    ROOM_TEMPERATURE_K,
    rt,
)


def _gap_delta_g(row):
    """Effective ΔΔG‡ for a row = electronic gap + differential (vib+tunnelling) delta.

    Returns None if the electronic gap is pending.
    """
    dd = row.get("electronic_ddg_kcal")
    if dd is None:
        return None
    delta = row.get("differential_correction_kcal") or 0.0
    return dd + delta


def classify(row):
    """RESOLVED / MARGINAL / UNRESOLVED / PENDING for one competitor row."""
    gap = _gap_delta_g(row)
    if gap is None:
        return {"status": "PENDING", "effective_ddg_kcal": None, "intended_wins": None}
    unc = row.get("uncertainty_kcal")
    if unc is None:
        return {"status": "PENDING_UNCERTAINTY", "effective_ddg_kcal": gap, "intended_wins": gap > 0}
    margin = abs(gap) - unc
    # Two meaningful states: the gap either clears its own uncertainty or it
    # does not. A gap smaller than its uncertainty is not a resolved ranking,
    # whatever its sign.
    status = "RESOLVED" if margin > 0 else "UNRESOLVED"
    return {"status": status, "effective_ddg_kcal": round(gap, 3),
            "uncertainty_kcal": unc, "intended_wins": gap > 0,
            "margin_kcal": round(margin, 3)}


def read_ledger(path, temperature_k=ROOM_TEMPERATURE_K, program_length=100,
                proofreading_stages=0):
    data = json.loads(Path(path).read_text())
    rows = data["competitors"]
    verdicts = {r["competitor"]: classify(r) for r in rows}

    resolved = {k: v for k, v in verdicts.items() if v["status"] == "RESOLVED"}
    losing = [k for k, v in resolved.items() if v["intended_wins"] is False]
    # The binding constraint is the smallest gap among competitors the intended
    # pathway actually beats; a competitor it does not beat is a hard failure.
    winning_gaps = [v["effective_ddg_kcal"] for v in resolved.values() if v["intended_wins"]]
    limiting = min(winning_gaps) if winning_gaps else None

    report = {
        "reactive_site": data.get("reactive_site", "unspecified"),
        "temperature_k": temperature_k,
        "rt_kcal": round(rt(temperature_k), 4),
        "competitors": verdicts,
        "all_resolved": all(v["status"] == "RESOLVED" for v in verdicts.values()),
        "intended_loses_to": losing,
        "limiting_discrimination_kcal": limiting,
    }

    if losing:
        report["assembler_reading"] = (
            "The intended pathway does not win against " + ", ".join(losing) +
            ": no positive fidelity is definable until the site or tool is changed.")
    elif limiting is None:
        report["assembler_reading"] = "No competitor resolved yet; ranking pending."
    else:
        asm = Assembler(per_step_discrimination_kcal=limiting,
                        proofreading_stages=proofreading_stages,
                        temperature_k=temperature_k)
        out = asm.run(AssemblyProgram(blocks=tuple(range(program_length))))
        report["assembler_reading"] = {
            "limiting_discrimination_kcal": limiting,
            "program_length": program_length,
            "proofreading_stages": proofreading_stages,
            "per_step_fidelity": out["per_step_fidelity"],
            "correct_full_length_yield": out["correct_full_length_yield"],
            "note": "Uses the smallest RESOLVED winning gap as the per-step "
                    "discrimination. Pending competitors could lower it.",
        }
    return report


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).parent / "discrimination-ledger.json")
    print(json.dumps(read_ledger(path), indent=1))
