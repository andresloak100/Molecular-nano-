"""Tests for the programmable-assembly fidelity model."""

import math

import pytest

from assembler_model import (
    Assembler,
    AssemblyProgram,
    fidelity_from_discrimination,
    product_yield,
    required_discrimination,
    required_fidelity,
    rt,
    ROOM_TEMPERATURE_K,
)


def test_rt_room_temperature():
    assert rt(ROOM_TEMPERATURE_K) == pytest.approx(0.59248, abs=1e-4)


def test_zero_discrimination_is_coin_flip():
    assert fidelity_from_discrimination(0.0) == pytest.approx(0.5)


def test_fidelity_increases_with_discrimination():
    f1 = fidelity_from_discrimination(1.0)
    f3 = fidelity_from_discrimination(3.0)
    assert 0.5 < f1 < f3 < 1.0


def test_proofreading_reduces_error_by_power():
    ddg = 2.0
    f0 = fidelity_from_discrimination(ddg, proofreading_stages=0)
    f1 = fidelity_from_discrimination(ddg, proofreading_stages=1)
    e0, e1 = 1 - f0, 1 - f1
    # one proofreading stage squares the error
    assert e1 == pytest.approx(e0 ** 2, rel=1e-9)


def test_yield_compounds():
    assert product_yield(100, 0.99) == pytest.approx(0.99 ** 100)
    # a 100-mer at 99% per step is mostly garbage
    assert product_yield(100, 0.99) < 0.4


def test_required_fidelity_roundtrip():
    f = required_fidelity(200, 0.9)
    assert product_yield(200, f) == pytest.approx(0.9, rel=1e-9)


def test_required_discrimination_roundtrip_no_proofreading():
    n, y = 100, 0.9
    ddg = required_discrimination(n, y)
    f = fidelity_from_discrimination(ddg)
    assert product_yield(n, f) == pytest.approx(y, rel=1e-6)


def test_required_discrimination_roundtrip_with_proofreading():
    n, y, stages = 100, 0.9, 1
    ddg = required_discrimination(n, y, proofreading_stages=stages)
    f = fidelity_from_discrimination(ddg, proofreading_stages=stages)
    assert product_yield(n, f) == pytest.approx(y, rel=1e-6)


def test_proofreading_lowers_required_discrimination():
    n, y = 100, 0.9
    ddg0 = required_discrimination(n, y, proofreading_stages=0)
    ddg1 = required_discrimination(n, y, proofreading_stages=1)
    assert ddg1 < ddg0  # proofreading buys accuracy, so chemistry can be less selective


def test_assembler_run_reports_defects():
    asm = Assembler(per_step_discrimination_kcal=4.0, proofreading_stages=0)
    prog = AssemblyProgram(blocks=tuple(range(50)))
    out = asm.run(prog)
    assert out["n_steps"] == 50
    assert 0.0 < out["correct_full_length_yield"] <= 1.0
    assert out["expected_defects_per_assembly"] == pytest.approx(50 * out["per_step_error"])


def test_a1_style_input_two_kcal_is_insufficient_for_long_assembly():
    # A site-selectivity DeltaDeltaG around 2 kcal/mol (the scale A1 probes)
    # gives decent single-step fidelity but poor long-chain yield without
    # proofreading -- the decision-relevant point this model exists to make.
    asm = Assembler(per_step_discrimination_kcal=2.0)
    assert asm.step_fidelity() > 0.9          # fine for one step
    assert product_yield(100, asm.step_fidelity()) < 0.1   # hopeless for 100


def test_invalid_inputs():
    with pytest.raises(ValueError):
        product_yield(-1, 0.5)
    with pytest.raises(ValueError):
        product_yield(10, 1.5)
    with pytest.raises(ValueError):
        required_fidelity(0, 0.9)
    with pytest.raises(ValueError):
        rt(0)


# --- discrimination ledger reader ---
import json as _json
import tempfile as _tempfile
import os as _os
from discrimination_ledger import read_ledger, classify


def _write(d):
    p = _tempfile.mktemp(suffix=".json")
    open(p, "w").write(_json.dumps(d))
    return p


def test_pending_row_is_pending():
    assert classify({"competitor": "x"})["status"] == "PENDING"


def test_resolved_row_when_gap_exceeds_uncertainty():
    v = classify({"electronic_ddg_kcal": 3.0, "uncertainty_kcal": 0.5})
    assert v["status"] == "RESOLVED" and v["intended_wins"] is True


def test_unresolved_when_uncertainty_swamps_gap():
    v = classify({"electronic_ddg_kcal": 0.3, "uncertainty_kcal": 2.0})
    assert v["status"] == "UNRESOLVED"


def test_differential_correction_folds_into_gap():
    v = classify({"electronic_ddg_kcal": 3.0, "differential_correction_kcal": -1.0,
                  "uncertainty_kcal": 0.5})
    assert v["effective_ddg_kcal"] == pytest.approx(2.0)


def test_ledger_reports_pending_when_empty():
    p = _write({"competitors": [{"competitor": "a"}]})
    r = read_ledger(p)
    assert r["all_resolved"] is False
    _os.unlink(p)


def test_ledger_flags_losing_competitor():
    p = _write({"competitors": [
        {"competitor": "loser", "electronic_ddg_kcal": -2.0, "uncertainty_kcal": 0.3}]})
    r = read_ledger(p)
    assert "loser" in r["intended_loses_to"]
    _os.unlink(p)


def test_ledger_limiting_gap_drives_yield():
    p = _write({"competitors": [
        {"competitor": "a", "electronic_ddg_kcal": 5.0, "uncertainty_kcal": 0.5},
        {"competitor": "b", "electronic_ddg_kcal": 3.0, "uncertainty_kcal": 0.5}]})
    r = read_ledger(p, program_length=50)
    assert r["limiting_discrimination_kcal"] == 3.0  # the smaller winning gap binds
    assert 0.0 < r["assembler_reading"]["correct_full_length_yield"] < 1.0
    _os.unlink(p)
