"""Independent acceptance checks for a defaulted DFT guess setting.

Synthetic integration fixtures only; no quantum calculations or science claims.
"""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import pytest
import nanodesign.campaign as campaign
from nanodesign.quantum import QuantumSettings
from compatibility_fixtures import fixture, file_hashes


@pytest.mark.parametrize("saved_status", [None, "completed", "failed", "not_converged"])
def test_legacy_campaign_report_preserves_evidence(tmp_path, saved_status):
    root = fixture(tmp_path, saved_status=saved_status)
    before = file_hashes(root)
    result = campaign.campaign_report(root)
    assert result["counts"] == {saved_status or "pending": 1}
    assert result["design_validated"] is False
    assert result["electronic_state_identity_verified"] is False
    assert before == file_hashes(root)


def test_completed_legacy_campaign_resume_does_not_recalculate(tmp_path, monkeypatch):
    root = fixture(tmp_path)
    before = file_hashes(root)
    monkeypatch.setattr(campaign, "run", lambda *args, **kwargs: pytest.fail("Unexpected quantum job"))
    assert campaign.run_campaign(root)["counts"] == {"completed": 1}
    after = file_hashes(root)
    assert all(after[path] == digest for path, digest in before.items())
    assert set(after) - set(before) <= {".campaign.lock"}


@pytest.mark.parametrize("missing", ["basis", "spin", "xc", "threads", "density_fit"])
def test_compatibility_does_not_fill_unrelated_missing_settings(tmp_path, missing):
    root = fixture(tmp_path, drop_setting=missing)
    with pytest.raises(ValueError, match="quantum settings"):
        campaign.campaign_report(root)


def test_unknown_saved_setting_is_rejected(tmp_path):
    root = fixture(tmp_path, extra_setting={"invented_setting": True})
    with pytest.raises(ValueError, match="quantum settings"):
        campaign.campaign_report(root)


@pytest.mark.parametrize("bad_guess", [None, False, 0, [], {}])
def test_explicit_invalid_guess_does_not_become_legacy_default(tmp_path, bad_guess):
    root = fixture(tmp_path, extra_setting={"scf_initial_guess": bad_guess})
    with pytest.raises(ValueError, match="quantum settings"):
        campaign.campaign_report(root)


@pytest.mark.skipif(not hasattr(QuantumSettings(), "scf_initial_guess"), reason="Feature not yet implemented")
@pytest.mark.parametrize("input_guess,result_guess,accepted", [
    (None, "minao", True), ("minao", None, True), ("atom", "atom", True),
    (None, "atom", False), ("atom", None, False), ("atom", "minao", False),
    ("minao", "invented", False), ("minao", "", False),
])
def test_missing_guess_only_matches_legacy_minao(tmp_path, input_guess, result_guess, accepted):
    root = fixture(tmp_path, input_guess=input_guess, result_guess=result_guess)
    before = file_hashes(root)
    if accepted:
        assert campaign.campaign_report(root)["counts"] == {"completed": 1}
    else:
        with pytest.raises(ValueError, match="quantum settings"):
            campaign.campaign_report(root)
    assert before == file_hashes(root)


def test_compatibility_does_not_bypass_result_hash_check(tmp_path):
    root = fixture(tmp_path)
    result_path = root / "runs/pose-0001/attempt-0001/result.json"
    result = json.loads(result_path.read_text())
    result["quantum_settings"]["scf_initial_guess"] = "minao"
    result_path.write_text(json.dumps(result))
    with pytest.raises(ValueError, match="changed"):
        campaign.campaign_report(root)
