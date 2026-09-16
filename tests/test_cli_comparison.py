import json

from nanodesign.cli import main
import nanodesign.method_comparison as comparison


def test_cli_passes_explicit_cc_guess_and_reports_actual_summary_path(tmp_path, monkeypatch, capsys):
    calls = []
    def run(output, directory, basis, *, cc_initial_guess):
        calls.append((output, directory, basis, cc_initial_guess))
        return {"status": "completed", "computed": {"example": 1}}
    monkeypatch.setattr(comparison, "run_method_comparison", run)
    output = tmp_path / "paired"
    assert main(["compare-methods", "--out", str(output), "--cc-initial-guess", "atom"]) == 0
    assert calls == [(str(output), None, "cc-pvdz", "atom")]
    result = json.loads(capsys.readouterr().out)
    assert result["result"] == str(output / "method_comparison.json")
    assert result["method_validated"] is False
