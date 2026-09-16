"""Audit committed JSON evidence for open-shell CC results with missing or
suspect SCF-solution provenance.

Motivation: open-shell HF in this repository has multiple converged solutions
that PySCF's stability analysis certifies as stable. For the ethynyl radical
at the reference geometry, the default `minao` guess converges to a solution
8.7 kcal/mol above the one reached from `atom`/`1e`, and CCSD(T) built on the
high solution is ~14 kcal/mol wrong. A number whose initial guess is
unrecorded cannot be audited by a reader.

Classification per file:
  OK_LABELED    every open-shell CC result carries an explicit scf_initial_guess
  OK_SCANNED    guesses were scanned per species and the selection is recorded
  UNLABELED     open-shell CC results present, no guess provenance anywhere
  ARTIFACT      an S^2 fingerprint matches a known wrong-solution value

Known UHF S^2 fingerprints at the published cc-pVDZ geometries (tolerance
+/- 0.005, chosen well inside the 0.4 gap between solutions):
  ethynyl radical: 0.7591 = high/artifact solution, 1.2212 = lowest known
  transition structure: 0.7595 = high solution, 1.2143 = lowest known

Read-only: this script never modifies evidence. Run from the repository root:
  python research/evidence-audit/audit_scf_provenance.py
Exit code 1 if any UNLABELED or ARTIFACT file is found.
"""

import json
import re
import subprocess
import sys

TOL = 0.005
ETHYNYL = {0.7591: "high/artifact solution", 1.2212: "lowest known solution"}
TRANSITION = {0.7595: "high solution", 1.2143: "lowest known solution"}


def leaves(obj, path=""):
    if isinstance(obj, dict):
        # Carry a species/source label into the path so list-indexed entries
        # (species[2]...) can still be fingerprinted by name.
        label = obj.get("species") or obj.get("source_label")
        if isinstance(label, str):
            path = f"{path}<{label}>" if path else f"<{label}>"
        for key, value in obj.items():
            sub = f"{path}.{key}" if path else key
            if isinstance(value, (dict, list)):
                yield from leaves(value, sub)
            else:
                yield sub, key, value
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            yield from leaves(value, f"{path}[{index}]")


def fingerprint(path, value):
    segments = re.split(r"[.\[\]<>]+", path.lower())
    table = None
    for segment in segments:
        if "ethynyl" in segment or segment == "cch":
            table = ETHYNYL
        if segment == "ts" or segment.endswith("_ts") or "transition" in segment:
            table = TRANSITION
    if table:
        for reference, label in table.items():
            if abs(value - reference) <= TOL:
                return label
    return None


def audit_file(path):
    with open(path) as handle:
        data = json.load(handle)

    has_cc = False
    guesses = []
    scanned = False
    open_shell_s2 = []
    for full, key, value in leaves(data):
        lowered = key.lower()
        if "ccsd" in lowered or ("cc" in lowered and "energy" in lowered):
            has_cc = True
        if lowered == "scf_initial_guess" or lowered == "selected_initial_guess":
            guesses.append((full, value))
        if "scan" in lowered and value is True:
            scanned = True
        if "scf_attempts" in full:
            scanned = True
        if (lowered.endswith("_s2") or lowered == "s2") and "expected" not in lowered \
                and "deviation" not in lowered and isinstance(value, (int, float)) \
                and value > 0.05:
            open_shell_s2.append((full, value))

    if not has_cc:
        return None

    notes = []
    for full, value in open_shell_s2:
        label = fingerprint(full, value)
        if label:
            notes.append(f"{full} = {value:.4f} -> {label}")

    artifact = any("artifact" in note or "high solution" in note for note in notes)
    recorded = [g for g in guesses if g[1]]
    if open_shell_s2 and not recorded and not scanned:
        verdict = "UNLABELED"
    elif artifact and not scanned:
        verdict = "ARTIFACT"
    elif scanned:
        verdict = "OK_SCANNED"
    else:
        verdict = "OK_LABELED"
    return verdict, recorded, notes


def main():
    tracked = subprocess.run(
        ["git", "ls-files", "*.json"], capture_output=True, text=True, check=True
    ).stdout.split()
    failures = 0
    for path in tracked:
        try:
            result = audit_file(path)
        except (OSError, json.JSONDecodeError) as error:
            print(f"UNREADABLE  {path}  ({error})")
            failures += 1
            continue
        if result is None:
            continue
        verdict, recorded, notes = result
        print(f"{verdict:<11} {path}")
        for full, value in recorded[:3]:
            print(f"            guess: {full} = {value}")
        for note in notes:
            print(f"            {note}")
        if verdict in ("UNLABELED", "ARTIFACT"):
            failures += 1
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
